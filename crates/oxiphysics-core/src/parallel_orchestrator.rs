// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parallel solver orchestration for multi-stage physics pipelines.
//!
//! This module provides a dependency-aware scheduler that organizes solver stages
//! into parallel waves via topological sorting. Stages within the same wave have
//! no mutual dependencies and can conceptually execute in parallel, while stages
//! in later waves depend on earlier ones.
//!
//! # Architecture
//!
//! The orchestrator works in three phases:
//! 1. **Registration** - stages and their dependencies are added via [`ParallelOrchestrator::add_stage`].
//! 2. **Scheduling** - [`ParallelOrchestrator::compute_schedule`] performs a topological sort
//!    to group independent stages into waves (Kahn's algorithm with cycle detection).
//! 3. **Execution** - [`ParallelOrchestrator::execute`] runs each wave sequentially,
//!    executing stages within each wave. Per-stage wall-clock timings are accumulated.
//!
//! # Future: rayon-based parallelism
//!
//! Currently stages within a wave run sequentially. When rayon is added as a
//! workspace dependency, intra-wave parallelism can be enabled behind a feature
//! flag (e.g. `parallel-rayon`) by replacing the sequential loop with
//! `rayon::scope` or `rayon::join`.

use std::fmt;
use std::time::Instant;

// ─── Error ──────────────────────────────────────────────────────────────────

/// Errors that can occur during orchestration scheduling or execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrchestratorError {
    /// A dependency cycle was detected among stages.
    CycleDetected {
        /// Human-readable description of the cycle.
        description: String,
    },
    /// A stage index referenced in a dependency is out of range.
    InvalidStageIndex {
        /// The invalid index that was referenced.
        index: usize,
        /// The total number of registered stages.
        total_stages: usize,
    },
    /// A stage index used during execution is out of range.
    ExecutionIndexOutOfRange {
        /// The invalid index.
        index: usize,
        /// The length of the stages slice provided.
        stages_len: usize,
    },
}

impl fmt::Display for OrchestratorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrchestratorError::CycleDetected { description } => {
                write!(f, "dependency cycle detected: {description}")
            }
            OrchestratorError::InvalidStageIndex {
                index,
                total_stages,
            } => {
                write!(
                    f,
                    "invalid stage index {index} (total stages: {total_stages})"
                )
            }
            OrchestratorError::ExecutionIndexOutOfRange { index, stages_len } => {
                write!(
                    f,
                    "execution index {index} out of range (stages slice length: {stages_len})"
                )
            }
        }
    }
}

impl std::error::Error for OrchestratorError {}

// ─── SolverStage trait ──────────────────────────────────────────────────────

/// A single stage in a physics solver pipeline.
///
/// Each stage performs a portion of the simulation step (e.g. broadphase collision,
/// constraint solving, integration). Stages declare an estimated cost for
/// load-balancing purposes and must be `Send + Sync` to support future
/// parallel execution.
pub trait SolverStage: Send + Sync {
    /// Returns the human-readable name of this stage.
    fn name(&self) -> &str;

    /// Advances the stage by `dt` seconds of simulation time.
    fn step(&mut self, dt: f64);

    /// Returns an estimated computational cost (arbitrary units) for load balancing.
    ///
    /// Higher values indicate more expensive stages. The orchestrator may use this
    /// to reorder stages within a wave for better load distribution.
    fn estimated_cost(&self) -> f64;
}

// ─── StageDependency ────────────────────────────────────────────────────────

/// Describes ordering constraints for a single stage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageDependency {
    /// The index of the stage this dependency record belongs to.
    pub stage_idx: usize,
    /// Indices of stages that must complete before this stage can run.
    pub depends_on: Vec<usize>,
}

// ─── PipelineSchedule ───────────────────────────────────────────────────────

/// A computed execution schedule produced by topological sorting.
///
/// Stages are grouped into *waves*. All stages within the same wave are
/// independent of each other and can conceptually run in parallel. Waves
/// themselves execute in order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PipelineSchedule {
    /// Groups of stage indices that can run in parallel within each wave.
    pub waves: Vec<Vec<usize>>,
}

impl PipelineSchedule {
    /// Returns the total number of waves.
    pub fn num_waves(&self) -> usize {
        self.waves.len()
    }

    /// Returns the total number of stages across all waves.
    pub fn num_stages(&self) -> usize {
        self.waves.iter().map(|w| w.len()).sum()
    }
}

// ─── ParallelOrchestrator ───────────────────────────────────────────────────

/// Orchestrates the execution of multiple solver stages respecting dependencies.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_core::parallel_orchestrator::{ParallelOrchestrator, SolverStage};
///
/// struct SimpleStage { name: String, cost: f64, step_count: u64 }
///
/// impl SolverStage for SimpleStage {
///     fn name(&self) -> &str { &self.name }
///     fn step(&mut self, _dt: f64) { self.step_count += 1; }
///     fn estimated_cost(&self) -> f64 { self.cost }
/// }
///
/// let mut orch = ParallelOrchestrator::new();
/// let a = orch.add_stage("broadphase", &[]);
/// let b = orch.add_stage("narrowphase", &[a]);
/// let c = orch.add_stage("solver", &[b]);
///
/// let mut stages: Vec<Box<dyn SolverStage>> = vec![
///     Box::new(SimpleStage { name: "broadphase".into(), cost: 1.0, step_count: 0 }),
///     Box::new(SimpleStage { name: "narrowphase".into(), cost: 2.0, step_count: 0 }),
///     Box::new(SimpleStage { name: "solver".into(), cost: 3.0, step_count: 0 }),
/// ];
///
/// orch.execute(&mut stages, 0.016).expect("execution failed");
/// ```
#[derive(Debug, Clone)]
pub struct ParallelOrchestrator {
    /// Registered stage names (index = stage id).
    stage_names: Vec<String>,
    /// Dependency records for each stage.
    dependencies: Vec<StageDependency>,
    /// Accumulated wall-clock timings per stage (seconds).
    timings: Vec<f64>,
}

impl ParallelOrchestrator {
    /// Creates a new empty orchestrator.
    pub fn new() -> Self {
        Self {
            stage_names: Vec::new(),
            dependencies: Vec::new(),
            timings: Vec::new(),
        }
    }

    /// Registers a new stage with the given name and dependency list.
    ///
    /// Returns the index of the newly added stage, which can be used as a
    /// dependency for later stages.
    ///
    /// # Arguments
    /// * `name` - Human-readable stage name.
    /// * `depends_on` - Indices of stages that must run before this one.
    pub fn add_stage(&mut self, name: &str, depends_on: &[usize]) -> usize {
        let idx = self.stage_names.len();
        self.stage_names.push(name.to_string());
        self.dependencies.push(StageDependency {
            stage_idx: idx,
            depends_on: depends_on.to_vec(),
        });
        self.timings.push(0.0);
        idx
    }

    /// Returns the number of registered stages.
    pub fn num_stages(&self) -> usize {
        self.stage_names.len()
    }

    /// Returns the registered stage names.
    pub fn stage_names(&self) -> &[String] {
        &self.stage_names
    }

    /// Computes a [`PipelineSchedule`] by topologically sorting stages into waves.
    ///
    /// Returns an error if the dependency graph contains a cycle or references
    /// an invalid stage index.
    pub fn compute_schedule(&self) -> Result<PipelineSchedule, OrchestratorError> {
        let waves = topological_sort(self.stage_names.len(), &self.dependencies)?;
        Ok(PipelineSchedule { waves })
    }

    /// Executes all registered stages in dependency order.
    ///
    /// Stages within the same wave are currently run sequentially. The timings
    /// for each stage are accumulated across repeated calls to `execute`.
    ///
    /// # Errors
    /// Returns an error if scheduling fails (cycle or invalid index) or if a
    /// scheduled stage index is out of range for the provided `stages` slice.
    pub fn execute(
        &mut self,
        stages: &mut [Box<dyn SolverStage>],
        dt: f64,
    ) -> Result<(), OrchestratorError> {
        let schedule = self.compute_schedule()?;

        for wave in &schedule.waves {
            // NOTE: future rayon parallelism would replace this sequential loop
            // with rayon::scope or par_iter over the wave indices, using unsafe
            // cell or split_at_mut to grant mutable access to disjoint elements.
            for &stage_idx in wave {
                if stage_idx >= stages.len() {
                    return Err(OrchestratorError::ExecutionIndexOutOfRange {
                        index: stage_idx,
                        stages_len: stages.len(),
                    });
                }
                let start = Instant::now();
                stages[stage_idx].step(dt);
                let elapsed = start.elapsed().as_secs_f64();
                if stage_idx < self.timings.len() {
                    self.timings[stage_idx] += elapsed;
                }
            }
        }

        Ok(())
    }

    /// Returns per-stage accumulated wall-clock timings in seconds.
    pub fn timings(&self) -> &[f64] {
        &self.timings
    }

    /// Returns the total accumulated wall-clock time across all stages.
    pub fn total_time(&self) -> f64 {
        self.timings.iter().sum()
    }

    /// Resets all accumulated timings to zero.
    pub fn reset_timings(&mut self) {
        for t in &mut self.timings {
            *t = 0.0;
        }
    }
}

impl Default for ParallelOrchestrator {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Topological sort ───────────────────────────────────────────────────────

/// Performs a topological sort of `n` stages using Kahn's algorithm.
///
/// Stages are grouped into waves: each wave contains stages whose dependencies
/// have all been satisfied by previous waves. This is a BFS-layered variant of
/// Kahn's algorithm.
///
/// # Errors
/// - [`OrchestratorError::InvalidStageIndex`] if any dependency index is >= `n`.
/// - [`OrchestratorError::CycleDetected`] if the graph contains a cycle.
pub fn topological_sort(
    n: usize,
    deps: &[StageDependency],
) -> Result<Vec<Vec<usize>>, OrchestratorError> {
    if n == 0 {
        return Ok(Vec::new());
    }

    // Validate all dependency indices.
    for dep in deps {
        if dep.stage_idx >= n {
            return Err(OrchestratorError::InvalidStageIndex {
                index: dep.stage_idx,
                total_stages: n,
            });
        }
        for &d in &dep.depends_on {
            if d >= n {
                return Err(OrchestratorError::InvalidStageIndex {
                    index: d,
                    total_stages: n,
                });
            }
        }
    }

    // Build adjacency list and in-degree count.
    // Edge: dependency -> stage (dependency must come first).
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut in_degree: Vec<usize> = vec![0; n];

    for dep in deps {
        for &d in &dep.depends_on {
            adjacency[d].push(dep.stage_idx);
            in_degree[dep.stage_idx] += 1;
        }
    }

    // BFS-layered Kahn's algorithm.
    let mut waves: Vec<Vec<usize>> = Vec::new();
    let mut current_wave: Vec<usize> = Vec::new();

    // Seed with all stages that have no dependencies.
    for (i, &deg) in in_degree.iter().enumerate() {
        if deg == 0 {
            current_wave.push(i);
        }
    }

    let mut processed = 0usize;

    while !current_wave.is_empty() {
        // Sort the wave for deterministic output.
        current_wave.sort_unstable();
        processed += current_wave.len();

        let mut next_wave: Vec<usize> = Vec::new();
        for &stage in &current_wave {
            for &neighbor in &adjacency[stage] {
                in_degree[neighbor] -= 1;
                if in_degree[neighbor] == 0 {
                    next_wave.push(neighbor);
                }
            }
        }

        waves.push(std::mem::take(&mut current_wave));
        current_wave = next_wave;
    }

    if processed != n {
        // Some stages were not processed => cycle exists.
        let remaining: Vec<usize> = (0..n).filter(|&i| in_degree[i] > 0).collect();
        let names: Vec<String> = remaining
            .iter()
            .filter_map(|&i| {
                deps.iter()
                    .find(|d| d.stage_idx == i)
                    .map(|_| format!("stage {i}"))
            })
            .collect();
        let description = if names.is_empty() {
            format!("cycle involving {remaining:?}")
        } else {
            format!("cycle involving: {}", names.join(", "))
        };
        return Err(OrchestratorError::CycleDetected { description });
    }

    Ok(waves)
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    /// A simple test stage that tracks how many times `step` was called
    /// and records the order of execution into a shared log.
    struct TestStage {
        stage_name: String,
        cost: f64,
        call_count: u64,
        execution_log: Arc<std::sync::Mutex<Vec<String>>>,
    }

    impl TestStage {
        fn new(name: &str, cost: f64, log: Arc<std::sync::Mutex<Vec<String>>>) -> Self {
            Self {
                stage_name: name.to_string(),
                cost,
                call_count: 0,
                execution_log: log,
            }
        }
    }

    impl SolverStage for TestStage {
        fn name(&self) -> &str {
            &self.stage_name
        }

        fn step(&mut self, _dt: f64) {
            self.call_count += 1;
            if let Ok(mut log) = self.execution_log.lock() {
                log.push(self.stage_name.clone());
            }
        }

        fn estimated_cost(&self) -> f64 {
            self.cost
        }
    }

    /// A stage that does a small busy-wait to produce measurable timing.
    struct TimedStage {
        stage_name: String,
        spin_iters: u64,
    }

    impl TimedStage {
        fn new(name: &str, spin_iters: u64) -> Self {
            Self {
                stage_name: name.to_string(),
                spin_iters,
            }
        }
    }

    impl SolverStage for TimedStage {
        fn name(&self) -> &str {
            &self.stage_name
        }

        fn step(&mut self, _dt: f64) {
            // Busy-spin to burn a deterministic amount of CPU work. `black_box`
            // on the loop variable forces the optimizer to actually iterate (it
            // cannot fold the sum into a closed form), so the measured wall-clock
            // cost scales with `spin_iters`.
            let mut acc = 0u64;
            for i in 0..self.spin_iters {
                acc = acc.wrapping_add(std::hint::black_box(i));
            }
            std::hint::black_box(acc);
        }

        fn estimated_cost(&self) -> f64 {
            self.spin_iters as f64
        }
    }

    // ── Linear pipeline: A -> B -> C ────────────────────────────────────

    #[test]
    fn test_linear_pipeline() {
        let mut orch = ParallelOrchestrator::new();
        let a = orch.add_stage("A", &[]);
        let b = orch.add_stage("B", &[a]);
        let _c = orch.add_stage("C", &[b]);

        let schedule = orch.compute_schedule().expect("scheduling should succeed");
        assert_eq!(schedule.waves.len(), 3, "linear pipeline needs 3 waves");
        assert_eq!(schedule.waves[0], vec![0]);
        assert_eq!(schedule.waves[1], vec![1]);
        assert_eq!(schedule.waves[2], vec![2]);

        // Verify execution order.
        let log = Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut stages: Vec<Box<dyn SolverStage>> = vec![
            Box::new(TestStage::new("A", 1.0, Arc::clone(&log))),
            Box::new(TestStage::new("B", 1.0, Arc::clone(&log))),
            Box::new(TestStage::new("C", 1.0, Arc::clone(&log))),
        ];

        orch.execute(&mut stages, 0.016)
            .expect("execute should succeed");

        let recorded = log.lock().expect("lock should not be poisoned");
        assert_eq!(&*recorded, &["A", "B", "C"]);
    }

    // ── Diamond dependency: A -> B, A -> C, B -> D, C -> D ──────────────

    #[test]
    fn test_diamond_dependency() {
        let mut orch = ParallelOrchestrator::new();
        let a = orch.add_stage("A", &[]);
        let b = orch.add_stage("B", &[a]);
        let c = orch.add_stage("C", &[a]);
        let _d = orch.add_stage("D", &[b, c]);

        let schedule = orch.compute_schedule().expect("scheduling should succeed");
        assert_eq!(schedule.waves.len(), 3, "diamond needs 3 waves");
        assert_eq!(schedule.waves[0], vec![0], "wave 0 has A");
        // B and C should be in the same wave (sorted).
        let mut wave1 = schedule.waves[1].clone();
        wave1.sort_unstable();
        assert_eq!(wave1, vec![1, 2], "wave 1 has B and C");
        assert_eq!(schedule.waves[2], vec![3], "wave 2 has D");
    }

    // ── Cycle detection ─────────────────────────────────────────────────

    #[test]
    fn test_cycle_detection() {
        // A -> B -> C -> A (cycle)
        let deps = vec![
            StageDependency {
                stage_idx: 0,
                depends_on: vec![2],
            },
            StageDependency {
                stage_idx: 1,
                depends_on: vec![0],
            },
            StageDependency {
                stage_idx: 2,
                depends_on: vec![1],
            },
        ];

        let result = topological_sort(3, &deps);
        assert!(result.is_err(), "cycle should produce an error");
        match result {
            Err(OrchestratorError::CycleDetected { description }) => {
                assert!(
                    description.contains("cycle"),
                    "error should mention cycle: {description}"
                );
            }
            other => panic!("expected CycleDetected, got {other:?}"),
        }
    }

    // ── Empty pipeline ──────────────────────────────────────────────────

    #[test]
    fn test_empty_pipeline() {
        let orch = ParallelOrchestrator::new();
        let schedule = orch
            .compute_schedule()
            .expect("empty schedule should succeed");
        assert!(schedule.waves.is_empty(), "empty pipeline has no waves");
        assert_eq!(schedule.num_waves(), 0);
        assert_eq!(schedule.num_stages(), 0);
    }

    // ── Single stage ────────────────────────────────────────────────────

    #[test]
    fn test_single_stage() {
        let mut orch = ParallelOrchestrator::new();
        orch.add_stage("only", &[]);

        let schedule = orch
            .compute_schedule()
            .expect("single stage should succeed");
        assert_eq!(schedule.waves.len(), 1);
        assert_eq!(schedule.waves[0], vec![0]);

        let log = Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut stages: Vec<Box<dyn SolverStage>> =
            vec![Box::new(TestStage::new("only", 5.0, Arc::clone(&log)))];

        orch.execute(&mut stages, 1.0)
            .expect("execute should succeed");

        let recorded = log.lock().expect("lock should not be poisoned");
        assert_eq!(&*recorded, &["only"]);
    }

    // ── Timing accumulation ─────────────────────────────────────────────

    #[test]
    fn test_timing_accumulation() {
        let mut orch = ParallelOrchestrator::new();
        orch.add_stage("fast", &[]);
        orch.add_stage("slow", &[]);

        // The 10× work gap is sized so the slow stage's wall time dominates
        // scheduler jitter even when the whole suite runs under heavy parallel
        // load, keeping `timings[1] > timings[0]` robust rather than flaky.
        let mut stages: Vec<Box<dyn SolverStage>> = vec![
            Box::new(TimedStage::new("fast", 5_000_000)),
            Box::new(TimedStage::new("slow", 50_000_000)),
        ];

        // Run multiple times to accumulate.
        for _ in 0..3 {
            orch.execute(&mut stages, 0.01)
                .expect("execute should succeed");
        }

        let timings = orch.timings();
        assert_eq!(timings.len(), 2);
        // Both stages should have recorded some positive time.
        assert!(
            timings[0] > 0.0,
            "fast stage should have positive timing: {}",
            timings[0]
        );
        assert!(
            timings[1] > 0.0,
            "slow stage should have positive timing: {}",
            timings[1]
        );
        // Total should equal sum.
        let total = orch.total_time();
        let sum = timings[0] + timings[1];
        assert!(
            (total - sum).abs() < 1e-15,
            "total {total} should equal sum {sum}"
        );
        // The slow stage should generally take more time than the fast one.
        // (Not a strict assertion since OS scheduling can vary, but 1000x
        // difference in iterations should be enough.)
        assert!(
            timings[1] > timings[0],
            "slow stage ({}) should take longer than fast stage ({})",
            timings[1],
            timings[0]
        );
    }

    // ── Invalid stage index ─────────────────────────────────────────────

    #[test]
    fn test_invalid_stage_index() {
        let deps = vec![StageDependency {
            stage_idx: 0,
            depends_on: vec![5], // 5 is invalid for n=2
        }];

        let result = topological_sort(2, &deps);
        assert!(result.is_err());
        match result {
            Err(OrchestratorError::InvalidStageIndex {
                index: 5,
                total_stages: 2,
            }) => {} // expected
            other => panic!("expected InvalidStageIndex, got {other:?}"),
        }
    }

    // ── Multiple independent stages in one wave ─────────────────────────

    #[test]
    fn test_all_independent() {
        let mut orch = ParallelOrchestrator::new();
        orch.add_stage("X", &[]);
        orch.add_stage("Y", &[]);
        orch.add_stage("Z", &[]);

        let schedule = orch.compute_schedule().expect("should succeed");
        assert_eq!(
            schedule.waves.len(),
            1,
            "all-independent stages fit in one wave"
        );
        assert_eq!(schedule.waves[0], vec![0, 1, 2]);
    }

    // ── Wide diamond (fan-out + fan-in) ─────────────────────────────────

    #[test]
    fn test_wide_fan_out_fan_in() {
        // Root -> 4 parallel stages -> sink
        let mut orch = ParallelOrchestrator::new();
        let root = orch.add_stage("root", &[]);
        let mid: Vec<usize> = (0..4)
            .map(|i| orch.add_stage(&format!("mid_{i}"), &[root]))
            .collect();
        let _sink = orch.add_stage("sink", &mid);

        let schedule = orch.compute_schedule().expect("should succeed");
        assert_eq!(schedule.waves.len(), 3);
        assert_eq!(schedule.waves[0], vec![0]); // root
        assert_eq!(schedule.waves[1], vec![1, 2, 3, 4]); // mid_0..mid_3
        assert_eq!(schedule.waves[2], vec![5]); // sink
    }

    // ── Execution with diamond verifies correct ordering ────────────────

    #[test]
    fn test_diamond_execution_order() {
        let mut orch = ParallelOrchestrator::new();
        let a = orch.add_stage("A", &[]);
        let b = orch.add_stage("B", &[a]);
        let c = orch.add_stage("C", &[a]);
        let _d = orch.add_stage("D", &[b, c]);

        let log = Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut stages: Vec<Box<dyn SolverStage>> = vec![
            Box::new(TestStage::new("A", 1.0, Arc::clone(&log))),
            Box::new(TestStage::new("B", 2.0, Arc::clone(&log))),
            Box::new(TestStage::new("C", 1.5, Arc::clone(&log))),
            Box::new(TestStage::new("D", 3.0, Arc::clone(&log))),
        ];

        orch.execute(&mut stages, 0.01)
            .expect("execute should succeed");

        let recorded = log.lock().expect("lock should not be poisoned");
        // A must come before B, C. D must come after both B and C.
        let pos_a = recorded
            .iter()
            .position(|s| s == "A")
            .expect("A should be in log");
        let pos_b = recorded
            .iter()
            .position(|s| s == "B")
            .expect("B should be in log");
        let pos_c = recorded
            .iter()
            .position(|s| s == "C")
            .expect("C should be in log");
        let pos_d = recorded
            .iter()
            .position(|s| s == "D")
            .expect("D should be in log");

        assert!(pos_a < pos_b, "A must run before B");
        assert!(pos_a < pos_c, "A must run before C");
        assert!(pos_b < pos_d, "B must run before D");
        assert!(pos_c < pos_d, "C must run before D");
    }

    // ── Reset timings ───────────────────────────────────────────────────

    #[test]
    fn test_reset_timings() {
        let mut orch = ParallelOrchestrator::new();
        orch.add_stage("A", &[]);

        let mut stages: Vec<Box<dyn SolverStage>> = vec![Box::new(TimedStage::new("A", 100_000))];

        orch.execute(&mut stages, 0.01)
            .expect("execute should succeed");
        assert!(orch.total_time() > 0.0);

        orch.reset_timings();
        assert!(
            orch.total_time().abs() < 1e-15,
            "timings should be zero after reset"
        );
    }

    // ── Self-dependency is a cycle ──────────────────────────────────────

    #[test]
    fn test_self_dependency_cycle() {
        let deps = vec![StageDependency {
            stage_idx: 0,
            depends_on: vec![0],
        }];

        let result = topological_sort(1, &deps);
        assert!(result.is_err(), "self-dependency should be a cycle");
        match result {
            Err(OrchestratorError::CycleDetected { .. }) => {} // expected
            other => panic!("expected CycleDetected, got {other:?}"),
        }
    }

    // ── PipelineSchedule helpers ────────────────────────────────────────

    #[test]
    fn test_pipeline_schedule_helpers() {
        let schedule = PipelineSchedule {
            waves: vec![vec![0, 1], vec![2], vec![3, 4, 5]],
        };
        assert_eq!(schedule.num_waves(), 3);
        assert_eq!(schedule.num_stages(), 6);
    }

    // ── OrchestratorError display ───────────────────────────────────────

    #[test]
    fn test_error_display() {
        let err = OrchestratorError::CycleDetected {
            description: "A -> B -> A".to_string(),
        };
        let msg = format!("{err}");
        assert!(msg.contains("cycle"));
        assert!(msg.contains("A -> B -> A"));

        let err2 = OrchestratorError::InvalidStageIndex {
            index: 10,
            total_stages: 3,
        };
        let msg2 = format!("{err2}");
        assert!(msg2.contains("10"));
        assert!(msg2.contains("3"));

        let err3 = OrchestratorError::ExecutionIndexOutOfRange {
            index: 5,
            stages_len: 2,
        };
        let msg3 = format!("{err3}");
        assert!(msg3.contains("5"));
        assert!(msg3.contains("2"));
    }

    // ── Default trait impl ──────────────────────────────────────────────

    #[test]
    fn test_default_orchestrator() {
        let orch = ParallelOrchestrator::default();
        assert_eq!(orch.num_stages(), 0);
        assert!(orch.timings().is_empty());
        assert!(orch.total_time().abs() < 1e-15);
    }
}
