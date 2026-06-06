// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Priority-based physics step budget allocation.
//!
//! In a large simulation, not every body or island needs the same computational
//! fidelity every frame.  This module assigns each simulation island a
//! `Priority` level and produces a `FrameSchedule` that tells the engine
//! how many substeps to give each island.
//!
//! ## Life Cycle
//!
//! ```rust,no_run
//! use oxiphysics::scheduler::{Scheduler, Priority};
//!
//! let mut sched = Scheduler::new(2, 4); // base = 2, max = 4
//!
//! let island_a = sched.add_island(Priority::Critical, 10);
//! let island_b = sched.add_island(Priority::Normal, 50);
//! let island_c = sched.add_island(Priority::Low, 200);
//!
//! // Simulate a few frames
//! for _ in 0..60 {
//!     let frame = sched.schedule(1.0 / 60.0);
//!     for alloc in &frame.allocations {
//!         // run `alloc.substeps` sub-steps of dt = `alloc.dt_per_substep`
//!         let _ = alloc.substeps;
//!         let _ = alloc.dt_per_substep;
//!     }
//!     // update velocity hints so the scheduler can auto-sleep idle islands
//!     sched.update_velocity(island_b, 0.001);
//!     sched.update_velocity(island_c, 0.0);
//!     sched.auto_sleep_step();
//! }
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Priority
// ---------------------------------------------------------------------------

/// Priority level assigned to a simulation island.
///
/// Per-priority island-ID groups returned by [`Scheduler::islands_by_priority`].
///
/// Tuple layout: `(critical, high, normal, low, sleeping)`.
pub type IslandsByPriority = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);

/// Higher-priority islands receive more substeps per frame and are never
/// skipped.  Lower-priority islands may receive fewer substeps or be suspended
/// entirely when sleeping.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Priority {
    /// No simulation; island is suspended.
    ///
    /// Use [`Scheduler::wake_island`] to restore an island from this state.
    Sleeping = 0,
    /// Reduced substeps; may be skipped under heavy load.
    Low = 1,
    /// Default priority — receives `base_substeps` per frame.
    Normal = 2,
    /// Elevated priority — receives `max_substeps` per frame.
    High = 3,
    /// Always receives `max_substeps`; never skipped or downgraded automatically.
    Critical = 4,
}

// ---------------------------------------------------------------------------
// IslandEntry
// ---------------------------------------------------------------------------

/// Metadata for a simulation island tracked by the scheduler.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IslandEntry {
    /// Unique island identifier (assigned by [`Scheduler::add_island`]).
    pub id: u32,
    /// Current simulation priority.
    pub priority: Priority,
    /// Number of rigid bodies in this island (informational).
    pub body_count: usize,
    /// Latest velocity magnitude hint (m/s).
    ///
    /// Updated by [`Scheduler::update_velocity`].  Used by
    /// [`Scheduler::auto_sleep_step`] to detect idle islands.
    pub velocity_magnitude: f64,
    /// Consecutive frames this island has been below `sleep_threshold`.
    pub steps_since_active: u32,
}

// ---------------------------------------------------------------------------
// FrameSchedule / StepAllocation
// ---------------------------------------------------------------------------

/// Substep allocation for a single simulation island within one frame.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StepAllocation {
    /// Island this allocation applies to.
    pub island_id: u32,
    /// Number of substeps to run.
    pub substeps: u32,
    /// Duration of each substep (s).
    pub dt_per_substep: f64,
}

/// The full schedule produced by [`Scheduler::schedule`] for one frame.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrameSchedule {
    /// Frame delta time (s).
    pub frame_dt: f64,
    /// Substep allocations for active islands (i.e., non-sleeping).
    pub allocations: Vec<StepAllocation>,
    /// Islands that were skipped this frame (sleeping or suspended).
    pub skipped_islands: Vec<u32>,
}

// ---------------------------------------------------------------------------
// Scheduler
// ---------------------------------------------------------------------------

/// Priority-based physics step scheduler.
///
/// Maintains a list of simulation islands and their priorities, then produces
/// per-frame substep allocations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Scheduler {
    islands: Vec<IslandEntry>,
    next_id: u32,
    /// Substep count for [`Priority::Normal`] islands.
    ///
    /// Must be ≥ 1.
    pub base_substeps: u32,
    /// Substep count for [`Priority::Critical`] and [`Priority::High`] islands.
    ///
    /// Must be ≥ `base_substeps`.
    pub max_substeps: u32,
    /// Velocity magnitude (m/s) below which a frame is counted toward the
    /// sleep delay.
    pub sleep_threshold: f64,
    /// Consecutive frames below `sleep_threshold` before auto-sleeping.
    pub sleep_delay: u32,
}

impl Scheduler {
    /// Create a scheduler.
    ///
    /// * `base_substeps` — substeps for `Normal` islands (clamped to `>= 1`).
    /// * `max_substeps`  — substeps for `High`/`Critical` islands
    ///   (clamped to `>= base_substeps`).
    pub fn new(base_substeps: u32, max_substeps: u32) -> Self {
        let base = base_substeps.max(1);
        let max = max_substeps.max(base);
        Self {
            islands: Vec::new(),
            next_id: 0,
            base_substeps: base,
            max_substeps: max,
            sleep_threshold: 0.01,
            sleep_delay: 60,
        }
    }

    // -----------------------------------------------------------------------
    // Island management
    // -----------------------------------------------------------------------

    /// Register a new simulation island and return its ID.
    pub fn add_island(&mut self, priority: Priority, body_count: usize) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.islands.push(IslandEntry {
            id,
            priority,
            body_count,
            velocity_magnitude: 0.0,
            steps_since_active: 0,
        });
        id
    }

    /// Unregister a simulation island.  No-op if `id` is not found.
    pub fn remove_island(&mut self, id: u32) {
        self.islands.retain(|e| e.id != id);
    }

    /// Immutable reference to an island entry.
    pub fn island(&self, id: u32) -> Option<&IslandEntry> {
        self.islands.iter().find(|e| e.id == id)
    }

    /// Mutable reference to an island entry.
    pub fn island_mut(&mut self, id: u32) -> Option<&mut IslandEntry> {
        self.islands.iter_mut().find(|e| e.id == id)
    }

    /// Change the priority of an island.
    pub fn set_priority(&mut self, id: u32, priority: Priority) {
        if let Some(e) = self.island_mut(id) {
            e.priority = priority;
        }
    }

    /// Update the velocity hint for an island.
    ///
    /// Should be called each frame with the maximum velocity magnitude of any
    /// body in the island.
    pub fn update_velocity(&mut self, id: u32, velocity: f64) {
        if let Some(e) = self.island_mut(id) {
            e.velocity_magnitude = velocity;
        }
    }

    /// Update the body count hint for an island.
    pub fn update_body_count(&mut self, id: u32, count: usize) {
        if let Some(e) = self.island_mut(id) {
            e.body_count = count;
        }
    }

    /// Wake a sleeping island by restoring it to [`Priority::Normal`].
    ///
    /// Has no effect on non-sleeping islands.
    pub fn wake_island(&mut self, id: u32) {
        if let Some(e) = self.island_mut(id)
            && e.priority == Priority::Sleeping
        {
            e.priority = Priority::Normal;
            e.steps_since_active = 0;
        }
    }

    // -----------------------------------------------------------------------
    // Scheduling
    // -----------------------------------------------------------------------

    /// Produce a [`FrameSchedule`] for a frame of duration `dt`.
    ///
    /// Substep allocations by priority:
    ///
    /// | Priority | Substeps |
    /// |----------|----------|
    /// | Sleeping | 0 (skipped) |
    /// | Low | `max(base_substeps / 2, 1)` |
    /// | Normal | `base_substeps` |
    /// | High | `max_substeps` |
    /// | Critical | `max_substeps` |
    pub fn schedule(&self, dt: f64) -> FrameSchedule {
        let mut allocations = Vec::new();
        let mut skipped = Vec::new();

        for island in &self.islands {
            let substeps = match island.priority {
                Priority::Sleeping => {
                    skipped.push(island.id);
                    continue;
                }
                Priority::Low => (self.base_substeps / 2).max(1),
                Priority::Normal => self.base_substeps,
                Priority::High | Priority::Critical => self.max_substeps,
            };
            allocations.push(StepAllocation {
                island_id: island.id,
                substeps,
                dt_per_substep: dt / substeps as f64,
            });
        }

        FrameSchedule {
            frame_dt: dt,
            allocations,
            skipped_islands: skipped,
        }
    }

    // -----------------------------------------------------------------------
    // Auto-sleep
    // -----------------------------------------------------------------------

    /// Advance the idle counter for all non-sleeping islands.
    ///
    /// Islands whose `velocity_magnitude` has been below `sleep_threshold` for
    /// `sleep_delay` consecutive frames are automatically set to
    /// [`Priority::Sleeping`].  [`Priority::Critical`] islands are exempt.
    ///
    /// Returns the IDs of islands that were newly put to sleep this call.
    pub fn auto_sleep_step(&mut self) -> Vec<u32> {
        let threshold = self.sleep_threshold;
        let delay = self.sleep_delay;
        let mut newly_sleeping = Vec::new();

        for e in &mut self.islands {
            if e.priority == Priority::Sleeping || e.priority == Priority::Critical {
                continue;
            }
            if e.velocity_magnitude < threshold {
                e.steps_since_active += 1;
                if e.steps_since_active >= delay {
                    e.priority = Priority::Sleeping;
                    newly_sleeping.push(e.id);
                }
            } else {
                e.steps_since_active = 0;
            }
        }

        newly_sleeping
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Total number of registered islands.
    pub fn island_count(&self) -> usize {
        self.islands.len()
    }

    /// Number of islands that are not sleeping.
    pub fn active_count(&self) -> usize {
        self.islands
            .iter()
            .filter(|e| e.priority != Priority::Sleeping)
            .count()
    }

    /// Iterate over all island entries.
    pub fn islands(&self) -> impl Iterator<Item = &IslandEntry> {
        self.islands.iter()
    }

    /// Collect island IDs grouped by priority tier.
    ///
    /// Returns `(critical, high, normal, low, sleeping)`.
    pub fn islands_by_priority(&self) -> IslandsByPriority {
        let mut critical = Vec::new();
        let mut high = Vec::new();
        let mut normal = Vec::new();
        let mut low = Vec::new();
        let mut sleeping = Vec::new();
        for e in &self.islands {
            match e.priority {
                Priority::Critical => critical.push(e.id),
                Priority::High => high.push(e.id),
                Priority::Normal => normal.push(e.id),
                Priority::Low => low.push(e.id),
                Priority::Sleeping => sleeping.push(e.id),
            }
        }
        (critical, high, normal, low, sleeping)
    }
}
