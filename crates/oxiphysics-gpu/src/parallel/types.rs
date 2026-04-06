//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// Result of load balancing: a set of ranges assigned to workers.
#[allow(dead_code)]
pub struct LoadBalancePlan {
    /// Ranges for each worker.
    pub ranges: Vec<std::ops::Range<usize>>,
    /// Total weight assigned to each worker (for Weighted strategy).
    pub weights: Vec<f64>,
}
impl LoadBalancePlan {
    /// Number of workers in this plan.
    #[allow(dead_code)]
    pub fn num_workers(&self) -> usize {
        self.ranges.len()
    }
    /// Maximum weight across workers (a measure of imbalance).
    #[allow(dead_code)]
    pub fn max_weight(&self) -> f64 {
        self.weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Imbalance ratio: max_weight / avg_weight. 1.0 is perfect balance.
    #[allow(dead_code)]
    pub fn imbalance_ratio(&self) -> f64 {
        if self.weights.is_empty() {
            return 1.0;
        }
        let total: f64 = self.weights.iter().sum();
        let avg = total / self.weights.len() as f64;
        if avg < 1e-15 {
            return 1.0;
        }
        self.max_weight() / avg
    }
}
/// A simple work-stealing queue for task-parallel scheduling.
///
/// Internally backed by a `Vec`T` with a front/back cursor pair.  "Stealing"
/// takes from the front (like a deque), while the owner pushes/pops from the
/// back.
#[allow(dead_code)]
pub struct WorkStealQueue<T> {
    pub(super) items: std::collections::VecDeque<T>,
}
#[allow(dead_code)]
impl<T: Send> WorkStealQueue<T> {
    /// Create an empty work-steal queue.
    pub fn new() -> Self {
        Self {
            items: std::collections::VecDeque::new(),
        }
    }
    /// Push a task onto the owner end (back).
    pub fn push(&mut self, task: T) {
        self.items.push_back(task);
    }
    /// Pop a task from the owner end (back).  Returns `None` if empty.
    pub fn pop(&mut self) -> Option<T> {
        self.items.pop_back()
    }
    /// Steal a task from the thief end (front).  Returns `None` if empty.
    pub fn steal(&mut self) -> Option<T> {
        self.items.pop_front()
    }
    /// Number of pending tasks.
    pub fn len(&self) -> usize {
        self.items.len()
    }
    /// Whether the queue is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}
/// Configuration for choosing optimal work group sizes.
///
/// Models GPU-like work group sizing where the total work is divided into
/// groups of a fixed size, potentially with padding in the last group.
#[allow(dead_code)]
pub struct WorkGroupConfig {
    /// Preferred work group size (e.g. 64, 128, 256).
    pub preferred_size: usize,
    /// Maximum work group size supported.
    pub max_size: usize,
    /// Minimum work group size (avoid groups too small for efficiency).
    pub min_size: usize,
}
impl WorkGroupConfig {
    /// Create a new config with preferred group size.
    #[allow(dead_code)]
    pub fn new(preferred_size: usize) -> Self {
        Self {
            preferred_size: preferred_size.max(1),
            max_size: 1024,
            min_size: 32,
        }
    }
    /// Create a default config suitable for CPU-side Rayon parallelism.
    #[allow(dead_code)]
    pub fn cpu_default() -> Self {
        let threads = rayon::current_num_threads().max(1);
        Self {
            preferred_size: 64,
            max_size: 1024,
            min_size: threads,
        }
    }
    /// Compute the optimal work group size for `total` items.
    ///
    /// Returns a size in `\[min_size, max_size\]` that balances occupancy.
    /// Prefers `preferred_size` but adjusts if `total` is small.
    #[allow(dead_code)]
    pub fn optimal_size(&self, total: usize) -> usize {
        if total == 0 {
            return self.min_size;
        }
        if total <= self.preferred_size {
            return total.max(self.min_size).min(self.max_size);
        }
        let preferred_groups = total.div_ceil(self.preferred_size);
        let preferred_waste = preferred_groups * self.preferred_size - total;
        let preferred_waste_ratio = preferred_waste as f64 / total as f64;
        if preferred_waste_ratio < 0.25 {
            return self.preferred_size;
        }
        let mut best_size = self.preferred_size;
        let mut best_waste = preferred_waste;
        for candidate in (self.min_size..=self.max_size).step_by(self.min_size) {
            let groups = total.div_ceil(candidate);
            let waste = groups * candidate - total;
            if waste < best_waste {
                best_waste = waste;
                best_size = candidate;
            }
        }
        best_size
    }
    /// Compute the number of work groups needed for `total` items.
    #[allow(dead_code)]
    pub fn num_groups(&self, total: usize) -> usize {
        let size = self.optimal_size(total);
        total.div_ceil(size)
    }
    /// Return ranges for each work group covering `0..total`.
    #[allow(dead_code)]
    pub fn group_ranges(&self, total: usize) -> Vec<std::ops::Range<usize>> {
        let size = self.optimal_size(total);
        (0..total)
            .step_by(size.max(1))
            .map(|start| start..(start + size).min(total))
            .collect()
    }
}
/// Load balancing strategy for distributing work across threads.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum LoadBalanceStrategy {
    /// Static: divide work evenly by index count.
    Static,
    /// Weighted: divide work so each thread gets roughly equal weight.
    Weighted,
    /// Guided: start with large chunks, decrease chunk size as work progresses.
    Guided,
}
/// Splits `n` particle indices into chunks suitable for parallel work.
///
/// `chunk_size` is chosen so that each Rayon worker thread gets at least one
/// chunk.
pub struct WorkChunker {
    /// Total number of items.
    pub n: usize,
    /// Size of each chunk.
    pub chunk_size: usize,
}
impl WorkChunker {
    /// Create a new `WorkChunker` for `n` items.
    ///
    /// `chunk_size` is set to `n / rayon::current_num_threads() + 1`.
    pub fn new(n: usize) -> Self {
        let threads = rayon::current_num_threads().max(1);
        let chunk_size = n / threads + 1;
        Self { n, chunk_size }
    }
    /// Return contiguous index ranges covering `0..n` without gaps or overlaps.
    pub fn chunks(&self) -> Vec<std::ops::Range<usize>> {
        let cs = self.chunk_size.max(1);
        (0..self.n)
            .step_by(cs)
            .map(|start| start..(start + cs).min(self.n))
            .collect()
    }
}
