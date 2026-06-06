//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::sync::{Arc, Mutex};
use std::thread;

/// Shared work queue type alias for boxed closures.
type WorkQueue = Arc<Mutex<Vec<Box<dyn FnOnce() + Send + 'static>>>>;

/// Statistics collected from a thread pool execution.
#[derive(Debug, Clone, Default)]
pub struct ThreadPoolStats {
    /// Total tasks submitted.
    pub tasks_submitted: usize,
    /// Total tasks completed.
    pub tasks_completed: usize,
    /// Number of worker threads.
    pub n_workers: usize,
    /// Maximum queue depth observed.
    pub max_queue_depth: usize,
}
impl ThreadPoolStats {
    /// Creates a new stats tracker.
    pub fn new(n_workers: usize) -> Self {
        Self {
            tasks_submitted: 0,
            tasks_completed: 0,
            n_workers,
            max_queue_depth: 0,
        }
    }
    /// Tasks per worker (average load).
    pub fn tasks_per_worker(&self) -> f64 {
        if self.n_workers == 0 {
            return 0.0;
        }
        self.tasks_completed as f64 / self.n_workers as f64
    }
    /// Utilization estimate: completed / submitted.
    pub fn completion_rate(&self) -> f64 {
        if self.tasks_submitted == 0 {
            return 1.0;
        }
        self.tasks_completed as f64 / self.tasks_submitted as f64
    }
}
/// Range-based parallel for loop with chunk splitting.
///
/// Executes the closure `f(index)` for every index in `start..end`,
/// splitting the range into `n_chunks` chunks and running each chunk on a
/// separate thread.
///
/// # Type constraints
/// - `f: Fn(usize) + Send + Sync` -- the closure must be thread-safe.
pub struct ParallelFor {
    /// Number of chunks to split the range into.
    pub n_chunks: usize,
}
impl ParallelFor {
    /// Create a new `ParallelFor` using the number of available hardware threads.
    pub fn new() -> Self {
        let n = thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(4);
        Self { n_chunks: n }
    }
    /// Create a new `ParallelFor` with an explicit number of chunks.
    pub fn with_chunks(n_chunks: usize) -> Self {
        Self {
            n_chunks: n_chunks.max(1),
        }
    }
    /// Execute `f(i)` for every `i` in `start..end` in parallel.
    pub fn run(&self, start: usize, end: usize, f: impl Fn(usize) + Send + Sync + 'static) {
        if start >= end {
            return;
        }
        let n = (end - start).min(self.n_chunks);
        let chunk_size = (end - start).div_ceil(n);
        let f = Arc::new(f);
        let mut handles = Vec::new();
        let mut cs = start;
        while cs < end {
            let ce = (cs + chunk_size).min(end);
            let f_clone = Arc::clone(&f);
            let h = thread::spawn(move || {
                for i in cs..ce {
                    f_clone(i);
                }
            });
            handles.push(h);
            cs = ce;
        }
        for h in handles {
            h.join().expect("ParallelFor: worker thread panicked");
        }
    }
    /// Execute `f(i)` for every `i` in `0..n` in parallel.
    pub fn run_n(&self, n: usize, f: impl Fn(usize) + Send + Sync + 'static) {
        self.run(0, n, f);
    }
}
/// A first-in-first-out (serial) work queue that simulates work-stealing
/// semantics without OS threads.
///
/// In a real work-stealing implementation each worker would have its own
/// `deque` and could *steal* from other workers' deques.  Here we provide the
/// API contract with a single shared queue for deterministic testing.
pub struct SerialWorkQueue<T> {
    pub(super) items: std::collections::VecDeque<T>,
}
impl<T> SerialWorkQueue<T> {
    /// Create an empty queue.
    pub fn new() -> Self {
        Self {
            items: std::collections::VecDeque::new(),
        }
    }
    /// Push an item to the back of the queue (producer end).
    pub fn push(&mut self, item: T) {
        self.items.push_back(item);
    }
    /// Pop an item from the front of the queue (consumer end).
    ///
    /// Returns `None` when the queue is empty.
    pub fn pop(&mut self) -> Option<T> {
        self.items.pop_front()
    }
    /// "Steal" an item from the back of the queue (work-stealing from tail).
    ///
    /// Returns `None` when the queue is empty.
    pub fn steal(&mut self) -> Option<T> {
        self.items.pop_back()
    }
    /// Number of items currently in the queue.
    pub fn len(&self) -> usize {
        self.items.len()
    }
    /// Returns `true` if the queue is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    /// Drain all items, executing `f` on each.
    pub fn drain_and_run(&mut self, mut f: impl FnMut(T)) {
        while let Some(item) = self.pop() {
            f(item);
        }
    }
}
/// Represents a range of work items split into chunks for parallel processing.
pub struct WorkRange {
    /// Index of the first work item (inclusive).
    pub start: usize,
    /// Index past the last work item (exclusive).
    pub end: usize,
    /// Preferred number of items per chunk.
    pub chunk_size: usize,
}
impl WorkRange {
    /// Creates a new `WorkRange`.
    pub fn new(start: usize, end: usize, chunk_size: usize) -> Self {
        Self {
            start,
            end,
            chunk_size,
        }
    }
    /// Returns `(chunk_start, chunk_end)` pairs covering the full range.
    pub fn chunks(&self) -> Vec<(usize, usize)> {
        if self.start >= self.end || self.chunk_size == 0 {
            return vec![];
        }
        let mut result = Vec::new();
        let mut pos = self.start;
        while pos < self.end {
            let next = (pos + self.chunk_size).min(self.end);
            result.push((pos, next));
            pos = next;
        }
        result
    }
    /// Returns the number of chunks.
    pub fn n_chunks(&self) -> usize {
        self.chunks().len()
    }
    /// Returns the total number of work items.
    pub fn total_work(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}
/// Thread pool metadata container (configuration only, no OS threads allocated).
///
/// For actual parallel execution see [`WorkStealingPool`].
pub struct ThreadPool {
    /// Number of worker threads.
    pub n_threads: usize,
}
impl ThreadPool {
    /// Creates a `ThreadPool` with `n_threads` workers.
    pub fn new(n_threads: usize) -> Self {
        Self { n_threads }
    }
    /// Returns the number of worker threads.
    pub fn n_threads(&self) -> usize {
        self.n_threads
    }
    /// Suggests a chunk size: `n_items / (n_threads * 4)`, minimum 1.
    pub fn suggested_chunk_size(&self, n_items: usize) -> usize {
        let divisor = self.n_threads.saturating_mul(4).max(1);
        (n_items / divisor).max(1)
    }
}
/// A simple thread pool using `std::thread` and channels.
///
/// Tasks are dispatched via a shared work queue (protected by a `Mutex`).
/// Each worker thread pops tasks and executes them until the pool is dropped.
///
/// # Example
/// ```no_run
/// # use oxiphysics_core::parallel::WorkStealingPool;
/// let pool = WorkStealingPool::new(4);
/// pool.join(); // wait for all tasks to complete
/// ```
pub struct WorkStealingPool {
    /// Shared queue of boxed closures.
    pub(super) queue: WorkQueue,
    /// Join handles for worker threads.
    pub(super) handles: Vec<thread::JoinHandle<()>>,
    /// Number of worker threads.
    pub n_threads: usize,
    /// Statistics tracker.
    pub(super) stats: Arc<Mutex<ThreadPoolStats>>,
}
impl WorkStealingPool {
    /// Create a new `WorkStealingPool` with `n_threads` worker threads.
    pub fn new(n_threads: usize) -> Self {
        let n = n_threads.max(1);
        let queue: WorkQueue = Arc::new(Mutex::new(Vec::new()));
        let stats = Arc::new(Mutex::new(ThreadPoolStats::new(n)));
        let handles = (0..n)
            .map(|_| {
                let q = Arc::clone(&queue);
                thread::spawn(move || {
                    loop {
                        let task = {
                            let mut locked = q.lock().unwrap_or_else(|e| e.into_inner());
                            if locked.is_empty() {
                                break;
                            }
                            locked.pop()
                        };
                        if let Some(f) = task {
                            f();
                        }
                    }
                })
            })
            .collect();
        Self {
            queue,
            handles,
            n_threads: n,
            stats,
        }
    }
    /// Submit a task to the pool.
    pub fn submit(&self, f: impl FnOnce() + Send + 'static) {
        let mut q = self.queue.lock().unwrap_or_else(|e| e.into_inner());
        q.push(Box::new(f));
        let depth = q.len();
        drop(q);
        let mut s = self.stats.lock().unwrap_or_else(|e| e.into_inner());
        s.tasks_submitted += 1;
        if depth > s.max_queue_depth {
            s.max_queue_depth = depth;
        }
    }
    /// Wait for all submitted tasks to complete.
    ///
    /// This consumes the pool. Any remaining tasks in the queue are executed
    /// by the calling thread after workers have exited.
    pub fn join(self) {
        let remaining: Vec<_> = {
            let mut q = self.queue.lock().unwrap_or_else(|e| e.into_inner());
            q.drain(..).collect()
        };
        let count = remaining.len();
        for task in remaining {
            task();
        }
        for h in self.handles {
            let _ = h.join();
        }
        let mut s = self.stats.lock().unwrap_or_else(|e| e.into_inner());
        s.tasks_completed = s.tasks_submitted;
        let _ = count;
    }
    /// Returns a snapshot of pool statistics.
    pub fn stats(&self) -> ThreadPoolStats {
        self.stats.lock().unwrap_or_else(|e| e.into_inner()).clone()
    }
}
/// Extended thread pool statistics with timing information.
#[derive(Debug, Clone, Default)]
pub struct ExtendedPoolStats {
    /// Base statistics.
    pub base: ThreadPoolStats,
    /// Total wall-clock nanoseconds consumed (simulated).
    pub total_ns: u64,
    /// Peak memory usage in bytes (simulated).
    pub peak_memory_bytes: usize,
}
impl ExtendedPoolStats {
    /// Create a new extended statistics tracker.
    pub fn new(n_workers: usize) -> Self {
        Self {
            base: ThreadPoolStats::new(n_workers),
            total_ns: 0,
            peak_memory_bytes: 0,
        }
    }
    /// Throughput: tasks per microsecond.
    pub fn throughput_tasks_per_us(&self) -> f64 {
        if self.total_ns == 0 {
            return 0.0;
        }
        self.base.tasks_completed as f64 / (self.total_ns as f64 / 1000.0)
    }
    /// Memory efficiency: tasks completed per KB of peak memory.
    pub fn memory_efficiency(&self) -> f64 {
        if self.peak_memory_bytes == 0 {
            return 0.0;
        }
        self.base.tasks_completed as f64 / (self.peak_memory_bytes as f64 / 1024.0)
    }
}
/// A double-ended work queue with steal-from-back semantics.
///
/// Worker threads pop from the front; thieves steal from the back.
/// This is the standard Chase-Lev deque API (serial simulation).
pub struct WorkStealingDeque<T> {
    pub(super) items: std::collections::VecDeque<T>,
    pub(super) steals: usize,
    pub(super) pops: usize,
}
impl<T> WorkStealingDeque<T> {
    /// Create an empty deque.
    pub fn new() -> Self {
        Self {
            items: std::collections::VecDeque::new(),
            steals: 0,
            pops: 0,
        }
    }
    /// Push a task to the bottom of the deque (local end).
    pub fn push_bottom(&mut self, item: T) {
        self.items.push_back(item);
    }
    /// Pop a task from the bottom (local worker).
    pub fn pop_bottom(&mut self) -> Option<T> {
        let v = self.items.pop_back();
        if v.is_some() {
            self.pops += 1;
        }
        v
    }
    /// Steal a task from the top (remote thief).
    pub fn steal_top(&mut self) -> Option<T> {
        let v = self.items.pop_front();
        if v.is_some() {
            self.steals += 1;
        }
        v
    }
    /// Number of items in the deque.
    pub fn len(&self) -> usize {
        self.items.len()
    }
    /// Returns `true` if the deque is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    /// Number of successful steal operations.
    pub fn steal_count(&self) -> usize {
        self.steals
    }
    /// Number of successful pop operations.
    pub fn pop_count(&self) -> usize {
        self.pops
    }
}
/// Structure-of-Arrays (SoA) layout for f64 data, SIMD-friendly.
///
/// Instead of `Vec<[f64; 4]>` (AoS), stores each component separately so
/// hardware SIMD units can operate on contiguous memory.
pub struct SoaVec3 {
    /// X components.
    pub xs: Vec<f64>,
    /// Y components.
    pub ys: Vec<f64>,
    /// Z components.
    pub zs: Vec<f64>,
}
impl SoaVec3 {
    /// Create an empty SoA vector.
    pub fn new() -> Self {
        Self {
            xs: Vec::new(),
            ys: Vec::new(),
            zs: Vec::new(),
        }
    }
    /// Create a SoA vector with pre-allocated capacity.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            xs: Vec::with_capacity(n),
            ys: Vec::with_capacity(n),
            zs: Vec::with_capacity(n),
        }
    }
    /// Push a 3-vector `(x, y, z)`.
    pub fn push(&mut self, x: f64, y: f64, z: f64) {
        self.xs.push(x);
        self.ys.push(y);
        self.zs.push(z);
    }
    /// Number of 3-vectors stored.
    pub fn len(&self) -> usize {
        self.xs.len()
    }
    /// Returns `true` if no 3-vectors have been pushed.
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }
    /// Get the i-th vector as `(x, y, z)`.
    pub fn get(&self, i: usize) -> (f64, f64, f64) {
        (self.xs[i], self.ys[i], self.zs[i])
    }
    /// Compute the dot product of all i-th vectors with a fixed vector `(ax, ay, az)`.
    ///
    /// Returns a Vec of dot products — this can be done with SIMD fma operations.
    pub fn dot_with(&self, ax: f64, ay: f64, az: f64) -> Vec<f64> {
        self.xs
            .iter()
            .zip(self.ys.iter())
            .zip(self.zs.iter())
            .map(|((&x, &y), &z)| x * ax + y * ay + z * az)
            .collect()
    }
    /// Compute squared norms of all vectors.
    pub fn norms_sq(&self) -> Vec<f64> {
        self.xs
            .iter()
            .zip(self.ys.iter())
            .zip(self.zs.iter())
            .map(|((&x, &y), &z)| x * x + y * y + z * z)
            .collect()
    }
}
