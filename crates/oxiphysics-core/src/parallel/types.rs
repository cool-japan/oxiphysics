//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::any::Any;
use std::collections::VecDeque;
use std::panic::{AssertUnwindSafe, catch_unwind, resume_unwind};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread;

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
    ///
    /// All chunks are joined before the call returns.  If a chunk panicked, the
    /// payload of the first panic (in chunk order) is re-raised on the calling thread.
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
        let mut first_panic = None;
        for h in handles {
            if let Err(payload) = h.join() {
                first_panic.get_or_insert(payload);
            }
        }
        if let Some(payload) = first_panic {
            resume_unwind(payload);
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
/// A boxed task executed by a [`WorkStealingPool`] worker.
type PoolTask = Box<dyn FnOnce() + Send + 'static>;

/// Mutable state of a [`WorkStealingPool`], protected by one mutex.
struct PoolState {
    /// Tasks waiting for a free worker, in submission order.
    queue: VecDeque<PoolTask>,
    /// Set when the pool is closed: workers exit once the queue is empty.
    shutdown: bool,
    /// Number of tasks currently being executed.
    running: usize,
    /// Counters reported by [`WorkStealingPool::stats`].
    stats: ThreadPoolStats,
    /// Payload of the first task that panicked, re-raised by `join`.
    first_panic: Option<Box<dyn Any + Send + 'static>>,
}

/// State shared between a [`WorkStealingPool`] and its worker threads.
struct PoolShared {
    state: Mutex<PoolState>,
    /// Signalled when a task is queued or the pool is closed.
    work_available: Condvar,
    /// Signalled when the last outstanding task has finished.
    idle: Condvar,
}

impl PoolShared {
    /// Locks the pool state, recovering from a poisoned mutex (the state is only
    /// modified by short, panic-free critical sections).
    fn lock(&self) -> MutexGuard<'_, PoolState> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }
    /// Runs a task that has already been counted in `running`, catching a panic
    /// so that it cannot kill the calling worker, and records the outcome.
    fn run_task(&self, task: PoolTask) {
        let outcome = catch_unwind(AssertUnwindSafe(task));
        let mut state = self.lock();
        state.running -= 1;
        match outcome {
            Ok(()) => state.stats.tasks_completed += 1,
            Err(payload) => {
                if state.first_panic.is_none() {
                    state.first_panic = Some(payload);
                }
            }
        }
        if state.queue.is_empty() && state.running == 0 {
            self.idle.notify_all();
        }
    }
    /// Main loop of a worker thread: execute queued tasks until the pool is
    /// closed and the queue is drained.
    fn worker_loop(&self) {
        loop {
            let task = {
                let mut state = self.lock();
                loop {
                    if let Some(task) = state.queue.pop_front() {
                        state.running += 1;
                        break Some(task);
                    }
                    if state.shutdown {
                        break None;
                    }
                    state = self
                        .work_available
                        .wait(state)
                        .unwrap_or_else(|e| e.into_inner());
                }
            };
            match task {
                Some(task) => self.run_task(task),
                None => return,
            }
        }
    }
}

/// A fixed-size thread pool with one shared first-in-first-out task queue.
///
/// `n_threads` worker threads are started by [`new`](Self::new) and sleep on a
/// condition variable while the queue is empty.  [`submit`](Self::submit) queues
/// a task and wakes a worker, so independent tasks run concurrently; tasks are
/// handed out in submission order.
///
/// * [`wait_idle`](Self::wait_idle) blocks until every task submitted so far has
///   finished, and [`stats`](Self::stats) reports how many tasks were really
///   submitted and completed.
/// * [`join`](Self::join) closes the pool, waits for all queued tasks to finish
///   and stops the workers.  Dropping the pool without calling `join` does the
///   same.
/// * A task that panics does not stop its worker.  The payload of the first
///   panic is re-raised by [`join`](Self::join) (it is discarded when the pool is
///   simply dropped); a panicked task is not counted as completed.
///
/// # Example
/// ```
/// use oxiphysics_core::parallel::WorkStealingPool;
/// use std::sync::Arc;
/// use std::sync::atomic::{AtomicUsize, Ordering};
///
/// let pool = WorkStealingPool::new(4);
/// let counter = Arc::new(AtomicUsize::new(0));
/// for _ in 0..16 {
///     let counter = Arc::clone(&counter);
///     pool.submit(move || {
///         counter.fetch_add(1, Ordering::Relaxed);
///     });
/// }
/// pool.join(); // wait for all tasks to complete
/// assert_eq!(counter.load(Ordering::Relaxed), 16);
/// ```
pub struct WorkStealingPool {
    /// State shared with the worker threads.
    shared: Arc<PoolShared>,
    /// Join handles of the worker threads.
    handles: Vec<thread::JoinHandle<()>>,
    /// Number of worker threads.
    pub n_threads: usize,
}
impl WorkStealingPool {
    /// Create a new `WorkStealingPool` with `n_threads` worker threads
    /// (at least one).
    ///
    /// If the operating system refuses to start some of the threads the pool
    /// works with the ones that did start (`n_threads` then holds the actual
    /// count); a pool without any worker runs every submitted task directly on
    /// the calling thread.
    pub fn new(n_threads: usize) -> Self {
        Self::with_spawner(n_threads, |index, shared| {
            thread::Builder::new()
                .name(format!("oxiphysics-pool-{index}"))
                .spawn(move || shared.worker_loop())
                .ok()
        })
    }
    /// Builds a pool whose worker threads are started by `spawn(index, shared)`;
    /// workers for which `spawn` returns `None` are simply missing.
    fn with_spawner(
        n_threads: usize,
        spawn: impl Fn(usize, Arc<PoolShared>) -> Option<thread::JoinHandle<()>>,
    ) -> Self {
        let requested = n_threads.max(1);
        let shared = Arc::new(PoolShared {
            state: Mutex::new(PoolState {
                queue: VecDeque::new(),
                shutdown: false,
                running: 0,
                stats: ThreadPoolStats::new(requested),
                first_panic: None,
            }),
            work_available: Condvar::new(),
            idle: Condvar::new(),
        });
        let handles: Vec<thread::JoinHandle<()>> = (0..requested)
            .filter_map(|index| spawn(index, Arc::clone(&shared)))
            .collect();
        let started = handles.len();
        shared.lock().stats.n_workers = started;
        Self {
            shared,
            handles,
            n_threads: started,
        }
    }
    /// Submit a task to the pool.
    ///
    /// The task is queued and picked up by the next free worker.  Tasks are
    /// started in the order in which they were submitted.
    pub fn submit(&self, f: impl FnOnce() + Send + 'static) {
        let task: PoolTask = Box::new(f);
        let mut state = self.shared.lock();
        state.stats.tasks_submitted += 1;
        if self.handles.is_empty() {
            state.running += 1;
            drop(state);
            self.shared.run_task(task);
            return;
        }
        state.queue.push_back(task);
        let depth = state.queue.len();
        if depth > state.stats.max_queue_depth {
            state.stats.max_queue_depth = depth;
        }
        drop(state);
        self.shared.work_available.notify_one();
    }
    /// Block until every task submitted so far has finished running.
    ///
    /// The pool stays usable afterwards: more tasks can be submitted.  The call
    /// must not be made from inside a task of the same pool, because that task
    /// would wait for itself.
    pub fn wait_idle(&self) {
        let mut state = self.shared.lock();
        while !(state.queue.is_empty() && state.running == 0) {
            state = self
                .shared
                .idle
                .wait(state)
                .unwrap_or_else(|e| e.into_inner());
        }
    }
    /// Closes the pool, lets the workers finish the queued tasks and joins them.
    ///
    /// A worker cannot join itself: when the last owner of the pool is dropped by
    /// one of the pool's own workers (a task holding the final `Arc`), that worker
    /// is skipped and ends by itself as soon as the queue is drained.
    ///
    /// Returns the payload of the first panic raised by a task or a worker.
    fn shutdown_and_join(&mut self) -> Option<Box<dyn Any + Send + 'static>> {
        self.shared.lock().shutdown = true;
        self.shared.work_available.notify_all();
        let handles = std::mem::take(&mut self.handles);
        let current = (!handles.is_empty()).then(|| thread::current().id());
        let mut worker_panic = None;
        for handle in handles {
            if Some(handle.thread().id()) == current {
                continue;
            }
            if let Err(payload) = handle.join() {
                worker_panic.get_or_insert(payload);
            }
        }
        let mut state = self.shared.lock();
        state.first_panic.take().or(worker_panic)
    }
    /// Wait for all submitted tasks to complete and stop the worker threads.
    ///
    /// This consumes the pool.  Workers keep taking tasks from the queue until it
    /// is empty, so every task submitted before the call has run when `join`
    /// returns.  If a task panicked, the payload of the first panic is re-raised
    /// on the calling thread once all workers have stopped.
    pub fn join(mut self) {
        if let Some(payload) = self.shutdown_and_join() {
            resume_unwind(payload);
        }
    }
    /// Returns a snapshot of pool statistics.
    ///
    /// `tasks_submitted` counts calls to [`submit`](Self::submit), `tasks_completed`
    /// the tasks that have run to completion (a task that panicked is not
    /// counted), and `max_queue_depth` the longest queue observed.
    pub fn stats(&self) -> ThreadPoolStats {
        self.shared.lock().stats.clone()
    }
}
impl Drop for WorkStealingPool {
    /// Waits for the queued tasks to finish and stops the workers; a task panic
    /// is discarded (use [`WorkStealingPool::join`] to observe it).
    fn drop(&mut self) {
        drop(self.shutdown_and_join());
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
#[cfg(test)]
mod degraded_pool_tests {
    use super::WorkStealingPool;
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::mpsc;
    use std::thread;
    use std::time::{Duration, Instant};

    #[test]
    fn pool_without_any_worker_runs_tasks_on_the_calling_thread() {
        let pool = WorkStealingPool::with_spawner(4, |_index, _shared| None);
        assert_eq!(pool.n_threads, 0);
        assert_eq!(pool.stats().n_workers, 0);
        let counter = Arc::new(AtomicUsize::new(0));
        for _ in 0..5 {
            let counter = Arc::clone(&counter);
            pool.submit(move || {
                counter.fetch_add(1, Ordering::SeqCst);
            });
            // Executed synchronously: the effect is visible as soon as `submit` returns.
        }
        assert_eq!(counter.load(Ordering::SeqCst), 5);
        pool.wait_idle();
        let stats = pool.stats();
        assert_eq!(stats.tasks_submitted, 5);
        assert_eq!(stats.tasks_completed, 5);
        pool.join();
    }
    #[test]
    fn inline_task_panic_is_reraised_by_join() {
        let pool = WorkStealingPool::with_spawner(2, |_index, _shared| None);
        pool.submit(|| panic!("inline failure"));
        assert_eq!(pool.stats().tasks_completed, 0);
        let outcome = catch_unwind(AssertUnwindSafe(move || pool.join()));
        let payload = outcome.expect_err("join must re-raise the task panic");
        assert_eq!(payload.downcast_ref::<&str>(), Some(&"inline failure"));
    }
    #[test]
    fn dropping_the_last_pool_handle_on_a_worker_thread_does_not_join_itself() {
        // The task owns the final `Arc` of the pool; dropping it on the worker thread
        // must not try to join that very thread (which deadlocks or panics).
        let pool = Arc::new(WorkStealingPool::new(2));
        let (release, gate) = mpsc::channel::<()>();
        let finished = Arc::new(AtomicUsize::new(0));
        let inner_pool = Arc::clone(&pool);
        let inner_finished = Arc::clone(&finished);
        pool.submit(move || {
            let _ = gate.recv();
            drop(inner_pool);
            inner_finished.fetch_add(1, Ordering::SeqCst);
        });
        drop(pool);
        release.send(()).expect("the task still owns the receiver");
        let deadline = Instant::now() + Duration::from_secs(10);
        while finished.load(Ordering::SeqCst) == 0 && Instant::now() < deadline {
            thread::sleep(Duration::from_millis(5));
        }
        assert_eq!(finished.load(Ordering::SeqCst), 1);
    }
    #[test]
    fn pool_keeps_the_workers_that_could_be_started() {
        let pool = WorkStealingPool::with_spawner(4, |index, shared| {
            index
                .is_multiple_of(2)
                .then(|| thread::spawn(move || shared.worker_loop()))
        });
        assert_eq!(pool.n_threads, 2);
        assert_eq!(pool.stats().n_workers, 2);
        let counter = Arc::new(AtomicUsize::new(0));
        for _ in 0..10 {
            let counter = Arc::clone(&counter);
            pool.submit(move || {
                counter.fetch_add(1, Ordering::SeqCst);
            });
        }
        pool.join();
        assert_eq!(counter.load(Ordering::SeqCst), 10);
    }
}
