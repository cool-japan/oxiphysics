//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::sync::Arc;
use std::thread;

/// Applies `f` to every element of `data` and returns the results.
///
/// Sequential fallback -- the API is parallel-ready (requires `Sync`).
pub fn parallel_map_f64(data: &[f64], f: impl Fn(f64) -> f64 + Sync) -> Vec<f64> {
    data.iter().map(|&x| f(x)).collect()
}
/// Reduces `data` to a single value using `f`, starting from `identity`.
///
/// Sequential fallback -- the API is parallel-ready (requires `Sync + Send`).
pub fn parallel_reduce_f64(
    data: &[f64],
    identity: f64,
    f: impl Fn(f64, f64) -> f64 + Sync + Send,
) -> f64 {
    data.iter().copied().fold(identity, f)
}
/// Splits `data` into `n_parts` roughly equal slices.
///
/// If `n_parts` is zero or `data` is empty the result is an empty `Vec`.
/// Extra items are distributed to the leading partitions.
pub fn scatter_gather<T: Clone + Default>(data: &[T], n_parts: usize) -> Vec<Vec<T>> {
    if n_parts == 0 || data.is_empty() {
        return vec![];
    }
    let len = data.len();
    let base = len / n_parts;
    let remainder = len % n_parts;
    let mut result = Vec::with_capacity(n_parts);
    let mut offset = 0;
    for i in 0..n_parts {
        let size = base + if i < remainder { 1 } else { 0 };
        result.push(data[offset..offset + size].to_vec());
        offset += size;
    }
    result
}
/// Concatenates a list of parts back into a single `Vec`.
pub fn gather<T: Clone>(parts: Vec<Vec<T>>) -> Vec<T> {
    parts.into_iter().flatten().collect()
}
/// Computes the dot product of `a` and `b`.
///
/// Panics if the slices have different lengths.
pub fn parallel_dot_product(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "parallel_dot_product: length mismatch");
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}
/// Multiplies each row of `matrix_rows` by `vector` and returns the result.
pub fn parallel_matrix_vec_multiply(matrix_rows: &[Vec<f64>], vector: &[f64]) -> Vec<f64> {
    matrix_rows
        .iter()
        .map(|row| {
            assert_eq!(
                row.len(),
                vector.len(),
                "parallel_matrix_vec_multiply: dimension mismatch"
            );
            row.iter().zip(vector.iter()).map(|(&a, &b)| a * b).sum()
        })
        .collect()
}
/// Computes the sum of all elements in `data`.
pub fn parallel_sum(data: &[f64]) -> f64 {
    data.iter().copied().sum()
}
/// Apply `f` to every item in `items` in parallel (no return value).
///
/// Like `parallel_map` but discards results, useful for side-effecting work.
pub fn parallel_for_each<T>(items: &[T], f: impl Fn(&T) + Send + Sync + 'static)
where
    T: Send + Sync + Clone + 'static,
{
    if items.is_empty() {
        return;
    }
    let n_threads = thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .min(items.len());
    let f = Arc::new(f);
    let items_arc: Arc<[T]> = items.into();
    let chunk_size = items.len().div_ceil(n_threads);
    let mut handles = Vec::new();
    let mut start = 0;
    while start < items.len() {
        let end = (start + chunk_size).min(items.len());
        let f_clone = Arc::clone(&f);
        let items_clone = Arc::clone(&items_arc);
        let h = thread::spawn(move || {
            for i in start..end {
                f_clone(&items_clone[i]);
            }
        });
        handles.push(h);
        start = end;
    }
    for h in handles {
        h.join().expect("parallel_for_each: worker thread panicked");
    }
}
/// Apply `f` to every item in `items` in parallel and return the results.
///
/// Uses `std::thread::spawn` to fan out work across chunks, then collects
/// results in the original order.
///
/// # Type constraints
/// - `T: Send + Sync` -- items must be shareable across threads.
/// - `R: Send` -- results must be sendable back to the main thread.
/// - `f: Fn(&T) -> R + Send + Sync` -- the closure must be thread-safe.
pub fn parallel_map<T, R>(items: &[T], f: impl Fn(&T) -> R + Send + Sync + 'static) -> Vec<R>
where
    T: Send + Sync + Clone + 'static,
    R: Send + 'static,
{
    if items.is_empty() {
        return Vec::new();
    }
    let n_threads = thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .min(items.len());
    let f = Arc::new(f);
    let items_arc: Arc<[T]> = items.into();
    let chunk_size = items.len().div_ceil(n_threads);
    let mut handles: Vec<thread::JoinHandle<Vec<R>>> = Vec::new();
    let mut start = 0;
    while start < items.len() {
        let end = (start + chunk_size).min(items.len());
        let f_clone = Arc::clone(&f);
        let items_clone = Arc::clone(&items_arc);
        let h = thread::spawn(move || (start..end).map(|i| f_clone(&items_clone[i])).collect());
        handles.push(h);
        start = end;
    }
    let mut result = Vec::with_capacity(items.len());
    for h in handles {
        result.extend(h.join().expect("parallel_map: worker thread panicked"));
    }
    result
}
/// Filter `items` in parallel, keeping only those for which `predicate` returns `true`.
///
/// Preserves original order.
pub fn parallel_filter<T>(
    items: &[T],
    predicate: impl Fn(&T) -> bool + Send + Sync + 'static,
) -> Vec<T>
where
    T: Send + Sync + Clone + 'static,
{
    if items.is_empty() {
        return Vec::new();
    }
    let n_threads = thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .min(items.len());
    let predicate = Arc::new(predicate);
    let items_arc: Arc<[T]> = items.into();
    let chunk_size = items.len().div_ceil(n_threads);
    let mut handles: Vec<thread::JoinHandle<Vec<T>>> = Vec::new();
    let mut start = 0;
    while start < items.len() {
        let end = (start + chunk_size).min(items.len());
        let pred_clone = Arc::clone(&predicate);
        let items_clone = Arc::clone(&items_arc);
        let h = thread::spawn(move || {
            (start..end)
                .filter(|&i| pred_clone(&items_clone[i]))
                .map(|i| items_clone[i].clone())
                .collect()
        });
        handles.push(h);
        start = end;
    }
    let mut result = Vec::new();
    for h in handles {
        result.extend(h.join().expect("parallel_filter: worker thread panicked"));
    }
    result
}
/// Reduce `items` to a single value using `f` with the given `identity` element.
///
/// Fans out the reduction into chunks (one per available hardware thread),
/// reduces each chunk sequentially, then combines the chunk results.
///
/// # Type constraints
/// - `T: Clone + Send + Sync` -- items and identity must be cloneable and thread-safe.
/// - `f: Fn(T, T) -> T + Send + Sync` -- the combining function must be thread-safe.
pub fn parallel_reduce<T>(
    items: &[T],
    f: impl Fn(T, T) -> T + Send + Sync + 'static,
    identity: T,
) -> T
where
    T: Clone + Send + Sync + 'static,
{
    if items.is_empty() {
        return identity;
    }
    let n_threads = thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .min(items.len());
    let f = Arc::new(f);
    let items_arc: Arc<[T]> = items.into();
    let chunk_size = items.len().div_ceil(n_threads);
    let mut handles: Vec<thread::JoinHandle<T>> = Vec::new();
    let mut start = 0;
    while start < items.len() {
        let end = (start + chunk_size).min(items.len());
        let f_clone = Arc::clone(&f);
        let items_clone = Arc::clone(&items_arc);
        let id = identity.clone();
        let h = thread::spawn(move || {
            let mut acc = id;
            for i in start..end {
                acc = f_clone(acc, items_clone[i].clone());
            }
            acc
        });
        handles.push(h);
        start = end;
    }
    let mut acc = identity;
    for h in handles {
        let chunk_result = h.join().expect("parallel_reduce: worker thread panicked");
        acc = f(acc, chunk_result);
    }
    acc
}
/// Reduce with a custom operator struct for more complex reductions.
///
/// The operator must provide `identity`, `combine`, and `finalize` methods.
pub trait ReduceOperator: Send + Sync + 'static {
    /// The accumulator type.
    type Acc: Clone + Send + Sync + 'static;
    /// The input element type.
    type Item: Clone + Send + Sync + 'static;
    /// The final result type.
    type Result;
    /// Identity element for the accumulator.
    fn identity(&self) -> Self::Acc;
    /// Fold a single item into the accumulator.
    fn fold(&self, acc: Self::Acc, item: Self::Item) -> Self::Acc;
    /// Combine two accumulators.
    fn combine(&self, left: Self::Acc, right: Self::Acc) -> Self::Acc;
    /// Convert the final accumulator to the result.
    fn finalize(&self, acc: Self::Acc) -> Self::Result;
}
/// Reduce with a custom operator.
pub fn parallel_reduce_with_op<Op>(items: &[Op::Item], op: &Op) -> Op::Result
where
    Op: ReduceOperator,
{
    if items.is_empty() {
        return op.finalize(op.identity());
    }
    let mut acc = op.identity();
    for item in items {
        acc = op.fold(acc, item.clone());
    }
    op.finalize(acc)
}
/// Parallel merge sort for slices of `f64`.
///
/// Splits the data into chunks, sorts each chunk on a separate thread,
/// then merges the sorted chunks pairwise.
pub fn parallel_merge_sort(data: &mut [f64]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    if n <= 64 {
        insertion_sort(data);
        return;
    }
    let n_threads = thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4)
        .min(n / 32)
        .max(1);
    let chunk_size = n.div_ceil(n_threads);
    let chunks: Vec<Vec<f64>> = data.chunks(chunk_size).map(|c| c.to_vec()).collect();
    let handles: Vec<thread::JoinHandle<Vec<f64>>> = chunks
        .into_iter()
        .map(|mut chunk| {
            thread::spawn(move || {
                chunk.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                chunk
            })
        })
        .collect();
    let mut sorted_chunks: Vec<Vec<f64>> = handles
        .into_iter()
        .map(|h| h.join().expect("parallel_merge_sort: worker panicked"))
        .collect();
    while sorted_chunks.len() > 1 {
        let mut merged = Vec::new();
        let mut i = 0;
        while i < sorted_chunks.len() {
            if i + 1 < sorted_chunks.len() {
                merged.push(merge_sorted(&sorted_chunks[i], &sorted_chunks[i + 1]));
                i += 2;
            } else {
                merged.push(sorted_chunks[i].clone());
                i += 1;
            }
        }
        sorted_chunks = merged;
    }
    data.copy_from_slice(&sorted_chunks[0]);
}
/// Merge two sorted slices into a single sorted vector.
pub fn merge_sorted(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] <= b[j] {
            result.push(a[i]);
            i += 1;
        } else {
            result.push(b[j]);
            j += 1;
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}
/// Simple insertion sort for small arrays.
pub fn insertion_sort(data: &mut [f64]) {
    for i in 1..data.len() {
        let key = data[i];
        let mut j = i;
        while j > 0 && data[j - 1] > key {
            data[j] = data[j - 1];
            j -= 1;
        }
        data[j] = key;
    }
}
/// Computes the inclusive prefix sum of `data`.
///
/// Returns a vector where element `i` is `sum(data[0..=i])`.
pub fn prefix_sum(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0_f64;
    for &x in data {
        acc += x;
        result.push(acc);
    }
    result
}
/// Computes the exclusive prefix sum of `data`.
///
/// Returns a vector where element `i` is `sum(data[0..i])`.
pub fn exclusive_prefix_sum(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut acc = 0.0_f64;
    for &x in data {
        result.push(acc);
        acc += x;
    }
    result
}
/// Find the minimum value in `data` in parallel.
///
/// Returns `f64::INFINITY` for empty input.
pub fn parallel_min(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::INFINITY;
    }
    parallel_reduce(data, |a: f64, b: f64| a.min(b), f64::INFINITY)
}
/// Find the maximum value in `data` in parallel.
///
/// Returns `f64::NEG_INFINITY` for empty input.
pub fn parallel_max(data: &[f64]) -> f64 {
    if data.is_empty() {
        return f64::NEG_INFINITY;
    }
    parallel_reduce(data, |a: f64, b: f64| a.max(b), f64::NEG_INFINITY)
}
/// Returns the number of logical CPU cores available to this process.
///
/// Falls back to 1 if the OS query fails.
pub fn available_threads() -> usize {
    thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}
/// Returns a suggested number of threads for CPU-bound work.
///
/// Uses half the available cores (minimum 1) to leave headroom for the OS and
/// other processes.
pub fn suggested_thread_count() -> usize {
    (available_threads() / 2).max(1)
}
/// Process `data` in fixed-size chunks, applying `f` to each chunk and
/// collecting the results into a flat `Vec`R`.
///
/// This pattern mirrors hand-written SIMD loops where a fixed "vector width"
/// is processed per iteration.  Here `chunk_size` plays the role of the SIMD
/// lane count.
pub fn chunk_process<T, R, F>(data: &[T], chunk_size: usize, f: F) -> Vec<R>
where
    T: Clone,
    F: Fn(&[T]) -> Vec<R>,
{
    if data.is_empty() || chunk_size == 0 {
        return Vec::new();
    }
    data.chunks(chunk_size).flat_map(f).collect()
}
/// Apply `f` element-wise on `a` and `b` in chunks of `chunk_size`, returning
/// a `Vec`f64`.
///
/// Panics if `a.len() != b.len()`.
pub fn chunk_zip_map(
    a: &[f64],
    b: &[f64],
    chunk_size: usize,
    f: impl Fn(f64, f64) -> f64,
) -> Vec<f64> {
    assert_eq!(
        a.len(),
        b.len(),
        "chunk_zip_map: slices must be same length"
    );
    if a.is_empty() || chunk_size == 0 {
        return Vec::new();
    }
    let n = a.len();
    let mut result = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        let end = (i + chunk_size).min(n);
        for k in i..end {
            result.push(f(a[k], b[k]));
        }
        i = end;
    }
    result
}
/// Compute the dot product of `a` and `b` using chunk-based accumulation.
///
/// Accumulates partial sums over chunks of `chunk_size`, then sums the
/// partials — mimics a SIMD reduction.
pub fn chunk_dot_product(a: &[f64], b: &[f64], chunk_size: usize) -> f64 {
    assert_eq!(a.len(), b.len(), "chunk_dot_product: length mismatch");
    if a.is_empty() {
        return 0.0;
    }
    let chunk_size = chunk_size.max(1);
    let mut total = 0.0;
    let mut i = 0;
    while i < a.len() {
        let end = (i + chunk_size).min(a.len());
        let mut partial = 0.0;
        for k in i..end {
            partial += a[k] * b[k];
        }
        total += partial;
        i = end;
    }
    total
}
/// In-place quicksort (serial reference for a parallel sort API).
///
/// Uses the median-of-three pivot selection and falls back to insertion sort
/// for subarrays with ≤ 16 elements.
pub fn parallel_sort(data: &mut [f64]) {
    quicksort_recursive(data);
}
/// Recursive quicksort implementation for f64 slices.
pub fn quicksort_recursive(data: &mut [f64]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    if n <= 16 {
        insertion_sort(data);
        return;
    }
    let pivot_idx = median_of_three(data);
    data.swap(pivot_idx, n - 1);
    let pivot = data[n - 1];
    let mut store = 0;
    for i in 0..n - 1 {
        if data[i] <= pivot {
            data.swap(i, store);
            store += 1;
        }
    }
    data.swap(store, n - 1);
    let (left, right) = data.split_at_mut(store);
    quicksort_recursive(left);
    quicksort_recursive(&mut right[1..]);
}
/// Returns the index of the median-of-three pivot (first, mid, last).
pub fn median_of_three(data: &[f64]) -> usize {
    let n = data.len();
    let mid = n / 2;
    let (a, b, c) = (data[0], data[mid], data[n - 1]);
    if (a <= b && b <= c) || (c <= b && b <= a) {
        mid
    } else if (b <= a && a <= c) || (c <= a && a <= b) {
        0
    } else {
        n - 1
    }
}
/// Sort a `Vec`T` using `parallel_sort` (f64 specialisation) and return sorted copy.
pub fn sorted_copy(data: &[f64]) -> Vec<f64> {
    let mut v = data.to_vec();
    parallel_sort(&mut v);
    v
}
/// Accumulate a histogram of `data` into `n_bins` equal-width bins in parallel.
///
/// Returns a `Vec`usize` of length `n_bins` with counts.
/// Empty input or `n_bins == 0` returns an empty `Vec`.
pub fn parallel_histogram(data: &[f64], n_bins: usize) -> Vec<usize> {
    if data.is_empty() || n_bins == 0 {
        return vec![];
    }
    let min = data.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let width = if (max - min).abs() < f64::EPSILON {
        1.0
    } else {
        (max - min) / n_bins as f64
    };
    let mut counts = vec![0usize; n_bins];
    for &v in data {
        let idx = ((v - min) / width) as usize;
        counts[idx.min(n_bins - 1)] += 1;
    }
    counts
}
