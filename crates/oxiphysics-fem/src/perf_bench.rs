// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Performance benchmark utilities for the OxiPhysics FEM crate.
//!
//! This module provides timing helpers, representative benchmark kernels, and
//! FLOP-rate estimators that can be driven by either `criterion` benchmarks or
//! simple wall-clock timing in examples.
//!
//! # Benchmark suite
//!
//! | Benchmark | Description | Dominant operation |
//! |---|---|---|
//! | [`bench_spmv`] | Sparse matrix–vector multiply | CSR SpMV, O(nnz) |
//! | [`bench_pcg`] | Parallel conjugate gradient | SpMV × iterations |
//! | [`bench_assembly`] | Parallel FE stiffness assembly | Ke scatter, O(n_elem) |
//! | [`bench_element_stiffness`] | Per-element Ke computation | Dense 12×12 |
//! | [`bench_gmres`] | Restarted GMRES | Arnoldi + backsolve |
//!
//! ## Quick usage
//!
//! ```
//! use oxiphysics_fem::perf_bench::{BenchHarness, bench_spmv, bench_pcg};
//!
//! let mut h = BenchHarness::new();
//!
//! // Benchmark SpMV with a 1000-DOF tridiagonal system
//! let report = h.run("spmv_n1000", || bench_spmv(1000));
//! println!("{}", report);
//!
//! // Benchmark PCG to 1e-6 tolerance on the same system
//! let report = h.run("pcg_n1000", || bench_pcg(1000, 100, 1e-6));
//! println!("{}", report);
//! ```

#![allow(dead_code)]

use crate::parallel_solver::{
    AssemblyTask, CsrMatrix, ParallelAssembler, ParallelGmresSolver, ParallelPcgSolver,
};
use std::time::{Duration, Instant};

// ── BenchReport ───────────────────────────────────────────────────────────────

/// Result of one benchmark run.
#[derive(Debug, Clone)]
pub struct BenchReport {
    /// Human-readable benchmark name.
    pub name: String,
    /// Number of iterations actually run (for amortised timing).
    pub iterations: u32,
    /// Total wall-clock time for all iterations.
    pub total_time: Duration,
    /// Mean time per iteration.
    pub mean_time: Duration,
    /// Estimated throughput in MFLOPs (if computed by the benchmark).
    pub mflops: Option<f64>,
    /// Problem size (N) as a diagnostic dimension.
    pub n: usize,
    /// Optional extra note (e.g. solver iterations to convergence).
    pub note: Option<String>,
}

impl std::fmt::Display for BenchReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{}] n={} iters={} mean={:.3}µs",
            self.name,
            self.n,
            self.iterations,
            self.mean_time.as_secs_f64() * 1e6
        )?;
        if let Some(mf) = self.mflops {
            write!(f, " {:.1} MFLOPs", mf)?;
        }
        if let Some(ref note) = self.note {
            write!(f, " ({})", note)?;
        }
        Ok(())
    }
}

// ── BenchHarness ─────────────────────────────────────────────────────────────

/// Timing harness that runs a closure multiple times and collects statistics.
pub struct BenchHarness {
    /// Number of warm-up iterations (not timed).
    pub warmup: u32,
    /// Number of timed iterations.
    pub iterations: u32,
    /// Accumulated reports.
    pub reports: Vec<BenchReport>,
}

impl BenchHarness {
    /// Create a harness with 3 warm-up and 10 timed iterations.
    pub fn new() -> Self {
        Self {
            warmup: 3,
            iterations: 10,
            reports: Vec::new(),
        }
    }

    /// Create a faster harness (1 warm-up, 5 timed).
    pub fn fast() -> Self {
        Self {
            warmup: 1,
            iterations: 5,
            reports: Vec::new(),
        }
    }

    /// Run `f` and return a `BenchReport`.
    pub fn run<F>(&mut self, name: &str, mut f: F) -> BenchReport
    where
        F: FnMut() -> BenchResult,
    {
        // Warm-up (result discarded)
        let mut last = BenchResult::default();
        for _ in 0..self.warmup {
            last = f();
        }

        // Timed runs
        let start = Instant::now();
        for _ in 0..self.iterations {
            last = f();
        }
        let total = start.elapsed();
        let mean = total / self.iterations;

        let report = BenchReport {
            name: name.to_string(),
            iterations: self.iterations,
            total_time: total,
            mean_time: mean,
            mflops: last.mflops,
            n: last.n,
            note: last.note,
        };
        self.reports.push(report.clone());
        report
    }

    /// Print a summary table of all collected reports.
    pub fn print_summary(&self) {
        println!("\n{:=<70}", "");
        println!(
            "{:<30} {:>8} {:>12} {:>10}",
            "Benchmark", "N", "Mean (µs)", "MFLOPs"
        );
        println!("{:=<70}", "");
        for r in &self.reports {
            let mf = r.mflops.map_or("—".to_string(), |m| format!("{:.1}", m));
            println!(
                "{:<30} {:>8} {:>12.3} {:>10}",
                r.name,
                r.n,
                r.mean_time.as_secs_f64() * 1e6,
                mf
            );
        }
        println!("{:=<70}", "");
    }
}

impl Default for BenchHarness {
    fn default() -> Self {
        Self::new()
    }
}

// ── BenchResult ───────────────────────────────────────────────────────────────

/// Internal return type from a benchmark closure.
#[derive(Debug, Default, Clone)]
pub struct BenchResult {
    /// Problem size.
    pub n: usize,
    /// Estimated MFLOPs.
    pub mflops: Option<f64>,
    /// Optional note (e.g. solver iterations).
    pub note: Option<String>,
}

// ── Matrix generators ─────────────────────────────────────────────────────────

/// Build an N×N tridiagonal CSR matrix (1D Poisson-like).
///
/// A = tridiag(−1, 2, −1).  Symmetric positive definite.  Used as the
/// canonical test problem for SpMV and iterative solver benchmarks.
pub fn tridiagonal_csr(n: usize) -> CsrMatrix {
    let nnz = if n == 1 { 1 } else { 3 * n - 2 };
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::with_capacity(nnz);
    let mut values = Vec::with_capacity(nnz);

    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }
        col_indices.push(i);
        values.push(2.0);
        if i < n - 1 {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        row_offsets[i + 1] = col_indices.len();
    }

    CsrMatrix {
        nrows: n,
        ncols: n,
        row_offsets,
        col_indices,
        values,
    }
}

/// Build a banded N×N CSR matrix with bandwidth `bw` (2bw+1 diagonals).
///
/// Useful for testing SpMV performance with varying sparsity patterns.
pub fn banded_csr(n: usize, bw: usize) -> CsrMatrix {
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    for i in 0..n {
        let lo = i.saturating_sub(bw);
        let hi = (i + bw).min(n - 1);
        for j in lo..=hi {
            col_indices.push(j);
            if j == i {
                values.push(2.0 * (hi - lo + 1) as f64);
            }
            // diagonal dominant
            else {
                values.push(-1.0);
            }
        }
        row_offsets[i + 1] = col_indices.len();
    }

    CsrMatrix {
        nrows: n,
        ncols: n,
        row_offsets,
        col_indices,
        values,
    }
}

/// Generate an assembly task for a regular 3-node triangle element.
///
/// The element stiffness Ke is a 6×6 identity-like matrix for simplicity.
/// In production this would be the actual FE stiffness from integration.
pub fn example_assembly_task(element_id: usize, n_dof_per_node: usize) -> AssemblyTask {
    let n_nodes = 3;
    let n_dof = n_nodes * n_dof_per_node;
    let global_start = element_id * n_dof_per_node;

    let global_dofs: Vec<usize> = (0..n_nodes)
        .flat_map(|node| (0..n_dof_per_node).map(move |d| global_start + node * n_dof_per_node + d))
        .collect();

    // Ke: identity scaled by element index (to avoid all-zero patterns)
    let scale = 1.0 + element_id as f64 * 0.01;
    let mut ke = vec![0.0; n_dof * n_dof];
    for d in 0..n_dof {
        ke[d * n_dof + d] = scale;
    }

    AssemblyTask { global_dofs, ke }
}

// ── Benchmark kernels ─────────────────────────────────────────────────────────

/// Benchmark: sparse matrix–vector multiplication on an N-DOF tridiagonal system.
///
/// Performs the multiply A × x where x = [1, 2, …, N] and returns the
/// MFLOPs estimate.
///
/// FLOPs = 2 × nnz  (one multiply + one add per non-zero).
pub fn bench_spmv(n: usize) -> BenchResult {
    let a = tridiagonal_csr(n);
    let x: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let mut y = vec![0.0; n];

    let t0 = Instant::now();
    a.spmv(&x, &mut y);
    let elapsed = t0.elapsed().as_secs_f64();

    let nnz = a.nnz();
    let flops = 2.0 * nnz as f64;
    let mflops = if elapsed > 0.0 {
        flops / elapsed / 1e6
    } else {
        f64::INFINITY
    };

    // Verify correctness of first/last element (2.0*1.0 - 2.0*1.0 = 0 intentional for the check)
    #[allow(clippy::eq_op)]
    {
        debug_assert!((y[0] - (2.0 * 1.0 - 2.0)).abs() < 1e-6 || n == 1);
    }

    BenchResult {
        n,
        mflops: Some(mflops),
        note: Some(format!("nnz={}", nnz)),
    }
}

/// Benchmark: parallel SpMV using Rayon.
pub fn bench_spmv_parallel(n: usize) -> BenchResult {
    let a = tridiagonal_csr(n);
    let x: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let mut y = vec![0.0; n];

    let t0 = Instant::now();
    a.spmv_par(&x, &mut y);
    let elapsed = t0.elapsed().as_secs_f64();

    let nnz = a.nnz();
    let mflops = 2.0 * nnz as f64 / elapsed / 1e6;
    BenchResult {
        n,
        mflops: Some(mflops),
        note: Some(format!("rayon nnz={}", nnz)),
    }
}

/// Benchmark: Parallel Conjugate Gradient to tolerance `tol` or `max_iter` steps.
///
/// Solves A x = b where A is the N-DOF tridiagonal system and b = ones.
/// Returns number of iterations in the note field.
pub fn bench_pcg(n: usize, max_iter: usize, tol: f64) -> BenchResult {
    let a = tridiagonal_csr(n);
    let b = vec![1.0; n];
    let mut x = vec![0.0; n];

    let solver = ParallelPcgSolver::new(max_iter, tol);

    let t0 = Instant::now();
    let stats = solver.solve(&a, &b, &mut x);
    let elapsed = t0.elapsed().as_secs_f64();
    let iters = stats.iterations;

    // FLOPs per iteration ≈ 4 × 2 × nnz (4 SpMVs + dot products)
    let nnz = a.nnz();
    let flops = iters as f64 * 6.0 * nnz as f64;
    let mflops = flops / elapsed / 1e6;

    BenchResult {
        n,
        mflops: Some(mflops),
        note: Some(format!("iters={} tol={:.0e}", iters, tol)),
    }
}

/// Benchmark: restarted GMRES on the tridiagonal system.
pub fn bench_gmres(n: usize, restart: usize, max_iter: usize, tol: f64) -> BenchResult {
    let a = tridiagonal_csr(n);
    let b = vec![1.0; n];
    let mut x = vec![0.0; n];

    let solver = ParallelGmresSolver {
        krylov_dim: restart,
        max_restarts: max_iter,
        tolerance: tol,
    };
    let t0 = Instant::now();
    let stats = solver.solve(&a, &b, &mut x);
    let elapsed = t0.elapsed().as_secs_f64();
    let iters = stats.iterations;

    let nnz = a.nnz();
    let flops = iters as f64 * restart as f64 * 2.0 * nnz as f64;
    let mflops = flops / elapsed.max(1e-9) / 1e6;

    BenchResult {
        n,
        mflops: Some(mflops),
        note: Some(format!(
            "restart={} iters={} tol={:.0e}",
            restart, iters, tol
        )),
    }
}

/// Benchmark: parallel FE stiffness assembly for `n_elements` triangle elements.
///
/// Uses [`ParallelAssembler`] to scatter element stiffness matrices
/// into a global CSR system.
pub fn bench_assembly(n_elements: usize, ndofs: usize) -> BenchResult {
    let dof_per_node = 2;
    let tasks: Vec<AssemblyTask> = (0..n_elements)
        .map(|e| example_assembly_task(e % (ndofs / (3 * dof_per_node)).max(1), dof_per_node))
        .collect();

    let assembler = ParallelAssembler { ndofs };
    let t0 = Instant::now();
    let _k = assembler.assemble(&tasks);
    let elapsed = t0.elapsed().as_secs_f64();

    // FLOPs: n_elements × n_dof² (one add per Ke entry)
    let n_dof_elem = 3 * dof_per_node;
    let flops = n_elements as f64 * (n_dof_elem * n_dof_elem) as f64;
    let mflops = flops / elapsed.max(1e-9) / 1e6;

    BenchResult {
        n: n_elements,
        mflops: Some(mflops),
        note: Some(format!("ndofs={}", ndofs)),
    }
}

/// Benchmark: per-element stiffness computation.
///
/// Generates `n_elements` stiffness matrices.  This measures the raw compute
/// throughput of building Ke before assembly.  For a linear tetrahedron Ke
/// is 12×12 and requires ≈ 2 × 12² = 288 FLOPs per element.
pub fn bench_element_stiffness(n_elements: usize) -> BenchResult {
    let n_dof = 12; // 4-node tet × 3 DOFs
    let t0 = Instant::now();
    let mut checksum = 0.0_f64;
    for e in 0..n_elements {
        let ke = compute_linear_tet_ke_stub(e, n_dof);
        checksum += ke[0]; // prevent dead-code elimination
    }
    let elapsed = t0.elapsed().as_secs_f64();
    let _ = checksum;

    let flops = n_elements as f64 * 2.0 * (n_dof * n_dof) as f64;
    let mflops = flops / elapsed.max(1e-9) / 1e6;

    BenchResult {
        n: n_elements,
        mflops: Some(mflops),
        note: Some(format!("n_dof_per_elem={}", n_dof)),
    }
}

/// Stub: compute a 12×12 linear tetrahedron stiffness matrix.
///
/// In production this would use the actual FE shape function gradients and
/// material constitutive matrix.  Here we use a synthetic SPD matrix to
/// measure pure throughput.
fn compute_linear_tet_ke_stub(element_id: usize, n_dof: usize) -> Vec<f64> {
    let scale = 1.0 + element_id as f64 * 1e-4;
    let mut ke = vec![0.0_f64; n_dof * n_dof];
    // Make it diagonally dominant
    for i in 0..n_dof {
        for j in 0..n_dof {
            ke[i * n_dof + j] = if i == j { 4.0 * scale } else { -scale };
        }
    }
    ke
}

// ── Comparative benchmark ─────────────────────────────────────────────────────

/// Configuration for the FEM benchmark suite.
///
/// Use [`SuiteConfig::full`] for production-scale benchmarking and
/// [`SuiteConfig::smoke`] for tiny smoke tests that exercise the codepaths in
/// well under a second.
#[derive(Debug, Clone)]
pub struct SuiteConfig {
    /// Warm-up iterations passed to [`BenchHarness`].
    pub warmup: u32,
    /// Timed iterations passed to [`BenchHarness`].
    pub iterations: u32,
    /// DOF counts for the SpMV (sequential and parallel) and PCG benchmarks.
    pub spmv_pcg_sizes: Vec<usize>,
    /// PCG solver limit applied to every PCG benchmark.
    pub pcg_max_iter: usize,
    /// PCG convergence tolerance applied to every PCG benchmark.
    pub pcg_tol: f64,
    /// `(n, restart, max_iter, tol)` for a single GMRES benchmark, or `None`
    /// to skip GMRES.
    pub gmres: Option<(usize, usize, usize, f64)>,
    /// `(n_elements, ndofs)` pairs for assembly benchmarks.
    pub assembly: Vec<(usize, usize)>,
    /// Element counts for the per-element stiffness benchmark.
    pub element_stiffness: Vec<usize>,
}

impl SuiteConfig {
    /// Production-scale configuration used by [`run_full_suite`].
    ///
    /// Total wall time is dominated by the GMRES restart and the largest
    /// PCG and assembly cases; expect several minutes in release mode.
    pub fn full() -> Self {
        Self {
            warmup: 1,
            iterations: 5,
            spmv_pcg_sizes: vec![100, 500, 1000, 5000],
            pcg_max_iter: 200,
            pcg_tol: 1e-8,
            gmres: Some((1000, 50, 500, 1e-8)),
            assembly: vec![(1000, 600), (10000, 6000)],
            element_stiffness: vec![1000, 10000],
        }
    }

    /// Tiny configuration for smoke tests; completes in milliseconds.
    pub fn smoke() -> Self {
        Self {
            warmup: 0,
            iterations: 1,
            spmv_pcg_sizes: vec![32],
            pcg_max_iter: 16,
            pcg_tol: 1e-4,
            gmres: Some((32, 8, 4, 1e-3)),
            assembly: vec![(16, 60)],
            element_stiffness: vec![16],
        }
    }
}

/// Run the FEM benchmark suite with the supplied configuration and return a
/// formatted report string.
///
/// ```
/// use oxiphysics_fem::perf_bench::{run_suite, SuiteConfig};
/// let report = run_suite(&SuiteConfig::smoke(), false);
/// assert!(!report.is_empty());
/// ```
pub fn run_suite(cfg: &SuiteConfig, verbose: bool) -> String {
    let mut h = BenchHarness {
        warmup: cfg.warmup,
        iterations: cfg.iterations.max(1),
        reports: Vec::new(),
    };

    for &n in &cfg.spmv_pcg_sizes {
        h.run(&format!("spmv_seq_n{}", n), || bench_spmv(n));
        h.run(&format!("spmv_par_n{}", n), || bench_spmv_parallel(n));
        h.run(&format!("pcg_n{}", n), || {
            bench_pcg(n, cfg.pcg_max_iter, cfg.pcg_tol)
        });
    }
    if let Some((n, restart, max_iter, tol)) = cfg.gmres {
        h.run(&format!("gmres_n{}_r{}", n, restart), || {
            bench_gmres(n, restart, max_iter, tol)
        });
    }
    for &(ne, nd) in &cfg.assembly {
        h.run(&format!("assembly_{}_elem", ne), || bench_assembly(ne, nd));
    }
    for &ne in &cfg.element_stiffness {
        h.run(&format!("ke_stiffness_{}", ne), || {
            bench_element_stiffness(ne)
        });
    }

    let mut out = String::new();
    if verbose {
        for r in &h.reports {
            out.push_str(&format!("{}\n", r));
        }
    } else {
        out = format!("{} benchmarks completed", h.reports.len());
    }
    out
}

/// Run the full FEM benchmark suite and return a formatted report string.
///
/// This is the function to call from a top-level bench binary when you want a
/// comprehensive timing snapshot.  See [`SuiteConfig::full`] for the exact
/// problem sizes used.
///
/// ```no_run
/// use oxiphysics_fem::perf_bench::run_full_suite;
/// let report = run_full_suite(false);
/// assert!(!report.is_empty());
/// ```
pub fn run_full_suite(verbose: bool) -> String {
    run_suite(&SuiteConfig::full(), verbose)
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tridiagonal_csr_shape() {
        let a = tridiagonal_csr(5);
        assert_eq!(a.nrows, 5);
        assert_eq!(a.ncols, 5);
        // 3*5 - 2 = 13 non-zeros
        assert_eq!(a.col_indices.len(), 13);
    }

    #[test]
    fn test_tridiagonal_csr_spmv() {
        // A = tridiag(-1,2,-1), x = [1,1,1,1,1]
        // Ax = [2-1, -1+2-1, -1+2-1, -1+2-1, -1+2] = [1,0,0,0,1]
        let a = tridiagonal_csr(5);
        let x = vec![1.0; 5];
        let mut y = vec![0.0; 5];
        a.spmv(&x, &mut y);
        assert!((y[0] - 1.0).abs() < 1e-10, "y[0]={}", y[0]);
        assert!((y[2]).abs() < 1e-10, "y[2]={}", y[2]);
        assert!((y[4] - 1.0).abs() < 1e-10, "y[4]={}", y[4]);
    }

    #[test]
    fn test_bench_spmv_runs() {
        let result = bench_spmv(100);
        assert_eq!(result.n, 100);
        assert!(result.mflops.is_some());
    }

    #[test]
    fn test_bench_pcg_converges() {
        let result = bench_pcg(50, 100, 1e-8);
        assert!(result.note.is_some());
        let note = result.note.unwrap();
        assert!(note.contains("iters="), "note: {}", note);
    }

    #[test]
    fn test_bench_gmres_runs() {
        let result = bench_gmres(50, 20, 100, 1e-6);
        assert_eq!(result.n, 50);
    }

    #[test]
    fn test_bench_assembly_runs() {
        let result = bench_assembly(100, 60);
        assert_eq!(result.n, 100);
        assert!(result.mflops.is_some());
    }

    #[test]
    fn test_bench_element_stiffness() {
        let result = bench_element_stiffness(100);
        assert_eq!(result.n, 100);
    }

    #[test]
    fn test_banded_csr() {
        let a = banded_csr(10, 2);
        assert_eq!(a.nrows, 10);
        // Each row has up to 5 entries; boundary rows have fewer
        assert!(a.col_indices.len() > 10);
    }

    #[test]
    fn test_bench_harness_warmup() {
        let mut h = BenchHarness {
            warmup: 1,
            iterations: 3,
            reports: Vec::new(),
        };
        let report = h.run("test", || bench_spmv(10));
        assert_eq!(report.iterations, 3);
        assert_eq!(h.reports.len(), 1);
    }

    #[test]
    fn test_run_full_suite() {
        // The full suite includes N=5000 PCG, N=1000 GMRES and 10k-element
        // assembly which together take several minutes in release mode.  For
        // CI we exercise the same code paths through the smoke configuration.
        let s = run_suite(&SuiteConfig::smoke(), false);
        assert!(s.contains("benchmarks"), "output: {}", s);
    }
}
