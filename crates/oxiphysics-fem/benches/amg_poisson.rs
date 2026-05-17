// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmark: AMG vs PCG-Jacobi on 2D Poisson.
//! Run with: cargo bench -p oxiphysics-fem --bench amg_poisson

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use oxiphysics_fem::{
    PcgWithPrecond,
    parallel_solver::ParallelPcgSolver,
    solvers::amg::{
        classical::{AmgClassical, make_2d_poisson},
        cycle::{CycleKind, amg_solve},
        preconditioner::AmgPreconditioner,
    },
};
use std::hint::black_box;

fn make_rhs(n: usize) -> Vec<f64> {
    vec![1.0; n]
}

fn bench_amg_vs_jacobi(c: &mut Criterion) {
    let mut group = c.benchmark_group("amg_poisson");
    // Only benchmark small sizes in CI (128³ would be 2M DOF — too slow for benchmarks that run in CI)
    // Use 16², 32² for the benchmarks
    for size in [16usize, 32] {
        let a = make_2d_poisson(size);
        let n = a.nrows;
        let b = make_rhs(n);

        // Benchmark PCG-Jacobi
        group.bench_with_input(
            BenchmarkId::new("pcg_jacobi", size * size),
            &size,
            |bench, _| {
                bench.iter(|| {
                    let mut x = vec![0.0f64; n];
                    let pcg = ParallelPcgSolver::new(500, 1e-6);
                    black_box(pcg.solve(&a, &b, &mut x))
                });
            },
        );

        // Benchmark PCG-AMG
        // Note: hierarchy.clone() is inside iter() because PcgWithPrecond takes ownership;
        // the clone cost (entire multilevel operator) is included in the measurement.
        group.bench_with_input(
            BenchmarkId::new("pcg_amg", size * size),
            &size,
            |bench, _| {
                let amg = AmgClassical::new();
                let hierarchy = amg.build(&a);
                let pcg_coarse = ParallelPcgSolver::new(500, 1e-10);
                bench.iter(|| {
                    let mut x = vec![0.0f64; n];
                    let precond = AmgPreconditioner {
                        hierarchy: hierarchy.clone(),
                        cycle_kind: CycleKind::V,
                        pcg: pcg_coarse.clone(),
                    };
                    let pcg = PcgWithPrecond::new(precond, 50, 1e-6);
                    black_box(pcg.solve(&a, &b, &mut x))
                });
            },
        );

        // Report iteration counts (informational eprintln, not part of Criterion measurement)
        let amg = AmgClassical::new();
        let hierarchy = amg.build(&a);
        let mut x = vec![0.0f64; n];
        let stats = amg_solve(&hierarchy, &b, &mut x, CycleKind::V, 50, 1e-6);
        eprintln!(
            "AMG {size}²: {} iterations, converged={}",
            stats.iterations, stats.converged
        );

        let mut x2 = vec![0.0f64; n];
        let pcg = ParallelPcgSolver::new(500, 1e-6);
        let stats2 = pcg.solve(&a, &b, &mut x2);
        eprintln!(
            "Jacobi-PCG {size}²: {} iterations, converged={}",
            stats2.iterations, stats2.converged
        );
    }
    group.finish();
}

criterion_group!(benches, bench_amg_vs_jacobi);
criterion_main!(benches);
