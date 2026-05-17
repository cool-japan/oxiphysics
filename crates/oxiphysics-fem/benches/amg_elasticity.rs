// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmark: AMG build and solve on an elasticity-proxy (2D Poisson) system.
//! Run with: cargo bench -p oxiphysics-fem --bench amg_elasticity

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use oxiphysics_fem::solvers::amg::{
    classical::{AmgClassical, make_2d_poisson},
    cycle::{CycleKind, amg_solve},
};
use std::hint::black_box;

fn bench_amg_setup_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("amg_elasticity_proxy");
    for size in [16usize, 32] {
        let a = make_2d_poisson(size);
        let n = a.nrows;
        let b = vec![1.0f64; n];

        group.bench_with_input(
            BenchmarkId::new("amg_build", size * size),
            &size,
            |bench, _| {
                let amg = AmgClassical::new();
                bench.iter(|| black_box(amg.build(&a)));
            },
        );

        group.bench_with_input(
            BenchmarkId::new("amg_solve", size * size),
            &size,
            |bench, _| {
                let amg = AmgClassical::new();
                let hierarchy = amg.build(&a);
                bench.iter(|| {
                    let mut x = vec![0.0f64; n];
                    black_box(amg_solve(&hierarchy, &b, &mut x, CycleKind::V, 20, 1e-6))
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_amg_setup_solve);
criterion_main!(benches);
