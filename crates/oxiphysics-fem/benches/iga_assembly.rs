// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmark: local sparse IGA assembly and the Jacobi-PCG solve of a
//! plane-stress cantilever on a bi-quadratic NURBS patch.
//! Run with: cargo bench -p oxiphysics-fem --bench iga_assembly
//!
//! Every knot-span element only touches its `(p + 1)(q + 1)` active control
//! points, so the assembly time grows linearly with the number of elements
//! (40 x 40 control points = 3200 DOFs).

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use oxiphysics_fem::isogeometric_fem::{
    IgaAssembler, IgaBoundaryCondition, IgaMaterial, NurbsSurface,
};
use std::hint::black_box;

/// Control points per direction.
const NET_SIZES: [usize; 3] = [10, 20, 40];

fn material() -> IgaMaterial {
    IgaMaterial::new(2.0e5, 0.3, 0.05)
}

fn bench_assembly(c: &mut Criterion) {
    let mut group = c.benchmark_group("iga_assembly");
    group.sample_size(10);
    for n in NET_SIZES {
        let surface = NurbsSurface::flat_rectangle(4.0, 1.0, n, n, 2);
        let assembler = IgaAssembler::new(n * n, material(), 3);
        group.bench_with_input(BenchmarkId::new("assemble", n * n), &n, |bench, _| {
            bench.iter(|| black_box(assembler.assemble(&surface)));
        });
    }
    group.finish();
}

fn bench_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("iga_solve");
    group.sample_size(10);
    for n in NET_SIZES {
        let surface = NurbsSurface::flat_rectangle(4.0, 1.0, n, n, 2);
        let assembler = IgaAssembler::new(n * n, material(), 3);
        let Ok(mut system) = assembler.assemble(&surface) else {
            continue;
        };
        if assembler
            .apply_body_force(&surface, &mut system, [0.0, -1.0])
            .is_err()
        {
            continue;
        }
        // clamp the edge u = 0: control points (0, j) of the row-major net
        let clamped: Vec<usize> = (0..n)
            .flat_map(|j| [2 * (j * n), 2 * (j * n) + 1])
            .collect();
        let clamp = IgaBoundaryCondition::fixed(clamped);
        group.bench_with_input(
            BenchmarkId::new("clamp_and_solve", n * n),
            &n,
            |bench, _| {
                bench.iter(|| {
                    let mut constrained = system.clone();
                    clamp.apply(&mut constrained, 1.0e8).ok()?;
                    constrained.solve(50_000, 1.0e-8).ok()
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_assembly, bench_solve);
criterion_main!(benches);
