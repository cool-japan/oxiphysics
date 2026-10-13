// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmark: periodic LBM streaming.
//!
//! Every group compares the production [`stream_2d`] / [`stream_3d`] (span
//! copies into the grid's persistent scratch buffer, followed by a buffer
//! swap) against the per-cell reference that was used before the scratch
//! buffer existed: it allocates and zero-fills a fresh `q × n` array on every
//! call and wraps every index with `rem_euclid`.
//!
//! Run with: `cargo bench -p oxiphysics-lbm --bench streaming`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use oxiphysics_lbm::grid::{LbmGrid2D, LbmGrid3D};
use oxiphysics_lbm::lattice::LatticeType;
use oxiphysics_lbm::streaming::{stream_2d, stream_3d};
use std::hint::black_box;

/// Allocating, `rem_euclid`-per-cell streaming of a 2-D grid (the old code).
fn reference_stream_2d(grid: &mut LbmGrid2D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    let n = nx * ny;
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];
    for (i, f_new_i) in f_new.iter_mut().enumerate() {
        let c = grid.lattice.velocity_2d(i);
        for y in 0..ny {
            for x in 0..nx {
                let src_x = (x as i64 - i64::from(c[0])).rem_euclid(nx as i64) as usize;
                let src_y = (y as i64 - i64::from(c[1])).rem_euclid(ny as i64) as usize;
                f_new_i[y * nx + x] = grid.f[i][src_y * nx + src_x];
            }
        }
    }
    grid.f = f_new;
}

/// Allocating, `rem_euclid`-per-cell streaming of a 3-D grid (the old code).
fn reference_stream_3d(grid: &mut LbmGrid3D) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let q = grid.lattice.q();
    let n = nx * ny * nz;
    let mut f_new: Vec<Vec<f64>> = vec![vec![0.0; n]; q];
    for (i, f_new_i) in f_new.iter_mut().enumerate() {
        let c = grid.lattice.velocity_3d(i);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let src_x = (x as i64 - i64::from(c[0])).rem_euclid(nx as i64) as usize;
                    let src_y = (y as i64 - i64::from(c[1])).rem_euclid(ny as i64) as usize;
                    let src_z = (z as i64 - i64::from(c[2])).rem_euclid(nz as i64) as usize;
                    f_new_i[z * ny * nx + y * nx + x] =
                        grid.f[i][src_z * ny * nx + src_y * nx + src_x];
                }
            }
        }
    }
    grid.f = f_new;
}

fn bench_stream_2d(c: &mut Criterion) {
    let mut group = c.benchmark_group("stream_2d_d2q9");
    for size in [128usize, 256] {
        let cells = (size * size) as u64;
        let mut grid = LbmGrid2D::new(size, size, LatticeType::D2Q9);
        group.throughput(criterion::Throughput::Elements(cells * 9));
        group.bench_with_input(
            BenchmarkId::new("reference_alloc_rem_euclid", size),
            &size,
            |b, _| b.iter(|| reference_stream_2d(black_box(&mut grid))),
        );
        let mut grid = LbmGrid2D::new(size, size, LatticeType::D2Q9);
        group.bench_with_input(BenchmarkId::new("span_copy_swap", size), &size, |b, _| {
            b.iter(|| stream_2d(black_box(&mut grid)))
        });
    }
    group.finish();
}

fn bench_stream_3d(c: &mut Criterion) {
    let mut group = c.benchmark_group("stream_3d_d3q19");
    group.sample_size(10);
    for size in [32usize, 64] {
        let cells = (size * size * size) as u64;
        let mut grid = LbmGrid3D::new(size, size, size, LatticeType::D3Q19);
        group.throughput(criterion::Throughput::Elements(cells * 19));
        group.bench_with_input(
            BenchmarkId::new("reference_alloc_rem_euclid", size),
            &size,
            |b, _| b.iter(|| reference_stream_3d(black_box(&mut grid))),
        );
        let mut grid = LbmGrid3D::new(size, size, size, LatticeType::D3Q19);
        group.bench_with_input(BenchmarkId::new("span_copy_swap", size), &size, |b, _| {
            b.iter(|| stream_3d(black_box(&mut grid)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_stream_2d, bench_stream_3d);
criterion_main!(benches);
