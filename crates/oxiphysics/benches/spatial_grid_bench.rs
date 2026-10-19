// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmarks for [`oxiphysics::spatial_grid::SpatialGrid`] with up to
//! 100 000 entries, next to naive linear scans over the same points to show the
//! gap the hash grid closes.
//!
//! The entries are spread uniformly at a constant density of one entry per unit
//! volume, so a query of radius `r` touches about `4.19 r³` entries however many
//! entries the grid holds.

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use oxiphysics::spatial_grid::SpatialGrid;
use std::hint::black_box;

/// Grid sizes of the scaling benchmarks.
const SIZES: [usize; 3] = [1_000, 10_000, 100_000];

/// Query radius (also the cell size of the radius/nearest benchmarks).
const RADIUS: f64 = 2.0;

/// Pair distance (also the cell size of the pair benchmark).
const PAIR_RADIUS: f64 = 1.0;

/// Deterministic xorshift64* generator.
struct Rng(u64);

impl Rng {
    fn next_unit(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let bits = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11;
        bits as f64 / (1_u64 << 53) as f64
    }

    fn point(&mut self, extent: f64) -> [f64; 3] {
        [
            self.next_unit() * extent,
            self.next_unit() * extent,
            self.next_unit() * extent,
        ]
    }
}

/// Edge length of the cube that holds `count` entries at unit density.
fn extent_for(count: usize) -> f64 {
    (count as f64).cbrt()
}

fn random_points(count: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = Rng(seed);
    let extent = extent_for(count);
    (0..count).map(|_| rng.point(extent)).collect()
}

/// Entry ID of the `index`-th point (the benchmark sizes all fit in `u32`).
fn entry_id(index: usize) -> u32 {
    u32::try_from(index).unwrap_or(u32::MAX)
}

fn build_grid(points: &[[f64; 3]], cell_size: f64) -> SpatialGrid {
    let mut grid = SpatialGrid::new(cell_size);
    for (index, &point) in points.iter().enumerate() {
        grid.insert(entry_id(index), point);
    }
    grid
}

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Query centres inside the point cloud, cycled through by the benchmarks.
fn probes(count: usize) -> Vec<[f64; 3]> {
    let mut rng = Rng(0xA11CE);
    let extent = extent_for(count);
    (0..64).map(|_| rng.point(extent)).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Building and updating
// ─────────────────────────────────────────────────────────────────────────────

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/build");
    for &n in &SIZES {
        let points = random_points(n, 1);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("insert_all", n), &points, |b, points| {
            b.iter(|| black_box(build_grid(black_box(points), RADIUS)));
        });
    }
    group.finish();
}

fn bench_update(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/update");
    let n = 100_000;
    let points = random_points(n, 2);
    let mut grid = build_grid(&points, RADIUS);
    let mut tick = 0_usize;
    group.throughput(Throughput::Elements(1_000));
    group.bench_function("move_1000_of_100000", |b| {
        b.iter(|| {
            tick += 1;
            let shift = if tick.is_multiple_of(2) { 0.75 } else { -0.75 };
            for index in (0..n).step_by(100) {
                let mut moved = points[index];
                moved[0] += shift;
                grid.update(entry_id(index), black_box(moved));
            }
        });
    });
    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Radius queries
// ─────────────────────────────────────────────────────────────────────────────

fn bench_query_radius(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/query_radius");
    for &n in &SIZES {
        let points = random_points(n, 3);
        let grid = build_grid(&points, RADIUS);
        let centers = probes(n);

        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("grid", n), &grid, |b, grid| {
            b.iter(|| {
                next += 1;
                let center = centers[next % centers.len()];
                black_box(grid.query_radius(black_box(center), RADIUS))
            });
        });

        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("linear_scan", n), &points, |b, points| {
            b.iter(|| {
                next += 1;
                let center = centers[next % centers.len()];
                let r2 = RADIUS * RADIUS;
                let hits: Vec<u32> = points
                    .iter()
                    .enumerate()
                    .filter(|&(_, &p)| dist_sq(p, black_box(center)) <= r2)
                    .map(|(index, _)| entry_id(index))
                    .collect();
                black_box(hits)
            });
        });
    }
    group.finish();
}

fn bench_query_aabb(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/query_aabb");
    for &n in &SIZES {
        let points = random_points(n, 4);
        let grid = build_grid(&points, RADIUS);
        let centers = probes(n);
        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("grid", n), &grid, |b, grid| {
            b.iter(|| {
                next += 1;
                let c = centers[next % centers.len()];
                let min = [c[0] - RADIUS, c[1] - RADIUS, c[2] - RADIUS];
                let max = [c[0] + RADIUS, c[1] + RADIUS, c[2] + RADIUS];
                black_box(grid.query_aabb(black_box(min), black_box(max)))
            });
        });
    }
    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Nearest neighbours
// ─────────────────────────────────────────────────────────────────────────────

fn bench_nearest(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/nearest");
    for &n in &SIZES {
        let points = random_points(n, 5);
        let grid = build_grid(&points, RADIUS);
        let centers = probes(n);

        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("grid", n), &grid, |b, grid| {
            b.iter(|| {
                next += 1;
                black_box(grid.nearest(black_box(centers[next % centers.len()])))
            });
        });

        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("linear_scan", n), &points, |b, points| {
            b.iter(|| {
                next += 1;
                let center = black_box(centers[next % centers.len()]);
                let best = points
                    .iter()
                    .enumerate()
                    .map(|(index, &p)| (dist_sq(p, center), index))
                    .min_by(|a, b| a.0.total_cmp(&b.0));
                black_box(best)
            });
        });
    }
    group.finish();
}

fn bench_k_nearest(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/k_nearest_8");
    for &n in &SIZES {
        let points = random_points(n, 6);
        let grid = build_grid(&points, RADIUS);
        let centers = probes(n);
        let mut next = 0_usize;
        group.bench_with_input(BenchmarkId::new("grid", n), &grid, |b, grid| {
            b.iter(|| {
                next += 1;
                black_box(grid.k_nearest(black_box(centers[next % centers.len()]), 8))
            });
        });
    }
    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Pairs (broad phase)
// ─────────────────────────────────────────────────────────────────────────────

fn bench_pairs(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial_grid/pairs_within_radius");
    group.sample_size(10);
    for &n in &SIZES {
        let points = random_points(n, 7);
        let grid = build_grid(&points, PAIR_RADIUS);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("grid", n), &grid, |b, grid| {
            b.iter(|| black_box(grid.pairs_within_radius(black_box(PAIR_RADIUS))));
        });

        // The quadratic scan is only practical for the smaller sizes (at 100 000
        // entries it needs 5 * 10^9 distance tests per iteration).
        if n <= 10_000 {
            group.bench_with_input(BenchmarkId::new("all_pairs", n), &points, |b, points| {
                b.iter(|| {
                    let r2 = PAIR_RADIUS * PAIR_RADIUS;
                    let mut pairs = Vec::new();
                    for (i, &a) in points.iter().enumerate() {
                        for (j, &p) in points.iter().enumerate().skip(i + 1) {
                            if dist_sq(a, p) <= r2 {
                                pairs.push((entry_id(i), entry_id(j)));
                            }
                        }
                    }
                    black_box(pairs)
                });
            });
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_build,
    bench_update,
    bench_query_radius,
    bench_query_aabb,
    bench_nearest,
    bench_k_nearest,
    bench_pairs,
);
criterion_main!(benches);
