// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics benchmark demo: compare broadphase algorithms and time
//! core physics operations.
//!
//! Times brute-force vs sweep-and-prune broadphase on random AABBs,
//! plus rigid body pipeline stepping, to give a rough performance
//! overview using `std::time::Instant`.

use oxiphysics::pipeline::PhysicsPipeline;
use oxiphysics_collision::BroadPhase;
use oxiphysics_collision::broadphase::{BruteForceBroadPhase, SweepAndPrune};
use oxiphysics_core::math::Vec3;
use oxiphysics_core::{Aabb, Transform};
use oxiphysics_geometry::{Shape, Sphere};
use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use std::sync::Arc;
use std::time::Instant;

/// Simple deterministic LCG for reproducible benchmarks (no rand crate).
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = (self.0 >> 33) as f64;
        bits / (u32::MAX as f64 + 1.0)
    }

    /// Random f64 in \[lo, hi).
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
    }
}

fn generate_aabbs(n: usize, rng: &mut Lcg) -> Vec<Aabb> {
    (0..n)
        .map(|_| {
            let cx = rng.range(-50.0, 50.0);
            let cy = rng.range(-50.0, 50.0);
            let cz = rng.range(-50.0, 50.0);
            let half = rng.range(0.1, 2.0);
            Aabb::new(
                Vec3::new(cx - half, cy - half, cz - half),
                Vec3::new(cx + half, cy + half, cz + half),
            )
        })
        .collect()
}

fn bench_broadphase<B: BroadPhase>(name: &str, bp: &B, aabbs: &[Aabb], iterations: usize) {
    // Warm-up.
    let _ = bp.find_pairs(aabbs);

    let start = Instant::now();
    let mut total_pairs = 0_usize;
    for _ in 0..iterations {
        let pairs = bp.find_pairs(aabbs);
        total_pairs += pairs.len();
    }
    let elapsed = start.elapsed();

    let avg_us = elapsed.as_micros() as f64 / iterations as f64;
    let avg_pairs = total_pairs / iterations;
    println!(
        "  {:<20} {:>8.1} us/iter  ({} pairs, {} iters)",
        name, avg_us, avg_pairs, iterations
    );
}

fn bench_pipeline_step(n_bodies: usize, n_steps: usize) {
    let mut pipeline = PhysicsPipeline::new();
    let mut bodies = RigidBodySet::new();
    let mut colliders = ColliderSet::new();

    // Static floor.
    let mut floor = RigidBody::new_static();
    floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    let floor_h = bodies.insert(floor);
    let floor_shape: Arc<dyn Shape> = Arc::new(oxiphysics_geometry::BoxShape::new(Vec3::new(
        100.0, 0.5, 100.0,
    )));
    colliders.insert(Collider::new(floor_shape).with_body(floor_h));

    // Dynamic spheres spread in a grid.
    let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.25));
    let side = (n_bodies as f64).sqrt().ceil() as usize;
    for i in 0..n_bodies {
        let row = i / side;
        let col = i % side;
        let mut b = RigidBody::new(1.0);
        b.transform = Transform::from_position(Vec3::new(
            col as f64 * 0.7 - (side as f64 * 0.35),
            2.0 + (row as f64) * 0.6,
            0.0,
        ));
        let h = bodies.insert(b);
        colliders.insert(
            Collider::new(Arc::clone(&sphere_shape))
                .with_body(h)
                .with_restitution(0.3),
        );
    }

    let dt = 1.0 / 120.0;

    // Warm-up.
    for _ in 0..5 {
        pipeline.step(dt, &mut bodies, &colliders);
    }

    let start = Instant::now();
    for _ in 0..n_steps {
        pipeline.step(dt, &mut bodies, &colliders);
    }
    let elapsed = start.elapsed();

    let avg_us = elapsed.as_micros() as f64 / n_steps as f64;
    println!(
        "  Pipeline ({} bodies) {:>8.1} us/step  ({} steps)",
        n_bodies, avg_us, n_steps
    );
}

fn main() {
    println!("=== benchmark_demo: physics performance comparison ===\n");

    let mut rng = Lcg::new(42);

    // ── Broadphase comparison ───────────────────────────────────────
    for &n in &[100, 500, 1000] {
        let aabbs = generate_aabbs(n, &mut rng);
        let iterations = if n <= 500 { 200 } else { 50 };

        println!("Broadphase with {} AABBs:", n);
        bench_broadphase("Brute Force", &BruteForceBroadPhase, &aabbs, iterations);
        bench_broadphase(
            "Sweep & Prune (X)",
            &SweepAndPrune::x_axis(),
            &aabbs,
            iterations,
        );
        println!();
    }

    // ── Pipeline stepping ───────────────────────────────────────────
    println!("Pipeline stepping (rigid body sim):");
    bench_pipeline_step(10, 500);
    bench_pipeline_step(50, 200);
    bench_pipeline_step(100, 100);
    println!();

    // ── Summary ─────────────────────────────────────────────────────
    println!("NOTE: These are wall-clock micro-benchmarks for illustration.");
    println!("For rigorous benchmarks, use `cargo bench` with Criterion.");
}
