// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Criterion benchmarks for the OxiPhysics engine.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use oxiphysics::pipeline::PhysicsPipeline;
use oxiphysics_collision::broadphase::{BruteForceBroadPhase, SweepAndPrune};
use oxiphysics_collision::gjk_epa::{GjkBox, GjkSphere, gjk_distance, gjk_intersect};
use oxiphysics_collision::{BroadPhase, CollisionPair, NarrowPhaseDispatcher};
use oxiphysics_core::math::Vec3;
use oxiphysics_core::{Aabb, Transform};
use oxiphysics_fem::assembly::assemble_stiffness;
use oxiphysics_fem::constitutive::LinearElasticMaterial;
use oxiphysics_fem::mesh::TetrahedralMesh;
use oxiphysics_geometry::{BoxShape, Shape, Sphere};
use oxiphysics_lbm::{
    LatticeType, LbmGrid2D, LbmGrid3D, bgk_collide_2d, bgk_collide_3d, stream_2d, stream_3d,
};
use oxiphysics_md::atom::AtomSet;
use oxiphysics_md::forcefield::{ForceField, PairForceField};
use oxiphysics_md::neighbor::PeriodicBox;
use oxiphysics_md::potential::LennardJones;
use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
use oxiphysics_sph::kernel::CubicSplineKernel;
use oxiphysics_sph::particle::{ParticleSet, SphParticle};
use oxiphysics_sph::wcsph::compute_density;
use std::hint::black_box;
use std::sync::Arc;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn make_aabb(cx: f64, cy: f64, cz: f64, half: f64) -> Aabb {
    Aabb::new(
        Vec3::new(cx - half, cy - half, cz - half),
        Vec3::new(cx + half, cy + half, cz + half),
    )
}

/// Build `n` AABBs spread on a regular grid so roughly half overlap.
fn build_aabbs(n: usize) -> Vec<Aabb> {
    let side = (n as f64).cbrt().ceil() as usize;
    let spacing = 1.5;
    let half = 0.6;
    let mut aabbs = Vec::with_capacity(n);
    'outer: for iz in 0..side {
        for iy in 0..side {
            for ix in 0..side {
                aabbs.push(make_aabb(
                    ix as f64 * spacing,
                    iy as f64 * spacing,
                    iz as f64 * spacing,
                    half,
                ));
                if aabbs.len() == n {
                    break 'outer;
                }
            }
        }
    }
    aabbs
}

// ─────────────────────────────────────────────────────────────────────────────
// Broadphase benchmarks
// ─────────────────────────────────────────────────────────────────────────────

fn bench_broadphase(c: &mut Criterion) {
    let mut group = c.benchmark_group("broadphase");

    for &n in &[100_usize, 1000_usize] {
        let aabbs = build_aabbs(n);

        group.bench_with_input(BenchmarkId::new("SweepAndPrune", n), &aabbs, |b, aabbs| {
            let sap = SweepAndPrune::default();
            b.iter(|| {
                let pairs = sap.find_pairs(black_box(aabbs));
                black_box(pairs)
            });
        });

        group.bench_with_input(BenchmarkId::new("BruteForce", n), &aabbs, |b, aabbs| {
            let bf = BruteForceBroadPhase;
            b.iter(|| {
                let pairs = bf.find_pairs(black_box(aabbs));
                black_box(pairs)
            });
        });
    }

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Broadphase SAP incremental benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn bench_broadphase_sap(c: &mut Criterion) {
    let mut group = c.benchmark_group("broadphase_sap");

    // Benchmark SAP on a spread of sizes
    for &n in &[200_usize, 500_usize, 2000_usize] {
        let aabbs = build_aabbs(n);

        group.bench_with_input(
            BenchmarkId::new("SweepAndPrune_pairs", n),
            &aabbs,
            |b, aabbs| {
                let sap = SweepAndPrune::default();
                b.iter(|| {
                    let pairs = sap.find_pairs(black_box(aabbs));
                    black_box(pairs)
                });
            },
        );

        // Benchmark just the find_pairs call on a freshly-built input each iteration.
        group.bench_with_input(BenchmarkId::new("SweepAndPrune_rebuild", n), &n, |b, &n| {
            b.iter(|| {
                let aabbs = build_aabbs(black_box(n));
                let sap = SweepAndPrune::default();
                let pairs = sap.find_pairs(&aabbs);
                black_box(pairs)
            });
        });
    }

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Narrowphase (GJK dispatcher) benchmarks
// ─────────────────────────────────────────────────────────────────────────────

fn bench_narrowphase(c: &mut Criterion) {
    let mut group = c.benchmark_group("narrowphase");

    let n_pairs = 100;

    // sphere-sphere
    let sphere_a = Sphere::new(0.5);
    let sphere_b = Sphere::new(0.5);
    let transforms_ss: Vec<(Transform, Transform)> = (0..n_pairs)
        .map(|i| {
            let offset = 0.8 + (i as f64) * 0.001;
            (
                Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
                Transform::from_position(Vec3::new(offset, 0.0, 0.0)),
            )
        })
        .collect();

    group.bench_function("sphere_sphere_100pairs", |b| {
        b.iter(|| {
            for (ta, tb) in black_box(&transforms_ss) {
                let result = NarrowPhaseDispatcher::generate_contacts(
                    &sphere_a,
                    ta,
                    &sphere_b,
                    tb,
                    CollisionPair::new(0, 1),
                );
                black_box(result);
            }
        });
    });

    // sphere-box
    let box_shape = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let transforms_sb: Vec<(Transform, Transform)> = (0..n_pairs)
        .map(|i| {
            let offset = 0.8 + (i as f64) * 0.001;
            (
                Transform::from_position(Vec3::new(0.0, offset, 0.0)),
                Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
            )
        })
        .collect();

    group.bench_function("sphere_box_100pairs", |b| {
        b.iter(|| {
            for (ta, tb) in black_box(&transforms_sb) {
                let result = NarrowPhaseDispatcher::generate_contacts(
                    &sphere_a,
                    ta,
                    &box_shape,
                    tb,
                    CollisionPair::new(0, 1),
                );
                black_box(result);
            }
        });
    });

    // box-box
    let box_a = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let box_b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let transforms_bb: Vec<(Transform, Transform)> = (0..n_pairs)
        .map(|i| {
            let offset = 0.8 + (i as f64) * 0.001;
            (
                Transform::from_position(Vec3::new(0.0, 0.0, 0.0)),
                Transform::from_position(Vec3::new(offset, 0.0, 0.0)),
            )
        })
        .collect();

    group.bench_function("box_box_100pairs", |b| {
        b.iter(|| {
            for (ta, tb) in black_box(&transforms_bb) {
                let result = NarrowPhaseDispatcher::generate_contacts(
                    &box_a,
                    ta,
                    &box_b,
                    tb,
                    CollisionPair::new(0, 1),
                );
                black_box(result);
            }
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// GJK distance / intersection benchmarks
// ─────────────────────────────────────────────────────────────────────────────

fn bench_gjk(c: &mut Criterion) {
    let mut group = c.benchmark_group("gjk");

    // Two separated spheres.
    let sa = GjkSphere {
        center: [0.0, 0.0, 0.0],
        radius: 0.5,
    };
    let sb_sep = GjkSphere {
        center: [2.0, 0.0, 0.0],
        radius: 0.5,
    };
    // Two overlapping spheres.
    let sb_ov = GjkSphere {
        center: [0.6, 0.0, 0.0],
        radius: 0.5,
    };

    group.bench_function("gjk_intersect_sphere_sphere_separated", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&sa), black_box(&sb_sep));
            black_box(result)
        });
    });

    group.bench_function("gjk_intersect_sphere_sphere_overlapping", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&sa), black_box(&sb_ov));
            black_box(result)
        });
    });

    group.bench_function("gjk_distance_sphere_sphere_separated", |b| {
        b.iter(|| {
            let result = gjk_distance(black_box(&sa), black_box(&sb_sep));
            black_box(result)
        });
    });

    // Box vs sphere
    let ba = GjkBox {
        center: [0.0, 0.0, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };
    let sphere_far = GjkSphere {
        center: [3.0, 0.0, 0.0],
        radius: 0.3,
    };
    let sphere_near = GjkSphere {
        center: [0.7, 0.0, 0.0],
        radius: 0.3,
    };

    group.bench_function("gjk_intersect_box_sphere_separated", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&ba), black_box(&sphere_far));
            black_box(result)
        });
    });

    group.bench_function("gjk_intersect_box_sphere_overlapping", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&ba), black_box(&sphere_near));
            black_box(result)
        });
    });

    group.bench_function("gjk_distance_box_sphere_separated", |b| {
        b.iter(|| {
            let result = gjk_distance(black_box(&ba), black_box(&sphere_far));
            black_box(result)
        });
    });

    // Box vs box
    let bb_sep = GjkBox {
        center: [3.0, 0.0, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };
    let bb_ov = GjkBox {
        center: [0.8, 0.0, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };

    group.bench_function("gjk_intersect_box_box_separated", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&ba), black_box(&bb_sep));
            black_box(result)
        });
    });

    group.bench_function("gjk_intersect_box_box_overlapping", |b| {
        b.iter(|| {
            let result = gjk_intersect(black_box(&ba), black_box(&bb_ov));
            black_box(result)
        });
    });

    // Batch GJK over 200 sphere pairs
    let pairs_200: Vec<(GjkSphere, GjkSphere)> = (0..200)
        .map(|i| {
            let x = i as f64 * 0.01;
            (
                GjkSphere {
                    center: [0.0, 0.0, 0.0],
                    radius: 0.5,
                },
                GjkSphere {
                    center: [x + 0.8, 0.0, 0.0],
                    radius: 0.3,
                },
            )
        })
        .collect();

    group.bench_function("gjk_intersect_200_sphere_pairs", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for (s1, s2) in black_box(&pairs_200) {
                if gjk_intersect(s1, s2) {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Collision detection pipeline benchmark (narrowphase dispatcher only)
// ─────────────────────────────────────────────────────────────────────────────

fn bench_collision_detection_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("collision_detection_pipeline");

    let n_pairs = 500;

    // Prepare overlapping sphere pairs at a fixed slight overlap.
    let sphere = Arc::new(Sphere::new(0.5));
    let pairs_sphere: Vec<(Transform, Transform)> = (0..n_pairs)
        .map(|i| {
            let z = i as f64 * 5.0;
            (
                Transform::from_position(Vec3::new(0.0, 0.0, z)),
                Transform::from_position(Vec3::new(0.82, 0.0, z)),
            )
        })
        .collect();

    let sphere_ref = sphere.as_ref();

    group.bench_function("dispatcher_sphere_sphere_500pairs", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for (ta, tb) in black_box(&pairs_sphere) {
                if NarrowPhaseDispatcher::generate_contacts(
                    sphere_ref,
                    ta,
                    sphere_ref,
                    tb,
                    CollisionPair::new(0, 1),
                )
                .is_some()
                {
                    hits += 1;
                }
            }
            black_box(hits)
        });
    });

    // Prepare overlapping box pairs.
    let box_shape = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
    let pairs_box: Vec<(Transform, Transform)> = (0..n_pairs)
        .map(|i| {
            let z = i as f64 * 5.0;
            (
                Transform::from_position(Vec3::new(0.0, 0.0, z)),
                Transform::from_position(Vec3::new(0.82, 0.0, z)),
            )
        })
        .collect();
    let box_ref = box_shape.as_ref();

    group.bench_function("dispatcher_box_box_500pairs", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for (ta, tb) in black_box(&pairs_box) {
                if NarrowPhaseDispatcher::generate_contacts(
                    box_ref,
                    ta,
                    box_ref,
                    tb,
                    CollisionPair::new(0, 1),
                )
                .is_some()
                {
                    hits += 1;
                }
            }
            black_box(hits)
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Rigid body pipeline benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn build_falling_sphere_scene(n: usize) -> (PhysicsPipeline, RigidBodySet, ColliderSet) {
    let mut pipeline = PhysicsPipeline::new();
    pipeline.config.gravity = Vec3::new(0.0, -9.81, 0.0);

    let mut bodies = RigidBodySet::new();
    let mut colliders = ColliderSet::new();

    // Static floor
    let mut floor = RigidBody::new_static();
    floor.transform = Transform::from_position(Vec3::new(0.0, -1.0, 0.0));
    let fh = bodies.insert(floor);
    let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(20.0, 0.5, 20.0)));
    colliders.insert(Collider::new(floor_shape).with_body(fh));

    // Dynamic spheres
    let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.3));
    let side = (n as f64).sqrt().ceil() as usize;
    let mut count = 0;
    'outer: for ix in 0..side {
        for iz in 0..side {
            let mut body = RigidBody::new(1.0);
            body.transform = Transform::from_position(Vec3::new(
                ix as f64 * 1.5 - side as f64 * 0.75,
                5.0 + (count as f64) * 0.1,
                iz as f64 * 1.5 - side as f64 * 0.75,
            ));
            let bh = bodies.insert(body);
            colliders.insert(Collider::new(Arc::clone(&sphere_shape)).with_body(bh));
            count += 1;
            if count == n {
                break 'outer;
            }
        }
    }

    (pipeline, bodies, colliders)
}

fn bench_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("pipeline");

    group.bench_function("step_50_dynamic_spheres", |b| {
        let (mut pipeline, mut bodies, colliders) = build_falling_sphere_scene(50);
        b.iter(|| {
            pipeline.step(black_box(1.0 / 60.0), &mut bodies, &colliders);
        });
    });

    group.bench_function("step_100_dynamic_spheres", |b| {
        let (mut pipeline, mut bodies, colliders) = build_falling_sphere_scene(100);
        b.iter(|| {
            pipeline.step(black_box(1.0 / 60.0), &mut bodies, &colliders);
        });
    });

    group.bench_function("step_200_dynamic_spheres", |b| {
        let (mut pipeline, mut bodies, colliders) = build_falling_sphere_scene(200);
        b.iter(|| {
            pipeline.step(black_box(1.0 / 60.0), &mut bodies, &colliders);
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Constraint solving benchmark
// ─────────────────────────────────────────────────────────────────────────────

/// Benchmark the sequential impulse constraint solver embedded in the pipeline.
///
/// We create a stack of N spheres on a floor and time a full step
/// (which exercises the PGS contact solver with N active constraints).
fn bench_constraint_solving(c: &mut Criterion) {
    let mut group = c.benchmark_group("constraint_solving");

    for &n in &[5_usize, 20_usize, 50_usize] {
        group.bench_with_input(
            BenchmarkId::new("pgs_contacts_n_spheres_on_floor", n),
            &n,
            |b, &n| {
                let mut pipeline = PhysicsPipeline::new();
                pipeline.config.gravity = Vec3::new(0.0, -9.81, 0.0);
                pipeline.config.solver_iterations = 10_u32;

                let mut bodies = RigidBodySet::new();
                let mut colliders = ColliderSet::new();

                // Static floor
                let mut floor = RigidBody::new_static();
                floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
                let fh = bodies.insert(floor);
                let floor_shape: Arc<dyn Shape> =
                    Arc::new(BoxShape::new(Vec3::new(20.0, 0.5, 20.0)));
                colliders.insert(Collider::new(floor_shape).with_body(fh));

                // Stack of spheres (each just touching the one below).
                let r = 0.4;
                let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(r));
                for i in 0..n {
                    let mut body = RigidBody::new(1.0);
                    body.transform = Transform::from_position(Vec3::new(
                        0.0,
                        r * 2.0 * i as f64 + r - 0.02,
                        0.0,
                    ));
                    let bh = bodies.insert(body);
                    colliders.insert(Collider::new(Arc::clone(&sphere_shape)).with_body(bh));
                }

                // Warm up the warmstarting cache.
                for _ in 0..5 {
                    pipeline.step(1.0 / 60.0, &mut bodies, &colliders);
                }

                b.iter(|| {
                    pipeline.step(black_box(1.0 / 60.0), &mut bodies, &colliders);
                });
            },
        );
    }

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// SPH density benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn build_sph_particles(n: usize) -> (ParticleSet, Vec<Vec<usize>>) {
    let h = 0.15;
    let spacing = 0.1;
    let side = (n as f64).cbrt().ceil() as usize;

    let mut particles = ParticleSet::with_capacity(n);
    let mut count = 0;
    'outer: for iz in 0..side {
        for iy in 0..side {
            for ix in 0..side {
                let p = SphParticle::new(
                    Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    ),
                    Vec3::zeros(),
                    0.001,
                );
                particles.add_particle(&p);
                count += 1;
                if count == n {
                    break 'outer;
                }
            }
        }
    }

    // Build naive neighbor lists (within 2h)
    let r2 = (2.0 * h) * (2.0 * h);
    let n_actual = particles.len();
    let mut neighbors = vec![Vec::new(); n_actual];
    for (i, neighbor_list) in neighbors.iter_mut().enumerate() {
        for j in 0..n_actual {
            if i != j {
                let dr = particles.positions[i] - particles.positions[j];
                if dr.norm_squared() < r2 {
                    neighbor_list.push(j);
                }
            }
        }
    }

    (particles, neighbors)
}

fn bench_sph(c: &mut Criterion) {
    let mut group = c.benchmark_group("sph");

    let (particles_proto, neighbors) = build_sph_particles(500);
    let kernel = CubicSplineKernel;
    let h = 0.15;

    group.bench_function("density_500_particles", |b| {
        b.iter(|| {
            let mut particles = particles_proto.clone();
            compute_density(
                black_box(&mut particles),
                black_box(&neighbors),
                black_box(&kernel),
                black_box(h),
            );
            black_box(&particles.densities[0]);
        });
    });

    // Larger: 1000 particles
    let (particles_1k, neighbors_1k) = build_sph_particles(1000);

    group.bench_function("density_1000_particles", |b| {
        b.iter(|| {
            let mut particles = particles_1k.clone();
            compute_density(
                black_box(&mut particles),
                black_box(&neighbors_1k),
                black_box(&kernel),
                black_box(h),
            );
            black_box(&particles.densities[0]);
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// LBM collision step benchmark (2D and 3D)
// ─────────────────────────────────────────────────────────────────────────────

fn bench_lbm(c: &mut Criterion) {
    let mut group = c.benchmark_group("lbm");

    // 2D BGK collision step
    let nx2 = 64;
    let ny2 = 64;
    let omega = 1.0;

    group.bench_function("bgk_collide_2d_64x64", |b| {
        let mut grid = LbmGrid2D::new(nx2, ny2, LatticeType::D2Q9);
        grid.set_equilibrium(32, 32, 1.05, 0.02, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            bgk_collide_2d(black_box(&mut grid), black_box(omega));
        });
    });

    group.bench_function("stream_2d_64x64", |b| {
        let mut grid = LbmGrid2D::new(nx2, ny2, LatticeType::D2Q9);
        grid.set_equilibrium(32, 32, 1.0, 0.01, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            stream_2d(black_box(&mut grid));
        });
    });

    group.bench_function("bgk_collide_then_stream_2d_64x64", |b| {
        let mut grid = LbmGrid2D::new(nx2, ny2, LatticeType::D2Q9);
        grid.set_equilibrium(32, 32, 1.05, 0.02, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            bgk_collide_2d(black_box(&mut grid), black_box(omega));
            stream_2d(black_box(&mut grid));
        });
    });

    // 2D – larger grid
    let nx2_large = 128;
    let ny2_large = 128;

    group.bench_function("bgk_collide_2d_128x128", |b| {
        let mut grid = LbmGrid2D::new(nx2_large, ny2_large, LatticeType::D2Q9);
        grid.set_equilibrium(64, 64, 1.0, 0.01, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            bgk_collide_2d(black_box(&mut grid), black_box(omega));
        });
    });

    // 3D BGK collision step
    let nx3 = 16;
    let ny3 = 16;
    let nz3 = 16;

    group.bench_function("bgk_collide_3d_16x16x16", |b| {
        let mut grid = LbmGrid3D::new(nx3, ny3, nz3, LatticeType::D3Q19);
        grid.set_equilibrium(8, 8, 8, 1.05, 0.01, 0.0, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            bgk_collide_3d(black_box(&mut grid), black_box(omega));
        });
    });

    group.bench_function("stream_3d_16x16x16", |b| {
        let mut grid = LbmGrid3D::new(nx3, ny3, nz3, LatticeType::D3Q19);
        grid.set_equilibrium(8, 8, 8, 1.0, 0.01, 0.0, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            stream_3d(black_box(&mut grid));
        });
    });

    group.bench_function("bgk_collide_then_stream_3d_16x16x16", |b| {
        let mut grid = LbmGrid3D::new(nx3, ny3, nz3, LatticeType::D3Q19);
        grid.set_equilibrium(8, 8, 8, 1.05, 0.01, 0.0, 0.0);
        grid.compute_macroscopic();
        b.iter(|| {
            bgk_collide_3d(black_box(&mut grid), black_box(omega));
            stream_3d(black_box(&mut grid));
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// FEM stiffness assembly benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn bench_fem(c: &mut Criterion) {
    let mut group = c.benchmark_group("fem");

    // 10x2x2 beam mesh
    let mesh = TetrahedralMesh::generate_beam(1.0, 0.2, 0.2, 10, 2, 2);
    let material = LinearElasticMaterial::new(200e9, 0.3); // steel

    group.bench_function("assemble_stiffness_10x2x2_beam", |b| {
        b.iter(|| {
            let k = assemble_stiffness(black_box(&mesh), black_box(&material));
            black_box(k)
        });
    });

    // Larger mesh: 20x2x2
    let mesh_large = TetrahedralMesh::generate_beam(2.0, 0.2, 0.2, 20, 2, 2);

    group.bench_function("assemble_stiffness_20x2x2_beam", |b| {
        b.iter(|| {
            let k = assemble_stiffness(black_box(&mesh_large), black_box(&material));
            black_box(k)
        });
    });

    // Wider mesh: 10x4x4
    let mesh_wide = TetrahedralMesh::generate_beam(1.0, 0.4, 0.4, 10, 4, 4);

    group.bench_function("assemble_stiffness_10x4x4_beam", |b| {
        b.iter(|| {
            let k = assemble_stiffness(black_box(&mesh_wide), black_box(&material));
            black_box(k)
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// MD LJ force computation benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn build_md_system(n: usize) -> (AtomSet, PeriodicBox, PairForceField) {
    let box_len = (n as f64).cbrt() * 3.5; // ~3.5 Å spacing, argon-like
    let pbox = PeriodicBox::cubic(box_len);

    let mut atoms = AtomSet::with_capacity(n);
    let side = (n as f64).cbrt().ceil() as usize;
    let spacing = box_len / side as f64;
    let mut count = 0;

    'outer: for iz in 0..side {
        for iy in 0..side {
            for ix in 0..side {
                atoms.add_atom(
                    Vec3::new(
                        ix as f64 * spacing + 0.01,
                        iy as f64 * spacing + 0.01,
                        iz as f64 * spacing + 0.01,
                    ),
                    Vec3::zeros(),
                    39.948, // argon mass
                    0.0,
                    0,
                );
                count += 1;
                if count == n {
                    break 'outer;
                }
            }
        }
    }

    let lj = LennardJones::new(0.0104, 3.4, 8.5); // argon LJ params (eV, Å)
    let mut ff = PairForceField::new();
    ff.add_interaction(0, 0, lj);

    (atoms, pbox, ff)
}

fn bench_md(c: &mut Criterion) {
    let mut group = c.benchmark_group("md");

    let (atoms_proto, pbox, ff) = build_md_system(200);

    group.bench_function("lj_forces_200_atoms", |b| {
        b.iter(|| {
            let mut atoms = atoms_proto.clone();
            atoms.clear_forces();
            let energy = ff.compute_forces(black_box(&mut atoms), black_box(&pbox));
            black_box(energy)
        });
    });

    let (atoms_proto_500, pbox_500, ff_500) = build_md_system(500);

    group.bench_function("lj_forces_500_atoms", |b| {
        b.iter(|| {
            let mut atoms = atoms_proto_500.clone();
            atoms.clear_forces();
            let energy = ff_500.compute_forces(black_box(&mut atoms), black_box(&pbox_500));
            black_box(energy)
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// AABB operations benchmark
// ─────────────────────────────────────────────────────────────────────────────

fn bench_aabb(c: &mut Criterion) {
    let mut group = c.benchmark_group("aabb");

    const N: usize = 10_000;

    let aabbs_a: Vec<Aabb> = (0..N)
        .map(|i| make_aabb(i as f64 * 2.0, 0.0, 0.0, 0.5))
        .collect();
    let aabbs_b: Vec<Aabb> = (0..N)
        .map(|i| make_aabb(i as f64 * 2.0 + 0.5, 0.0, 0.0, 0.5))
        .collect();
    let points: Vec<Vec3> = (0..N)
        .map(|i| Vec3::new(i as f64 * 2.0 + 0.25, 0.0, 0.0))
        .collect();

    group.bench_function("merge_10000", |b| {
        b.iter(|| {
            let mut result = Aabb::new(Vec3::zeros(), Vec3::zeros());
            for (a, b_aabb) in black_box(&aabbs_a).iter().zip(black_box(&aabbs_b).iter()) {
                result = a.merge(b_aabb);
            }
            black_box(result)
        });
    });

    group.bench_function("contains_point_10000", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for (aabb, pt) in black_box(&aabbs_a).iter().zip(black_box(&points).iter()) {
                if aabb.contains_point(pt) {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.bench_function("intersects_10000", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for (a, b_aabb) in black_box(&aabbs_a).iter().zip(black_box(&aabbs_b).iter()) {
                if a.intersects(b_aabb) {
                    count += 1;
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

// ─────────────────────────────────────────────────────────────────────────────
// Criterion entry points
// ─────────────────────────────────────────────────────────────────────────────

criterion_group!(broadphase, bench_broadphase);
criterion_group!(broadphase_sap, bench_broadphase_sap);
criterion_group!(narrowphase, bench_narrowphase);
criterion_group!(gjk, bench_gjk);
criterion_group!(collision_pipeline, bench_collision_detection_pipeline);
criterion_group!(pipeline, bench_pipeline);
criterion_group!(constraint_solving, bench_constraint_solving);
criterion_group!(sph, bench_sph);
criterion_group!(lbm, bench_lbm);
criterion_group!(fem, bench_fem);
criterion_group!(md, bench_md);
criterion_group!(aabb, bench_aabb);

criterion_main!(
    broadphase,
    broadphase_sap,
    narrowphase,
    gjk,
    collision_pipeline,
    pipeline,
    constraint_solving,
    sph,
    lbm,
    fem,
    md,
    aabb
);
