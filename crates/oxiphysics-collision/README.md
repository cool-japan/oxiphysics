# oxiphysics-collision

**[Stable]** High-performance collision detection pipeline for the OxiPhysics engine.

[![Tests](https://img.shields.io/badge/tests-2443-brightgreen)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-collision)](https://docs.rs/oxiphysics_collision)
[![version](https://img.shields.io/badge/version-0.1.3-blue)](https://crates.io/crates/oxiphysics-collision)

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-collision` implements the full collision detection stack: broad-phase pruning,
narrow-phase exactalgorithms, continuous collision detection, contact generation and caching,
deformable/soft-body support, and parallel batch processing. 34+ modules, 2,625 public items,
2,443 passing tests, 0 stubs.

## Features

- **Broad phase** — Sweep-and-Prune (`SweepAndPrune`), Dynamic BVH (`DynamicBvh`), SAP and DBVT variants
- **Narrow phase** — GJK+EPA (`Gjk`, `Epa`), enhanced/extended GJK, SAT, `NarrowPhaseDispatcher`, `BatchNarrowPhase`
- **Continuous collision** — `CcdPipeline`, shape casting (`shape_cast`)
- **Contact management** — `ContactManifold`, `ContactGraph`, manifold cache, contact generation
- **Ray queries** — `ray_aabb`, `ray_sphere`, `ray_triangle`; full spatial query API
- **Deformable & soft-body** — `deformable_collision`, `soft_body_collision`
- **Mesh & voxel** — `mesh_collision`, `voxel_collision`, `terrain_collision`
- **Filtering** — `CollisionFilter`, proximity detection, compound shape dispatch
- **Parallel** — `parallel_collision` for multi-threaded broad/narrow phase
- **k-d tree** — `kdtree_collision` for static scene queries

## Modules (34+)

`broadphase` · `ccd` · `compound_shapes` · `contact_generation` · `contact_graph` ·
`contact_manifold` · `dbvt` · `deformable_collision` · `gjk_epa` · `gjk_enhanced` ·
`gjk_extended` · `kdtree_collision` · `manifold_cache` · `mesh_collision` · `narrowphase` ·
`parallel_collision` · `proximity` · `ray_casting` · `sap` · `sat_collision` · `shape_cast` ·
`soft_body_collision` · `spatial_queries` · `sweep` · `terrain_collision` · `types` ·
`voxel_collision` · …

## Quick Example

The Rust block below is compiled and run by `cargo test --doc`, so it cannot
drift away from the real API.

```rust
use oxiphysics_collision::narrowphase::ShapeType;
use oxiphysics_collision::{BroadPhase, CollisionPair, Gjk, NarrowPhaseDispatcher, SweepAndPrune};
use oxiphysics_core::math::Vec3;
use oxiphysics_core::{Aabb, Transform};
use oxiphysics_geometry::{BoxShape, Sphere};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // A unit sphere at the origin and a 1 m cube whose centre is 1.2 m away on x.
    let sphere = Sphere::new(1.0);
    let cube = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let sphere_pose = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
    let cube_pose = Transform::from_position(Vec3::new(1.2, 0.0, 0.0));

    // Broad phase: sweep-and-prune over the world-space AABBs yields candidate pairs.
    let aabbs = [
        Aabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0)),
        Aabb::new(Vec3::new(0.7, -0.5, -0.5), Vec3::new(1.7, 0.5, 0.5)),
    ];
    let pairs = SweepAndPrune::x_axis().find_pairs(&aabbs);
    assert_eq!(pairs, vec![CollisionPair::new(0, 1)]);

    // Narrow phase, yes/no: GJK on the posed shapes.
    assert!(Gjk::intersect(&sphere, &sphere_pose, &cube, &cube_pose));

    // Narrow phase, contacts: the dispatcher picks an analytic sphere/box test or
    // falls back to GJK + EPA, and returns a `ContactManifold` when the shapes touch.
    let dispatcher = NarrowPhaseDispatcher::default();
    let result = dispatcher.dispatch(
        &sphere,
        ShapeType::Sphere,
        &sphere_pose,
        &cube,
        ShapeType::Box,
        &cube_pose,
        pairs[0],
    );
    let manifold = result.manifold.ok_or("overlapping shapes must produce a manifold")?;
    println!(
        "{} contact point(s), deepest penetration {:.4}",
        manifold.len(),
        manifold.max_depth()
    );
    assert!(manifold.max_depth() > 0.0);
    Ok(())
}
```

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team KitaSan)
