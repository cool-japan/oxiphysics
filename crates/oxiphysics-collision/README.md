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

```rust
use oxiphysics_collision::{BroadPhase, SweepAndPrune, Gjk, ContactManifold};
use oxiphysics_geometry::{Sphere, BoxShape};

let mut broad = SweepAndPrune::new();
broad.insert(0, sphere_aabb);
broad.insert(1, box_aabb);

let pairs = broad.overlapping_pairs();
for (a, b) in pairs {
    let mut gjk = Gjk::new();
    if let Some(contact) = gjk.intersect(&shapes[a], &shapes[b]) {
        let manifold = ContactManifold::from_contact(contact);
        println!("penetration depth: {:.4}", manifold.depth());
    }
}
```

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
