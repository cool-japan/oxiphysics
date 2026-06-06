# oxiphysics-collision TODO

**Version:** 0.1.2 | **Updated:** 2026-06-06 | **Status:** ✅ COMPLETE

## Milestone 1: Types & Foundation — ✅ COMPLETE
- [x] Core collision types (`types`)
- [x] Basic error handling
- [x] Unit tests (2,443 tests passing)

## Milestone 2: Broad Phase — ✅ COMPLETE
- [x] Sweep-and-Prune (`sap`, `sweep`)
- [x] Dynamic BVH (`dbvt`)
- [x] k-d tree collision (`kdtree_collision`)
- [x] Broad-phase dispatcher (`broadphase`)
- [x] Parallel broad phase (`parallel_collision`)

## Milestone 3: Narrow Phase — ✅ COMPLETE
- [x] GJK + EPA (`gjk_epa`)
- [x] Enhanced GJK (`gjk_enhanced`)
- [x] Extended GJK (`gjk_extended`)
- [x] SAT collision (`sat_collision`)
- [x] Narrow-phase dispatcher (`narrowphase`, `NarrowPhaseDispatcher`)
- [x] Batch narrow phase (`BatchNarrowPhase`)
- [x] Compound shape dispatch (`compound_shapes`)

## Milestone 4: Contact Management — ✅ COMPLETE
- [x] Contact generation (`contact_generation`)
- [x] Contact manifold (`contact_manifold`, `ContactManifold`)
- [x] Contact graph (`contact_graph`, `ContactGraph`)
- [x] Manifold cache (`manifold_cache`)
- [x] Proximity detection (`proximity`)
- [x] Collision filtering (`CollisionFilter`)

## Milestone 5: CCD & Ray Casting — ✅ COMPLETE
- [x] Continuous collision detection (`ccd`, `CcdPipeline`)
- [x] Shape casting (`shape_cast`)
- [x] Ray casting (`ray_casting`): `ray_aabb`, `ray_sphere`, `ray_triangle`
- [x] Spatial queries (`spatial_queries`)

## Milestone 6: Deformable, Mesh & Voxel — ✅ COMPLETE
- [x] Deformable collision (`deformable_collision`)
- [x] Soft-body collision (`soft_body_collision`)
- [x] Mesh collision (`mesh_collision`)
- [x] Voxel collision (`voxel_collision`)
- [x] Terrain collision (`terrain_collision`)

## Milestone 7: Polish — ✅ COMPLETE
- [x] Documentation (rustdoc, 0 warnings)
- [x] Examples
- [x] Benchmarks (criterion)
- [x] 0 stubs
