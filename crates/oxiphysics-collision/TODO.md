# oxiphysics-collision TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 67,817 SLoC | 2,443 tests

## Completed (v0.1.0 – v0.1.3)

All v0.1.x history is preserved below; new work tracks under the version sections that follow.

### Milestone 1: Types & foundation
- [x] Core collision types (`types`)
- [x] Basic error handling
- [x] Unit tests (2,443 tests passing)

### Milestone 2: Broad phase
- [x] Sweep-and-Prune (`sap`, `sweep`)
- [x] Dynamic BVH (`dbvt`)
- [x] k-d tree collision (`kdtree_collision`)
- [x] Broad-phase dispatcher (`broadphase`)
- [x] Parallel broad phase (`parallel_collision`)

### Milestone 3: Narrow phase
- [x] GJK + EPA (`gjk_epa`)
- [x] Enhanced GJK (`gjk_enhanced`)
- [x] Extended GJK (`gjk_extended`)
- [x] SAT collision (`sat_collision`)
- [x] Narrow-phase dispatcher (`narrowphase`, `NarrowPhaseDispatcher`)
- [x] Batch narrow phase (`BatchNarrowPhase`)
- [x] Compound shape dispatch (`compound_shapes`)

### Milestone 4: Contact management
- [x] Contact generation (`contact_generation`)
- [x] Contact manifold (`contact_manifold`, `ContactManifold`)
- [x] Contact graph (`contact_graph`, `ContactGraph`)
- [x] Manifold cache (`manifold_cache`)
- [x] Proximity detection (`proximity`)
- [x] Collision filtering (`CollisionFilter`)

### Milestone 5: CCD & ray casting
- [x] Continuous collision detection (`ccd`, `CcdPipeline`)
- [x] Shape casting (`shape_cast`)
- [x] Ray casting (`ray_casting`): `ray_aabb`, `ray_sphere`, `ray_triangle`
- [x] Spatial queries (`spatial_queries`)

### Milestone 6: Deformable, mesh & voxel
- [x] Deformable collision (`deformable_collision`)
- [x] Soft-body collision (`soft_body_collision`)
- [x] Mesh collision (`mesh_collision`)
- [x] Voxel collision (`voxel_collision`)
- [x] Terrain collision (`terrain_collision`)

### Milestone 7: Polish
- [x] Documentation (rustdoc, 0 warnings)
- [x] Examples
- [x] Benchmarks (criterion)
- [x] 0 stubs

### Roadmap foundations verified in code (2026-06-11)
- [x] Incremental GJK warm-started support caching already exists (`narrowphase/gjk_cache/` — warmstartedgjk, simplexcache, supportcache, eviction policies); default-on dispatch wiring is the open v0.2.0 item
- [x] Manifold clipping primitives exist (`clip_polygon_by_plane`, `find_incident_face`, Sutherland-Hodgman at `contact_generation.rs:196-340`) and 4-point reduction exists (`reduce_to_4`, `contact_manifold.rs:352`) — dispatch wiring is the open v0.2.0 item
- [x] Discrete VF/EE proximity exists (`deformable_collision/functions.rs:145,225`) — the continuous time-of-impact path is the open v0.2.0 item
- [x] Shape-ordinal canonicalization helper exists (`narrowphase/dispatch/functions.rs:28`) — but it guards the unsafe downcasts only via `debug_assert!`, so release-build soundness is the PRIORITY v0.2.0 item below

## v0.2.0 — Narrowphase soundness and manifold quality

Theme: make the narrowphase dispatch memory-sound first, then upgrade contact quality — full one-shot manifolds, warm-started GJK by default, MPR fallback, and continuous VF/EE queries for soft bodies.
Exit gate: Miri-clean dispatch (local script), 4-point one-shot manifolds resting without rocking, and zero missed collisions on the cloth-vs-cloth CCD benchmark.

### Priority: dispatch soundness

- [x] **PRIORITY** — Fix narrowphase dispatch raw-pointer downcast UB (planned 2026-06-11)
  - **Goal:** zero `unsafe` downcasts in dispatch; Miri-clean; release builds safe for unordered shape pairs.
  - **Design:** verified real at `narrowphase/dispatch/functions.rs:67-76` — `&*(shape_a as *const dyn Shape as *const Sphere)` guarded only by `debug_assert!` (shape-ordinal canonicalization helper exists at line 28 but release builds skip the check → latent UB). Refactor `dispatch` to canonicalize/swap pairs first, then downcast via `Any::downcast_ref` or a closed `ShapeRef` enum; flip-normal on swapped output.
  - **Files:** `narrowphase/dispatch/functions.rs`
  - **Tests:** swapped-pair property test (manifold equals flipped-normal of the unswapped result); Miri run over dispatch paths via local script
  - **Risk:** dispatch is the hot path — benchmark before/after; a closed `ShapeRef` enum fixes the shape set, so document the extension story for user shapes
  - **Files:** oxiphysics-geometry/src/shape.rs (add `as_any` to Shape trait + all impls), collision narrowphase/dispatch/{functions,types}.rs, all `impl Shape` sites workspace-wide
  - **Tests:** reversed-order dispatch tests for every registered pair (manifold + normal sign both arg orders); degenerate-shape safe-fallback test; full 2,443 existing suite; targeted Miri if installed (skip gracefully)
  - **Risk:** Shape trait change ripples to every impl — grep `impl Shape for` first; plain `dispatch()` behavior change (canonicalization) — grep callers first

### Manifold quality & continuous queries

- [x] One-shot full-manifold generation in dispatch — box-box done (completed 2026-06-12)
  - **Goal:** box-box and convex-convex produce 4-point manifolds in a single query (no multi-frame accumulation); 1-box-on-ground rests without rocking in 1 solver iteration.
  - **Design:** clipping primitives already exist (`clip_polygon_by_plane`, `find_incident_face`, Sutherland-Hodgman at `contact_generation.rs:196-340`) and reduction exists (`reduce_to_4`, `contact_manifold.rs:352`) — wire them into the dispatch paths replacing single-point GJK/EPA results; feature-ID contact matching for warm-start continuity.
  - **Files:** `narrowphase/dispatch/`, `contact_generation.rs`, `contact_manifold.rs`
  - **Tests:** box-on-ground single-iteration rest; manifold point-count and feature-ID stability across frames
  - **Files:** narrowphase/dispatch/functions.rs (box-box arm), contact_generation.rs, contact_manifold.rs
  - **Tests:** box-on-ground 4-point manifold rests in 1 iter (no rocking); tilted-edge 2-point; feature-ID stable; normal sign both dispatch orders
  - **Risk:** SAT axis/face mapping must be exact — unit-test face selection first

- [x] Convex-convex one-shot full manifolds (general polyhedra) (completed 2026-06-12)
  - **Goal:** extend the box-box face-clip manifold (shipped 2026-06-12 in `narrowphase/dispatch/box_manifold.rs`) to arbitrary convex hulls — Sutherland-Hodgman reference/incident face clipping for ConvexHull vs ConvexHull and ConvexHull vs Box, replacing single-point GJK/EPA results.
  - **Done:** new `narrowphase/dispatch/convex_manifold.rs` exposes `polyhedron_manifold_from_normal` + `convex_pair_manifold`; box-box face path now delegates to the shared helper (regression-tested identical). ConvexHull×ConvexHull and Box×ConvexHull dispatch arms registered, producing full 4-point manifolds from the EPA normal, with an edge-edge witness fallback (`FACE_ALIGN_MIN` heuristic) and curved shapes kept single-point. Face data is extracted robustly via direct supporting-plane enumeration (`extract_convex_faces`) rather than the triangle-soup `ConvexHull3D` builder, which produced malformed hulls (dropped vertices, duplicate triangles) for prisms.
  - **Files:** `narrowphase/dispatch/convex_manifold.rs` (new), `narrowphase/dispatch/box_manifold.rs` (face path delegates), `narrowphase/dispatch/functions.rs` (`convex_manifold_dispatch`, `world_hull_vertices`), `narrowphase/dispatch/narrowphasedispatcher_traits.rs` (registration)
  - **Tests:** `tests/convex_manifold.rs` (box-box regression, tetra-on-box 3-point + COM-stationary, hull-hull face-to-face, edge-edge, both dispatch orders) + source unit tests (hex-prism 8-face extraction, cube merge, witness consistency). Full suite 2472/2472 green; clippy clean.

- [x] Enable the GJK warm-start cache by default in the dispatch pipeline (planned 2026-06-11)
  - **Goal:** >90% warm-start hit rate on coherent scenes, measured by a new hit-rate benchmark.
  - **Design:** the cache itself already exists (`narrowphase/gjk_cache/` — warmstartedgjk, simplexcache, supportcache, eviction policies); this item is the default-on wiring in dispatch plus a criterion hit-rate benchmark.
  - **Files:** `narrowphase/gjk_cache/`, `narrowphase/dispatch/`
  - **Tests:** hit-rate benchmark on a coherent tumbling-bodies scene; correctness parity cold-vs-warm
  - **Files:** narrowphase/dispatch/types.rs (DispatchConfig + dispatcher GJK fallback path), gjk_cache integration
  - **Tests:** coherent-motion 60-frame scene asserting hit-ratio > 50%; cache invalidation on body removal
  - **Risk:** dispatcher may need &mut/state for the registry — keep API additive

- [x] MPR (Minkowski Portal Refinement) (done 2026-06-14)
  - **Goal:** penetration depth within 1% of EPA on 1e5 random convex pairs at >2x EPA speed; used as EPA fallback on simplex degeneracy.
  - **Design:** Snethen XenoCollide portal discovery + EPA minimum-penetration refinement.
  - **Tests:** 1e5 random convex-pair depth comparison vs EPA; degenerate-simplex fallback trigger cases
    - **Files:** NEW narrowphase/gjk/mpr.rs (full XenoCollide); MODIFY narrowphase/gjk/functions.rs (remove/relocate simplified mpr_query), narrowphase/gjk/types.rs, narrowphase/dispatch/functions.rs (EPA-fallback wiring), narrowphase/gjk/mod.rs, narrowphase/mod.rs, lib.rs; NEW tests/mpr_epa_parity.rs
    - **Tests:** 1e5 random convex pairs depth within 1% of EPA; degenerate-simplex fallback triggers; smoke case parity; criterion bench >2x (env-gated)
    - **Risk:** Portal sub-region replacement signs — validate discovery phase on sphere/sphere first
    - **Done:** `mpr.rs` (631 lines): XenoCollide portal discovery to bracket the origin, then a self-contained, bounded minimum-penetration refinement seeded from a GJK-termination boundary tetrahedron (the interior-center seed was the bug — it corrupted the EPA polytope for curved/deep overlaps). `mpr_full`/`mpr_contact` shim via `functions.rs::mpr_query`; EPA-fallback wired in `dispatch::gjk_epa` (EPA `None` → `mpr_contact`). Parity verified: env-gated 1e5 sweep max relative depth error = 0.04% (gate 1%); 1k always-on subset and rotated-box cases pass. 10/10 in `tests/mpr_epa_parity.rs`, 11 in-crate mpr unit tests, full crate 2486 pass, clippy clean. NOTE: the env-gated criterion `>2x EPA speed` bench was not added (no criterion dev-dep / benches dir); the 1% depth-parity goal is met and is the load-bearing acceptance criterion.

- [x] Vertex-face/edge-edge continuous CCD with cubic root solver
  - **Goal:** zero missed collisions on cloth-vs-cloth benchmark at 5 m/s, dt=1/60.
  - **Design:** coplanarity cubic (Provot 1997 / Bridson 2002), conservative Bernstein root isolation + Newton polish; discrete VF/EE proximity already exists (`deformable_collision/functions.rs:145,225`) — add the time-of-impact path for oxiphysics-softbody.
  - **Files:** `deformable_collision/functions.rs`, new TOI path under `ccd/`
  - **Note:** primary consumer is oxiphysics-softbody's v0.3.0 IPC item (CCD-filtered line search)

## v0.3.0 — Hierarchy and scale

Theme: wide SIMD hierarchies, two-level instancing, GPU builds, compressed meshes, and large-world robustness.
Exit gate: 1M-triangle query speedups, <2 ms 10k-instance midphase, and origin-shift jitter parity at 100 km.

- [ ] QBVH: 4/8-wide SIMD BVH
  - **Goal:** ray-cast and AABB query ≥3x over binary DBVT on 1M-triangle mesh.
  - **Design:** SoA quantized child AABBs, `core::simd` traversal (pure Rust), Embree/Parry-style layout; refit path for dynamic bodies.
  - **Tests:** 1M-triangle ray-cast/AABB-query criterion benches vs `dbvt`; refit correctness under motion

- [ ] TLAS/BLAS two-level hierarchy
  - **Goal:** 10k instanced meshes broadphase+midphase <2 ms/frame; per-instance transform updates without BLAS rebuild.
  - **Design:** BLAS per shape (QBVH), TLAS over instances with refit-vs-rebuild SAH heuristic.
  - **Tests:** 10k-instance scene frame budget; SAH refit-vs-rebuild decision unit tests

- [ ] GPU LBVH radix build
  - **Goal:** 100k AABB BVH built on GPU <1 ms; bit-identical topology to CPU reference.
  - **Design:** Karras 2012 morton-code LBVH in WGSL; reuse oxiphysics-gpu radix sort primitives (`gpu_sort.rs`, `parallel_sort.rs`); cudarc variant feature-gated non-default (see deferred track).
  - **Tests:** GPU-vs-CPU topology bit-equality on randomized AABB sets

- [ ] Quantized/compressed trimesh midphase
  - **Goal:** ≤6 bytes/triangle BVH memory (vs ~32), <10% query slowdown.
  - **Design:** 16-bit quantized AABBs in parent frame, triangle-strip leaf packing; optional oxiarc-zstd offline asset compression via oxiphysics-io.
  - **Note:** compression goes through oxiarc-* only, per COOLJAPAN compression policy

- [ ] Large-world origin shifting
  - **Goal:** simulation 100 km from origin shows identical contact jitter to origin-local (f64 broadphase, f32-safe export).
  - **Design:** per-island local frames + explicit `shift_origin(delta)` API across broadphase/DBVT/manifold caches.
  - **Tests:** jitter-parity scenario at 0 km vs 100 km; cache-consistency checks across `shift_origin`

- [ ] Trimesh contact clustering
  - **Goal:** body sliding over 10k-triangle terrain produces ≤8 stable clustered contacts (no per-triangle popping).
  - **Design:** normal-bucketed clustering + deepest-point retention on top of existing per-manifold `reduce_to_4`.
  - **Tests:** 10k-triangle terrain slide with contact-count and popping metrics

## v1.0 — Production & validation

Theme: an unsafe-free, fuzz-hardened narrowphase with analytic accuracy gates and a documented no-tunneling contract.

- [ ] Unsafe-free + fuzzed narrowphase
  - **Goal:** `#![deny(unsafe_code)]` or documented justification per block; proptest 1e6 random shape pairs (degenerate: zero-size, near-touching, deep overlap) — zero panics, GJK/EPA distance error <1e-6.
  - **Design:** proptest strategies per shape; Miri job in CI (local script).
  - **Note:** builds directly on the v0.2.0 PRIORITY dispatch-soundness fix

- [ ] Accuracy conformance vs analytic baselines
  - **Goal:** distance/TOI for sphere-sphere, box-box, capsule-tri match closed-form to 1e-9.
  - **Design:** golden-value test suite; criterion gates (SAP 100k AABBs < 5 ms).

- [ ] CCD no-tunneling guarantee
  - **Goal:** stress suite (bullets 100 m/s vs 1 cm walls, 10k trials) zero tunnels.
  - **Design:** conservative advancement + speculative margin contract documented as API guarantee.
  - **Note:** pairs with oxiphysics-rigid's v0.2.0 solver-integrated speculative contacts

## Deferred / research track
- [~] cudarc (CUDA) variant of the GPU LBVH radix build (feature-gated, non-default per Pure Rust policy) — Ready when: local NVIDIA RTX hardware is available to validate cudarc kernels; the wgpu/WGSL path in v0.3.0 ships independently of this.

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 225 (critical 1, high 85, medium 107, low 32); 20 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (1)
- [x] **F0198** (critical, bug, medium/small) `narrowphase/specialized/functions.rs:177` — specialized::box_box_sat: A→B normal, centre-based points, unnormalized edge depth

### High (85)
- [x] **F0059** (high, fabrication, medium/medium) `broadphase/bvhbroadphase_traits.rs:25` — BvhBroadphase::find_pairs empty unless prebuilt; O(N^2), tree unused
- [x] **F0066** (high, bug, medium/small) `broadphase/types.rs:660` — BvhBroadphase find_pairs is O(n²), and empty when the tree is unbuilt
- [ ] **F0069** (high, bug, medium/small) `ccd.rs:133` — time_of_impact angular bound is a constant; slow approach exhausts iterations
- [ ] **F0070** (high, bug, medium/small) `ccd.rs:430` — gjk_distance_approx returns an upper bound; conservative advancement overshoots
- [x] **F0073** (high, fabrication, hard/large) `compound_shapes.rs:540` — ConvexDecomposition concavity is always zero, so it never splits
- [x] **F0076** (high, fabrication, easy/small) `compound_shapes.rs:1332` — ImportedCollisionMesh convex hull keeps only 6 extremal points
- [x] **F0078** (high, bug, medium/medium) `contact_generation.rs:518` — BoxBoxManifold one-sided SAT, deepest-axis choice, summed-extent depth
- [ ] **F0079** (high, bug, easy/small) `contact_generation.rs:592` — CapsuleCapsuleManifold does not re-solve s after clamping t
- [x] **F0082** (high, bug, easy/small) `contact_generation.rs:944` — PersistentManifold::update discards ages; stale points never expire
- [x] **F0092** (high, perf, medium/small) `dbvt/types.rs:419` — DynamicBvh::query_overlapping_pairs brute-forces all leaf pairs
- [ ] **F0097** (high, bug, medium/medium) `deformable_collision/functions.rs:284` — B-vertex and edge-edge contacts stored as mesh-A vertex indices
- [ ] **F0098** (high, bug, medium/medium) `deformable_collision/functions.rs:815` — ccd_deformable_deformable skips B→A, edge-edge and swept broadphase
- [ ] **F0099** (high, bug, medium/small) `deformable_collision/functions.rs:898` — generate_constraints targets inconsistent; overshoot for sphere and plane contacts
- [ ] **F0104** (high, bug, easy/small) `fluid_collision.rs:736` — SPH-rigid coupling: wrong cutoff, attractive sign, acceleration units
- [x] **F0106** (high, bug, medium/small) `gjk_enhanced/functions.rs:199` — gjk_enhanced exits on separating-axis test, so distances are wrong
- [x] **F0107** (high, bug, medium/medium) `gjk_enhanced/functions.rs:457` — epa_enhanced removes only one face per expansion (non-convex polytope)
- [x] **F0108** (high, bug, medium/small) `gjk_enhanced/functions.rs:611` — gjk_time_of_impact checks only t=1 and misses tunnelling
- [x] **F0109** (high, bug, medium/small) `gjk_enhanced/types.rs:235` — Warm-start seeds stale simplex; cached contact yields intersecting=false
- [x] **F0113** (high, fabrication, medium/small) `gjk_epa.rs:289` — gjk_distance discards min_dist and returns wrong distance
- [x] **F0114** (high, bug, medium/medium) `gjk_epa.rs:317` — gjk_distance returns support-vertex gap, not GJK distance
- [x] **F0115** (high, fabrication, medium/medium) `gjk_epa.rs:409` — EPA contact points query support at Minkowski vertices used as directions
- [x] **F0116** (high, fabrication, medium/small) `gjk_epa.rs:694` — gjk_distance_with_witnesses walks support vertices instead of running GJK
- [x] **F0119** (high, bug, medium/small) `gjk_extended.rs:461` — ExpandingPolytope::expand never finds horizon edges; polytope destroyed
- [x] **F0120** (high, fabrication, medium/small) `gjk_extended.rs:565` — EpaExpander fabricates normal +Z, depth 0 and zero contact points
- [x] **F0121** (high, bug, medium/medium) `gjk_extended.rs:637` — GjkSolver::query reports distance 0 for any non-enclosing tetrahedron
- [x] **F0122** (high, fabrication, medium/small) `gjk_extended.rs:716` — Gjk2D::distance never updates direction; returns 0 for separated polygons
- [x] **F0124** (high, fabrication, hard/medium) `gjk_extended.rs:783` — GjkRaycast::cast hit test uses an unrelated +Z support point
- [x] **F0125** (high, fabrication, medium/small) `gjk_extended.rs:861` — GjkTimeOfImpact::compute ignores velocity and never returns a TOI
- [ ] **F0128** (high, fabrication, medium/medium) `joint_collision.rs:874` — JointLimitSolver::solve fakes PGS: no velocity update, impulse scales with iterations
- [ ] **F0135** (high, bug, medium/small) `mesh_collision.rs:365` — triangle_triangle_intersect reports every coplanar pair as intersecting
- [x] **F0136** (high, panic, easy/trivial) `mesh_collision.rs:446` — TriMeshBvh::build overflows the stack on an empty triangle list
- [x] **F0137** (high, perf, medium/small) `mesh_collision.rs:800` — cloth_self_collision is O(V²·T): triangle scan per vertex pair
- [ ] **F0143** (high, bug, medium/small) `mesh_collision.rs:1064` — swept_triangle_vs_mesh tests only moving vertices against static faces
- [x] **F0153** (high, fabrication, medium/medium) `narrowphase/dispatch/functions.rs:1348` — dispatch_concave_vs_convex ignores convex shape, fabricates depth and pair ids
- [ ] **F0158** (high, perf, medium/medium) `narrowphase/dispatch/types.rs:408` — generate_contacts formats Debug per pair, bypasses registered manifold table
- [ ] **F0160** (high, bug, medium/medium) `narrowphase/epa/functions/functions_impl.rs:107` — EPA/MPR contact normal is A→B, contradicting docs and solver
- [x] **F0165** (high, fabrication, medium/medium) `narrowphase/functions.rs:260` — gjk_fallback is 3-axis world SAT with unsigned normal
- [x] **F0166** (high, fabrication, medium/medium) `narrowphase/functions.rs:285` — gjk_fallback is a 3-axis world SAT with unsigned normal
- [x] **F0167** (high, fabrication, medium/medium) `narrowphase/functions.rs:393` — convex_vs_mesh replaces each triangle with its bounding sphere
- [x] **F0174** (high, bug, easy/small) `narrowphase/gjk/functions.rs:549` — gjk_ray_cast reports hits for rays that laterally miss
- [x] **F0178** (high, fabrication, medium/small) `narrowphase/gjk/types.rs:527` — Gjk::query returns body origins as witnesses and unconverged distance
- [x] **F0180** (high, bug, easy/small) `narrowphase/gjk_cache/functions.rs:596` — gjk_cache::gjk_distance termination inverted; returns first-support distance
- [x] **F0181** (high, bug, medium/small) `narrowphase/gjk_cache/types.rs:932` — closest_points_cached reports overlapping shapes as separated; NaN otherwise
- [x] **F0183** (high, bug, medium/small) `narrowphase/manifold.rs:272` — build_box_manifold: tangent equals normal axis; every depth equals half extent
- [x] **F0190** (high, fabrication, medium/small) `narrowphase/simple.rs:306` — simple::box_box_sat: centre-based points, unnormalized edge depths, A→B normal
- [x] **F0195** (high, bug, easy/small) `narrowphase/specialized/functions.rs:14` — sphere_sphere/capsule_capsule normals A→B; solver expects B→A
- [x] **F0196** (high, bug, easy/small) `narrowphase/specialized/functions.rs:70` — sphere_box: centre inside box gives depth = radius; witnesses coincide
- [x] **F0197** (high, fabrication, medium/small) `narrowphase/specialized/functions.rs:159` — box_box_sat: unnormalized edge depth, centre-based contact points
- [x] **F0205** (high, bug, easy/trivial) `narrowphase/specialized/functions.rs:736` — torus_plane tests the ring point farthest from the plane
- [x] **F0207** (high, fabrication, medium/medium) `narrowphase/specialized/functions.rs:1096` — convex_hull_vs_convex_hull fakes EPA with centre normal, ignores rotation
- [x] **F0208** (high, bug, easy/small) `narrowphase/types.rs:376` — ShapeKind Plane support constant; box/capsule/convex vs plane never collide
- [x] **F0216** (high, bug, medium/small) `parallel_collision/types.rs:301` — AabbTreeInner::insert orphans every leaf after the second insert
- [x] **F0225** (high, fabrication, easy/small) `proximity.rs:732` — MeshMeshDist uses a centroid heuristic, not true triangle distance
- [x] **F0228** (high, bug, easy/trivial) `proximity_query.rs:901` — point_in_triangle_barycentric divides by |n| instead of |n|²
- [x] **F0230** (high, bug, easy/trivial) `query/functions.rs:656` — compute_time_of_impact_spheres returns the exit time when spheres already overlap
- [x] **F0232** (high, fabrication, medium/medium) `query/functions.rs:756` — project_point_convex fans vertices around the centroid; wrong for 3-D hulls
- [x] **F0239** (high, fabrication, hard/large) `recast/contour.rs:97` — Recast contours are convex hulls of region cells
- [x] **F0240** (high, stub, hard/large) `recast/mod.rs:22` — RecastConfig: four dead fields, no validation, overflow on cell_size
- [x] **F0245** (high, perf, hard/medium) `recast/region.rs:58` — Region build is O(N·W·H), has no watershed, and overflows u16 ids
- [x] **F0246** (high, bug, medium/small) `recast/walkable.rs:65` — Ledge filter compares floor top with the neighbour column's lowest span bottom
- [x] **F0248** (high, bug, easy/small) `sap/types.rs:604` — SweepAndPrune::reorder_axes breaks the overlap filter; full sort on every add/update
- [x] **F0253** (high, bug, easy/small) `sap/types.rs:886` — IncrementalSap::reorder_axes overwrites the X endpoint list with another axis
- [x] **F0254** (high, bug, easy/small) `sat_collision/functions.rs:193` — closest_points_on_segments omits Ericson re-clamp; shape_cast seg_seg_dist too
- [x] **F0255** (high, bug, medium/medium) `sat_collision/types.rs:230` — SAT contact points use the wrong support side and fake edges
- [x] **F0257** (high, bug, medium/small) `sat_collision/types.rs:985` — ObbCapsuleCollision::test uses the capsule point nearest the OBB centre
- [ ] **F0260** (high, bug, medium/small) `shape_cast.rs:479` — Conservative advancement with |v| step misses oblique contacts
- [ ] **F0261** (high, fabrication, medium/small) `shape_cast.rs:532` — CapsuleCast::vs_box tests endpoints only and hardcodes normal +Y
- [ ] **F0262** (high, bug, medium/small) `shape_cast.rs:631` — BoxCast::vs_triangle measures only triangle vertices; tunnels through faces
- [ ] **F0263** (high, fabrication, medium/medium) `shape_cast.rs:647` — BoxCast::vs_triangle measures triangle vertices only; tunnels through faces
- [ ] **F0264** (high, fabrication, medium/small) `shape_cast.rs:684` — ConvexCast 'GJK' is the closest vertex pair; overshoots and panics on empty input
- [ ] **F0269** (high, fabrication, medium/medium) `soft_body_collision.rs:397` — edge_edge_ccd is 64 fixed samples; batch skips cross-mesh pairs
- [ ] **F0270** (high, bug, easy/small) `soft_body_collision.rs:710` — SoftBodyBvh splits unpartitioned leaves, giving a degenerate chain and stack overflow
- [ ] **F0273** (high, fabrication, medium/medium) `softbody_collision.rs:807` — CcdDeformable::edge_edge_ccd reduces edges to their midpoints
- [ ] **F0274** (high, fabrication, easy/small) `softbody_collision.rs:917` — PrimitiveVsDeformable::query_box fabricates depth 1e-6 and +Y normal
- [ ] **F0275** (high, bug, easy/small) `softbody_collision.rs:1056` — ContactResponse uses friction as thickness; impulses applied to separating contacts
- [ ] **F0276** (high, bug, medium/medium) `spatial_queries.rs:255` — AabbTree::insert never refits ancestors, so queries miss leaves
- [ ] **F0277** (high, bug, easy/small) `spatial_queries.rs:351` — AabbTree::update tests overlap instead of fat-AABB containment
- [ ] **F0284** (high, bug, easy/trivial) `sweep.rs:981` — conservative_advancement_rotating step has wrong units and overshoots the TOI
- [ ] **F0285** (high, bug, medium/small) `sweep.rs:1075` — swept_sphere_triangle_mesh detects only face hits and misses edges and vertices
- [ ] **F0286** (high, fabrication, easy/small) `sweep.rs:1153` — swept_sphere_capsule replaces the capsule with one static sphere
- [ ] **F0291** (high, fabrication, medium/medium) `terrain_collision/types.rs:549` — VoxelTerrain::extract_isosurface emits fixed corner triangles, not marching cubes
- [ ] **F0295** (high, bug, medium/small) `terrain_collision/types.rs:1326` — TerrainRaycast loops in place on horizontal rays and approximates cells as flat
- [ ] **F0301** (high, bug, medium/medium) `voxel_collision.rs:383` — MeshVoxelizer sign from nearest face plane; O(voxels·T); index panic
- [ ] **F0302** (high, bug, easy/small) `voxel_collision.rs:509` — FastSweeping Eikonal update formulas are wrong; distances inflated
- [ ] **F0303** (high, bug, medium/small) `voxel_collision.rs:696` — VoxelRayQuery::cast aborts at the grid boundary and never hits from outside

### Medium (107)
- [ ] **F0062** (medium, perf, easy/small) `broadphase/sweepandprune_traits.rs:26` — Pipeline SAP fixed to X axis; vertical stacks degrade to O(N²)
- [ ] **F0063** (medium, panic, easy/medium) `broadphase/sweepandprune_traits.rs:33` — NaN-unsafe partial_cmp sort comparators can panic broadphase
- [ ] **F0064** (medium, perf, easy/small) `broadphase/types.rs:130` — BroadphaseSceneGraph rebuilds a full BVH on every query
- [ ] **F0065** (medium, fabrication, medium/medium) `broadphase/types.rs:306` — DynamicAabbTree claims incremental updates but does full rebuilds
- [ ] **F0068** (medium, perf, easy/small) `broadphase/types.rs:951` — Warm-start merge is O(n·m) and duplicates a full SAP pass
- [ ] **F0071** (medium, api, medium/small) `ccd.rs:1077` — CcdPipelineConfig iteration/tolerance dead; detect is O(N^2)
- [ ] **F0072** (medium, perf, easy/small) `ccd.rs:1121` — CcdPipeline::detect is O(n²) and recomputes swept AABBs per pair
- [ ] **F0074** (medium, bug, medium/small) `compound_shapes.rs:752` — point_inside ray parity double-counts shared edges; O(T) scans
- [ ] **F0075** (medium, fabrication, easy/small) `compound_shapes.rs:974` — MinkowskiSum::origin_in_md uses 6 axes, giving false-positive overlaps
- [ ] **F0077** (medium, fabrication, medium/medium) `compound_shapes.rs:1444` — shape_volume containment test is a support-function heuristic, not containment
- [ ] **F0081** (medium, bug, easy/small) `contact_generation.rs:840` — ManifoldReduction::reduce step 4 degenerate on coplanar points; picks duplicates
- [ ] **F0083** (medium, bug, easy/small) `contact_generation.rs:1263` — ContactEventSystem emits events in random HashSet order
- [ ] **F0084** (medium, bug, easy/small) `contact_graph.rs:228` — speculative_impulse ignores separation and bounces before contact
- [ ] **F0085** (medium, perf, easy/small) `contact_graph.rs:254` — ContactCache uses linear search and Vec::remove(0); mislabeled LRU
- [ ] **F0086** (medium, panic, easy/small) `contact_graph.rs:921` — ContactCacheInvalidator::update_and_check panics when body count grows
- [ ] **F0089** (medium, bug, easy/small) `contact_manifold.rs:380` — reduce_to_4 returns duplicate points for collinear or degenerate contacts
- [ ] **F0090** (medium, bug, easy/trivial) `contact_manifold.rs:500` — Sphere-AABB depth wrong when sphere centre is inside box
- [ ] **F0091** (medium, bug, medium/medium) `dbvt/types.rs:363` — DynamicBvh raw usize handles: stale remove silently deletes another proxy
- [ ] **F0093** (medium, perf, medium/small) `dbvt/types.rs:445` — DynamicBvh::query_overlapping_pairs is brute-force O(n²) over leaves
- [ ] **F0096** (medium, bug, easy/small) `deformable_collision/functions.rs:236` — Zero-distance vertex-face and edge-edge contacts are silently dropped
- [ ] **F0100** (medium, api, medium/medium) `error.rs:16` — Crate Error enum is never produced by any API
- [ ] **F0101** (medium, bug, easy/trivial) `fluid_collision.rs:230` — buoyancy_sphere centre of buoyancy wrong except at h=r
- [ ] **F0102** (medium, bug, easy/trivial) `fluid_collision.rs:342` — lift_force points down for positive C_l and spuriously along +z
- [ ] **F0105** (medium, fabrication, easy/trivial) `fluid_collision.rs:955` — underwater_contact_friction viscous coefficient has no effect on output
- [ ] **F0118** (medium, panic, easy/small) `gjk_extended.rs:381` — ExpandingPolytope::from_simplex asserts on simplices GJK legitimately returns
- [ ] **F0123** (medium, bug, medium/small) `gjk_extended.rs:744` — Gjk2D::distance never updates dir and returns nearest vertex
- [ ] **F0126** (medium, bug, easy/small) `joint_collision.rs:110` — Joint penetration doc says signed but returns magnitude; impulse units wrong
- [ ] **F0127** (medium, bug, easy/trivial) `joint_collision.rs:647` — max_friction_impulse returns a force, ignores dt
- [ ] **F0133** (medium, bug, easy/small) `manifold_cache.rs:631` — ManifoldCompressor::compress ignores max_count of 2, 3 and above 4
- [ ] **F0134** (medium, perf, easy/small) `manifold_cache.rs:880` — build_contact_islands is O(M·B) and non-deterministic across runs
- [ ] **F0138** (medium, perf, easy/small) `mesh_collision.rs:810` — cloth_self_collision is O(V²·T) with a triangle scan per pair
- [ ] **F0139** (medium, fabrication, easy/small) `mesh_collision.rs:856` — vertex_triangle_ccd claims a cubic solve but samples 64 steps
- [ ] **F0140** (medium, bug, medium/small) `mesh_collision.rs:893` — vertex_triangle_ccd samples 64 steps despite a 'cubic root' doc
- [ ] **F0141** (medium, bug, easy/trivial) `mesh_collision.rs:988` — closest_point_on_triangle clamps barycentrics instead of projecting
- [ ] **F0142** (medium, perf, medium/medium) `mesh_collision.rs:1000` — build_proximity_field brute-forces every grid cell against every triangle
- [ ] **F0144** (medium, fabrication, medium/small) `mesh_collision.rs:1097` — penetration_depth_estimate returns unsigned distance on a sample line
- [ ] **F0145** (medium, fabrication, medium/medium) `mesh_collision.rs:1183` — Mesh contact manifold: difference normal, centroid point, fake depth
- [ ] **F0148** (medium, perf, medium/medium) `narrowphase/dispatch/convex_manifold.rs:190` — extract_convex_faces is O(n^4) per hull per contact
- [ ] **F0149** (medium, perf, easy/small) `narrowphase/dispatch/functions.rs:325` — type_name formats the Debug output of both shapes on every generate_contacts call
- [ ] **F0150** (medium, bug, easy/small) `narrowphase/dispatch/functions.rs:389` — Compound 'AABB pruning' is hardcoded child-centre distance > 5
- [ ] **F0151** (medium, fabrication, easy/trivial) `narrowphase/dispatch/functions.rs:1236` — box_box_manifold_clip fakes a 4-point manifold by jittering box A
- [ ] **F0152** (medium, bug, easy/small) `narrowphase/dispatch/functions.rs:1293` — Speculative contact depth positive when separated; overlap gives margin+1
- [ ] **F0154** (medium, fabrication, medium/medium) `narrowphase/dispatch/functions.rs:1420` — Heightfield generic dispatch corrupts pair ids and over-reports per-vertex depth
- [ ] **F0155** (medium, api, medium/small) `narrowphase/dispatch/types.rs:42` — DispatchConfig fallback/iteration/tolerance fields are dead
- [ ] **F0162** (medium, fabrication, medium/medium) `narrowphase/epa/functions_2.rs:38` — Raw-EPA witness points are Minkowski-difference points, not shape points
- [ ] **F0164** (medium, fabrication, medium/small) `narrowphase/epa/types.rs:389` — EPA fallback_contact invents a contact (+Y normal) and blocks the MPR fallback
- [ ] **F0168** (medium, fabrication, medium/small) `narrowphase/functions.rs:591` — point_query on Convex returns nearest vertex, never inside
- [ ] **F0169** (medium, fabrication, medium/small) `narrowphase/functions.rs:622` — ray_cast on Convex uses a bogus bounding sphere
- [ ] **F0170** (medium, bug, easy/trivial) `narrowphase/functions.rs:708` — ray_vs_aabb gives non-unit wrong normals for −y/−z rays
- [ ] **F0171** (medium, perf, medium/medium) `narrowphase/gjk/functions.rs:66` — GJK simplex allocates a new Vec at every reduction step
- [ ] **F0176** (medium, bug, medium/small) `narrowphase/gjk/functions.rs:748` — gjk_toi_angular_sweep never accumulates rotation
- [ ] **F0177** (medium, bug, medium/small) `narrowphase/gjk/functions.rs:947` — gjk_hull_distance returns a search-direction length as the distance
- [ ] **F0185** (medium, fabrication, medium/medium) `narrowphase/obb_sat/functions.rs:121` — compute_contact_points claims face clipping but returns 1-2 guessed points
- [ ] **F0186** (medium, bug, easy/small) `narrowphase/obb_sat/types.rs:257` — OBB SAT normals A→B though documented B→A
- [ ] **F0189** (medium, bug, medium/small) `narrowphase/primitives.rs:720` — sphere_cylinder computes capsule distance, ignoring flat caps
- [ ] **F0191** (medium, bug, medium/small) `narrowphase/simple.rs:881` — obb_capsule closest point checks only segment endpoints and midpoint
- [ ] **F0192** (medium, bug, medium/small) `narrowphase/simple.rs:925` — aabb_capsule tests only the segment point nearest the box centre
- [ ] **F0194** (medium, fabrication, easy/small) `narrowphase/simple.rs:1049` — heightfield_box samples 5 points despite claiming to iterate grid cells
- [ ] **F0201** (medium, bug, medium/small) `narrowphase/specialized/functions.rs:305` — cone_sphere takes the closest point at the sphere's height, not on the slant
- [ ] **F0203** (medium, bug, medium/small) `narrowphase/specialized/functions.rs:631` — capsule_triangle_mesh samples 6 points per triangle and misses contacts
- [ ] **F0204** (medium, bug, easy/small) `narrowphase/specialized/functions.rs:696` — Heightfield sphere and box contacts hardcode a +Y normal and ignore slope
- [ ] **F0209** (medium, fabrication, easy/small) `parallel_collision/functions.rs:158` — par_sphere_box and par_box_box invent normal, depth and contact point
- [ ] **F0211** (medium, bug, easy/trivial) `parallel_collision/functions.rs:235` — ccd_sphere_sphere reports exit time for initially overlapping pairs
- [ ] **F0212** (medium, bug, easy/small) `parallel_collision/functions.rs:284` — Contact reduction and island order depend on HashMap randomness
- [ ] **F0213** (medium, fabrication, medium/medium) `parallel_collision/functions.rs:333` — run_parallel_collision is serial O(N²) and discards its SoA buffer
- [ ] **F0214** (medium, fabrication, easy/small) `parallel_collision/functions.rs:339` — run_parallel_collision discards SoA buffer, ignores margin, serial O(N^2)
- [ ] **F0215** (medium, bug, medium/small) `parallel_collision/types.rs:49` — ParallelUnionFind::union TOCTOU race loses merges under concurrency
- [ ] **F0220** (medium, missing-feature, medium/medium) `parallel_collision/types.rs:726` — parallel_collision module performs no parallel work
- [ ] **F0222** (medium, bug, easy/trivial) `proximity.rs:434` — CapsuleCapsuleDist witnesses pushed to far side of capsules
- [ ] **F0223** (medium, perf, easy/small) `proximity.rs:462` — KnnQuery claims kd-tree but full-sorts all points per query
- [ ] **F0224** (medium, perf, medium/small) `proximity.rs:674` — PointMeshDist 'BVH-accelerated' is brute force; assert! panics
- [ ] **F0226** (medium, fabrication, easy/trivial) `proximity.rs:772` — SweptVolumeDist samples 9 times; treats capsules as points
- [ ] **F0227** (medium, bug, easy/small) `proximity_query.rs:768` — triangle_triangle_distance returns >0 for intersecting triangles
- [ ] **F0229** (medium, bug, medium/small) `query/functions.rs:572` — sphere_cast_scene/sphere_sweep_vs_box use square Minkowski box
- [ ] **F0231** (medium, bug, medium/medium) `query/functions.rs:729` — project_point_convex fans unordered 3-D vertices; wrong result
- [ ] **F0234** (medium, bug, easy/small) `ray_casting/functions.rs:354` — ray_heightfield marches past t_max while outside grid
- [ ] **F0237** (medium, perf, easy/small) `ray_casting/types.rs:315` — Bvh::intersect lacks closest-hit pruning and near-first order
- [ ] **F0238** (medium, panic, easy/small) `ray_casting/types.rs:647` — Heightfield::height_at underflows/out-of-bounds; constructor unvalidated
- [ ] **F0241** (medium, bug, easy/small) `recast/polymesh.rs:51` — PolyMesh indices silently truncated by `as u16`
- [ ] **F0242** (medium, bug, medium/medium) `recast/polymesh.rs:75` — 'ear_clip_triangulate' is a fan; u16 index wrap; polygons never merged
- [ ] **F0243** (medium, bug, easy/small) `recast/rasterize.rs:139` — Span merge ORs walkable flag; abs(ny) marks undersides walkable
- [ ] **F0247** (medium, bug, easy/small) `sap/types.rs:106` — GridBroadphase insert hangs on huge AABB or zero cell_size
- [ ] **F0249** (medium, perf, medium/medium) `sap/types.rs:692` — IncrementalSap update scans all endpoints; query builds O(n²) sets
- [ ] **F0250** (medium, perf, medium/medium) `sap/types.rs:756` — IncrementalSap::query_pairs fully sorts 3 axes and builds O(n²) per-axis sets
- [ ] **F0251** (medium, bug, medium/medium) `sap/types.rs:774` — IncrementalSap drops touching pairs by insertion order; reorder_axes corrupts X
- [ ] **F0252** (medium, perf, easy/small) `sap/types.rs:817` — insert_aabb/update_aabb run a full query_pairs per body
- [ ] **F0256** (medium, perf, medium/medium) `sat_collision/types.rs:866` — ObbTree nodes are AABBs yet each visit runs full 15-axis SAT
- [ ] **F0258** (medium, bug, easy/trivial) `sat_collision/types.rs:1037` — ConvexPolyhedraSat returns unoriented contact normal
- [ ] **F0265** (medium, panic, easy/trivial) `shape_cast.rs:823` — NonLinearCast::cast loops forever when dt <= 0
- [ ] **F0266** (medium, perf, medium/medium) `shape_cast.rs:855` — MeshBvhCast/PointMeshDist claim BVH but are brute force
- [ ] **F0267** (medium, perf, medium/small) `shape_cast.rs:864` — MeshBvhCast has no BVH and asserts on malformed indices
- [ ] **F0268** (medium, perf, easy/small) `soft_body_collision.rs:160` — detect_vertex_triangle_contacts brute-force O(V·T) despite SoftBodyBvh
- [ ] **F0271** (medium, perf, medium/small) `softbody_collision.rs:339` — DeformableVsBvh is a brute-force O(V·T) scan with a wrong-side normal
- [ ] **F0272** (medium, bug, easy/small) `softbody_collision.rs:461` — Self-collision hash ignores threshold > cell_size, missing pairs
- [ ] **F0278** (medium, missing-feature, easy/trivial) `spatial_queries.rs:744` — ContactPairFilter never reads groups despite documented group check
- [ ] **F0279** (medium, bug, easy/trivial) `sweep.rs:169` — swept_sphere_sphere impact-branch normal violates documented B→A convention
- [ ] **F0281** (medium, bug, easy/trivial) `sweep.rs:548` — swept_capsule_plane uses radius+half_height for every plane normal
- [ ] **F0283** (medium, fabrication, medium/small) `sweep.rs:793` — swept_capsule_aabb sweeps midpoint sphere of radius r+half_len
- [ ] **F0287** (medium, bug, easy/small) `sweep.rs:1172` — swept_sphere_capsule replaces capsule by one sphere at midpoint projection
- [ ] **F0288** (medium, bug, easy/small) `terrain_collision/types.rs:155` — TerrainShadowMap marches 20 fixed samples over half the terrain
- [ ] **F0290** (medium, bug, easy/trivial) `terrain_collision/types.rs:447` — Sphere buoyancy uses linear submerged fraction, not spherical cap
- [ ] **F0292** (medium, panic, easy/small) `terrain_collision/types.rs:758` — HeightfieldShape::new panics; height_at underflows for nx<2
- [ ] **F0293** (medium, panic, easy/small) `terrain_collision/types.rs:776` — HeightfieldShape::height_at underflows for nx<2; extrapolates past edge
- [ ] **F0294** (medium, bug, medium/small) `terrain_collision/types.rs:1131` — Sphere-terrain depth measured vertically, misses slope contacts
- [ ] **F0296** (medium, bug, medium/small) `terrain_collision/types.rs:1519` — TerrainCcd is discrete sampling; steps=0 yields NaN and silent miss
- [ ] **F0298** (medium, duplicate, hard/large) `types.rs:224` — Duplicate public ContactManifold, ContactPoint, SweepAndPrune, ConvexShape types
- [ ] **F0299** (medium, bug, easy/small) `types.rs:369` — warm_start_from matches Unknown features, duplicating impulses

### Low (32) — tracked in the external ledger
F0060, F0067, F0087, F0094, F0095, F0103, F0110, F0111, F0112, F0117, F0129, F0130, F0146, F0147, F0156, F0157, F0161, F0173, F0179, F0184, F0193, F0200, F0202, F0210, F0233, F0235, F0236, F0244, F0259, F0280, F0282, F0300

<!-- AUDIT-2026-10-12 END -->
