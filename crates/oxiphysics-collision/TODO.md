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

- [ ] MPR (Minkowski Portal Refinement)
  - **Goal:** penetration depth within 1% of EPA on 1e5 random convex pairs at >2x EPA speed; used as EPA fallback on simplex degeneracy.
  - **Design:** Snethen XenoCollide portal refinement.
  - **Tests:** 1e5 random convex-pair depth comparison vs EPA; degenerate-simplex fallback trigger cases

- [ ] Vertex-face/edge-edge continuous CCD with cubic root solver
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
