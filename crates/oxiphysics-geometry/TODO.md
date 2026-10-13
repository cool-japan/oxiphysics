# oxiphysics-geometry TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 85,006 SLoC | 3,120 tests

Geometry kernel of the OxiPhysics workspace: primitive shapes, triangle meshes
and the full mesh-processing pipeline, computational geometry, parametric and
implicit surfaces, procedural generation. The forward roadmap below is
code-verified (2026-06-11) and focuses on exact-arithmetic robustness, 3D
meshing, advanced surface operations, and analysis-suitable (IGA-ready)
geometry — all pure Rust.

**Cross-crate contract** — dependencies and consumers of this roadmap:
- Exact-arithmetic robust booleans (v0.2.0) consume `oxiphysics-core`'s Shewchuk adaptive predicates (core v0.2.0); the same core predicates also serve `oxiphysics-collision` narrow-phase robustness. Sequence the boolean rebuild after core lands them.
- Convex decomposition (VHACD shipped; CoACD-style v2 at v1.0) supplies collision proxies consumed by `oxiphysics-collision`.
- T-spline / subdivision-limit evaluation (v1.0) feeds `oxiphysics-fem`'s isogeometric-analysis upgrade (fem v0.3.0).

## Completed (v0.1.0 – v0.1.3)

### Milestone 1: Primitive Shapes
- [x] Define shape trait hierarchy (`shape`)
- [x] Implement `Sphere`, `BoxShape`, `Capsule`, `Cylinder`, `Cone`, `Torus`
- [x] `Compound` shape aggregation
- [x] Unit tests (3,120 tests passing)

### Milestone 2: Meshes & Computational Geometry
- [x] Triangle mesh (`triangle_mesh`)
- [x] Convex hull (`convex_hull`, `quickhull`)
- [x] Height field (`heightfield`)
- [x] Voronoi diagrams (`voronoi`)
- [x] Medial axis transform (`medial_axis`)
- [x] CSG operations (`csg`, `mesh_boolean`)
- [x] Point cloud processing (`point_cloud`)
- [x] Spatial hashing (`spatial_hash`)

### Milestone 3: Mesh Processing Pipeline
- [x] Mesh operations (`mesh_ops`)
- [x] Mesh repair (`mesh_repair`)
- [x] Mesh quality analysis (`mesh_quality`)
- [x] Mesh simplification / decimation (`mesh_simplification`, `decimation`)
- [x] Remeshing (`remesh`)
- [x] Subdivision surfaces (`subdivision`)

### Milestone 4: Parametric & Implicit Surfaces
- [x] NURBS geometry (`nurbs_geometry`)
- [x] B-splines (`bspline`)
- [x] Parametric surfaces (`parametric`)
- [x] Implicit surfaces (`implicit_surfaces`)
- [x] Level sets (`level_set`)
- [x] Signed distance fields (`signed_distance_field`)
- [x] Offset surfaces (`offset_surface`)
- [x] Swept volumes (`swept`)

### Milestone 5: Procedural & Topology
- [x] Fractal geometry (`fractal_geometry`)
- [x] Geodesic domes (`geodesic`)
- [x] Origami folding geometry (`origami`)
- [x] Terrain processing (`terrain_processing`)
- [x] Procedural geometry (`procedural_geometry`)
- [x] Robot geometry primitives (`robot_geometry`)
- [x] Topology analysis (`topology`, `topology_geometry`)

### Milestone 6: Polish
- [x] Documentation (rustdoc, 0 warnings)
- [x] Examples
- [x] Benchmarks (criterion)
- [x] 0 stubs

### Verified baseline for the forward roadmap (code-checked 2026-06-11)
- [x] Heat-method geodesics (`geodesic_geometry` `HeatMethodParams`) + Dijkstra graph geodesics (`geodesic/`) — heat method is NOT re-opened; only exact MMP is open below
- [x] FMM signed-distance fields + octree acceleration (`signed_distance_field/fmm.rs`, `octree.rs`) — FMM is NOT re-opened; JFA + Poisson reconstruction are the open parts
- [x] VHACD approximate convex decomposition (`vhacd.rs`) — comparison baseline for the CoACD-style v2 gate
- [x] 2D Delaunay + Voronoi-from-Delaunay (`computational_geometry`: `delaunay_2d`, `voronoi_from_delaunay`) — 3D/CDT/α-shapes are the open parts
- [x] QEM simplification/decimation (`mesh_simplification`, `decimation`) — attribute-preserving v2 is the open part
- [x] Simplified LSCM + Tutte parameterization (`mesh_param.rs`) — ABF++ and a proper sparse LSCM are the open parts
- [x] Medial axis transform (`medial_axis.rs`) — distinct structure from the planned straight skeleton
- [x] Float-predicate mesh booleans/CSG (`mesh_boolean`, `csg`) — functional today; exact-arithmetic rebuild is the open part

## v0.2.0 — Robustness & meshing

- [x] Exact-arithmetic robust booleans — degenerate-band delivered 2026-06-12
  - **Files:** mesh_boolean.rs (degenerate-band sidedness routed through `oxiphysics_core::exact_predicates::orient3d`)
  - **Tests:** coincident-cube Union/Difference watertight + Euler V−E+F=2; 1000 near-degenerate booleans (cube∩cube rotated + cube∩octahedron) zero panic/NaN, finite tri counts; volume identity |A∪B|=|A|+|B|−|A∩B| on a clean overlap
  - **Risk (resolved):** the winding-number core did resist full exact threading, so a focused exact-classify helper handles the degenerate band while keeping the `SimpleMesh`/`mesh_boolean` API stable.
  - **Goal (met for coincident faces):** boolean of two coincident-face cubes yields a watertight manifold (Union & Difference, Euler=2); 1000 near-degenerate booleans, 0 failures.
  - **Design:** consumes core Shewchuk `orient3d` + Simulation-of-Simplicity (parity-of-lowest-index tie-break) in `mesh_boolean`. Shewchuk 1997; Edelsbrunner-Mücke SoS.
  - Delivered: exact ray-parity point-in-mesh (`exact_classify_point_vs_mesh`), exact ray–triangle crossing (`exact_ray_crosses_triangle` via 5 `orient3d` tests), SoS tie-break (`sos_sign`), coincident-face router (`TriClass`/`classify_triangle_exact`/`keep_triangle`); `mesh_boolean` now classifies via the exact, coincident-face-aware path. File 1952 lines (<2000); clippy clean (-D warnings); no `#[allow]`; no production unwrap.
  - [x] Replace float orientation/intersection sign tests with core `orient3d`-family predicates (centroid sidedness + coincident-face plane tests now exact)
  - [x] Simulation-of-simplicity (SoS) symbolic perturbation to eliminate degenerate branches (`sos_sign`: degenerate `orient3d` → parity of lowest global vertex index, order-independent → consistent virtual general position)
  - [x] Exact face classification (inside/outside/on-boundary) via `classify_triangle_exact` → `TriClass{Inside,Outside,CoplanarSame,CoplanarOpposite}`
  - [x] Watertight coincident-face result with manifoldness post-conditions (Union & Difference: `is_watertight` + Euler=2 asserted)
  - [x] Stress corpus: 1000 near-degenerate booleans (coincident faces, near-coincident rotated/translated), 0 failures; coincident-cube watertight gate
  - [ ] full-BSP exact booleans (intersection-curve re-triangulation; the degenerate-band classifier shipped now makes coincident-FACE booleans watertight via orient3d+SoS, but arbitrary mid-triangle intersections still need a BSP/CDT rebuild for guaranteed manifoldness and an exact triangle-triangle intersection segment). Boundary: whole-triangle centroid selection cannot re-triangulate where a face is pierced through its interior, so the volume identity holds only to a whole-triangle-selection tolerance there and Intersection of face-touching solids is the zero-volume degenerate case (correctly removed).

- [x] 3D Delaunay + constrained Delaunay + alpha shapes (delivered 2026-06-12; CDT facet recovery deferred — see sub-item)
  - **Goal:** tetrahedralize 10⁴ points with empty-sphere property holding for all tets; α-shape recovers known surface.
  - **Files:** `computational_geometry/delaunay_3d.rs` (new, ~1440 L < 2000); wired through `computational_geometry/mod.rs` + `lib.rs` (`pub use computational_geometry::delaunay_3d`).
  - **Design:** incremental Bowyer-Watson 3D + star-shaped-cavity repair + α-complex. Edelsbrunner-Mücke 1994; Bowyer/Watson 1981; predicates Shewchuk 1997; SoS Edelsbrunner-Mücke 1990.
  - Verified: `computational_geometry` has 2D `delaunay_2d` + `voronoi_from_delaunay` (see Completed); 3D + α now delivered, CDT still open.
  - Depends on: core exact `insphere` predicate for robust empty-sphere tests (consumed via `oxiphysics_core::exact_predicates::{insphere, orient3d}`).
  - Delivered: exact Bowyer-Watson driven by core `orient3d`/`insphere`; every tet kept positively oriented; all degeneracies (cospherical/coplanar/grid) resolved by parity-of-lowest-index SoS (matching `mesh_boolean::sos_sign`); walk-based point location with brute-force fallback; circumsphere flood + **star-shaped-cavity repair** (face-visibility test that absorbs non-visible neighbours, guaranteeing the refill is gap/overlap-free — this was the key to correct volumes under cospherical degeneracy); full face-adjacency (`neigh`) rebuild; super-tetra strip + compaction. `Tetrahedralization::{circumradius_sq, alpha_shape, verify_empty_sphere}`. Tests: single-tet, bipyramid 2-tet known case, cube (8 cospherical corners, volume==1), random 400/1000-pt empty-sphere audits (default-run) + env-gated 10⁴-pt audit, α-shape closed-surface (sphere cloud + convex hull, every edge twice), α=0 empty, duplicate tolerance. clippy `-D warnings` clean; no `#[allow]`; no production `unwrap()`.
  - [x] Incremental Bowyer-Watson 3D with walk-based point location
  - [x] Cavity retriangulation + degenerate-input handling via exact predicates
  - [ ] Constrained Delaunay: boundary facet recovery — follow-on; unconstrained Delaunay + α-shapes shipped now, but enforcing prescribed boundary facets (Steiner-point insertion / facet recovery) is a distinct effort filed for a later pass.
  - [x] α-complex filtration + α-shape extraction with known-surface recovery test
  - [x] Empty-sphere audit over every tet of a 10⁴-point tetrahedralization (`verify_empty_sphere`; `random_cloud_10k_empty_sphere_audit`, env-gated `OXIPHYSICS_DELAUNAY_10K`)

- [ ] ABF++ parameterization + LSCM v2
  - **Goal:** disk-topology mesh UV area distortion < 5% vs current simplified LSCM.
  - **Design:** angle-based flattening with proper sparse solve. Sheffer 2005.
  - Verified: `mesh_param.rs` has simplified LSCM + Tutte (see Completed); ABF is missing.
  - [ ] ABF++ angle-space optimization (Newton on the Lagrangian with sparse linear solves)
  - [ ] Angle-to-UV reconstruction (least-squares layout)
  - [ ] LSCM v2 with a proper sparse least-squares solve replacing the simplified path
  - [ ] Distortion metrics (area/angle) + disk-topology benchmark meshes; < 5% area-distortion gate

## v0.3.0 — Advanced operations

- [ ] Straight skeleton + 2D Minkowski sums
  - **Goal:** skeleton of a non-convex polygon topologically correct; Minkowski sum area matches analytic for convex pairs to 1e-9.
  - **Design:** wavefront-propagation straight skeleton + convolution Minkowski sum. Aichholzer 1996.
  - Verified: `medial_axis.rs` is a distinct structure (see Completed); straight skeleton + general Minkowski are missing.
  - [ ] Wavefront propagation with edge and split events (priority-queue kinetic simulation)
  - [ ] Skeleton arc/face extraction + non-convex topology tests
  - [ ] Convex-convex Minkowski sum (edge merge) + convolution method for general polygons
  - [ ] Analytic-area validation for convex pairs to 1e-9; polygon-offsetting application test

- [ ] Exact geodesics (MMP)
  - **Goal:** MMP distance matches analytic on unfolded cube to 1e-9.
  - **Design:** Mitchell-Mount-Papadimitriou exact polyhedral geodesics. Mitchell 1987.
  - Verified: `geodesic/` Dijkstra + `geodesic_geometry` heat method already ship (see Completed) — exact MMP is the only open part.
  - [ ] Edge-window data structure + window propagation queue
  - [ ] Window trimming/merging with exact unfolding arithmetic
  - [ ] Distance-field extraction + geodesic path backtracing
  - [ ] Unfolded-cube analytic validation to 1e-9; cross-check against the heat-method approximation

- [ ] Quad-dominant remeshing (cross-field)
  - **Goal:** quad-dominant remesh with >85% quads, singularities only at field defects.
  - **Design:** instant-meshes-style cross-field + position-field optimization. Jakob 2015.
  - Verified: `remesh.rs` is triangle remeshing; the quad path is missing.
  - [ ] 4-RoSy cross-field smoothing (hierarchical/extrinsic)
  - [ ] Position-field optimization aligned to the cross-field
  - [ ] Quad extraction + T-vertex cleanup
  - [ ] Metrics gate: >85% quad ratio, singularities only at field defects

- [~] (planned 2026-06-14) Jump-flooding SDF + screened-Poisson reconstruction
  - **Goal:** JFA SDF within 1 voxel of exact; Poisson recon reproduces sphere to <1% radius error.
  - **Design:** JFA voxel SDF + Kazhdan screened-Poisson surface reconstruction from oriented points. Rong-Tan 2006; Kazhdan-Hoppe 2013.
  - Verified: `signed_distance_field/fmm.rs` + `octree.rs` already ship (see Completed) — JFA + Poisson reconstruction are the open parts.
  - [ ] JFA seed injection + 1+JFA flood passes on a voxel grid
  - [ ] Sign resolution + accuracy audit (≤ 1 voxel vs exact SDF)
  - [ ] Screened-Poisson: octree-discretized Laplacian, normal-field divergence RHS, screening term
  - [ ] Iso-surface extraction + sphere-reconstruction gate (<1% radius error)

## v1.0 — Production & validation

### Collision-grade decomposition & simplification quality
- [ ] CoACD-style approximate convex decomposition v2
  - **Goal:** ≥20% fewer hulls than VHACD at equal concavity tolerance on a benchmark mesh set.
  - **Design:** MCTS concavity-aware cutting. Wei 2022.
  - Verified: `vhacd.rs` ships (see Completed) and is the comparison baseline.
  - Consumers: collision proxies for `oxiphysics-collision`.
  - [ ] Concavity metric (Hausdorff-based) + cutting-plane candidate generation
  - [ ] MCTS search over cut sequences with concavity-aware rollouts
  - [ ] Benchmark harness: hull count at equal concavity tolerance vs existing VHACD (≥20% reduction gate)

- [ ] Quadric simplification v2 with attribute preservation
  - **Goal:** 90% decimation keeps Hausdorff < 0.5% bbox and UV-seam error < 1px-equiv.
  - **Design:** QEM extended to preserve UV/normals/seams + mesh-quality validation harness. Garland-Heckbert 1998 (attributes).
  - Verified: `mesh_simplification`/`decimation` QEM ship (see Completed); the attribute-aware variant is missing.
  - [ ] Extended quadrics over position + UV + normal space
  - [ ] Seam/border constraint handling during edge collapse
  - [ ] Validation harness: Hausdorff (< 0.5% bbox) + UV-seam error (< 1px-equiv) at 90% decimation

### IGA-ready surfaces
- [ ] T-spline / subdivision-limit evaluation toward IGA
  - **Goal:** T-spline surface watertight across T-junctions; basis partition-of-unity to 1e-12.
  - **Design:** analysis-suitable T-splines + watertight subdivision limit evaluation feeding FEM IGA. Sederberg 2003; Scott 2012.
  - Verified: B-spline/NURBS/subdivision ship (see Completed); no T-junction support — real T-splines are missing.
  - Consumers: `oxiphysics-fem` isogeometric-analysis upgrade (fem v0.3.0).
  - [ ] T-mesh data structure with analysis-suitability (T-junction validity) checks
  - [ ] Blending-function extraction; partition-of-unity verification to 1e-12
  - [ ] Watertight subdivision-limit evaluation (exact limit surface, including extraordinary vertices)
  - [ ] Bézier extraction toward FEM IGA consumption

## Deferred / research track
- [~] GPU ports of flood/extraction kernels (JFA passes, marching cubes) — Ready when: `oxiphysics-gpu` exposes a general compute-kernel API beyond its current SPH/MD/LBM pipelines.
- [~] Quality tetrahedral volume-meshing pipeline for FEM (sizing fields + sliver removal on top of CDT) — Ready when: 3D Delaunay + constrained Delaunay (v0.2.0 above) have landed.

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 263 (critical 3, high 86, medium 150, low 24); 42 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (3)
- [ ] **F1546** (critical, fabrication, easy/small) `implicit_geometry/functions.rs:202` — extract_isosurface uses garbage edge table and fan triangulation
- [ ] **F1555** (critical, stub, easy/small) `level_set/functions.rs:123` — Level-set MC tri table truncated (223/256 cases empty); edge table corrupt
- [ ] **F1723** (critical, bug, medium/small) `torus.rs:620` — Torus quartic solver mixes Ferrari conventions; ray casts wrong/miss

### High (86)
- [ ] **F1440** (high, bug, easy/trivial) `architectural_geometry.rs:562` — Cable-net form finding applies T/L force instead of T
- [ ] **F1443** (high, bug, easy/trivial) `architectural_geometry.rs:1053` — Solar azimuth computed from north but used as from south
- [ ] **F1447** (high, bug, easy/small) `architectural_geometry.rs:1190` — fit_plane returns zero normal for any x=c or y=c plane
- [x] **F1449** (high, bug, medium/small) `boolean_ops/functions.rs:400` — polygon2d_difference intersects outer half-planes; returns empty
- [x] **F1456** (high, fabrication, easy/small) `bspline/functions.rs:20` — solve_linear_system_row ignores matrix, returns data points as control points
- [x] **F1457** (high, panic, easy/small) `bspline/functions.rs:33` — solve_3x_systems fixed [f64;7] rows panic for n_ctrl >= 5
- [ ] **F1464** (high, bug, easy/trivial) `capsule.rs:62` — Capsule hemisphere inertia over-counts by 9/64 m r²
- [x] **F1470** (high, bug, easy/trivial) `compound/compound_traits.rs:88` — Compound::inertia_tensor ignores child rotation
- [ ] **F1471** (high, bug, easy/trivial) `compound/compound_traits.rs:97` — Compound::inertia_tensor ignores child rotations
- [ ] **F1472** (high, bug, easy/small) `compound/types.rs:243` — CompoundShapeEx::inertia_tensor ignores child LocalTransform rotation
- [ ] **F1473** (high, bug, easy/trivial) `compound/types.rs:539` — CompoundShape capsule child inertia uses wrong hemisphere offset
- [ ] **F1479** (high, bug, easy/small) `computational_geometry/functions.rs:330` — minkowski_sum loops forever pushing 1e18 sentinel vertices
- [x] **F1485** (high, bug, easy/trivial) `cone.rs:663` — Cone inertia tensor is about base diameter, not center of mass
- [ ] **F1486** (high, bug, easy/trivial) `cone.rs:670` — Cone inertia tensor is about base centre, not centre of mass
- [x] **F1490** (high, fabrication, hard/medium) `convex_decomposition.rs:481` — ApproxConvexDecomp concavity always 0, so decompose never splits
- [x] **F1492** (high, bug, medium/medium) `convex_decomposition.rs:688` — GjkDistance::distance returns bogus search-direction norm
- [ ] **F1493** (high, fabrication, medium/medium) `convex_decomposition.rs:731` — GjkDistance returns wrong distances; simplex handling is bogus
- [ ] **F1494** (high, bug, easy/small) `convex_decomposition.rs:926` — InertiaFromHull doubles cross terms, ignores COM, masks errors with abs
- [x] **F1495** (high, fabrication, medium/medium) `convex_hull/convexhull_traits.rs:67` — ConvexHull Shape impl: ray_cast None, AABB inertia, bogus volume/COM
- [x] **F1498** (high, bug, easy/small) `convex_hull/types.rs:22` — ConvexHull volume/COM sums tets over vertex pairs; cube volume 0
- [ ] **F1499** (high, fabrication, easy/small) `convex_hull/types.rs:30` — ConvexHull Shape volume/COM sum vertex pairs; inertia is AABB box
- [ ] **F1504** (high, fabrication, easy/small) `csg/functions.rs:261` — extract_mesh_simple emits flat z-plane triangles, not marching cubes
- [ ] **F1511** (high, perf, medium/medium) `decimation/types.rs:88` — QEM decimation recomputes all quadrics/edges per collapse: O(F²)
- [ ] **F1514** (high, bug, easy/small) `discrete_geometry.rs:218` — voronoi_area pairs cotangents with wrong edges; O(V·F) per call
- [ ] **F1515** (high, bug, medium/medium) `discrete_geometry.rs:542` — Heat-method geodesics: explicit Euler without mass matrix, O(n·nnz) Gauss-Seidel
- [ ] **F1516** (high, perf, medium/medium) `discrete_geometry.rs:599` — HeatMethodSolver Poisson Gauss-Seidel scans all triplets per vertex
- [ ] **F1518** (high, fabrication, hard/medium) `discrete_geometry.rs:957` — Hodge decomposition co-exact part is a normal vector, not N×∇g
- [ ] **F1521** (high, fabrication, medium/medium) `discrete_geometry.rs:1310` — ricci_flow_step ignores conformal factors; curvature never changes
- [ ] **F1527** (high, fabrication, medium/medium) `geodesic/functions.rs:317` — heat_geodesic discards heat/divergence, returns BFS-tree path lengths
- [ ] **F1536** (high, bug, medium/medium) `geodesic_geometry/functions.rs:408` — Heat-method divergence misses half the terms; Jacobi on singular Laplacian
- [ ] **F1541** (high, fabrication, medium/medium) `implicit.rs:474` — SdfReinitialize is not fast marching: one 1-D min pass
- [ ] **F1543** (high, bug, medium/medium) `implicit.rs:1003` — DualContouring: degenerate quads, no QEF, misses active cells
- [ ] **F1544** (high, bug, medium/small) `implicit.rs:1068` — MarchingTetrahedra: 5-tet split covers 2/3 of cube; case table wrong
- [ ] **F1545** (high, bug, medium/small) `implicit.rs:1154` — MarchingTetrahedra: degenerate tet leaves holes, normals hardcoded +Y
- [ ] **F1547** (high, bug, easy/small) `implicit_geometry/sdfcone_traits.rs:16` — SdfCone reports points above the base as inside
- [ ] **F1548** (high, bug, easy/trivial) `implicit_geometry/sdfhexagonalprism_traits.rs:20` — SdfHexagonalPrism reflection clamp inverted; outside points read inside
- [ ] **F1549** (high, bug, medium/small) `implicit_geometry/sdfpyramid_traits.rs:34` — SdfPyramid has wrong sign and magnitude (above apex reads inside)
- [ ] **F1550** (high, bug, easy/small) `implicit_geometry/sdftriangularprism_traits.rs:27` — SdfTriangularPrism never negative inside; points below prism read 0
- [ ] **F1551** (high, fabrication, easy/small) `implicit_surfaces.rs:532` — marching_cubes fan-triangulates, has a corrupt edge table, and evaluates nodes 8x
- [ ] **F1559** (high, bug, easy/small) `level_set/types.rs:572` — LevelSetField::mean_curvature: bogus φxy term, φxz/φyz hard-coded zero
- [ ] **F1560** (high, fabrication, medium/medium) `level_set/types.rs:789` — fast_sweep discards its Godunov work and only propagates along x
- [ ] **F1562** (high, fabrication, hard/large) `medial_axis.rs:65` — MedialAxis::from_points is a 5-point straight-line stub
- [ ] **F1570** (high, bug, easy/small) `medial_axis.rs:689` — distance_transform_exact lower envelope always collapses
- [ ] **F1575** (high, fabrication, easy/small) `medical_geometry.rs:925` — Medical marching_cubes fans active edges instead of using tri-table
- [ ] **F1577** (high, missing-feature, hard/large) `mesh_boolean.rs:846` — mesh_boolean keeps/drops whole triangles; never splits along intersections
- [ ] **F1580** (high, fabrication, medium/small) `mesh_param.rs:287` — LscmParameterization is fake LSCM; every v stays 0
- [ ] **F1581** (high, bug, medium/medium) `mesh_param.rs:466` — Tutte/Harmonic boundary placed in vertex-id order, not loop order
- [ ] **F1586** (high, bug, medium/medium) `mesh_processing/functions.rs:339` — qem_decimate collapses self-edges of dead faces, corrupting indices
- [ ] **F1587** (high, bug, easy/small) `mesh_processing/functions.rs:588` — fix_non_manifold_vertices is always a no-op
- [ ] **F1591** (high, bug, easy/trivial) `mesh_quality.rs:181` — Tet dihedral normals inconsistent: regular tet reports 109.47° max
- [ ] **F1593** (high, bug, easy/small) `mesh_quality.rs:409` — flip_edges_2d never fixes non-Delaunay edges, breaks Delaunay ones
- [ ] **F1594** (high, bug, easy/trivial) `mesh_quality.rs:512` — Scaled Jacobian mis-normalised: regular tet scores 0.5, not 1
- [ ] **F1598** (high, bug, easy/small) `mesh_repair/functions.rs:365` — Hole filling emits triangles with inverted winding
- [ ] **F1599** (high, bug, easy/small) `mesh_repair/functions.rs:449` — Orientation check/repair uses normal dot product instead of edge direction
- [ ] **F1601** (high, bug, medium/small) `mesh_repair_ext.rs:552` — RemeshingUniform::split_long_edges hangs for L<=0 and cracks mesh
- [ ] **F1607** (high, bug, hard/medium) `mesh_simplification/functions.rs:254` — qem_simplify applies stale heap candidates; O(T) scan per collapse
- [ ] **F1610** (high, stub, medium/medium) `mesh_simplification/functions.rs:560` — check_normal_deviation stub; QemConfig normal/seam/boundary options never read
- [ ] **F1614** (high, bug, easy/small) `nurbs_geometry.rs:877` — Bézier degree elevation blending coefficients swapped
- [ ] **F1615** (high, fabrication, medium/medium) `nurbs_geometry.rs:969` — NURBS degree elevation fabricated: sampled points, swapped coefficients, wrong knots
- [ ] **F1616** (high, stub, medium/small) `nurbs_geometry.rs:1150` — NurbsTrimming::is_inside always returns true; trimming never applied
- [ ] **F1635** (high, bug, easy/trivial) `point_cloud/types.rs:359` — Smallest-eigenvector power iteration fails for axis-aligned planar neighbourhoods
- [ ] **F1642** (high, bug, medium/small) `quickhull.rs:742` — hull_inertia_tensor drops tetra cross terms; uses vertex mean not COM
- [ ] **F1646** (high, bug, medium/medium) `remesh.rs:304` — Edge flips use stale edge map, invert winding, skip validity checks
- [ ] **F1648** (high, bug, easy/trivial) `robot_geometry.rs:413` — jacobian_6dof_numerical angular rows negated and in body frame
- [ ] **F1653** (high, fabrication, medium/small) `signed_distance_field/fmm.rs:231` — FMM is Manhattan BFS: eikonal quadratic dead, no signed propagation
- [ ] **F1666** (high, bug, medium/small) `sphere.rs:824` — Welzl MES: depth-n recursion, no shuffle, sign-flipped circumsphere_3 centre
- [ ] **F1667** (high, bug, easy/trivial) `sphere.rs:863` — circumsphere_3 returns mirrored centre; minimum_enclosing_sphere wrong
- [ ] **F1668** (high, bug, medium/small) `sphere.rs:891` — circumsphere_4 fallback not enclosing; Welzl recursion depth n, unshuffled
- [ ] **F1675** (high, bug, easy/small) `sphere_packing/types.rs:752` — HCP lattice uses cubic FCC cell, generated spheres overlap
- [ ] **F1676** (high, bug, easy/small) `sphere_packing/types.rs:776` — HCP lattice uses cubic cell; generated spheres overlap
- [ ] **F1682** (high, bug, medium/small) `subdivision.rs:349` — Butterfly stencil weights sum to 7/8, not affine invariant
- [ ] **F1684** (high, bug, medium/small) `subdivision.rs:541` — Catmull-Clark applies interior rule to boundary/corner vertices
- [ ] **F1687** (high, bug, easy/trivial) `swept/functions.rs:155` — Swept sphere TOI returns exit time/None when initially overlapping
- [ ] **F1689** (high, bug, easy/trivial) `swept/functions.rs:440` — capsule_capsule_toi bisects non-monotone distance, tunnels through contacts
- [ ] **F1690** (high, bug, medium/small) `swept/functions.rs:878` — ccd_sphere_vs_capsule samples only t=0.5/1, misses contacts
- [ ] **F1694** (high, bug, easy/trivial) `terrain_processing.rs:145` — compute_aspect mirrors north and south
- [ ] **F1696** (high, bug, medium/small) `terrain_processing.rs:339` — D-infinity flow direction is not Tarboton's algorithm
- [ ] **F1706** (high, fabrication, medium/medium) `topology/functions.rs:135` — betti_numbers_approx uses invented heuristic formulas
- [ ] **F1709** (high, fabrication, medium/medium) `topology/functions.rs:297` — Medial-axis 'circumcenter' is always midpoint of b and c
- [ ] **F1715** (high, stub, hard/medium) `topology/types.rs:1071` — TopologicalSurgery ops are no-ops/corrupting yet return success
- [ ] **F1719** (high, fabrication, hard/medium) `topology_geometry.rs:542` — ReebGraph::compute drops saddles and chains extrema by value
- [ ] **F1720** (high, fabrication, medium/medium) `topology_geometry.rs:661` — MorseComplex labels every regular vertex a saddle
- [ ] **F1727** (high, fabrication, medium/small) `triangle_mesh.rs:843` — TriangleMesh::inertia_tensor is AABB box inertia, not mesh inertia
- [ ] **F1730** (high, fabrication, hard/medium) `vhacd.rs:195` — V-HACD emits axis-aligned boxes, concavity is box fill ratio
- [ ] **F1732** (high, bug, medium/small) `voronoi/functions.rs:151` — Hull-site Voronoi cells closed from circumcentres only
- [ ] **F1733** (high, fabrication, easy/small) `voronoi/functions.rs:416` — voronoi_cell_volume_3d returns mean distance, not a volume

### Medium (150)
- [ ] **F1438** (medium, perf, medium/small) `crates/oxiphysics-geometry/Cargo.toml:16` — geometry depends on FEM for one LSCM solve; discards convergence
- [ ] **F1441** (medium, fabrication, easy/small) `architectural_geometry.rs:614` — prestress_energy uses fabricated rest length t/1e3
- [ ] **F1442** (medium, fabrication, medium/small) `architectural_geometry.rs:925` — Glass plate stress/deflection coefficients hard-coded to square-plate values
- [ ] **F1445** (medium, bug, medium/small) `architectural_geometry.rs:1115` — annual_solar_exposure sums 12 days, mislabelled kWh/m²/year
- [ ] **F1448** (medium, fabrication, hard/medium) `boolean_ops/functions.rs:390` — polygon2d_union returns convex hull of both polygons
- [ ] **F1451** (medium, missing-feature, hard/large) `boolean_ops/types.rs:378` — MeshBoolean keeps or drops whole triangles; no intersection splitting
- [ ] **F1453** (medium, bug, medium/small) `boolean_ops/types.rs:949` — RobustMeshBoolean tests distance to infinite planes, not triangles
- [ ] **F1454** (medium, fabrication, easy/small) `boolean_ops/types.rs:972` — ManifoldBoolean::execute computes manifold checks then discards them
- [ ] **F1455** (medium, bug, easy/trivial) `box_shape.rs:435` — obb_inertia_tensor_rotated uses R·I·Rᵀ with row-basis convention
- [ ] **F1460** (medium, panic, easy/trivial) `bspline/types.rs:748` — eval_deriv with k > degree indexes ders[k] out of bounds
- [ ] **F1461** (medium, panic, easy/small) `bspline/types.rs:909` — Knot/fit/NURBS constructors panic on invalid public input
- [ ] **F1462** (medium, fabrication, medium/medium) `bspline/types_2.rs:134` — NurbsSurface::fit_to_grid copies nearest samples, not least squares
- [ ] **F1463** (medium, bug, medium/small) `bspline/types_2.rs:397` — BsplineSurface::refine_knots uses the stale original knot vector for each insertion
- [ ] **F1465** (medium, bug, easy/small) `capsule.rs:951` — Capsule hemisphere inertia adds a spurious 9mr²/64 term
- [ ] **F1467** (medium, bug, medium/small) `cell_complex.rs:547` — rank_of_matrix division-free Gauss-Jordan overflows i32
- [ ] **F1468** (medium, fabrication, medium/medium) `cell_complex.rs:920` — CellularApproximation::map_degree fabricated from Euler-characteristic ratio
- [ ] **F1474** (medium, bug, easy/small) `compound/types.rs:947` — penetration_depth_sphere misses spheres whose centre is deep inside a child
- [ ] **F1475** (medium, perf, easy/small) `computational_geometry/delaunay_3d.rs:531` — is_duplicate linear scan makes delaunay_3d O(n²)
- [ ] **F1476** (medium, perf, medium/small) `computational_geometry/delaunay_3d.rs:817` — delaunay_3d allocates O(T) cavity mark array per insertion
- [ ] **F1480** (medium, bug, medium/small) `computational_geometry/functions.rs:392` — Visibility test lets segments pass through obstacle interiors
- [ ] **F1481** (medium, perf, medium/medium) `computational_geometry/functions.rs:528` — delaunay_2d is O(n²) Bowyer-Watson with a scale-dependent epsilon
- [ ] **F1482** (medium, bug, easy/small) `computational_geometry/types.rs:355` — SlabPointLocator::locate returns the topmost segment; face labels ignored
- [ ] **F1483** (medium, bug, medium/small) `cone.rs:101` — Cone::closest_point projects horizontally instead of onto slant
- [ ] **F1488** (medium, panic, medium/small) `convex_decomposition.rs:122` — ConvexHull3d::build panics on fewer than 4 points; coplanar seed gives a degenerate hull
- [ ] **F1491** (medium, fabrication, medium/small) `convex_decomposition.rs:634` — MinkowskiSum::contains_point tests only 10 fixed directions
- [ ] **F1500** (medium, perf, hard/large) `convex_hull/types.rs:136` — ConvexHull3D::build is naive incremental, not Quickhull as documented
- [ ] **F1503** (medium, bug, medium/small) `csg/functions.rs:136` — sphere_march_relaxed has no overstep fallback and tunnels thin features
- [ ] **F1506** (medium, bug, easy/small) `csg/functions.rs:476` — csg_surface_area staircase estimator returns about 1.5x the true area
- [ ] **F1507** (medium, bug, easy/small) `csg/functions.rs:955` — marching_cubes_sample uses the min step, truncating non-cubic bounds
- [ ] **F1508** (medium, bug, medium/small) `csg/sdfcone_traits.rs:39` — SdfCone::sdf overestimates distance below the apex; dead code
- [ ] **F1509** (medium, bug, medium/medium) `decimation/functions.rs:298` — QEM decimation: full recompute per collapse, no link check, feature vertices moved
- [ ] **F1510** (medium, bug, medium/medium) `decimation/functions.rs:337` — QEM decimators move feature vertices, skip link check, O(T^2)
- [ ] **F1513** (medium, bug, medium/small) `decimation/types.rs:511` — preserve_boundary moves or deletes boundary vertices
- [ ] **F1517** (medium, perf, medium/small) `discrete_geometry.rs:864` — hodge_decompose: two O(V·nnz) Gauss-Seidel loops, O(V·F) normals
- [ ] **F1519** (medium, fabrication, medium/medium) `discrete_geometry.rs:1096` — laplacian_eigenvalues_power returns k iterates, not k eigenvalues
- [ ] **F1520** (medium, fabrication, medium/medium) `discrete_geometry.rs:1112` — laplacian_eigenvalues_power returns k iterates of one power sequence
- [ ] **F1523** (medium, api, medium/medium) `error.rs:339` — Validation helpers are dead code; TriangleMesh::new unvalidated
- [ ] **F1524** (medium, panic, easy/medium) `fractal_geometry/types.rs:480` — diamond_square panics for roughness<=0; generators use unseeded RNG
- [ ] **F1525** (medium, perf, medium/small) `fractal_geometry/types.rs:745` — correlation_dimension stores all n²/2 pair distances; fake Richardson blend
- [ ] **F1529** (medium, fabrication, medium/medium) `geodesic/functions.rs:644` — intrinsic_delaunay ignores intrinsic lengths and uses 3-D chords
- [ ] **F1530** (medium, bug, medium/medium) `geodesic/functions.rs:654` — intrinsic_delaunay uses extrinsic 3D chords and rebuilds adjacency per flip
- [ ] **F1531** (medium, fabrication, easy/small) `geodesic/functions.rs:774` — parallel_transport discards angle and flips the along-edge component
- [ ] **F1532** (medium, fabrication, medium/small) `geodesic/functions.rs:820` — varadhan_distances: discarded heat_geodesic call, 5-hop saturated heat
- [ ] **F1533** (medium, fabrication, medium/small) `geodesic/functions.rs:826` — varadhan_distances discards the heat solve and uses unscaled 5-step diffusion
- [ ] **F1534** (medium, bug, medium/small) `geodesic_geometry/functions.rs:126` — FMM eikonal update solves an unfounded quadratic
- [ ] **F1535** (medium, panic, easy/trivial) `geodesic_geometry/functions.rs:212` — fast_marching_geodesic expect() panics on faces with repeated indices
- [ ] **F1537** (medium, perf, easy/small) `geodesic_geometry/types.rs:194` — voronoi_area/vertex_normal rebuild the vertex-to-face map on every call
- [ ] **F1539** (medium, panic, easy/small) `heightfield/types.rs:393` — HeightField deserialize/new accept degenerate sizes and panic downstream
- [ ] **F1540** (medium, perf, medium/medium) `heightfield/types.rs:851` — HeightField::ray_cast tessellates the whole terrain on every ray
- [ ] **F1542** (medium, bug, easy/small) `implicit.rs:909` — MarchingCubes vertex cache is keyed per cell, so vertices are never shared
- [ ] **F1552** (medium, bug, medium/small) `implicit_surfaces.rs:587` — RBF fit: CG on indefinite 2-D TPS matrix, no polynomial term
- [ ] **F1553** (medium, bug, medium/small) `implicit_surfaces.rs:600` — RBFImplicit::fit runs CG on an indefinite, possibly non-symmetric TPS system
- [ ] **F1556** (medium, perf, easy/small) `level_set/types.rs:59` — SignedDistanceTransform is brute force O(grid×points); promised fast sweep missing
- [ ] **F1557** (medium, fabrication, medium/small) `level_set/types.rs:90` — compute_grid_signed picks the sign with a sphere around the centroid
- [ ] **F1558** (medium, bug, easy/small) `level_set/types.rs:251` — normal_velocity_extension uses non-upwind central |grad phi|
- [ ] **F1563** (medium, fabrication, easy/small) `medial_axis.rs:172` — Voronoi2dSkeleton returns every triplet circumcentre, not Voronoi vertices
- [ ] **F1564** (medium, bug, easy/trivial) `medial_axis.rs:222` — StraightSkeleton offset miter uses sin(θ/2); collinear vertex scaled 1e9
- [ ] **F1565** (medium, fabrication, hard/large) `medial_axis.rs:238` — StraightSkeleton::skeleton_edges returns vertex-to-centroid spokes
- [ ] **F1566** (medium, fabrication, hard/medium) `medial_axis.rs:249` — StraightSkeleton::skeleton_edges just links every vertex to the centroid
- [ ] **F1567** (medium, fabrication, medium/small) `medial_axis.rs:445` — shape_diameter returns constant 2*radius, ignoring point and normal
- [ ] **F1568** (medium, stub, medium/small) `medial_axis.rs:451` — shape_diameter ignores point and normal and always returns 2*radius
- [ ] **F1569** (medium, bug, easy/trivial) `medial_axis.rs:590` — distance_transform is documented as Euclidean but computes L1
- [ ] **F1573** (medium, bug, easy/small) `medical_geometry.rs:491` — extract_centerline_approx loops forever if voxel_size[2] <= 0
- [ ] **F1576** (medium, perf, easy/small) `mesh_boolean.rs:733` — weld_vertices is an O(n²) linear scan run on every mesh_boolean
- [ ] **F1578** (medium, perf, medium/medium) `mesh_boolean.rs:1168` — Boolean classification and contour extraction are O(Ta·Tb) with no spatial index
- [ ] **F1579** (medium, perf, medium/small) `mesh_ops.rs:656` — quadric_decimate rebuilds all edges per collapse and places vertices at midpoints
- [ ] **F1582** (medium, bug, easy/trivial) `mesh_param.rs:827` — fill_holes fan triangles have inverted winding
- [ ] **F1583** (medium, perf, medium/small) `mesh_param.rs:900` — HarmonicMap Gauss-Seidel rescans all triangles per neighbour per sweep
- [ ] **F1584** (medium, bug, medium/small) `mesh_processing/functions.rs:231` — catmull_clark_subdivide applies wrong CC rules; vertices invert
- [ ] **F1588** (medium, missing-feature, easy/small) `mesh_processing/functions.rs:667` — isotropic_remesh claims an edge-flip step but never flips edges
- [ ] **F1589** (medium, bug, medium/small) `mesh_processing/functions.rs:1060` — LSCM output is clamped to [0,1], collapsing the conformal map
- [ ] **F1590** (medium, fabrication, medium/medium) `mesh_processing/functions.rs:1070` — generate_texture_atlas ignores atlas_size and returns one patch
- [ ] **F1592** (medium, bug, easy/small) `mesh_quality.rs:353` — smooth_laplacian 'tet_count >= 4' boundary test moves surface vertices
- [ ] **F1595** (medium, perf, easy/small) `mesh_quality.rs:892` — quality_guided_smoothing scans all tets per vertex and accepts inverted tets
- [ ] **F1596** (medium, perf, easy/small) `mesh_repair/functions.rs:105` — weld_vertices / merge_duplicate_vertices are O(n²)
- [ ] **F1597** (medium, bug, medium/small) `mesh_repair/functions.rs:207` — repair_t_junctions is O(T·V) and fixes only one edge per face
- [ ] **F1600** (medium, bug, medium/small) `mesh_repair/types.rs:240` — FlipMesh::flip_pass flips stale triangles without convexity test
- [ ] **F1602** (medium, bug, medium/small) `mesh_repair_ext.rs:575` — RemeshingUniform::split_long_edges cracks the mesh
- [ ] **F1603** (medium, stub, medium/medium) `mesh_repair_ext.rs:673` — IsotropicRemeshing omits edge flips and surface projection
- [ ] **F1604** (medium, perf, medium/small) `mesh_repair_ext.rs:788` — MeshDecimation::decimate is O(F²) brute force with midpoint target
- [ ] **F1605** (medium, fabrication, medium/small) `mesh_repair_ext.rs:882` — Loop subdivision uses plain midpoints; no boundary rule (two copies)
- [ ] **F1606** (medium, fabrication, medium/medium) `mesh_repair_ext.rs:1022` — MeshBoolean claims CSG; only concatenates and counts bbox overlaps
- [ ] **F1609** (medium, panic, easy/trivial) `mesh_simplification/functions.rs:531` — hausdorff_one_sided panics on n_samples=0
- [ ] **F1611** (medium, stub, medium/medium) `mesh_simplification/types.rs:237` — QemConfig normal/seam/boundary options silently ignored
- [ ] **F1612** (medium, stub, medium/medium) `mesh_simplification/types.rs:310` — ProgressiveMesh::build never records collapses
- [ ] **F1613** (medium, bug, easy/small) `nurbs_geometry.rs:710` — NURBS mean_curvature uses wrong shape-operator trace
- [ ] **F1617** (medium, bug, medium/small) `nurbs_geometry.rs:1214` — IgaElement stores only u-direction non-rational basis
- [ ] **F1618** (medium, fabrication, medium/small) `nurbs_geometry.rs:1294` — intersect_newton is fixed-step gradient descent, not Newton
- [ ] **F1621** (medium, duplicate, hard/large) `offset_surface.rs:61` — Duplicate public types: 4 NurbsCurve, 5 HalfEdgeMesh, 2 Sdf traits
- [ ] **F1622** (medium, bug, easy/small) `offset_surface.rs:415` — VoxelSdf trilinear sampling off by half a voxel; empty grid panics
- [ ] **F1623** (medium, bug, easy/small) `offset_surface.rs:847` — Convex/concave edge classification ignores the dihedral angle
- [ ] **F1625** (medium, fabrication, hard/medium) `origami.rs:336` — RigidOrigami::update_positions is not rigid folding
- [ ] **F1626** (medium, bug, medium/small) `origami.rs:417` — Miura-ori projected dimensions and Poisson ratio are wrong
- [ ] **F1629** (medium, bug, easy/trivial) `parametric/functions.rs:90` — bezier_arc_approximation middle control point is wrong
- [ ] **F1631** (medium, bug, easy/small) `parametric/functions.rs:159` — bezier_torsion third derivative doubled
- [ ] **F1632** (medium, perf, medium/small) `point_cloud/functions.rs:268` — fpfh_feature rebuilds a KD-tree per query; not real FPFH
- [ ] **F1636** (medium, fabrication, medium/small) `point_cloud/types.rs:614` — Principal curvatures returned as covariance eigenvalues
- [ ] **F1637** (medium, panic, easy/trivial) `point_cloud/types.rs:669` — ICP panics on empty target cloud
- [ ] **F1639** (medium, perf, hard/medium) `quickhull.rs:61` — 'QuickHull' is a naive O(n·F) incremental hull
- [ ] **F1640** (medium, panic, easy/small) `quickhull.rs:514` — NaN coordinates panic hull/voronoi builders via partial_cmp().expect
- [ ] **F1641** (medium, perf, hard/medium) `quickhull.rs:552` — IncrementalConvexHull rebuilds full hull on every insert
- [ ] **F1644** (medium, bug, medium/small) `remesh.rs:180` — isotropic_remesh split_long_edges creates T-junctions
- [ ] **F1645** (medium, bug, medium/small) `remesh.rs:227` — Edge collapse nondeterministic, O(E·V), no link condition
- [ ] **F1647** (medium, perf, medium/medium) `remesh.rs:365` — project_to_surface brute-force O(V·T) every iteration
- [ ] **F1649** (medium, panic, easy/small) `robot_geometry.rs:447` — jacobian_analytical panics/corrupts when q.len() < DOF
- [ ] **F1651** (medium, fabrication, medium/medium) `robot_geometry.rs:1102` — URDF joint origins mapped to fake DH parameters
- [ ] **F1652** (medium, fabrication, easy/small) `robot_geometry.rs:1201` — jacobian_condition_number is a fabricated proxy, not σmax/σmin
- [ ] **F1655** (medium, fabrication, medium/small) `signed_distance_field/marching_cubes.rs:199` — Marching cubes: hardcoded +Z normals, unwelded triangle soup
- [ ] **F1656** (medium, bug, easy/trivial) `signed_distance_field/marching_cubes.rs:444` — MC_TABLE row 203 omits crossing edge 1, leaving mesh holes
- [ ] **F1657** (medium, bug, easy/small) `signed_distance_field/octree.rs:104` — Octree SDF refinement misses surfaces in large cells
- [ ] **F1658** (medium, panic, easy/small) `spatial_hash.rs:82` — SpatialHash cell_size unvalidated; query_radius loop can hang/overflow
- [ ] **F1660** (medium, perf, medium/small) `spatial_hash.rs:262` — SpatialHashPos3D::k_nearest brute-forces all items despite doc
- [ ] **F1665** (medium, bug, easy/trivial) `sphere.rs:768` — Real SH Y_2^-2 and Y_3^-2 normalisation 2x too large
- [ ] **F1669** (medium, bug, easy/small) `sphere_packing/functions.rs:70` — structure_factor evaluates a single k=q(1,1,1), not isotropic S(q)
- [ ] **F1672** (medium, perf, easy/small) `sphere_packing/types.rs:234` — Void voxelization tests every voxel against every sphere
- [ ] **F1673** (medium, fabrication, medium/small) `sphere_packing/types.rs:334` — Void connectivity and percolation threshold are hardcoded
- [ ] **F1677** (medium, perf, medium/small) `sphere_packing/types.rs:957` — PackingOptimizer::step runs two O(n^2) all-pairs passes
- [ ] **F1679** (medium, perf, medium/small) `spline_geometry/functions.rs:198` — B-spline LSQ fit: dense Gauss-Jordan, per-row alloc, panics
- [ ] **F1680** (medium, bug, medium/small) `spline_geometry/types.rs:41` — SweptSurface frame flips 180° when |tangent.x| crosses 0.9
- [ ] **F1681** (medium, bug, easy/small) `subdivision.rs:255` — Loop subdivision boundary vertex and edge rules are wrong
- [ ] **F1685** (medium, stub, medium/small) `subdivision.rs:623` — AdaptiveCriterion::NearVertices is a stub returning false
- [ ] **F1686** (medium, fabrication, medium/medium) `subdivision.rs:675` — doo_sabin_step emits only shrunken face triangles
- [ ] **F1688** (medium, bug, easy/trivial) `swept/functions.rs:186` — swept_sphere_vs_plane misses initial penetration
- [ ] **F1691** (medium, bug, easy/trivial) `swept/types.rs:182` — SweptSphere::ray_intersect never hits a stationary sphere
- [ ] **F1692** (medium, bug, medium/small) `swept/types.rs:380` — SweptCapsule::toi_vs_sphere uses the capsule's bounding sphere
- [ ] **F1693** (medium, bug, easy/trivial) `swept/types.rs:570` — LinearExtrusion volume/surface_area use |sweep| instead of normal height
- [ ] **F1695** (medium, bug, easy/trivial) `terrain_processing.rs:193` — Curvature cross term f sign-inconsistent with g; profile denominator wrong
- [ ] **F1703** (medium, bug, easy/small) `terrain_processing.rs:949` — fill_sinks only raises single-cell pits
- [ ] **F1704** (medium, bug, medium/medium) `terrain_processing.rs:1024` — Marching squares returns an unordered point cloud and misses the last row/column
- [ ] **F1705** (medium, bug, easy/small) `topology/functions.rs:109` — compute_0d_persistence loses elder birth after union; one essential pair only
- [ ] **F1707** (medium, fabrication, medium/medium) `topology/functions.rs:141` — betti_numbers_approx uses invented formulas, no homology
- [ ] **F1708** (medium, fabrication, hard/large) `topology/functions.rs:257` — extract_curve_skeleton returns smoothed full mesh graph, not skeleton
- [ ] **F1710** (medium, bug, easy/trivial) `topology/types.rs:53` — Polygon offset miter length wrong by sqrt(2) factor
- [ ] **F1711** (medium, bug, easy/trivial) `topology/types.rs:201` — HalfEdgeMesh::num_edges counts by twin index position
- [ ] **F1712** (medium, bug, medium/small) `topology/types.rs:347` — HalfEdgeMesh::vertex_neighbors drops neighbours of boundary vertices
- [ ] **F1714** (medium, fabrication, medium/medium) `topology/types.rs:657` — PersistenceDiagram::bottleneck_distance is a one-sided nearest match
- [ ] **F1716** (medium, bug, easy/trivial) `topology_geometry.rs:209` — minimal_torus triangle list is non-manifold
- [ ] **F1717** (medium, perf, easy/small) `topology_geometry.rs:282` — compute_0d materialises and sorts all n²/2 pairwise edges
- [ ] **F1721** (medium, fabrication, hard/medium) `topology_geometry.rs:769` — AlphaComplex includes all short triples, not Delaunay subset
- [ ] **F1722** (medium, bug, easy/small) `torus.rs:448` — Torus::intersects_aabb false negatives from corner-only SDF test
- [ ] **F1724** (medium, bug, easy/trivial) `triangle_mesh.rs:319` — laplacian_smooth moves boundary vertices despite doc
- [ ] **F1725** (medium, fabrication, hard/medium) `triangle_mesh.rs:679` — Heat Kernel Signature is a per-vertex degree proxy, not HKS
- [ ] **F1726** (medium, perf, medium/small) `triangle_mesh.rs:732` — smooth_laplacian scans all triangles per vertex: O(V·F)
- [ ] **F1729** (medium, bug, easy/small) `vhacd.rs:109` — VHacdConfig unvalidated: resolution=0 leads to inf grid and likely panic
- [ ] **F1731** (medium, bug, medium/small) `voronoi/functions.rs:137` — delaunay_to_voronoi truncates hull cells; Lloyd/from_sites inherit the error
- [ ] **F1734** (medium, perf, easy/small) `voronoi/functions.rs:496` — voronoi_edges O(T^2) and delaunay_to_voronoi O(n·T)
- [ ] **F1736** (medium, bug, easy/small) `voronoi/functions.rs:544` — power_adjacent tests only the midpoint and ignores weights
- [ ] **F1737** (medium, fabrication, hard/medium) `voronoi/functions.rs:570` — fortune_delaunay claims O(n log n) sweep, is Bowyer-Watson O(n^2)
- [ ] **F1739** (medium, fabrication, medium/small) `voronoi/functions.rs:1572` — power_diagram cells are convex hulls of 60x60 grid samples
- [ ] **F1740** (medium, bug, medium/medium) `voronoi/types.rs:148` — Voronoi cells of convex-hull sites truncated (unbounded region dropped)

### Low (24) — tracked in the external ledger
F1450, F1452, F1466, F1487, F1497, F1501, F1522, F1528, F1538, F1561, F1572, F1574, F1624, F1634, F1643, F1659, F1663, F1670, F1674, F1683, F1700, F1718, F1728, F1742

<!-- AUDIT-2026-10-12 END -->
