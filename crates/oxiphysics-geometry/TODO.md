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

- [ ] 3D Delaunay + constrained Delaunay + alpha shapes
  - **Goal:** tetrahedralize 10⁴ points with empty-sphere property holding for all tets; α-shape recovers known surface.
  - **Design:** incremental Bowyer-Watson 3D + CDT + α-complex. Edelsbrunner-Mücke 1994.
  - Verified: `computational_geometry` has 2D `delaunay_2d` + `voronoi_from_delaunay` (see Completed); 3D/CDT/alpha are the open parts.
  - Depends on: core exact `insphere` predicate for robust empty-sphere tests.
  - [ ] Incremental Bowyer-Watson 3D with walk-based point location
  - [ ] Cavity retriangulation + degenerate-input handling via exact predicates
  - [ ] Constrained Delaunay: boundary facet recovery
  - [ ] α-complex filtration + α-shape extraction with known-surface recovery test
  - [ ] Empty-sphere audit over every tet of a 10⁴-point tetrahedralization

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

- [ ] Jump-flooding SDF + screened-Poisson reconstruction
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
