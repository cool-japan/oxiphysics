# oxiphysics-geometry

**[Alpha]** Geometric shapes, mesh processing, and spatial algorithms for the OxiPhysics engine.

[![Tests](https://img.shields.io/badge/tests-3120-brightgreen)](https://github.com/cool-japan/oxiphysics)
[![docs.rs](https://img.shields.io/docsrs/oxiphysics-geometry)](https://docs.rs/oxiphysics_geometry)
[![version](https://img.shields.io/badge/version-0.1.3-blue)](https://crates.io/crates/oxiphysics-geometry)

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-geometry` delivers a comprehensive geometry library covering primitive shapes,
advanced mesh processing, parametric curves and surfaces, computational geometry, and
procedural generation. 57+ modules, 3,245 public items, 3,120 passing tests, 0 stubs.

## Features

- **Primitive shapes** — `Sphere`, `BoxShape`, `Capsule`, `Cylinder`, `Cone`, `Torus`, `Compound`
- **Meshes** — `TriangleMesh`, `ConvexHull`, `HeightField`; mesh ops, boolean ops, repair, quality metrics
- **Parametric & implicit** — NURBS, B-splines, implicit surfaces, level sets, signed distance fields
- **Computational geometry** — Voronoi diagrams, Quickhull, convex hull, medial axis
- **Subdivision & remeshing** — Loop/Catmull-Clark subdivision, decimation, remeshing, simplification
- **Point clouds** — point cloud processing, spatial hashing
- **Procedural** — fractal geometry, terrain processing, origami, geodesic domes
- **Topology** — topology analysis, swept volumes, offset surfaces
- **Robot geometry** — kinematic chain geometry primitives
- **CSG** — constructive solid geometry (union, intersection, difference)

## Modules (57+)

`bspline` · `capsule` · `convex_hull` · `csg` · `cylinder` · `decimation` ·
`fractal_geometry` · `geodesic` · `heightfield` · `implicit_surfaces` · `level_set` ·
`medial_axis` · `mesh_boolean` · `mesh_ops` · `mesh_processing` · `mesh_quality` ·
`mesh_repair` · `mesh_simplification` · `nurbs_geometry` · `offset_surface` · `origami` ·
`parametric` · `point_cloud` · `procedural_geometry` · `quickhull` · `remesh` ·
`robot_geometry` · `shape` · `signed_distance_field` · `spatial_hash` · `sphere` ·
`subdivision` · `swept` · `terrain_processing` · `topology` · `topology_geometry` ·
`torus` · `triangle_mesh` · `voronoi` · …

## Quick Example

```rust
use oxiphysics_geometry::{TriangleMesh, ConvexHull, Sphere, SignedDistanceField};

// Build a triangle mesh and compute its convex hull
let mesh = TriangleMesh::from_obj("model.obj").expect("load mesh");
let hull = ConvexHull::from_mesh(&mesh);
println!("hull vertices: {}", hull.vertex_count());

// Query SDF at a point
let sphere = Sphere::new(1.0);
let sdf = SignedDistanceField::from_shape(&sphere, 0.05);
let dist = sdf.query([0.5, 0.0, 0.0].into());
println!("distance to surface: {dist:.4}");
```

## License

Apache-2.0 — Copyright 2026 COOLJAPAN OU (Team Kitasan)
