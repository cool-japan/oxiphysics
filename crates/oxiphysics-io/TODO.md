# oxiphysics-io TODO

Last updated: 2026-06-06 / v0.1.3

File-format I/O subcrate. All listed items are production code, tested.

## Phase 1: Foundation
- [x] Core types (trait-based reader/writer abstractions, common error type)
- [x] Basic error handling
- [x] Unit tests

## Phase 2: Mesh / geometry formats
- [x] VTK legacy / XML (VTU, PVTU if present)
- [x] glTF 2.0 (scene export for visualization)
- [x] OBJ (Wavefront)
- [x] STL (binary + ASCII)
- [x] Gmsh .msh

## Phase 3: Solver-specific formats
- [x] CalculiX (.inp / .frd) — FEM solver compat
- [x] LAMMPS (data / dump / restart) — MD solver compat
- [x] OpenFOAM case directory
- [x] Abaqus .inp subset

## Phase 4: Scientific / trajectory
- [x] HDF5 containers
- [x] NetCDF
- [x] XDMF (goes with VTK HDF5 backend)
- [x] CSV / JSON time-series

## Phase 5: Scene / state
- [x] Scene JSON round-trip (via oxiphysics scene module)
- [x] Snapshot binary format (fast save / restore)

## v0.1.2 fixes
- [x] (2026-06-01) Re-exported `particle_formats` public types at crate root (`DcdWriter`, `DcdReader`, `XyzWriter`, `XyzReader`, `ParticleFrame`, `ParticleTrajectory`, `TrajectoryStats`, `BinaryFrameReader`, `BinaryFrameWriter`, `GroReader`, `GroWriter`, `DcdHeader`). Fixes failing `DcdWriter` doctest.

## Outstanding (v0.2.0)
No additional items at this level. Format additions will be driven by user demand.
