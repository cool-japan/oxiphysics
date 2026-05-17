# oxiphysics-io TODO

Last updated: 2026-05-17 / v0.1.1

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

## Outstanding (v0.2.0)
No local v0.2.0 items at this level. Format additions will be driven by user demand; log feature requests in root `TODO.md` Deferred section.
