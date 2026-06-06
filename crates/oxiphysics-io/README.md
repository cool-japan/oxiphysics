# oxiphysics-io

**Status: Alpha** | Version 0.1.3 | 2026-06-06

Multi-format file I/O and serialization for the OxiPhysics engine. Pure Rust with zero unsafe code.

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-io` provides readers and writers for 80+ scientific and engineering file formats spanning
molecular simulation, computational fluid dynamics, finite element analysis, medical imaging, geophysics,
and more.

## Domains Covered

| Domain | Formats |
|---|---|
| Molecular / MD | PDB, LAMMPS dump, GROMACS (GRO), Amber |
| CFD / FEM | OpenFOAM, Exodus, CGNS, Abaqus, CalculiX, Fluent, SU2, Plot3D |
| Mesh / Geometry | STL, OBJ, GMSH, VTK, VTU, XDMF, GLTF |
| Scientific Data | HDF5, netCDF/xarray, NumPy-compatible, CSV, JSON, binary |
| Trajectory | Generic trajectory, XDMF temporal, VTK time series |
| Specialized | Tecplot, Seismic, DICOM-style medical, Weather, Quantum chemistry, Robotics |

## Key APIs

```rust
use oxiphysics_io::*;

// Molecular
let atoms = PdbReader::open("protein.pdb")?.read_model()?;
PdbWriter::new("out.pdb").write(&atoms)?;
let traj = LammpsDumpReader::open("dump.lammpstrj")?.frames().collect::<Vec<_>>();

// CFD / Mesh
VtkWriter::new("field.vtk").write_unstructured(&mesh, &fields)?;
VtuWriter::new("field.vtu").write_parallel(&mesh, &fields)?;
write_xdmf_particles("particles.xdmf", &positions, &attributes)?;
write_xdmf_temporal("time.xdmf", &frames)?;

// General
let records = CsvReader::open("data.csv")?.into_records()?;
let obj = ObjReader::open("model.obj")?.read()?;
```

## Statistics

- **5,978** public items
- **4,879** tests
- **0** stubs — fully implemented

