# oxiphysics-lbm

**Status: Alpha** | Version 0.1.1 | 2026-05-17

Lattice Boltzmann method (LBM) simulation for the OxiPhysics engine. Pure Rust.

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-lbm` provides a comprehensive LBM framework covering 40+ modules for single-/multi-phase
flow, thermal, turbulent, reactive, biological, geophysical, and quantum-inspired lattice-Boltzmann
simulations in 2D and 3D.

## Domains Covered

| Domain | Modules |
|---|---|
| Core lattice | D2Q9, D3Q19, D3Q27 lattices; BGK, MRT, TRT, Entropic collision |
| Multiphase | Shan-Chen, phase field, solidification |
| Thermal / Reactive | Thermal LBM, reactive flow, combustion |
| Turbulence | Smagorinsky SGS, turbulence model, wall model |
| Biofluid | Hemodynamics, non-Newtonian flow |
| Specialized | Acoustics, aeroacoustics, MHD, electrokinetic, porous media, particle-laden |
| Applied | Climate, geophysical, traffic flow, soft-matter, quantum LBM |
| Numerics | Immersed boundary, Zou-He BC, forcing, streaming, initialization |

## Key APIs

```rust
use oxiphysics_lbm::*;

// 2D BGK simulation
let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
bgk_collide_2d(&mut grid, omega)?;
stream_2d(&mut grid)?;

// 3D D3Q27
let mut lattice = D3Q27Lattice::new(nx, ny, nz);
bgk_collide_3d(&mut lattice, omega)?;
stream_3d(&mut lattice)?;

// Multi-relaxation-time
let mrt = MrtCollision2D::new(relaxation_rates)?;
mrt.collide(&mut grid)?;

// Multiphase Shan-Chen
let mut model = ShanChenModel::new(g_coupling);
model.step(&mut grid)?;

// Full simulation driver
let mut sim = LbmSimulation2D::new(config)?;
sim.run(n_steps)?;
```

## Statistics

- **5,649** public items
- **5,313** tests — **largest test suite in the workspace**
- **0** stubs — fully implemented

