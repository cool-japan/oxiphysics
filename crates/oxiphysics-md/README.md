# oxiphysics-md

**Status: Stable** | Version 0.1.3 | 2026-06-06

Molecular dynamics simulation for the OxiPhysics engine. Pure Rust.

Part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) project.

## Overview

`oxiphysics-md` is the largest crate in the workspace by public API surface (6,343 items). It covers
75+ modules spanning classical MD, enhanced sampling, QM/MM, coarse-grained models, free energy
calculations, and domain-specific simulations from protein folding to battery electrochemistry.

## Domains Covered

| Domain | Modules |
|---|---|
| Core MD | atom, integrator, neighbor_list, potential, simulation, trajectory |
| Force fields | Amber, CHARMM, coarse-grained, polarizable, water models |
| Thermodynamics | thermostat, barostat, NHC (Nosé-Hoover chain), REMD |
| Enhanced sampling | metadynamics, steered MD, path integral, rare event, REMD, FEP |
| Free energy | free_energy, FEP, thermodynamic integration |
| QM/MM | qmmm, ab initio MD, quantum MD, nn_potential, ml_potential |
| Biomolecular | protein, DNA mechanics, DNA simulation, lipid membrane, membrane sim, biomolecular |
| Materials | crystal growth, glass MD, metal alloy MD, solid state MD, nanoparticle MD |
| Electrostatics | Ewald, electrostatics, electrochemistry, ionic liquid MD |
| Applied | battery MD, polymer MD, zeolite MD, tribochemistry MD, adsorption |
| Quantum | photochemistry MD, quantum chemistry, quantum MD, path integral MD |

## Key APIs

```rust
use oxiphysics_md::*;
// Re-exports coarse_grained::*, metadynamics::*, qmmm::*

// Classical simulation
let mut sim = MdSimulation::new(config)?;
sim.add_atoms(atoms);
sim.run(n_steps)?;

// Thermostat / Barostat
let nhc = NoseHooverChain::new(target_temp, tau, n_chains)?;
let barostat = ParrinelloRahmanBarostat::new(target_pressure, tau)?;

// Free energy
let alchemical = FepSampler::new(lambda_schedule)?;
let ddg = alchemical.run(&mut sim)?;

// Metadynamics (re-exported)
let meta = MetadynamicsDriver::new(cvs, gaussian_height, sigma)?;
meta.run(&mut sim, n_steps)?;

// QM/MM (re-exported)
let qmmm = QmMmSystem::new(qm_region, mm_region)?;
```

## Statistics

- **6,343** public items — highest in the workspace
- **5,171** tests
- **0** stubs — fully implemented

