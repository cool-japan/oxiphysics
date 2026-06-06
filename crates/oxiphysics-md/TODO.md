# oxiphysics-md TODO

Last updated: 2026-06-06 / v0.1.2

Molecular dynamics subcrate. All listed items are production code, tested. This document enumerates what exists; root `TODO.md` Phase 21.5 / 21.6 cover v0.2.0 validation targets.

## Phase 1: Foundation
- [x] Core types (Atom, AtomSet, Topology, Bond, Angle, Dihedral, Improper)
- [x] Neighbor lists (Verlet + cell list)
- [x] Basic error handling
- [x] Unit tests

## Phase 2: Potentials
- [x] Lennard-Jones (with cutoff + shift)
- [x] Coulomb / direct electrostatics
- [x] Morse
- [x] Harmonic bond / angle
- [x] Cosine dihedral
- [x] Improper (harmonic + Fourier)

## Phase 3: Force fields
- [x] AMBER — Coulomb, non-bonded 1-2 / 1-3 / 1-4 exclusions, harmonic bonded terms (`amber/`)
- [x] CHARMM — Urey-Bradley, n-fold dihedrals, CMAP-style support (`charmm.rs`)
- [x] OPLS-AA
- [x] ReaxFF — bond-order reactive potential

## Phase 4: Long-range electrostatics
- [x] Ewald summation (direct erfc real-space + k-space reciprocal + self-energy correction) (`ewald/`)
- [x] Particle-Mesh Ewald (B-spline grid spreading, FFT-backed reciprocal) (`ewald/`)
- [x] Structure factors and virial/pressure computation

## Phase 5: Integration & thermostats / barostats
- [x] Velocity-Verlet
- [x] Leapfrog / position-Verlet
- [x] Berendsen thermostat
- [x] Nosé-Hoover thermostat / chain
- [x] Langevin thermostat
- [x] Parrinello-Rahman / Berendsen barostat

## Phase 6: Enhanced sampling & advanced
- [x] Replica exchange (REMD)
- [x] Umbrella sampling
- [x] Metadynamics
- [x] QM/MM coupling
- [x] Lipid bilayer helpers
- [x] Protein-folding templates

## Outstanding (v0.2.0)
See root `TODO.md`:
- Phase 21.5 — Lennard-Jones liquid RDF at ρ* = 0.85, T* = 0.71 vs Verlet (1967)
- Phase 21.6 — NaCl Madelung constant + dipole-array potential vs analytical
No local v0.2.0 items at this level.
