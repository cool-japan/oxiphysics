# oxiphysics-lbm TODO

Last updated: 2026-06-01 / v0.1.2

Lattice Boltzmann fluid-dynamics subcrate. All listed items are production code
with unit tests in-tree. This document enumerates what actually ships; it is not
a roadmap. See the repository root `TODO.md` for the v0.2.0 roadmap (in
particular Phase 21.3 — lid-driven cavity centerline validation against
Ghia, Ghia & Shin 1982).

## Phase 1: Foundation
- [x] Core types & traits (`Lattice`, `LatticeDimensions`, `LatticeType`,
      distribution arrays, velocity sets) in `lattice/`
- [x] Structured error type with domain-specific variants (`error.rs`)
- [x] Grid containers (`LbmGrid2D`, `LbmGrid3D`, D2Q9 specialised grid) in
      `grid/` and `lattice/`
- [x] Initialisation helpers (`initialization.rs`): equilibrium init, shear,
      Taylor–Green vortex, noise perturbation
- [x] Unit tests co-located with every module

## Phase 2: Lattices
- [x] D2Q9 (2-D nine-velocity) — `lattice/types.rs`
- [x] D3Q19 (standard isothermal 3-D) — `lattice/types.rs`, `d3q19_full.rs`
- [x] D3Q27 (higher moment isotropy) — `lattice/types.rs`, `d3q27.rs`

Note: D3Q15 is *not* implemented (kept off the list intentionally).

## Phase 3: Streaming (`streaming.rs`)
- [x] Standard push streaming (2-D and 3-D)
- [x] Pull-scheme streaming (2-D and 3-D)
- [x] Push-scheme streaming (2-D and 3-D) with mass-conservation tests
- [x] AA-pattern in-place streaming (cache-friendly, no second buffer)
- [x] Swap streaming
- [x] Streaming with body-force correction
- [x] Half-way bounce-back streaming, including moving-wall variant
- [x] Periodic shift helpers, bulk bounce-back helpers
- [x] Snapshot / restore helpers (`snapshot_2d`, `snapshot_3d`,
      `restore_snapshot_*`, `max_diff_snapshots`)

## Phase 4: Collision operators (`collision/types.rs`)
- [x] `BgkCollision` — single-relaxation-time BGK
- [x] `BgkOverrelaxation` — over-relaxed BGK variant
- [x] `TrtCollision` — two-relaxation-time, magic parameter
      `Λ = 3/16` for slip-free no-slip walls
- [x] `MrtCollision` / `mrt.rs` / `mrt3d.rs` — multi-relaxation-time
      (moment-space collision) for 2-D and 3-D
- [x] `CumulantCollision` — full cumulant space, Geier et al. 2015
      (high-Reynolds stability)
- [x] `RegularizedCollision` — classic regularised LBM
- [x] `RegularizedCollisionFull` — Hermite-reconstruction of Π⁽¹⁾
      non-equilibrium stress
- [x] `RecursiveRegularized` — recursive regularised LBM
- [x] `HybridRecursiveRegularized` — HRR with numerical dissipation blending
- [x] `EntropicCollision` — H-theorem-enforcing ELBM
- [x] `KbcCollision` — Karlin–Bösch–Chikatamarla entropic stabiliser
- [x] `HybridCollision` — BGK ↔ entropic mode switching
      (`CollisionMode` enum)
- [x] `CentralMomentCollision` — central-moment space
- [x] `CascadedCollision` — cascaded LBM
- [x] `RawMomentCollision` — raw-moment space
- [x] Dedicated entropic module (`entropic.rs`) with Newton-Raphson α-finder,
      H-function, KL / symmetric-KL divergence, entropy over-relaxation,
      neq-entropy indicator, KBC-D2Q9 diagnostics

## Phase 5: Forcing schemes (`forcing.rs`)
- [x] Guo forcing (`GuoForcing`, `GuoForcingScheme`, `apply_guo_forcing_d3q19`)
- [x] Guo body force in `lattice/types.rs`
- [x] He–Luo forcing
- [x] Shan–Chen forcing scheme
- [x] Exact Difference Scheme (EDS)
- [x] Body-force helpers (`BodyForce`, `BodyForceType`)
- [x] Gravitational forcing, Boussinesq forcing, oscillating force,
      rotating-frame force, Coriolis force
- [x] Smagorinsky forcing, MHD Lorentz force

## Phase 6: Boundary conditions (`boundary.rs`, `curved_boundary.rs`, `zou_he/`)
- [x] `BoundaryType::Periodic`
- [x] `BoundaryType::Velocity` (equilibrium inflow)
- [x] `BoundaryType::Pressure` (non-equilibrium extrapolation)
- [x] `BoundaryType::ZouHeVelocity` and `ZouHePressure` (full Zou–He set)
- [x] `BoundaryType::ConvectiveOutflow` (`u_conv`-based)
- [x] `BoundaryType::ExtrapolationOutflow` (zero-gradient Neumann)
- [x] `BoundaryType::MovingWall` (lid-driven cavity)
- [x] `BoundaryType::InterpolatedBounceBack` (Bouzidi curved-wall)
- [x] Bulk bounce-back, half-way bounce-back (in `streaming.rs`)
- [x] Wall-model boundary (`wall_model.rs`) for near-wall turbulence
- [x] Immersed-boundary method (`immersed_boundary/`, `immersed_boundary_lbm.rs`)

## Phase 7: Thermal / heat transfer (`thermal/`, `thermal_lbm.rs`,
`heat_transfer_lbm.rs`, `conjugate_heat.rs`)
- [x] `ThermalD2Q9`, `ThermalLbm` (double-distribution momentum+temperature)
- [x] Boussinesq coupling (`BoussinesqCoupling`)
- [x] Canonical setups: De Vahl Davis natural convection,
      Rayleigh–Bénard, natural-convection general setup
- [x] Nusselt-number computation
- [x] Conjugate heat transfer across solid/fluid interfaces
- [x] Standalone `heat_transfer_lbm.rs` application module

## Phase 8: Multiphase & phase-field (`multiphase/`, `phase_field/`,
`cahn_hilliard*.rs`, `phase_separation.rs`, `droplet_dynamics*`)
- [x] Shan–Chen pseudo-potential (`ShanChenModel`, plus multi-component
      `MultiComponentSC`) with multiple `PsiType` pseudo-potential forms
- [x] Free-energy model (`FreeEnergyModel`)
- [x] Cahn–Hilliard (`cahn_hilliard.rs`, `cahn_hilliard_lbm.rs`)
- [x] Phase-field module (`phase_field/`, `phase_field_lbm.rs`)
- [x] Spinodal decomposition parameters and bubble-state tracking
- [x] Phase-separation standalone module
- [x] Droplet-dynamics modules (`droplet_dynamics.rs`,
      `droplet_dynamics_lbm/`)

## Phase 9: Turbulence (`turbulence/`, `turbulence_model/`, `turbulent_*`)
- [x] Static Smagorinsky (`SmagorinskyModel`)
- [x] Dynamic Smagorinsky (`DynamicSmagorinsky`)
- [x] k-ε model (`KEpsilonState`)
- [x] k-ω SST (`KOmegaSst`)
- [x] LES-to-DNS transition helper (`LesToDns`)
- [x] Turbulent Prandtl presets (`TurbPrandtlPreset`)
- [x] Turbulent-channel canonical setup (`turbulent_channel/`)
- [x] Turbulent dispersion (`turbulent_dispersion_lbm.rs`)

## Phase 10: Non-Newtonian rheology (`non_newtonian/`)
- [x] Power-law fluid (`power_law.rs`)
- [x] Bingham plastic (`bingham.rs`)
- [x] Herschel–Bulkley (`herschel_bulkley.rs`)
- [x] Carreau (`carreau.rs`)
- [x] Viscoelastic (`viscoelastic.rs`, `viscoelastic_lbm.rs`)
- [x] LBM integration layer and dedicated test suite

## Phase 11: Particle coupling
- [x] Lagrangian particle coupling (`particle_coupling.rs`)
- [x] Particle-laden flow (`particle_laden_lbm.rs`)
- [x] LBM particles helper (`lbm_particles.rs`)
- [x] Suspension flow (`suspension_lbm.rs`)
- [x] Sedimentation (`sedimentation/`, `sediment_transport.rs`)
- [x] Granular LBM (`granular_lbm.rs`)

## Phase 12: Porous media & multiscale
- [x] Porous-media LBM (`porous/`, `porous_media/`, `porous_media_lbm.rs`)
- [x] Multiscale LBM (`multiscale/`, `multiscale_lbm.rs`)

## Phase 13: Reactive flow & combustion
- [x] Species transport (`reactive/species.rs`, `concentration.rs`)
- [x] Catalytic reactions (`reactive/catalytic.rs`)
- [x] Combustion (`reactive/combustion.rs`, `combustion_lbm.rs`)
- [x] Reactive LBM integration (`reactive/lbm_reactive.rs`,
      `reactive_flow.rs`)
- [x] Extended reactive kernels and test suites

## Phase 14: Electromagnetic / charged-fluid LBM
- [x] Electrokinetic (`electrokinetic/`, `electrokinetic_lbm.rs`,
      `electrokinetics.rs`, `electrokinetics_lbm.rs`)
- [x] Electro-osmotic (`electroosmotic_lbm.rs`)
- [x] Electrostatic (`electrostatic_lbm.rs`)
- [x] Magnetohydrodynamics (`magnetohydrodynamics_lbm.rs`)
- [x] Ferrofluid (`ferrofluid_lbm.rs`)
- [x] Plasma (`plasma_lbm.rs`)

## Phase 15: Acoustics & compressible
- [x] Acoustic LBM (`acoustic_lbm.rs`, `acoustics_lbm.rs`)
- [x] Acoustic streaming (`acoustic_streaming_lbm.rs`)
- [x] Aeroacoustics (`aeroacoustics/`, `aeroacoustics_lbm.rs`)
      including A-weighting filter and near-to-far-field traits
- [x] Compressible LBM (`compressible.rs`)

## Phase 16: Bio- & soft-matter applications
- [x] Hemodynamics (`hemodynamics_lbm.rs`)
- [x] Biofluids (`biofluid_lbm/`, `biofluids_lbm.rs`) with
      red-blood-cell membrane traits
- [x] Biofilm (`biofilm_lbm.rs`)
- [x] Soft matter (`soft_matter_lbm.rs`)
- [x] Polymer (`polymer_lbm.rs`)
- [x] Microfluidics (`microfluidics.rs`, `microfluidics_lbm.rs`)
- [x] Mixing (`mixing_lbm.rs`)
- [x] Solidification (`solidification_lbm.rs`)

## Phase 17: Large-scale / geophysical / specialty
- [x] Geophysical LBM (`geophysical_lbm.rs`)
- [x] Climate LBM (`climate_lbm/` — energy balance, greenhouse gas, ENSO,
      Hadley-cell, monsoon, carbon-cycle, cloud-feedback,
      ice-albedo, ocean-heat-uptake, permafrost, sea-level-rise traits)
- [x] Quantum LBM (`quantum_lbm.rs`)
- [x] Neural LBM (`neural_lbm.rs`)
- [x] Traffic-flow LBM (`traffic_flow_lbm.rs`)

## Phase 18: Simulation driver & I/O
- [x] Simulation orchestration (`simulation/` with statistics traits)
- [x] Hybrid-LBM dispatcher (`hybrid_lbm.rs`)
- [x] LBM optimisation helpers (`lbm_optimization.rs`)
- [x] Diffusion module (`diffusion.rs`)
- [x] Snapshot / restore helpers (in `streaming.rs`) used by restart paths
- [x] Unit and integration tests co-located with each module

Note: no VTK/Paraview writer is implemented in this crate; visualisation
currently goes through consumer crates.

## v0.1.2 correctness fixes (2026-06-01)
- [x] `LbmGrid3D::step(omega: f64)` — implemented real BGK collide-and-stream step (equilibrium collision + pull-scheme periodic streaming + `compute_macroscopic`). Removed dead empty `step_placeholder` (zero call sites). 3 new tests: mass conservation, equilibrium fixed-point, aggressive-omega.

## Outstanding (v0.2.0)

No subcrate-local v0.2.0 items justified by the current audit. The only
outstanding validation item that touches this crate lives in the root
`TODO.md` under **Phase 21.3**: centerline-velocity validation of the
lid-driven cavity at Re = 100, 400, 1000 against Ghia, Ghia & Shin (1982).
