# oxiphysics-lbm TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 104,521 SLoC | 5,320 tests

Lattice Boltzmann fluid-dynamics subcrate. Everything under **Completed** is
production code with unit tests in-tree — Phases 1–18 below enumerate what
actually ships, preserved as release history. The v0.2.0 / v0.3.0 / v1.0
sections are the forward roadmap, code-verified against the tree on
2026-06-11: where an item builds on infrastructure that already ships, the
shipped piece is recorded under *Code-audit confirmations* and only the
missing delta remains open. The former root-TODO Phase 21.3 pointer
(lid-driven cavity vs Ghia, Ghia & Shin 1982) is resolved here: Re=100
shipped with a regression baseline; the Re=400/1000 extension lives under
v1.0 below. All work is pure Rust; GPU paths target wgpu by default, with the
cudarc (CUDA) backend strictly feature-gated and non-default.

How to read this file: `[x]` shipped, `[ ]` open roadmap item, `[~]` deferred
until its stated readiness condition holds. Every open item carries a
measurable **Goal** (acceptance tolerance) and a **Design** sketch with its
primary reference.

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Core types & traits (`Lattice`, `LatticeDimensions`, `LatticeType`,
      distribution arrays, velocity sets) in `lattice/`
- [x] Structured error type with domain-specific variants (`error.rs`)
- [x] Grid containers (`LbmGrid2D`, `LbmGrid3D`, D2Q9 specialised grid) in
      `grid/` and `lattice/`
- [x] Initialisation helpers (`initialization.rs`): equilibrium init, shear,
      Taylor–Green vortex, noise perturbation
- [x] Unit tests co-located with every module

### Phase 2: Lattices
- [x] D2Q9 (2-D nine-velocity) — `lattice/types.rs`
- [x] D3Q19 (standard isothermal 3-D) — `lattice/types.rs`, `d3q19_full.rs`
- [x] D3Q27 (higher moment isotropy) — `lattice/types.rs`, `d3q27.rs`

Note: D3Q15 is *not* implemented (kept off the list intentionally — see
Deferred / research track).

### Phase 3: Streaming (`streaming.rs`)
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

### Phase 4: Collision operators (`collision/types.rs`)
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

### Phase 5: Forcing schemes (`forcing.rs`)
- [x] Guo forcing (`GuoForcing`, `GuoForcingScheme`, `apply_guo_forcing_d3q19`)
- [x] Guo body force in `lattice/types.rs`
- [x] He–Luo forcing
- [x] Shan–Chen forcing scheme
- [x] Exact Difference Scheme (EDS)
- [x] Body-force helpers (`BodyForce`, `BodyForceType`)
- [x] Gravitational forcing, Boussinesq forcing, oscillating force,
      rotating-frame force, Coriolis force
- [x] Smagorinsky forcing, MHD Lorentz force

### Phase 6: Boundary conditions (`boundary.rs`, `curved_boundary.rs`, `zou_he/`)
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

### Phase 7: Thermal / heat transfer (`thermal/`, `thermal_lbm.rs`, `heat_transfer_lbm.rs`, `conjugate_heat.rs`)
- [x] `ThermalD2Q9`, `ThermalLbm` (double-distribution momentum+temperature)
- [x] Boussinesq coupling (`BoussinesqCoupling`)
- [x] Canonical setups: De Vahl Davis natural convection,
      Rayleigh–Bénard, natural-convection general setup
- [x] Nusselt-number computation
- [x] Conjugate heat transfer across solid/fluid interfaces
- [x] Standalone `heat_transfer_lbm.rs` application module

### Phase 8: Multiphase & phase-field (`multiphase/`, `phase_field/`, `cahn_hilliard*.rs`, `phase_separation.rs`, `droplet_dynamics*`)
- [x] Shan–Chen pseudo-potential (`ShanChenModel`, plus multi-component
      `MultiComponentSC`) with multiple `PsiType` pseudo-potential forms
- [x] Free-energy model (`FreeEnergyModel`)
- [x] Cahn–Hilliard (`cahn_hilliard.rs`, `cahn_hilliard_lbm.rs`)
- [x] Phase-field module (`phase_field/`, `phase_field_lbm.rs`)
- [x] Spinodal decomposition parameters and bubble-state tracking
- [x] Phase-separation standalone module
- [x] Droplet-dynamics modules (`droplet_dynamics.rs`,
      `droplet_dynamics_lbm/`)

### Phase 9: Turbulence (`turbulence/`, `turbulence_model/`, `turbulent_*`)
- [x] Static Smagorinsky (`SmagorinskyModel`)
- [x] Dynamic Smagorinsky (`DynamicSmagorinsky`)
- [x] k-ε model (`KEpsilonState`)
- [x] k-ω SST (`KOmegaSst`)
- [x] LES-to-DNS transition helper (`LesToDns`)
- [x] Turbulent Prandtl presets (`TurbPrandtlPreset`)
- [x] Turbulent-channel canonical setup (`turbulent_channel/`)
- [x] Turbulent dispersion (`turbulent_dispersion_lbm.rs`)

### Phase 10: Non-Newtonian rheology (`non_newtonian/`)
- [x] Power-law fluid (`power_law.rs`)
- [x] Bingham plastic (`bingham.rs`)
- [x] Herschel–Bulkley (`herschel_bulkley.rs`)
- [x] Carreau (`carreau.rs`)
- [x] Viscoelastic (`viscoelastic.rs`, `viscoelastic_lbm.rs`)
- [x] LBM integration layer and dedicated test suite

### Phase 11: Particle coupling
- [x] Lagrangian particle coupling (`particle_coupling.rs`)
- [x] Particle-laden flow (`particle_laden_lbm.rs`)
- [x] LBM particles helper (`lbm_particles.rs`)
- [x] Suspension flow (`suspension_lbm.rs`)
- [x] Sedimentation (`sedimentation/`, `sediment_transport.rs`)
- [x] Granular LBM (`granular_lbm.rs`)

### Phase 12: Porous media & multiscale
- [x] Porous-media LBM (`porous/`, `porous_media/`, `porous_media_lbm.rs`)
- [x] Multiscale LBM (`multiscale/`, `multiscale_lbm.rs`)

### Phase 13: Reactive flow & combustion
- [x] Species transport (`reactive/species.rs`, `concentration.rs`)
- [x] Catalytic reactions (`reactive/catalytic.rs`)
- [x] Combustion (`reactive/combustion.rs`, `combustion_lbm.rs`)
- [x] Reactive LBM integration (`reactive/lbm_reactive.rs`,
      `reactive_flow.rs`)
- [x] Extended reactive kernels and test suites

### Phase 14: Electromagnetic / charged-fluid LBM
- [x] Electrokinetic (`electrokinetic/`, `electrokinetic_lbm.rs`,
      `electrokinetics.rs`, `electrokinetics_lbm.rs`)
- [x] Electro-osmotic (`electroosmotic_lbm.rs`)
- [x] Electrostatic (`electrostatic_lbm.rs`)
- [x] Magnetohydrodynamics (`magnetohydrodynamics_lbm.rs`)
- [x] Ferrofluid (`ferrofluid_lbm.rs`)
- [x] Plasma (`plasma_lbm.rs`)

### Phase 15: Acoustics & compressible
- [x] Acoustic LBM (`acoustic_lbm.rs`, `acoustics_lbm.rs`)
- [x] Acoustic streaming (`acoustic_streaming_lbm.rs`)
- [x] Aeroacoustics (`aeroacoustics/`, `aeroacoustics_lbm.rs`)
      including A-weighting filter and near-to-far-field traits
- [x] Compressible LBM (`compressible.rs`)

### Phase 16: Bio- & soft-matter applications
- [x] Hemodynamics (`hemodynamics_lbm.rs`)
- [x] Biofluids (`biofluid_lbm/`, `biofluids_lbm.rs`) with
      red-blood-cell membrane traits
- [x] Biofilm (`biofilm_lbm.rs`)
- [x] Soft matter (`soft_matter_lbm.rs`)
- [x] Polymer (`polymer_lbm.rs`)
- [x] Microfluidics (`microfluidics.rs`, `microfluidics_lbm.rs`)
- [x] Mixing (`mixing_lbm.rs`)
- [x] Solidification (`solidification_lbm.rs`)

### Phase 17: Large-scale / geophysical / specialty
- [x] Geophysical LBM (`geophysical_lbm.rs`)
- [x] Climate LBM (`climate_lbm/` — energy balance, greenhouse gas, ENSO,
      Hadley-cell, monsoon, carbon-cycle, cloud-feedback,
      ice-albedo, ocean-heat-uptake, permafrost, sea-level-rise traits)
- [x] Quantum LBM (`quantum_lbm.rs`)
- [x] Neural LBM (`neural_lbm.rs`)
- [x] Traffic-flow LBM (`traffic_flow_lbm.rs`)

### Phase 18: Simulation driver & I/O
- [x] Simulation orchestration (`simulation/` with statistics traits)
- [x] Hybrid-LBM dispatcher (`hybrid_lbm.rs`)
- [x] LBM optimisation helpers (`lbm_optimization.rs`)
- [x] Diffusion module (`diffusion.rs`)
- [x] Snapshot / restore helpers (in `streaming.rs`) used by restart paths
- [x] Unit and integration tests co-located with each module

Note: no VTK/Paraview writer is implemented in this crate; visualisation
currently goes through consumer crates (see Deferred / research track).

### v0.1.2 correctness fixes (2026-06-01)
- [x] `LbmGrid3D::step(omega: f64)` — implemented real BGK collide-and-stream
      step (equilibrium collision + pull-scheme periodic streaming +
      `compute_macroscopic`). Removed dead empty `step_placeholder` (zero call
      sites). 3 new tests: mass conservation, equilibrium fixed-point,
      aggressive-omega.

### Code-audit confirmations (2026-06-11) — shipped; must not be re-opened
- [x] Multiscale rescale primitives: `coarse_to_fine_pop`,
      `fine_to_coarse_pop`, `coarse_to_fine_rescale` (`multiscale/`) —
      baseline for the nested 2:1 refinement driver
- [x] Spalding and Werner–Wengle wall models (`wall_model.rs`) — baseline for
      the WMLES coupling; the wall models themselves are not new work
- [x] FW-H far-field aeroacoustics (`aeroacoustics/`: `FwhPanel`,
      `NearToFarField`) — complete; deliberately not re-listed under v1.0
      validation
- [x] Hemodynamics (`hemodynamics_lbm.rs`) and Lagrangian particles
      (`lbm_particles.rs`) — baseline for the Womersley / sphere-drag cases
- [x] Lid-driven cavity Ghia Re=100 centerline validation + stored regression
      baseline (root TODO Phase 21.3) — only the Re=400/1000 extension
      remains open, moved to v1.0 below
- [x] AA-pattern in-place CPU streaming (`streaming.rs`) — algorithmic
      baseline for the GPU esoteric-twist port
- [x] GPU D3Q19 BGK step in `oxiphysics-gpu` — baseline for in-place GPU
      streaming
- [x] `is_interface_cell` helper — seed for the free-surface VOF cell flags

## v0.2.0 — Geometry & near-wall fidelity

Theme: sharp interfaces, local resolution, and memory-efficient domains —
the geometry capabilities that diffuse-interface multiphase and dense grids
cannot cover.

- [ ] Free-surface LBM (VOF mass tracking) — **Goal:** 2D dam break front
      within 8% of experiment; mass conserved to 1e-4. **Design:** per-cell
      fill level φ, interface reconstruction, mass-exchange +
      fluid/interface/gas flag reinit. Körner 2005; Thürey. *(Shipped:
      multiphase is diffuse-interface only (Shan–Chen / Cahn–Hilliard) and an
      `is_interface_cell` helper exists; open delta: the mass-tracking VOF
      machinery.)*
- [ ] Grid refinement (2:1 nested) finalization — **Goal:** refined cylinder
      wake matches uniform-fine Strouhal within 3% at ¼ the cells.
      **Design:** ν-consistent population rescaling + temporal interpolation
      at refinement interfaces, driven from existing rescale primitives.
      Lagrava 2012. *(Shipped: `multiscale/` `coarse_to_fine_pop` /
      `fine_to_coarse_pop` / `coarse_to_fine_rescale`; open delta: the
      nested-domain driver + temporal interpolation.)*
- [ ] Sparse/indirect addressing — **Goal:** vascular/porous domain with 10%
      fluid fraction uses ≤1.3× fluid-cell memory vs dense grid. **Design:**
      list-based fluid-cell indexing + neighbor index arrays.
      waLBerla/FluidX3D sparse lists. *(Code-verified missing 2026-06-11.)*

**Exit criteria:** all three acceptance targets pass as tolerance-checked
tests; no regression across the existing 5,320-test suite.

## v0.3.0 — GPU & advanced kinetics

Theme: device-resident throughput and kinetic schemes beyond the
unit-CFL on-lattice baseline.

- [ ] GPU in-place streaming (esoteric-twist/AA) — **Goal:** ≥500 MLUPS D3Q19
      on mid-GPU; bit-stable vs CPU TRT to 1e-6. **Design:**
      single-population in-place stream+collide on wgpu, no second buffer.
      Geier–Schönherr 2017; Bailey 2009. *(Shipped: GPU D3Q19 BGK step in
      `oxiphysics-gpu` and AA-pattern CPU streaming in `streaming.rs`; open
      delta: the GPU esoteric-twist kernel.)*
      **Depends on:** `oxiphysics-gpu` v0.2.0 device-side reduction/sort work
      for fully resident state and on-device diagnostics. wgpu is the default
      backend (pure Rust); cudarc remains feature-gated, non-default.
- [ ] Semi-Lagrangian LBM — **Goal:** stable at CFL>1 (≈4) on Taylor-Green
      with controlled diffusion. **Design:** off-lattice departure-point
      interpolation. Krämer et al. 2017. *(Code-verified missing
      2026-06-11.)*
- [ ] Wall-modeled LES coupling — **Goal:** turbulent channel Re_τ=180 mean
      u⁺ within 5% of DNS in the log layer. **Design:** wire existing
      Spalding/Werner–Wengle wall stress to Smagorinsky/k-ω-SST as a WMLES
      BC. Malaspinas–Sagaut 2014. *(Shipped: `wall_model.rs` (Spalding,
      Werner–Wengle) and the LES models; open delta: only the WMLES coupling
      between them.)*
- [ ] Non-uniform / stretched grids — **Goal:** near-wall y⁺<1 with 40% fewer
      nodes vs uniform at equal accuracy. **Design:** grid-stretching with
      interpolated streaming. *(Code-verified missing 2026-06-11.)*

**Exit criteria:** GPU streaming hits the MLUPS target with CPU-equivalence
checks; WMLES channel and semi-Lagrangian TGV meet their tolerances.

### Sequencing notes
- The 2:1 nested refinement driver (v0.2.0) is a prerequisite for combining
  stretched/refined grids with WMLES in wall-bounded cases.
- The GPU esoteric-twist kernel gates on `oxiphysics-gpu` v0.2.0
  reduction/sort groundwork — track that crate's TODO before starting.

## v1.0 — Production & validation

Validation ladder with hard tolerances. Each case lands as a reproducible
example plus a tolerance-checked regression test with a stored baseline,
extending the pattern established by the Ghia Re=100 cavity baseline.

- [ ] Lid-driven cavity extension to Re=400/1000 — **Goal:** centerline u/v
      within 5% of Ghia 1982 at Re=1000. **Design:** extend existing Re=100
      baseline with MRT. *(Merged from the former "Outstanding (v0.2.0)" root
      Phase 21.3 pointer; Re=100 validation + regression baseline already
      shipped — only Re=400/1000 remain open.)*
- [ ] Taylor-Green decay order verification — **Goal:** energy-decay
      convergence slope ≥1.9; decay rate within 2% of analytic. **Design:**
      periodic TGV refinement study.
- [ ] Womersley pulsatile + sphere-drag DNS — **Goal:** Womersley velocity
      within 3% of Bessel solution; sphere C_D within 8% of Schiller–Naumann
      at Re=50–200. **Design:** oscillatory pressure-driven pipe +
      particle-resolved sphere. *(Builds on shipped `hemodynamics_lbm` +
      `lbm_particles`; FW-H far-field already present via `aeroacoustics/`
      `FwhPanel`/`NearToFarField` — not new work.)*

## Deferred / research track

- [~] D3Q15 lattice — Ready when: a concrete consumer needs the
      lower-isotropy stencil; D3Q19/D3Q27 cover all current use (intentional
      omission since v0.1.x).
- [~] In-crate VTK/ParaView writer — Ready when: a shared writer lands in
      `oxiphysics-io`; visualisation deliberately stays in consumer crates
      until then.
- [~] cudarc (CUDA) parity for the GPU in-place streaming kernel — Ready
      when: the wgpu esoteric-twist path (v0.3.0) has landed and
      `oxiphysics-gpu`'s cudarc feature has stabilised; stays feature-gated
      and non-default per the Pure Rust default policy.
