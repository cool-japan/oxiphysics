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

- [x] Free-surface LBM (VOF mass tracking) — **Goal:** 2D dam break front
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
- [x] Wall-modeled LES coupling (done 2026-06-14) — **Goal:** turbulent channel Re_τ=180 mean
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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 330 (critical 4, high 124, medium 151, low 51); 48 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (4)
- [x] **F2480** (critical, bug, easy/small) `d3q19_full.rs:49` — D3Q19_OPP maps edge directions to non-opposites
- [ ] **F2601** (critical, fabrication, hard/medium) `magnetohydrodynamics_lbm.rs:1064` — MhdLbm shrinks B threefold per step and has no Lorentz coupling
- [ ] **F2612** (critical, fabrication, hard/medium) `mrt.rs:1339` — CascadedCollision2D::collide reverses momentum; fake back-transform
- [ ] **F2785** (critical, bug, easy/small) `wall_model.rs:387` — u_tau_newton Jacobian wrong: u_tau collapses to 1e-15 in log layer

### High (124)
- [ ] **F2434** (high, fabrication, medium/medium) `biofluid_lbm/types.rs:264` — MicrofluidicsLbm::step never streams; ln(psi) recomputed 16x per cell
- [x] **F2435** (high, fabrication, easy/small) `biofluid_lbm/types.rs:291` — MicrofluidicsLbm::step never streams populations
- [ ] **F2440** (high, bug, medium/small) `boundary_lbm.rs:229` — InterpolatedBc uses wrong Bouzidi coefficients and wrong samples
- [ ] **F2448** (high, bug, easy/small) `cahn_hilliard_lbm.rs:668` — Variable mobility applied as M(C)∇²μ; Degenerate mobility clamped wrong
- [ ] **F2451** (high, bug, easy/small) `collision/functions.rs:959` — Regularized f1 lacks −cs²δ Hermite term; mass not conserved
- [ ] **F2453** (high, bug, medium/medium) `collision/functions.rs:1043` — RawMomentCollision relaxes normal stresses at bulk rate; ghosts discarded
- [ ] **F2455** (high, bug, medium/small) `collision/functions.rs:1213` — entropic_alpha_d3q19 bisection inverted; ELBM degenerates to feq or ω=2
- [x] **F2456** (high, fabrication, hard/medium) `collision/types.rs:182` — CumulantCollision is not cumulant; does not conserve mass
- [ ] **F2457** (high, fabrication, hard/medium) `collision/types.rs:318` — KbcCollision is not KBC: gamma is either omega_s or 2
- [x] **F2460** (high, fabrication, medium/medium) `collision/types.rs:714` — CentralMomentCollision maps moment k onto population k
- [ ] **F2462** (high, fabrication, medium/large) `combustion_lbm.rs:946` — CombustionLbm never advects species or temperature with the flow
- [ ] **F2464** (high, bug, medium/small) `compressible.rs:67` — Compressible equilibrium cubic term uses cs⁴; cs²=γ/3 breaks moments
- [ ] **F2465** (high, fabrication, hard/medium) `compressible.rs:92` — CompressibleLbm cs²=γ/3 on D2Q9 breaks mass/momentum conservation
- [ ] **F2467** (high, fabrication, medium/small) `compressible.rs:268` — EulerLbm claims inviscid ω=2 but sets f=feq (ω=1)
- [ ] **F2468** (high, bug, medium/medium) `compressible.rs:553` — NsCompressibleLbm temperature shrinks by factor cs² per relaxation
- [ ] **F2469** (high, fabrication, hard/large) `compressible.rs:832` — RayleighTaylorLbm::step only bumps u_y; rho/phase never evolve
- [ ] **F2471** (high, bug, easy/trivial) `compressible.rs:1076` — cj_pressure_ratio always returns exactly 1.0
- [ ] **F2472** (high, fabrication, medium/medium) `compressible.rs:1108` — DetonationLbm::step density update is identically zero
- [ ] **F2474** (high, fabrication, medium/medium) `conjugate_heat.rs:812` — Conjugate solver: interface flux always zero; LBM ignores solid mask
- [ ] **F2475** (high, fabrication, medium/medium) `curved_boundary.rs:129` — bouzidi_bc/interpolated_bounce_back use arbitrary array slots, not lattice links
- [ ] **F2476** (high, fabrication, medium/medium) `curved_boundary.rs:136` — bouzidi_bc/interpolated_bounce_back/moving_wall_bc use wrong indices
- [ ] **F2478** (high, fabrication, hard/large) `curved_boundary.rs:331` — CurvedBoundarySet::apply is a no-op/fabrication for all methods
- [ ] **F2479** (high, bug, easy/small) `curved_boundary.rs:521` — RotatingCylinder::boundary_set uses inverted SDF and zero wall velocity
- [ ] **F2482** (high, bug, medium/small) `d3q19_full.rs:284` — Guo forcing injects only (1-ω/2)F momentum per step
- [x] **F2483** (high, fabrication, easy/small) `d3q19_full.rs:327` — collide_mrt is per-direction BGK presented as MRT
- [x] **F2484** (high, bug, easy/small) `d3q19_full.rs:818` — TRT collision updates populations in place, reading stale opposites
- [x] **F2485** (high, bug, medium/small) `d3q19_full.rs:928` — D3q19Simulation Zou-He density/velocity formulas wrong, outlet forced inflow
- [ ] **F2490** (high, bug, easy/trivial) `d3q27.rs:389` — D3Q27 collide_mrt (TRT) reads already-overwritten populations in place
- [x] **F2491** (high, fabrication, medium/medium) `d3q27.rs:617` — D3Q27 collide_mrt_full is BGK with s[9]; matrix ignored
- [ ] **F2502** (high, bug, easy/trivial) `electrokinetic/types.rs:53` — Henry function formula wrong: 1.83–3.5 instead of 1–1.5
- [ ] **F2504** (high, bug, easy/trivial) `electrokinetic/types.rs:350` — PoissonSolver::solve_poisson_2d source sign inverted
- [ ] **F2509** (high, fabrication, easy/small) `electrokinetic_lbm.rs:436` — collide_stream discards applied E_y via multiplication by 0.0
- [ ] **F2510** (high, bug, easy/trivial) `electrokinetics.rs:349` — Nernst-Planck electromigration term has wrong sign
- [ ] **F2515** (high, fabrication, hard/large) `electrokinetics_lbm.rs:246` — ElectrokineticLBM never solves Poisson; tau from Pa·s viscosity
- [ ] **F2519** (high, bug, easy/trivial) `electroosmotic_lbm.rs:473` — Nernst-Planck migration velocity omits diffusivity, double-counts z
- [ ] **F2522** (high, bug, medium/medium) `electrostatic_lbm.rs:132` — PoissonLbm solves wrong equation and drops Dirichlet values
- [ ] **F2525** (high, stub, hard/large) `electrostatic_lbm.rs:517` — SemiconductorLbm::step only steps Poisson; carriers never evolve
- [ ] **F2526** (high, bug, medium/medium) `electrostatic_lbm.rs:576` — SpaceChargeLbm f never initialised; first step zeroes concentrations
- [ ] **F2528** (high, bug, medium/medium) `electrostatic_lbm.rs:826` — PlasmaLbm drives electrons and ions same direction; Ey dropped
- [ ] **F2529** (high, fabrication, medium/small) `entropic.rs:162` — find_entropic_alpha always returns 2; root search is dead code
- [ ] **F2530** (high, fabrication, medium/small) `entropic.rs:364` — Entropic/KBC/EOR collisions ignore nu; always inviscid omega=2
- [ ] **F2531** (high, fabrication, hard/large) `entropic.rs:664` — KbcD2Q9 ignores nu; alpha is almost always 2, giving an inviscid mirror
- [ ] **F2537** (high, bug, easy/trivial) `forcing.rs:597` — Boussinesq buoyancy sign wrong; test enshrines downward force
- [ ] **F2540** (high, bug, medium/small) `free_surface.rs:427` — Gas-side reconstruction overwrites the wrong population direction
- [ ] **F2541** (high, bug, medium/small) `free_surface.rs:501` — VOF mass exchange applies pre-stream formula to post-stream data
- [ ] **F2545** (high, fabrication, medium/medium) `geophysical_lbm.rs:335` — ShallowWaterLbm ignores gravity; uses Navier-Stokes equilibrium
- [ ] **F2547** (high, fabrication, medium/medium) `geophysical_lbm.rs:607` — Mantle/thermohaline scalars use own moment, never advected by flow
- [ ] **F2548** (high, fabrication, medium/small) `geophysical_lbm.rs:913` — ABL geostrophic forcing adds mass not momentum; wall stress discarded
- [ ] **F2569** (high, bug, medium/small) `hemodynamics_lbm.rs:117` — Rheology tau uses SI dynamic viscosity as lattice viscosity
- [ ] **F2571** (high, bug, easy/small) `hemodynamics_lbm.rs:745` — Push-stream from Wall cells overwrites bounce-back populations
- [ ] **F2581** (high, bug, easy/small) `immersed_boundary_lbm.rs:26` — delta_phi/delta_peskin4 violate partition of unity and are discontinuous
- [ ] **F2582** (high, bug, easy/small) `immersed_boundary_lbm.rs:300` — IBM velocity interpolation wrongly weighted by marker arc length ds
- [ ] **F2583** (high, bug, medium/small) `immersed_boundary_lbm.rs:765` — FlexibleIbmFiber bending force is anti-restoring; energy formula wrong
- [ ] **F2591** (high, bug, easy/trivial) `lattice/types.rs:239` — TrtCollision::collide reads already-overwritten opposite populations
- [ ] **F2593** (high, bug, easy/trivial) `lbm_optimization.rs:100` — optimal_relaxation uses Λ=τ⁺τ⁻, returns unstable omega_minus
- [ ] **F2598** (high, bug, easy/small) `magnetohydrodynamics_lbm.rs:202` — LorentzForce::compute has inverted Jz and drops in-plane Jx,Jy
- [ ] **F2602** (high, fabrication, medium/medium) `microfluidics_lbm.rs:425` — ElectroosmosisLBM::step never streams; EO profile cannot form
- [ ] **F2606** (high, fabrication, medium/large) `mixing_lbm/types.rs:581` — MixingLBM never streams, wrong weights, scalar mass lost
- [ ] **F2607** (high, fabrication, medium/medium) `mixing_lbm/types.rs:590` — Mixing LBMs: no streaming, wrong D2Q9 weights, scalar mass lost
- [ ] **F2608** (high, bug, easy/small) `mrt.rs:284` — D3Q19 mrt_transform_matrix is singular; public inverse returns garbage
- [ ] **F2609** (high, fabrication, easy/small) `mrt.rs:613` — Root-exported MrtD3Q19 collide_cell is plain BGK
- [ ] **F2611** (high, bug, easy/trivial) `mrt.rs:1139` — TrtCollision2D magic parameter Λ=τ⁺τ⁻ gives unstable s⁻>2
- [ ] **F2625** (high, bug, easy/trivial) `multiphase/types.rs:289` — Rayleigh-Plesset acceleration divides 1.5Ṙ² inertia term by density
- [ ] **F2626** (high, bug, easy/small) `multiphase/types.rs:403` — MultiComponentSC A–A self-interaction uses psi_B; g_bb never used
- [ ] **F2628** (high, fabrication, medium/small) `multiphase_lbm.rs:277` — phase_diagram claims Maxwell construction but returns spinodal heuristics
- [ ] **F2630** (high, bug, easy/trivial) `multiscale/types.rs:311` — TwoWayCoupling subtracts Guo force and omits −u·F term
- [ ] **F2632** (high, bug, medium/small) `multiscale/types.rs:610` — Grid-refinement tau_fine and non-equilibrium rescaling are inverted
- [ ] **F2636** (high, bug, easy/trivial) `multiscale/types.rs:1350` — PatchedGrid::restrict sums r² fine cells instead of averaging
- [ ] **F2639** (high, fabrication, medium/medium) `multiscale_lbm.rs:312` — AdaptiveLbm::step relaxes toward zero equilibrium; fine grids never advance
- [ ] **F2640** (high, fabrication, medium/large) `multiscale_lbm.rs:420` — Multiscale 'LBM' types relax toward zero; fake macro projection
- [ ] **F2643** (high, fabrication, hard/medium) `neural_lbm.rs:280` — PhysicsInformedLbm: feq=ρ/q, no streaming, untrained net, fake residual
- [ ] **F2645** (high, fabrication, hard/medium) `neural_lbm.rs:496` — DMD is fake: normalised snapshots as modes, Rayleigh quotients
- [ ] **F2646** (high, bug, easy/small) `non_newtonian/bingham.rs:50` — Bingham yield test compares viscous stress only, freezing yielded cells
- [ ] **F2647** (high, fabrication, medium/small) `non_newtonian/lbm_integration.rs:779` — apply_non_newtonian_collision uses |u| as the shear rate
- [ ] **F2648** (high, bug, medium/medium) `non_newtonian/lbm_integration.rs:824` — Non-Newtonian LBM shear rate is unscaled Pi^neq or velocity magnitude
- [ ] **F2657** (high, bug, easy/small) `particle_coupling.rs:215` — Two-way coupling returns particle gravity to the fluid
- [ ] **F2660** (high, bug, easy/trivial) `particle_laden_lbm.rs:294` — Schiller-Naumann drag is missing one factor of slip speed
- [ ] **F2662** (high, bug, easy/small) `particle_laden_lbm.rs:385` — buoyancy_force sign inverted: heavy particles accelerate upward
- [ ] **F2663** (high, bug, easy/small) `particle_laden_lbm.rs:455` — Soft-sphere DEM damping and friction have inverted signs
- [ ] **F2664** (high, bug, medium/medium) `particle_laden_lbm.rs:467` — DEM contact: reversed damping and friction signs, overlap-based tangential force
- [ ] **F2673** (high, fabrication, hard/medium) `phase_field_lbm.rs:254` — CahnHilliardLbm has no chemical-potential diffusion or surface tension
- [ ] **F2684** (high, fabrication, hard/medium) `plasma_lbm.rs:420` — PlasmaLbm::step fakes the Poisson E-field and skips the continuity update
- [ ] **F2685** (high, bug, medium/small) `polymer_lbm.rs:348` — Oldroyd-B/FENE-P stretching uses transposed velocity gradient; solver passes zeros
- [ ] **F2687** (high, fabrication, medium/medium) `polymer_lbm.rs:892` — PolymerLbmGrid2D evolves conformation with zero velocity gradients
- [ ] **F2696** (high, bug, medium/small) `porous_media/types.rs:714` — PorousMediaDriver Guo forcing lacks the −u term; Forchheimer force ignored
- [ ] **F2701** (high, bug, easy/trivial) `porous_media_lbm.rs:568` — PorousLBM Guo forcing term does not sum to zero; creates mass
- [ ] **F2703** (high, fabrication, hard/medium) `quantum_lbm.rs:317` — GrossPitaevskiiLbm::step does not solve the Gross-Pitaevskii equation
- [ ] **F2704** (high, fabrication, hard/medium) `quantum_lbm.rs:350` — GrossPitaevskiiLbm: inconsistent lattice and fabricated phase update
- [ ] **F2705** (high, bug, medium/medium) `quantum_lbm.rs:598` — SuperfluidLbm mutual friction has wrong sign; slip grows exponentially
- [ ] **F2712** (high, bug, easy/small) `reactive/lbm_reactive.rs:120` — Reaction rate multiplies all species concentrations, including products
- [ ] **F2714** (high, bug, easy/small) `reactive_flow.rs:269` — step_chemistry dimensionally wrong: no 1/(rho·cp), ignores fuel_stoich
- [ ] **F2722** (high, bug, medium/medium) `sedimentation/types.rs:867` — LbmDemSolver feeds SI density/gravity into lattice; one-way coupling
- [ ] **F2724** (high, fabrication, easy/small) `sedimentation/types.rs:1355` — Sedimentation3D::step not a flux divergence; phi grows everywhere
- [ ] **F2725** (high, bug, easy/trivial) `sedimentation/types.rs:1627` — SediLbm::mixture_omega lowers viscosity as particle loading rises
- [ ] **F2727** (high, bug, medium/medium) `simulation/types.rs:225` — LbmSimulation::step ignores dt; 3D drops body_force/Smagorinsky
- [ ] **F2728** (high, bug, easy/small) `simulation/types.rs:275` — 2D Guo forcing omits -u/cs² term, injects mass every step
- [ ] **F2729** (high, missing-feature, medium/medium) `simulation/types.rs:290` — 3D LbmSimulation silently ignores body force, Smagorinsky, boundaries
- [ ] **F2731** (high, fabrication, easy/small) `simulation/types.rs:641` — EnsembleRunner perturbation overwritten; all members identical
- [ ] **F2737** (high, bug, easy/trivial) `solidification_lbm.rs:341` — PhaseChangeLbm::apply_mould_wall bounce-back is a no-op double swap
- [x] **F2741** (high, perf, medium/medium) `streaming.rs:37` — stream_2d/stream_3d allocate full q×n buffer every step
- [ ] **F2744** (high, bug, medium/small) `streaming.rs:326` — stream_swap_2d streams wrong populations; AA pattern semantics wrong
- [x] **F2745** (high, bug, medium/small) `streaming.rs:336` — stream_swap_2d streams half the populations the wrong way
- [ ] **F2746** (high, bug, easy/trivial) `streaming.rs:672` — Half-way moving-wall bounce-back momentum term has wrong sign
- [x] **F2749** (high, bug, easy/trivial) `thermal/functions.rs:81` — D3Q7 thermal equilibrium advects scalar at half the fluid speed
- [ ] **F2751** (high, bug, easy/trivial) `thermal/types.rs:663` — BoussinesqCoupling buoyancy sign makes hot fluid sink
- [ ] **F2755** (high, bug, easy/small) `traffic_flow_lbm.rs:1051` — TrafficLbm feeds physical speed (m/s) into lattice equilibrium
- [ ] **F2756** (high, bug, easy/small) `turbulence/functions.rs:48` — Smagorinsky tau_eff divides by cs^4 twice, ~6x eddy viscosity
- [x] **F2764** (high, bug, easy/small) `turbulence_model/functions.rs:775` — eigenvalues_3x3_symmetric wrong cubic constant; Sigma SGS ν_t ≡ 0
- [ ] **F2769** (high, bug, medium/small) `turbulence_model/types.rs:576` — Dynamic Smagorinsky M sign inconsistent; Cs^2 clamps to zero
- [x] **F2771** (high, fabrication, easy/small) `turbulence_model/types.rs:1002` — AmdModel::eddy_viscosity always returns zero
- [x] **F2776** (high, fabrication, medium/medium) `turbulent_channel/types.rs:933` — WallFunctionLbm wall function overwrites only derived velocity (no-op)
- [x] **F2778** (high, bug, easy/trivial) `turbulent_dispersion_lbm.rs:264` — ScalarTransport bounce-back uses wrong opposite index
- [x] **F2779** (high, bug, easy/small) `turbulent_dispersion_lbm.rs:737` — DispersionLBM::advance accumulates diffusivity and nullifies collision
- [x] **F2782** (high, fabrication, hard/medium) `viscoelastic_lbm.rs:380` — ViscoelasticLbm uses |u| as shear rate; stress term is a mass source
- [ ] **F2787** (high, fabrication, easy/small) `wall_model.rs:901` — Werner-Wengle wall stress formula fabricated with wrong dimensions
- [ ] **F2789** (high, bug, medium/small) `zou_he/functions.rs:88` — D3Q19 Zou-He omits edge populations in transverse correction
- [ ] **F2790** (high, bug, medium/medium) `zou_he/functions.rs:106` — D3Q19 Zou-He faces use D2Q9 coefficients and omit edge terms
- [ ] **F2791** (high, bug, medium/medium) `zou_he/functions.rs:291` — D3Q19 Zou-He faces pair mirror populations instead of opposites
- [ ] **F2793** (high, bug, easy/small) `zou_he/functions.rs:615` — zou_he_pressure_left/right_d3q19 triple-count populations, mirror pairs
- [ ] **F2794** (high, bug, easy/small) `zou_he/functions.rs:630` — D3Q19 pressure_left/right triple-count density and inject 4/3 momentum
- [ ] **F2795** (high, bug, medium/small) `zou_he/functions.rs:768` — D3Q27 Zou-He inlet/outlet inject zero momentum, wrong density
- [ ] **F2796** (high, bug, medium/small) `zou_he/functions.rs:775` — D3Q27 Zou-He reduces to plain bounce-back with wrong density
- [ ] **F2798** (high, bug, easy/small) `zou_he/functions_2.rs:39` — D3Q19 moving-wall bounce-back quadratic in u_wall, wrong pairs
- [ ] **F2799** (high, bug, easy/small) `zou_he/functions_2.rs:46` — Moving-wall bounce-back squares wall velocity; wrong pair (16,17)

### Medium (151)
- [ ] **F2422** (medium, missing-feature, medium/medium) `crates/oxiphysics-lbm/Cargo.toml:13` — No parallel or serde features for the LBM solver
- [ ] **F2424** (medium, bug, easy/trivial) `acoustic_lbm.rs:45` — Linearized acoustic equilibrium loses momentum when rho0≠1
- [ ] **F2426** (medium, bug, easy/small) `acoustic_lbm.rs:274` — PML damps full populations, draining background density
- [ ] **F2427** (medium, pure-rust, easy/small) `acoustic_lbm.rs:578` — FrequencyResponse::from_signal is an O(N^2) hand-rolled DFT
- [ ] **F2428** (medium, fabrication, easy/small) `acoustic_lbm.rs:780` — FW-H integrate_surface hardcodes n·r̂ = 1
- [ ] **F2430** (medium, fabrication, easy/trivial) `acoustics_lbm.rs:725` — guo_forcing_streaming computes the Guo term then discards it
- [ ] **F2432** (medium, fabrication, medium/medium) `aeroacoustics/functions.rs:169` — Far-field FW-H uses only steady 1/r² near-field term
- [ ] **F2433** (medium, missing-feature, medium/medium) `aeroacoustics_lbm.rs:100` — AeroacousticLBM has no time stepping; record_pressure panics on mismatch
- [ ] **F2436** (medium, fabrication, medium/small) `biofluids_lbm.rs:228` — BloodFlowLbm::step is not LBM; wrong non-Newtonian profile
- [ ] **F2437** (medium, fabrication, medium/medium) `biofluids_lbm.rs:304` — Womersley profile fabricated, discontinuous, wrong units
- [ ] **F2439** (medium, bug, medium/small) `boundary.rs:596` — Interpolated bounce-back drops q<0.5 Bouzidi branch; unstable extrapolation
- [ ] **F2441** (medium, fabrication, easy/small) `boundary_lbm.rs:355` — apply_boundary silently turns 4 BC types into equilibrium overwrite
- [ ] **F2442** (medium, api, easy/small) `cahn_hilliard.rs:119` — Unseeded thread RNG in public init APIs blocks reproducibility
- [ ] **F2446** (medium, pure-rust, easy/small) `cahn_hilliard_lbm.rs:273` — structure_factor_1d computes 2D DFT by brute force O((nx*ny)^2)
- [ ] **F2449** (medium, bug, easy/small) `climate_lbm/types.rs:227` — °C temperature used as LBM density; velocity blows up near 0 °C
- [ ] **F2450** (medium, bug, easy/trivial) `climate_lbm/types.rs:243` — Radiative forcing applied to temp only; scaled by ω
- [ ] **F2452** (medium, fabrication, hard/medium) `collision/functions.rs:966` — compute_q1_recursive is a fabricated recursive-regularization term
- [ ] **F2458** (medium, fabrication, medium/small) `collision/types.rs:413` — Entropic collision never solves for mirror-state alpha; bisection inverted
- [ ] **F2459** (medium, fabrication, easy/small) `collision/types.rs:508` — MrtCollision ignores 15 of 19 rates; it is two-rate regularized BGK
- [ ] **F2466** (medium, perf, easy/medium) `compressible.rs:143` — Compressible solvers allocate nx*ny heap Vecs per stream step
- [ ] **F2470** (medium, fabrication, hard/large) `compressible.rs:840` — RayleighTaylorLbm::step never transports density or phase
- [ ] **F2477** (medium, bug, easy/small) `curved_boundary.rs:238` — moving_wall_bc uses ux+uy as e_i·u and len−1−dir opposite
- [ ] **F2481** (medium, perf, easy/small) `d3q19_full.rs:171` — D3Q19/D3Q27 streaming clones whole lattice and uses rem_euclid
- [ ] **F2487** (medium, fabrication, easy/small) `d3q19_full.rs:1061` — setup_poiseuille_channel discards its body-force argument
- [ ] **F2488** (medium, bug, medium/small) `d3q19_full.rs:1170` — total_force_on_solid sums all solid populations, not link momentum exchange
- [ ] **F2495** (medium, bug, easy/small) `droplet_dynamics.rs:343` — Taylor deformation uses Weber number instead of capillary number
- [ ] **F2496** (medium, stub, medium/medium) `droplet_dynamics.rs:449` — Phase-field DropletLBM has no dynamics; f_dist never initialised
- [ ] **F2497** (medium, stub, hard/large) `droplet_dynamics.rs:455` — DropletLBM has populations but no collide/stream; f_dist stays zero
- [ ] **F2500** (medium, bug, easy/small) `droplet_dynamics_lbm/types.rs:812` — Saffman lift is dimensionally wrong (extra |ω|) and has reversed direction
- [ ] **F2503** (medium, bug, easy/small) `electrokinetic/types.rs:211` — Surface charge from zeta has the wrong sign in two APIs
- [ ] **F2505** (medium, bug, easy/small) `electrokinetic/types.rs:673` — Transient EOF start-up decays 4x too fast
- [ ] **F2511** (medium, bug, easy/trivial) `electrokinetics.rs:466` — Double-layer charge density has the wrong sign; integrated charge is hardcoded
- [ ] **F2514** (medium, bug, medium/medium) `electrokinetics_lbm.rs:227` — tau built from dynamic viscosity in Pa·s; Poisson never solved
- [ ] **F2516** (medium, bug, easy/small) `electrokinetics_lbm.rs:285` — Code labelled 'Guo forcing' lacks the (1−1/2τ) factor and half-force velocity
- [ ] **F2520** (medium, fabrication, hard/large) `electroosmotic_lbm.rs:646` — ElectroosmoticFlow claims coupling but only advects ions
- [ ] **F2521** (medium, bug, easy/trivial) `electroosmotic_lbm.rs:707` — ElectrophoresisModel mobility missing the 2/3 Henry prefactor
- [ ] **F2524** (medium, bug, medium/medium) `electrostatic_lbm.rs:399` — Carriers lack a charge sign; Einstein relation claimed but not applied
- [ ] **F2527** (medium, bug, easy/small) `electrostatic_lbm.rs:619` — SpaceChargeLbm drift divides mobility by V_T, has no y drift, and asserts in new()
- [ ] **F2532** (medium, api, medium/large) `error.rs:21` — Error hierarchy unused: no fallible public API in crate
- [ ] **F2534** (medium, fabrication, medium/medium) `ferrofluid_lbm.rs:416` — Langevin alpha fabricated; magnetoviscosity never reaches the collision
- [ ] **F2539** (medium, bug, easy/trivial) `forcing.rs:1253` — CoriolisForce centrifugal term points inward (sign error)
- [ ] **F2542** (medium, bug, medium/small) `free_surface.rs:700` — Interface cells promoted from Gas keep stale distributions
- [ ] **F2544** (medium, bug, easy/trivial) `geophysical_lbm.rs:261` — enstrophy computes ∂v/∂y − ∂u/∂x instead of vorticity
- [ ] **F2546** (medium, fabrication, medium/small) `geophysical_lbm.rs:510` — Rossby 'streamfunction' stores wrong-component pseudo-vorticity
- [ ] **F2549** (medium, fabrication, medium/small) `geophysical_lbm.rs:923` — ABL log-law wall stress discarded; geostrophic nudging adds mass
- [ ] **F2550** (medium, fabrication, hard/medium) `geophysical_lbm.rs:1359` — Two-layer PlanetaryWave layers never coupled; g_prime unused
- [ ] **F2551** (medium, fabrication, medium/medium) `geophysical_lbm.rs:1416` — PlanetaryWave layers uncoupled and L_R wrong; tidal depth unused
- [ ] **F2552** (medium, bug, easy/trivial) `granular_lbm.rs:224` — mu(I) effective viscosity divides by shear rate squared (wrong units)
- [ ] **F2553** (medium, bug, medium/small) `granular_lbm.rs:327` — Granular shear rate ignores normal strain; global pressure; SI nu used as lattice nu
- [ ] **F2558** (medium, perf, medium/medium) `grid/grid_extended.rs:335` — Grid streaming clones the full population array every step
- [ ] **F2559** (medium, bug, easy/small) `grid/grid_extended.rs:841` — 3-D outlet zero-gradient copies from the y−1 neighbour, not x−1
- [ ] **F2560** (medium, bug, easy/small) `grid/types.rs:42` — CellGrid3D obstacle cells are transparent (no bounce-back)
- [ ] **F2563** (medium, stub, medium/medium) `grid/types.rs:419` — NodeFlag::Symmetry treated as plain fluid (no-op)
- [ ] **F2565** (medium, fabrication, medium/small) `heat_transfer_lbm.rs:201` — Rayleigh-Benard init promises perturbation never applied; no walls
- [ ] **F2566** (medium, fabrication, medium/medium) `heat_transfer_lbm.rs:214` — Rayleigh–Bénard setup has no walls, no perturbation, and an unused Ra
- [ ] **F2572** (medium, bug, easy/small) `hemodynamics_lbm.rs:861` — WSS returned unsigned, so compute_osi is always zero
- [ ] **F2574** (medium, bug, easy/small) `hybrid_lbm.rs:232` — mass_flux_interface projects velocity with (nx + ny)
- [ ] **F2575** (medium, bug, easy/small) `immersed_boundary/functions.rs:91` — IB interpolation scaled by dx^d; 2D spreading missing 1/dx²
- [ ] **F2576** (medium, fabrication, easy/small) `immersed_boundary/functions.rs:314` — ib_forcing_term hardcodes τ=1 prefactor and ignores the force
- [ ] **F2577** (medium, perf, easy/small) `immersed_boundary/functions.rs:340` — 3D IBM loops 5³ stencil with 375 cos() per marker
- [ ] **F2578** (medium, bug, easy/small) `immersed_boundary/functions.rs:398` — IBM interpolation multiplies by dx^d; 2D spreading lacks 1/dx²
- [ ] **F2584** (medium, bug, easy/small) `initialization.rs:348` — Taylor-Green decay rate doubled; init not divergence-free for nx≠ny
- [ ] **F2585** (medium, panic, easy/small) `initialization.rs:600` — Checkpoint::from_bytes overflow lets crafted header panic
- [ ] **F2586** (medium, panic, easy/small) `initialization.rs:618` — Checkpoint::from_bytes panics on crafted header (overflow, huge q)
- [ ] **F2587** (medium, fabrication, easy/small) `lattice/functions.rs:862` — equilibrium_from_hermite discards second-order coefficient a2
- [ ] **F2588** (medium, bug, easy/trivial) `lattice/functions.rs:941` — strain_rate_magnitude_d2q9 computes Σf_neq²|c|², not Π:Π
- [ ] **F2589** (medium, perf, easy/small) `lattice/streaming.rs:16` — Periodic streaming allocates a full lattice copy every step
- [ ] **F2590** (medium, bug, easy/trivial) `lattice/types.rs:33` — GuoBodyForce forcing term omits density (acceleration used as force)
- [ ] **F2592** (medium, perf, easy/medium) `lattice/types.rs:490` — LBM collision/stream loops all serial; crate has no rayon feature
- [ ] **F2594** (medium, bug, easy/trivial) `lbm_particles.rs:201` — Particle mass uses 2D disc area but drag is 3D Stokes
- [ ] **F2599** (medium, missing-feature, medium/medium) `magnetohydrodynamics_lbm.rs:331` — InductionEquation advertises ∇×(u×B) advection but only diffuses
- [ ] **F2600** (medium, fabrication, easy/small) `magnetohydrodynamics_lbm.rs:455` — slow_magnetosonic_speed uses invented |vA−cs|/√2 formula
- [ ] **F2604** (medium, bug, easy/trivial) `microfluidics_lbm.rs:1161` — MicrofluidicLBMGrid Guo term missing (c·u)(c·F)/cs⁴; stress_xy factor wrong
- [ ] **F2605** (medium, bug, easy/trivial) `mixing_lbm/types.rs:405` — Taylor-Aris channel 1/210 coefficient applied to half-width (4x low)
- [ ] **F2610** (medium, perf, easy/small) `mrt.rs:639` — Solver stream() clones the entire distribution array each step
- [ ] **F2614** (medium, bug, easy/small) `mrt3d.rs:583` — MRT D3Q19 stress reconstruction divides 3pxx moment by 6
- [ ] **F2615** (medium, perf, easy/trivial) `mrt3d.rs:685` — TRT collide_grid heap-allocates two Vecs per cell per step
- [ ] **F2617** (medium, bug, easy/trivial) `mrt3d.rs:848` — mrt_equilibrium_moments ε_eq has wrong sign and coefficient
- [ ] **F2618** (medium, bug, easy/trivial) `multiphase/functions.rs:276` — spinodal_growth_rate uses 2κ, inconsistent with module CH convention
- [ ] **F2619** (medium, bug, easy/trivial) `multiphase/functions.rs:408` — 2D Young-Laplace pressure returns 2σ/R instead of σ/R
- [ ] **F2620** (medium, bug, easy/trivial) `multiphase/functions.rs:450` — Hadamard–Rybczynski terminal velocity factor wrong (3x Stokes)
- [ ] **F2621** (medium, bug, easy/trivial) `multiphase/functions.rs:466` — Hadamard-Rybczynski drag has extra /3; terminal velocity factor also wrong
- [ ] **F2622** (medium, fabrication, easy/small) `multiphase/functions.rs:555` — Hertz-Knudsen: p_sat_ref=1 hard-coded, no latent heat; Leidenfrost dimensionless nonsense
- [ ] **F2623** (medium, bug, easy/small) `multiphase/types.rs:51` — Free-energy pressure uses B/4 instead of 3B/4; gradient-term signs disagree
- [ ] **F2624** (medium, api, easy/small) `multiphase/types.rs:184` — Public type names misspelled, one contains Cyrillic letters
- [ ] **F2629** (medium, fabrication, easy/small) `multiscale/functions.rs:27` — coarse_to_fine_pop hard-codes factor 2, ignores u_fine, clamps negatives
- [ ] **F2631** (medium, fabrication, medium/small) `multiscale/types.rs:373` — LbmDsmcCoupling discards DSMC temperature; flux_match is a 50/50 blend
- [ ] **F2633** (medium, bug, medium/small) `multiscale/types.rs:770` — Chapman-Enskog f1 off by 1/cs² versus its sigma_xy; f2 invented
- [ ] **F2634** (medium, perf, medium/small) `multiscale/types.rs:991` — LBM-MD velocity averaging is all-pairs O(nodes×particles); can panic
- [ ] **F2635** (medium, fabrication, medium/medium) `multiscale/types.rs:1088` — HMM micro_to_macro_viscosity runs no micro solve; shear_rate unused
- [ ] **F2637** (medium, bug, easy/trivial) `multiscale/types.rs:1506` — LbmNavierStokesBridge viscous stress has wrong sign and extra ρ
- [ ] **F2642** (medium, fabrication, easy/small) `multiscale_lbm.rs:622` — Burnett coefficients fabricated with invented Kn dependence
- [ ] **F2644** (medium, fabrication, medium/medium) `neural_lbm.rs:480` — DataDrivenLbm::dmd_modes is not DMD; rom_predict has no amplitudes
- [ ] **F2650** (medium, bug, easy/trivial) `non_newtonian/viscoelastic.rs:249` — MaxwellFluid::creep_compliance omits instantaneous elastic 1/G
- [ ] **F2651** (medium, fabrication, medium/small) `non_newtonian/viscoelastic.rs:353` — PTT steady-shear viscosity formula is fabricated (no thinning at ξ=0)
- [ ] **F2652** (medium, fabrication, medium/small) `non_newtonian/viscoelastic.rs:405` — Giesekus steady shear viscosity is an invented closed form
- [ ] **F2653** (medium, bug, easy/trivial) `non_newtonian/viscoelastic.rs:493` — von_mises_2d shear coefficient 1 instead of 3
- [ ] **F2654** (medium, bug, easy/trivial) `non_newtonian/viscoelastic.rs:630` — Johnson–Segalman steady shear uses ξ(1−ξ) instead of ξ(2−ξ)
- [ ] **F2656** (medium, panic, easy/trivial) `particle_coupling.rs:168` — interpolate_fluid_velocity panics for nx==1 or ny==1
- [ ] **F2658** (medium, fabrication, easy/small) `particle_coupling.rs:899` — LiftForce::saffman_lift is not the Saffman formula
- [ ] **F2659** (medium, bug, easy/trivial) `particle_laden_lbm.rs:179` — Bilinear interpolation extrapolates with negative weights outside domain
- [ ] **F2661** (medium, bug, easy/trivial) `particle_laden_lbm.rs:329` — Saffman lift missing a diameter factor and shear sign
- [ ] **F2665** (medium, perf, medium/medium) `particle_laden_lbm.rs:486` — compute_dem_collisions is O(n^2) all-pairs
- [ ] **F2668** (medium, bug, medium/small) `particle_laden_lbm.rs:994` — Explicit Euler on Stokes drag diverges for dt > 2τ_p
- [ ] **F2672** (medium, bug, easy/trivial) `phase_field/types.rs:651` — AllenCahn::step ignores dx in Laplacian
- [ ] **F2677** (medium, bug, easy/trivial) `phase_field_lbm.rs:753` — Spinodal growth_rate and fastest_growing_k have wrong signs
- [ ] **F2678** (medium, stub, hard/large) `phase_field_lbm.rs:971` — ThreePhaseFlow has no collision/step; phases never evolve
- [ ] **F2686** (medium, bug, easy/trivial) `polymer_lbm.rs:528` — prandtl_karman_friction iteration diverges; returns ~2 instead of 0.008
- [ ] **F2690** (medium, perf, medium/medium) `porous/types.rs:188` — LBM collide/stream loops serial; full-lattice clone each stream
- [ ] **F2691** (medium, perf, easy/small) `porous/types.rs:225` — Porous D2Q9 solvers clone whole lattice every stream
- [ ] **F2692** (medium, bug, medium/medium) `porous/types.rs:520` — Porous nodes use u·φ heuristic; DBF solver double-counts drag
- [ ] **F2695** (medium, bug, easy/small) `porous_media/functions.rs:712` — Guo forcing omits −u/cs² term, injecting mass every step
- [ ] **F2698** (medium, bug, easy/trivial) `porous_media/types.rs:886` — RevAveraging::average_porosity bounds y loop by x1
- [ ] **F2702** (medium, perf, easy/small) `porous_media_lbm.rs:606` — PorousLBM::step allocates a full lattice each step
- [ ] **F2706** (medium, bug, easy/small) `quantum_lbm.rs:708` — Vinen equation production and decay terms swapped
- [ ] **F2707** (medium, bug, medium/small) `quantum_lbm.rs:840` — SchroedingerLbm continuity coefficient wrong; phase unwrapped, updated in place
- [ ] **F2710** (medium, bug, easy/trivial) `reactive/combustion.rs:1128` — ZeroDReactor::ignition_delay hangs for dt <= 0
- [ ] **F2711** (medium, bug, easy/trivial) `reactive/extended.rs:687` — mixture_fraction nonzero in pure oxidiser stream
- [ ] **F2716** (medium, fabrication, easy/small) `reactive_flow.rs:353` — total_heat_release counts inert diluent as burned and is not a rate
- [ ] **F2717** (medium, bug, easy/trivial) `reactive_flow.rs:1278` — Chapman–Jouguet detonation velocity formula is incorrect
- [ ] **F2719** (medium, fabrication, easy/small) `sediment_transport.rs:237` — VanRijn transport rate is a made-up power of Shields θ
- [ ] **F2721** (medium, perf, easy/medium) `sedimentation/types.rs:164` — Darcy/Sedi LBM allocate nx*ny inner Vecs every stream step
- [ ] **F2723** (medium, fabrication, medium/medium) `sedimentation/types.rs:1178` — FluidizationLbm epsilon_mf nonsense formula; update_phi fake decay
- [ ] **F2726** (medium, bug, easy/small) `simulation/types.rs:188` — LbmSimulation::new accepts ν<=0 / zero dims without validation
- [ ] **F2730** (medium, missing-feature, easy/small) `simulation/types.rs:378` — LbmSimulation offers no mutable grid access or initial-condition setter
- [ ] **F2733** (medium, fabrication, medium/medium) `simulation/types.rs:1064` — Adaptive dt controller has no effect: step() ignores dt
- [ ] **F2734** (medium, fabrication, medium/small) `soft_matter_lbm.rs:345` — Beris-Edwards Q11 co-rotation and flow-alignment zeroed; corot12 sign wrong
- [ ] **F2736** (medium, fabrication, hard/medium) `soft_matter_lbm.rs:1132` — SoftMatterGrid::step has no fluid–LC coupling despite docs
- [ ] **F2738** (medium, bug, hard/medium) `solidification_lbm.rs:377` — PhaseChangeLbm enthalpy not conserved; inconsistent h-T mapping
- [ ] **F2739** (medium, bug, medium/small) `solidification_lbm.rs:562` — Dendritic phase-field driving force and latent heat signs inverted
- [ ] **F2740** (medium, fabrication, medium/medium) `solidification_lbm.rs:839` — Casting porosity marks every solid cell; eutectic composition no-op
- [ ] **F2743** (medium, bug, medium/medium) `streaming.rs:239` — stream_aa_2d even+odd performs two streams, contradicting docs
- [ ] **F2747** (medium, missing-feature, hard/large) `suspension_lbm.rs:151` — SuspensionLbm has no lattice; CouplingStrategy never used
- [ ] **F2748** (medium, bug, easy/trivial) `suspension_lbm.rs:788` — Flocculation::volume_fraction divides by current monomer count
- [ ] **F2750** (medium, bug, easy/trivial) `thermal/types.rs:83` — De Vahl Davis insulating walls applied only at bottom
- [ ] **F2754** (medium, bug, easy/trivial) `thermal_lbm.rs:351` — Boussinesq force makes hot fluid sink with documented negative gravity
- [ ] **F2757** (medium, bug, easy/trivial) `turbulence/functions.rs:94` — compute_strain_rate_d3q19 omits relaxation time τ
- [ ] **F2761** (medium, fabrication, medium/medium) `turbulence/types.rs:254` — DynamicSmagorinsky Cs² from fabricated scalar 'Germano' formula
- [ ] **F2763** (medium, bug, easy/small) `turbulence_model/functions.rs:170` — SST F1/F2 compute tanh(arg)^n instead of tanh(arg^n)
- [ ] **F2766** (medium, fabrication, medium/small) `turbulence_model/types.rs:69` — SyntheticTurbulenceInflow claims divergence-free Smirnov modes but isn't
- [ ] **F2767** (medium, fabrication, medium/small) `turbulence_model/types.rs:76` — Synthetic inflow not divergence-free; meaningless cos(ky*kz) factors
- [ ] **F2768** (medium, bug, easy/trivial) `turbulence_model/types.rs:207` — Vreman ν_t ignores grid spacing dx
- [ ] **F2770** (medium, bug, easy/trivial) `turbulence_model/types.rs:867` — KOmegaField::effective_tau returns omega, not tau
- [ ] **F2772** (medium, bug, easy/trivial) `turbulent_channel/types.rs:181` — ChannelFlow::no_slip_walls overwrites instead of swapping populations
- [ ] **F2773** (medium, bug, medium/small) `turbulent_channel/types.rs:205` — Velocity-shift forcing omits tau; injected force is F/tau
- [ ] **F2775** (medium, fabrication, medium/small) `turbulent_channel/types.rs:795` — computed_re_tau/friction_velocity echo the target; forcing off by 1/τ
- [ ] **F2777** (medium, fabrication, easy/small) `turbulent_channel/types.rs:1020` — KMM/MKM 'DNS reference data' is synthetic log-law plus Gaussian
- [ ] **F2780** (medium, bug, easy/small) `viscoelastic_lbm.rs:163` — Oldroyd-B evolve_stress lacks upper-convected terms; N1 always zero
- [ ] **F2783** (medium, bug, medium/small) `wall_model.rs:125` — Spalart-Allmaras production ignores molecular viscosity, wrong fv2
- [ ] **F2786** (medium, fabrication, medium/medium) `wall_model.rs:628` — WallModelSwitcher NEQWM branch is identical to EWM
- [ ] **F2792** (medium, bug, easy/small) `zou_he/functions.rs:554` — zou_he_pressure_outlet_left_d3q19 overwrites known −x populations

### Low (51) — tracked in the external ledger
F2423, F2429, F2431, F2443, F2445, F2454, F2461, F2463, F2493, F2501, F2506, F2508, F2513, F2517, F2518, F2523, F2543, F2554, F2557, F2567, F2570, F2580, F2595, F2596, F2597, F2613, F2616, F2627, F2638, F2641, F2671, F2676, F2679, F2680, F2681, F2683, F2694, F2699, F2700, F2713, F2715, F2720, F2732, F2742, F2752, F2753, F2759, F2760, F2774, F2781, F2784

<!-- AUDIT-2026-10-12 END -->
