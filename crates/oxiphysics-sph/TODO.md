# oxiphysics-sph TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 97,190 SLoC | 4,370 tests

Smoothed-particle hydrodynamics subcrate. Everything under **Completed** is
production code with co-located tests. The v0.2.0 / v0.3.0 / v1.0 sections are
the forward roadmap, code-verified against the tree on 2026-06-11: where an
item builds on infrastructure that already ships, the shipped piece is
recorded under *Code-audit confirmations* below and only the missing delta
remains open. All work is pure Rust; GPU paths target wgpu by default, with
the cudarc (CUDA) backend strictly feature-gated and non-default.

How to read this file: `[x]` shipped, `[ ]` open roadmap item, `[~]` deferred
until its stated readiness condition holds. Every open item carries a
measurable **Goal** (acceptance tolerance) and a **Design** sketch with its
primary reference.

## Completed (v0.1.0 – v0.1.3)

### Foundation, core implementation, polish
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests
- [x] Implement primary algorithms
- [x] Add integration tests
- [x] Performance benchmarks
- [x] Documentation
- [x] Examples
- [x] Optimization

### Implemented modules
- [x] Pressure solvers: IISPH, PCISPH, DFSPH, DFSPH (full), WCSPH
- [x] Kernel functions and neighbor search
- [x] Adaptive smoothing length (`adaptive_h`, `adaptive_sph`)
- [x] Boundary SPH and open boundary conditions
- [x] Free-surface tracking
- [x] Multiphase / immiscible fluid support
- [x] Surface tension (CSF model)
- [x] Turbulence and turbulence SPH models
- [x] Granular flow SPH
- [x] Viscosity models
- [x] Coupling layer (SPH↔rigid, SPH↔FEM, SPH↔DEM)
- [x] CFL-based adaptive timestepping
- [x] Simulation orchestration module
- [x] 4,370 tests passing, 0 stubs

### Code-audit confirmations (2026-06-11) — shipped; must not be re-opened
- [x] δ-SPH density diffusion: Molteni–Colagrossi `DeltaSph` plus
      `fourtakas_density_diffusion` (`wcsph_solver.rs`, `wcsph/functions.rs`)
- [x] Riemann-SPH (Godunov-flux particle interactions)
- [x] Morton Z-order spatial sorting (`ZOrderCurve`, CPU-side)
- [x] Particle split/merge adaptivity primitives (`ParticleSplit`,
      `can_merge` in `adaptive_refinement.rs`)
- [x] ISPH pressure projection with Jacobi PPE sweeps
      (`pressure_solvers::IsphPressure`) — baseline for the PCG/AMG upgrade
- [x] Elastic SPH carrying per-particle `deformation_gradient` and
      Green–Lagrange strain (`elastic_sph`) — baseline for formal TLSPH
- [x] Minimal 16-case marching-cubes surface extraction with color field
      (`surface.rs::MarchingCubesSph`) — baseline for reconstruction v2
- [x] Open-boundary Riemann-invariant BC (`open_boundary::RiemannInvariantBC`)
      — baseline for the active-absorption wave tank
- [x] GPU SPH kernels in `oxiphysics-gpu` (`gpu_sph_*`) — functional but
      readback-bound (`compute_density_gpu` returns `Vec<f32>` every step);
      baseline for the v0.3.0 resident pipeline

## v0.2.0 — Accuracy & boundaries

Theme: tighten pressure-field fidelity and wall treatment before scaling out.

### Density & pressure accuracy
- [ ] δ⁺-SPH (shifting-coupled density diffusion) — **Goal:** hydrostatic tank
      pressure noise std/mean < 1% over 10 s. **Design:** extend existing
      diffusion with free-surface-aware δ⁺ term + particle-shifting coupling.
      Sun et al. 2017. *(Shipped: `DeltaSph` (Molteni–Colagrossi) and
      `fourtakas_density_diffusion`; open delta: the δ⁺ variant only.)*
- [x] ISPH pressure-Poisson via PCG (2026-06-12) — **Goal:** ‖∇·v‖ reduced 100× vs
      50 Jacobi sweeps at equal cost; 3D dam break stable. **Design:** replace
      Jacobi PPE with matrix-free PCG using core AMG. Cummins–Rudman 1999.
      *(Shipped: `isph.rs` — self-contained `CsrMatrix` + Jacobi-PCG, symmetric
      Brookshaw PPE with an iterated (Richardson) pressure projection for a
      consistent divergence-free correction; enclosed-box hydrostatic, ≥90%
      divergence reduction, PCG-convergence, dam-break-stability, and 2D
      Taylor–Green tests. AMG left as a future drop-in preconditioner.)*
      *(Shipped: Jacobi-only `IsphPressure`; open delta: the Krylov/AMG path,
      reusing `oxiphysics-core` AMG — no new external dependency.)*

### Viscosity & wall treatment
- [x] Implicit/semi-implicit viscosity solver — **Goal:** stable lid-driven
      cavity at Re=1 with dt 100× the explicit viscous limit. **Design:**
      matrix-free CG on implicit velocity-Laplacian operator. Takahashi 2015.
      *(Implemented 2026-06-14: `viscosity_implicit.rs` — Morris symmetric Laplacian assembled as an SPD M-matrix `(I − dt·ν·L)`, solved with the in-crate Jacobi-PCG; stability gate at 100× the explicit viscous limit, explicit/implicit steady-state agreement, SPD, and CG-convergence tests.)*
- [ ] Semi-analytic boundary integrals — **Goal:** hydrostatic pressure on flat
      wall within 1% of analytic, no boundary particles. **Design:** wall
      renormalization factor γ + boundary flux (USAW). Ferrand et al. 2013.
      *(Code-verified missing 2026-06-11.)*

**Exit criteria:** all four acceptance targets above pass as tolerance-checked
tests; no regression across the existing 4,370-test suite.

## v0.3.0 — GPU & solids

Theme: device-resident throughput, plus solid mechanics and spatial
adaptivity on top of the v0.2.0 accuracy base.

### GPU-resident execution
- [ ] GPU-resident pipeline (zero readback) — **Goal:** ≥20× CPU at 10⁶
      particles, ≤1 readback per N steps. **Design:** port
      neighbor(Morton)→density→force→integrate to wgpu compute; state stays on
      device. Goswami 2010. *(Shipped: CPU `ZOrderCurve` Morton sort and
      `oxiphysics-gpu` `gpu_sph_*` kernels, currently readback-bound; open
      delta: the device-resident loop.)*
      **Depends on:** `oxiphysics-gpu` v0.2.0 device-side reduction/sort work
      (radix sort + parallel reductions). wgpu is the default backend (pure
      Rust); cudarc remains feature-gated, non-default.
- [ ] Surface reconstruction v2 + GPU marching cubes — **Goal:** splash surface
      Hausdorff error < 0.5h vs reference. **Design:** full 256-case MC +
      anisotropic kernels. Yu–Turk 2013. *(Shipped: 16-case minimal
      `MarchingCubesSph` + color field; open delta: full case table,
      anisotropic kernels, and the wgpu compute port via `oxiphysics-gpu`.)*

### Solid mechanics & resolution adaptivity
- [ ] Total-Lagrangian SPH for solids — **Goal:** cantilever tip deflection
      within 3% of analytic; no rank-deficiency instability. **Design:**
      reference-config kernels + kernel-correction matrix + hourglass control.
      Vignjevic 2006; Ganzenmüller 2015. *(Shipped: `elastic_sph` already
      carries `deformation_gradient` + Green–Lagrange strain — confirm TL vs
      updated formulation; open delta: formal TLSPH + hourglass control.)*
- [ ] Multi-resolution APR cross-resolution coupling — **Goal:** variable-res
      dam break conserves mass to 1e-6, surge front within 2% of uniform-res.
      **Design:** conservative interaction across adjacent resolution bands on
      top of split/merge. Vacondio 2013. *(Shipped: `ParticleSplit` /
      `can_merge` primitives; open delta: the conservative cross-resolution
      interaction operator.)*

**Exit criteria:** 10⁶-particle dam break runs device-resident at the stated
speedup; TLSPH cantilever and variable-resolution dam break meet tolerances.

### Sequencing notes
- The ISPH PCG/AMG upgrade (v0.2.0) reuses `oxiphysics-core` AMG; schedule it
  independently of the GPU track.
- The GPU-resident pipeline gates on `oxiphysics-gpu` v0.2.0 sort/reduction
  primitives — track that crate's TODO before starting.

## v1.0 — Production & validation

Validation ladder with hard tolerances. Each case lands as a reproducible
example plus a tolerance-checked regression test with stored reference
baselines (same pattern as the LBM Ghia Re=100 cavity baseline).

### Validation cases
- [ ] Kleefsman 3D dam break — **Goal:** peak pressure at sensor P2 within
      15%, front arrival within 10%. **Design:** box-with-obstacle, pressure
      sensors P1–P4 vs Kleefsman 2005 experiment.
- [ ] Taylor–Green vortex convergence — **Goal:** WCSPH spatial order ≥1.8 on
      log-log L2 velocity error; KE decay matches analytic. **Design:**
      periodic TGV, refinement study.
- [ ] Wave tank with active absorption — **Goal:** reflection coefficient < 5%
      for regular waves. **Design:** piston wavemaker + active absorption BC.
      *(Builds on shipped `open_boundary::RiemannInvariantBC`.)*

### Production hardening
- [ ] Solver benchmark matrix — criterion benches covering
      WCSPH/PCISPH/DFSPH/ISPH(PCG) on the dam-break and TGV cases, tracked
      across releases to catch performance regressions.
- [ ] Deterministic-mode accumulation — adopt `oxiphysics-core`
      compensated/ordered parallel reductions in density and force sums so
      multi-thread runs are bit-reproducible (cross-link: core v0.2.0
      "Numerical foundations & determinism" roadmap).
- [ ] API stabilisation pass for new solver configs — `#[non_exhaustive]`
      where variants may grow; rustdoc + runnable example for every new
      public type added in v0.2.0/v0.3.0.

## Deferred / research track

- [~] cudarc (CUDA) backend parity for the GPU-resident pipeline — Ready
      when: the wgpu-resident pipeline (v0.3.0) has landed and
      `oxiphysics-gpu`'s cudarc feature has stabilised; stays feature-gated
      and non-default per the Pure Rust default policy.
- [~] Device-side PPE solve (PCG on GPU for ISPH) — Ready when: the
      wgpu-resident pipeline lands and `oxiphysics-gpu` ships device
      reduction primitives (prerequisite for CG dot products on-device).
- [~] Device-side f64 accumulation for GPU density/PPE reductions — Ready
      when: wgpu/WGSL expose shader-f64 on target adapters, or
      `oxiphysics-gpu` ships compensated two-float accumulation primitives.

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 333 (critical 11, high 116, medium 169, low 37); 33 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (11)
- [x] **F4428** (critical, bug, medium/small) `dfsph.rs:198` — DFSPH kappa divides by alpha instead of multiplying
- [x] **F4436** (critical, bug, medium/small) `dfsph_solver.rs:684` — DFSPH stiffness divides by alpha instead of multiplying, all DFSPH solvers
- [x] **F4496** (critical, bug, hard/medium) `iisph.rs:379` — IISPH Jacobi uses positive a_ii: compressed fluid gets zero pressure
- [ ] **F4574** (critical, bug, easy/trivial) `particle/functions.rs:16` — Cubic-spline gradient 1.5× too large versus its kernel
- [ ] **F4580** (critical, bug, easy/trivial) `particle/types.rs:1293` — WCSPH pressure force double-counts every fluid-fluid pair
- [ ] **F4582** (critical, bug, easy/trivial) `particle/types.rs:1339` — Morris viscosity uses dynamic μ in kinematic formula (ρ₀× too viscous)
- [x] **F4698** (critical, bug, easy/small) `transport_sph/types.rs:852` — DiffusionSph Brookshaw Laplacian sign inverted: anti-diffusive
- [x] **F4700** (critical, bug, easy/trivial) `transport_sph/types.rs:1215` — ThermalConduction::compute_dtdt anti-diffusive: heat flows cold to hot
- [x] **F4715** (critical, bug, medium/small) `turbulence_sph.rs:904` — TurbulentViscosityForce is anti-viscous and fakes mass as ρh³
- [ ] **F4726** (critical, bug, easy/trivial) `viscosity.rs:1218` — delta_sph_diffusion_contribution anti-diffusive; test only checks finiteness
- [ ] **F4739** (critical, bug, easy/small) `wcsph/functions.rs:1158` — fourtakas_density_diffusion is anti-diffusive; test checks unaffected particles

### High (116)
- [ ] **F4383** (high, bug, medium/medium) `adaptive.rs:312` — adaptive_step merges distant particles, can panic, and runs in O(n^3)
- [ ] **F4384** (high, bug, medium/small) `adaptive.rs:341` — Merge pass joins distant particles and uses stale swap_remove indices
- [x] **F4385** (high, bug, easy/small) `adaptive.rs:502` — Duplicate ∂W/∂h and Ω sign errors in adaptive.rs
- [x] **F4386** (high, bug, easy/small) `adaptive_h.rs:134` — ∂W/∂h sign error and Ω sign error give Ω≈3
- [ ] **F4389** (high, fabrication, medium/small) `adaptive_refinement.rs:640` — Grad-h Omega is a hard-coded constant 1+1/d
- [x] **F4390** (high, bug, easy/trivial) `adaptive_sph.rs:92` — VariableHKernel Wendland C2 normalisation 8x too large
- [x] **F4394** (high, fabrication, medium/medium) `astrophysical_sph/functions.rs:158` — leapfrog_step ignores viscosity and computes no hydrodynamic forces
- [x] **F4397** (high, bug, easy/trivial) `astrophysics_sph.rs:472` — apply_pair reaction force wrong sign; momentum not conserved
- [ ] **F4400** (high, bug, medium/small) `astrophysics_sph.rs:601` — AdiabticSph: viscosity heats entropy but is missing from momentum
- [x] **F4409** (high, bug, easy/small) `climate_sph.rs:462` — Weiss CO2 solubility formula wrong: K0 ≈ 4e45
- [ ] **F4410** (high, fabrication, medium/medium) `climate_sph.rs:480` — ocean_pco2 ignores alkalinity and hardcodes pH; K0 off by 1e46
- [x] **F4411** (high, bug, easy/small) `coastal_sph.rs:431` — TidalFlat Manning flux flows uphill (sign inverted)
- [x] **F4413** (high, fabrication, medium/medium) `compressible_sph.rs:205` — CompressibleSph::step ignores dt and never moves particles
- [x] **F4415** (high, bug, easy/small) `coupling/functions.rs:161` — Coupling cubic_kernel uses 8/(πh³) with 2h support: 8x too large
- [ ] **F4431** (high, bug, easy/small) `dfsph_full.rs:211` — Divergence-free kappa missing rho factor, correction about 1/rho too weak
- [x] **F4434** (high, perf, easy/small) `dfsph_full.rs:509` — DfsphSimulation rebuilds neighbours with O(n²) brute force each step
- [x] **F4435** (high, bug, easy/small) `dfsph_solver.rs:70` — dfsph_solver Wendland C2/C4 kernels normalised 8x too large
- [x] **F4438** (high, bug, easy/trivial) `dfsph_solver.rs:901` — Morris viscosity is anti-dissipative and missing a density factor
- [x] **F4441** (high, perf, easy/small) `dfsph_solver.rs:1127` — DFSPH neighbour search probes 125 cells and clones lists per particle
- [x] **F4442** (high, bug, medium/small) `dfsph_solver.rs:1269` — Divergence solve re-applies accumulated pressure every iteration
- [x] **F4443** (high, perf, easy/small) `dfsph_solver.rs:1277` — Whole neighbor_lists cloned per particle per solver iteration
- [x] **F4444** (high, bug, medium/medium) `dfsph_solver.rs:1315` — solve_constant_density uses frozen density, accumulates and reapplies pressure
- [x] **F4445** (high, bug, medium/small) `dfsph_solver.rs:1325` — solve_constant_density uses frozen density and accumulates pressure unboundedly
- [x] **F4448** (high, bug, easy/trivial) `elastic_sph/functions.rs:236` — stress_divergence sign inverted: compressed particles attract
- [x] **F4449** (high, bug, easy/trivial) `elastic_sph/functions.rs:281` — hourglass_control injects energy and has wrong units
- [x] **F4451** (high, fabrication, medium/small) `elastic_sph/functions.rs:758` — deformation_gradient_sph returns an unnormalised kernel moment, not F
- [x] **F4452** (high, bug, medium/small) `elastic_sph/functions.rs:793` — deformation_gradient_sph is unnormalised; undeformed F ≠ I
- [x] **F4455** (high, bug, easy/trivial) `elastic_sph/types.rs:630` — Jaumann co-rotational term has the wrong sign
- [x] **F4462** (high, bug, easy/small) `explosion_sph.rs:104` — Kinney–Graham Ps/P0 ratio treated as kPa: overpressure ~100x low
- [x] **F4464** (high, bug, medium/small) `explosion_sph.rs:630` — Sod exact Riemann solver has reversed contact velocity
- [x] **F4468** (high, bug, easy/small) `fracture_sph/functions.rs:125` — max_eigenvector_3x3 converges to the minimum-eigenvalue direction
- [x] **F4469** (high, bug, easy/trivial) `fracture_sph/functions.rs:149` — Fracture Wendland C2 kernel normalisation 8x too large
- [x] **F4472** (high, fabrication, hard/medium) `fracture_sph/types.rs:714` — Fracture step never integrates momentum; stress does not affect motion
- [x] **F4473** (high, bug, easy/small) `free_surface.rs:133` — SurfaceTension normals/curvature hardcode h=1.0, ignore smoothing_length
- [x] **F4482** (high, bug, easy/trivial) `geothermal_sph.rs:240` — water_viscosity Vogel formula misapplied: ~3.6e4x too small
- [x] **F4483** (high, bug, easy/small) `geothermal_sph.rs:1036` — Thermal diffusion missing 1/rho_i and Brookshaw factor 2
- [ ] **F4486** (high, fabrication, medium/medium) `granular/types.rs:264` — GranularSim ignores friction parameters, so rolling resistance is a no-op
- [ ] **F4489** (high, bug, medium/medium) `granular/types.rs:731` — DEM tangential force accelerates sliding and has no spring history
- [x] **F4494** (high, fabrication, hard/medium) `iisph.rs:317` — IisphSolver dimensionally inconsistent, not Ihmsen IISPH
- [ ] **F4499** (high, bug, easy/small) `immersed_boundary_sph.rs:596` — Prescribed body motion compounds because centre is overwritten
- [ ] **F4502** (high, bug, medium/medium) `incompressible_sph.rs:343` — ISPH pressure Poisson Jacobi is not a Laplacian; compression gives zero pressure
- [ ] **F4503** (high, bug, easy/small) `incompressible_sph.rs:430` — ISPH viscosity and SGS terms are anti-diffusive and wrongly divided by ρ
- [x] **F4507** (high, bug, easy/trivial) `kernel.rs:1051` — QuinticWendlandC6Kernel uses 2D normalisation in 3D
- [x] **F4509** (high, bug, easy/trivial) `kernel.rs:1341` — sph_gradient_3d returns −∇f (uses x_j−x_i with dW/dr)
- [x] **F4510** (high, fabrication, medium/small) `lagrangian_particles.rs:176` — FTLE identically 1.0; deformation gradient and MSD use invented reference
- [ ] **F4518** (high, bug, easy/small) `magnetohydro.rs:229` — MHD curl_b and induction RHS return negated curls
- [x] **F4519** (high, fabrication, medium/medium) `magnetohydro.rs:328` — dedner_cleaning uniformly decays B instead of cleaning ∇·B
- [x] **F4524** (high, fabrication, medium/small) `mls_sph.rs:121` — MLS gradient formula is wrong and fails linear consistency
- [x] **F4529** (high, fabrication, medium/small) `multiphase/types.rs:783` — Multiphase curvature set to |∇c|; particle mass hardcoded 1.0
- [x] **F4532** (high, fabrication, medium/small) `multiphase_sph/functions.rs:125` — compute_curvature differences colour over particle index order
- [ ] **F4533** (high, fabrication, medium/small) `multiphase_sph/functions.rs:132` — compute_curvature takes a finite difference over particle index order
- [ ] **F4537** (high, fabrication, medium/medium) `multiphase_sph/types.rs:416` — MultiphaseSphSimulation::step reads curvature that is never computed
- [ ] **F4538** (high, fabrication, medium/small) `multiphase_sph/types.rs:423` — MultiphaseSphSimulation never computes curvature, so CSF force is always zero
- [ ] **F4545** (high, bug, easy/small) `nanoscale_sph.rs:1026` — DPD forces are attractive and anti-dissipative; integrate applies a half kick
- [ ] **F4546** (high, bug, easy/small) `neighbor/types.rs:168` — PbcLinkedCellList returns duplicate neighbours when an axis has <3 cells
- [ ] **F4548** (high, bug, easy/trivial) `neighbor/types.rs:666` — AdaptiveNeighborSearch scans only h_i and misses neighbours covered by h_j
- [ ] **F4549** (high, panic, easy/small) `neighbor/types.rs:979` — Octree build recurses forever on 9+ coincident or NaN positions
- [ ] **F4557** (high, fabrication, medium/small) `neural_sph.rs:248` — Neural divergence correction never changes velocities; pressure sign inverted
- [ ] **F4558** (high, fabrication, medium/small) `neural_sph.rs:335` — Divergence-free correction never changes velocity and has the wrong sign
- [ ] **F4560** (high, bug, easy/trivial) `nuclear_sph/functions.rs:138` — Pauli exclusion force is attractive; its test is a tautology
- [ ] **F4566** (high, bug, easy/trivial) `ocean_sph.rs:696` — wendland_c2_grad uses 7/(4π) normalisation instead of 21/(16π)
- [ ] **F4568** (high, perf, medium/medium) `ocean_sph.rs:978` — OceanSimulation::step builds an n-element Vec per particle, twice
- [ ] **F4578** (high, perf, medium/medium) `particle/types.rs:1024` — WCSPH density and force loops are O(n²) all-pairs
- [ ] **F4581** (high, bug, easy/trivial) `particle/types.rs:1325` — accumulate_pressure_acceleration double-counts every fluid-fluid pair
- [ ] **F4583** (high, bug, easy/trivial) `particle/types.rs:1371` — Morris viscosity treats dynamic mu (Pa·s) as kinematic nu
- [ ] **F4590** (high, bug, medium/small) `plasma_sph.rs:309` — induction_rhs returns about −½∇×(v×B)
- [ ] **F4591** (high, bug, easy/trivial) `plasma_sph.rs:333` — MhdSph::lorentz_force computes −∇×B, reversing the Lorentz force
- [ ] **F4592** (high, fabrication, medium/medium) `plasma_sph.rs:344` — Dedner divergence cleaning is a no-op; psi never sourced
- [ ] **F4594** (high, bug, medium/medium) `plasma_sph.rs:391` — MhdSph::euler_step: no gas pressure, and the curl operators are sign-inverted
- [ ] **F4597** (high, fabrication, hard/medium) `pressure_solvers/types_advanced.rs:365` — DFSPH and PCISPH iterative solvers never update state, so cannot converge
- [ ] **F4598** (high, fabrication, medium/medium) `pressure_solvers/types_advanced.rs:684` — PcisphSolver::predict_correct_step never applies pressure; densities frozen
- [x] **F4599** (high, fabrication, medium/medium) `pressure_solvers/types_core.rs:144` — IISPH, PCISPH corrector and projector use fake 1/r kernel gradients
- [x] **F4601** (high, bug, easy/small) `pressure_solvers/types_core.rs:561` — Continuity drho_dt and ISPH divergence have inverted signs
- [x] **F4605** (high, bug, easy/small) `quantum_sph.rs:259` — Gaussian Laplacian sign flipped and Bohm force off by m² factor
- [x] **F4614** (high, bug, easy/trivial) `sediment/types.rs:427` — Rubey settling uses d^2 instead of d^3 in viscous term
- [x] **F4615** (high, bug, easy/trivial) `sediment/types.rs:796` — SphSedimentAdvection::diffusion_rate is anti-diffusive (extra minus)
- [x] **F4617** (high, bug, easy/trivial) `sediment_sph.rs:142` — Ferguson-Church settling uses sqrt(R g d^4) instead of d^3
- [x] **F4618** (high, bug, easy/trivial) `sediment_sph.rs:457` — sph_diffusion_pair has anti-diffusive sign
- [x] **F4619** (high, fabrication, medium/medium) `sediment_sph.rs:988` — morphological_step fakes Exner update with placeholders and magic constants
- [ ] **F4623** (high, fabrication, hard/large) `sediment_transport.rs:734` — TurbidityCurrentSph has no particle interaction at all
- [x] **F4628** (high, bug, easy/small) `shifting.rs:114` — fickian_shift moves particles toward clusters and has m³ units
- [x] **F4629** (high, perf, medium/medium) `simulation/types.rs:139` — SphSim density/pressure/viscosity loops are O(N²) all-pairs
- [ ] **F4636** (high, bug, easy/trivial) `simulation/types_sim.rs:294` — SphSimulation Tait EOS divides stiffness by gamma twice
- [ ] **F4639** (high, panic, easy/trivial) `simulation/types_sim.rs:484` — apply_periodic_bc loops forever on infinite or huge coordinates
- [ ] **F4641** (high, bug, medium/medium) `simulation/types_sim.rs:725` — PeriodicSphSim wraps positions but has no periodic interactions
- [ ] **F4642** (high, bug, easy/small) `soil_sph.rs:391` — Drucker-Prager return map volumetric correction is 3x too small
- [ ] **F4645** (high, bug, medium/small) `soil_sph.rs:503` — Pore-pressure diffusion uses kernel W instead of a Laplacian
- [ ] **F4646** (high, fabrication, easy/trivial) `soil_sph.rs:542` — bishop_simplified ignores unit_weight, hardcodes 18
- [ ] **F4647** (high, perf, easy/small) `soil_sph.rs:931` — SoilSphSolver O(n²) neighbour build plus 3 deep clones per step
- [ ] **F4654** (high, bug, easy/trivial) `sph_analysis.rs:979` — sph_curl computes negated curl (cross product order reversed)
- [ ] **F4655** (high, bug, easy/trivial) `sph_analysis.rs:1019` — sph_laplacian Brookshaw operator returns negated Laplacian
- [ ] **F4656** (high, bug, easy/trivial) `sph_boundary.rs:182` — Lennard-Jones wall force attractive for r0 < r < 2r0
- [x] **F4658** (high, bug, easy/small) `sph_coupling.rs:142` — transfer_stress multiplies traction by volume, giving N·m not N
- [x] **F4660** (high, bug, easy/small) `sph_coupling.rs:335` — penalty_coupling force is attractive, unsigned and zero at contact
- [x] **F4665** (high, bug, easy/small) `sph_thermal.rs:191` — compute_heat_flux is missing the r factor and has an extra 2x
- [x] **F4666** (high, fabrication, medium/small) `sph_thermal.rs:266` — phase_change() latent heat never affects temperature evolution
- [x] **F4671** (high, bug, medium/small) `surface_tension/functions.rs:387` — Marangoni force has the wrong sign and is applied in the bulk, not tangentially
- [x] **F4674** (high, fabrication, medium/small) `surface_tension/types.rs:328` — Akinci curvature term replaced by magic 0.001 density-difference term
- [x] **F4679** (high, fabrication, medium/small) `thermal_sph/types.rs:249` — ThermalSphSimulation::step silently ignores configured phase change
- [x] **F4690** (high, test, easy/small) `transport_sph/functions.rs:168` — Vacuous transport tests hide the anti-diffusion sign bugs
- [x] **F4691** (high, fabrication, easy/small) `transport_sph/functions.rs:491` — Kernel-value-only 'Laplacians' are placeholders, not consistent SPH operators
- [x] **F4692** (high, fabrication, medium/small) `transport_sph/types.rs:64` — AdvectionScheme Upwind anti-diffusive; WenoSph limiter identically one
- [x] **F4695** (high, bug, easy/small) `transport_sph/types.rs:422` — HeatTransferSPH conduction: inverted sign and doubled conductivity
- [x] **F4699** (high, bug, medium/small) `transport_sph/types.rs:936` — AdvectionSPH anti-diffusive stabiliser; semi_lagrangian flag is fake
- [x] **F4701** (high, bug, medium/small) `transport_sph/types.rs:1315` — FickianDiffusion Laplacian anti-diffusive; anisotropic tensor silently ignored
- [x] **F4703** (high, bug, easy/small) `turbulence/types.rs:164` — SPS stress: C_I squared, factor 2 missing, density ignored
- [x] **F4704** (high, bug, easy/trivial) `turbulence/types.rs:334` — Smagorinsky SGS acceleration anti-dissipative: +∇·τ with τ = −2ρν_tS
- [x] **F4705** (high, fabrication, hard/medium) `turbulence/types.rs:383` — Dynamic Smagorinsky Germano coefficient is a constant: L = M/2
- [x] **F4707** (high, bug, easy/small) `turbulence/types.rs:672` — turbulence::SigmaModel singular values wrong: cubic q has wrong sign
- [x] **F4712** (high, fabrication, medium/small) `turbulence_sph.rs:610` — Root-exported SigmaModel::eddy_viscosity is fabricated Smagorinsky×0.1
- [ ] **F4719** (high, panic, easy/small) `variable_resolution.rs:364` — Octree insert recurses without bound on coincident points or leaf size 0
- [x] **F4720** (high, bug, easy/small) `viscoelastic_sph.rs:239` — Oldroyd-B/FENE-P/Giesekus polymer source term missing 1/λ factor
- [ ] **F4727** (high, bug, medium/small) `viscosity_implicit.rs:183` — Implicit viscosity multiplies a μ/ρ Laplacian by ν again; matrix asymmetric
- [ ] **F4734** (high, bug, easy/small) `wcsph/functions.rs:132` — WCSPH viscosity not scaled by ρ; cubic Laplacian sign-indefinite
- [ ] **F4740** (high, bug, easy/trivial) `wcsph/functions.rs:1461` — Riemann interface pressure dissipation term has anti-dissipative sign
- [ ] **F4741** (high, bug, easy/trivial) `wcsph/functions.rs:1541` — dynamic_boundary_pressure hydrostatic term sign inverted; test tautological

### Medium (169)
- [ ] **F4380** (medium, missing-feature, medium/medium) `crates/oxiphysics-sph/Cargo.toml:13` — No serde/parallel features or benches in the SPH crate
- [ ] **F4381** (medium, perf, easy/small) `acoustic_sph.rs:161` — LinearAcousticSph::step computes all-pairs O(N²) right-hand sides
- [ ] **F4382** (medium, bug, easy/small) `acoustic_sph.rs:212` — 3D acoustic SPH uses the 1D kernel constant and halves the pressure gradient
- [ ] **F4387** (medium, bug, medium/medium) `adaptive_h.rs:610` — coarsen_fine_particles merges particles by array index, not by location
- [ ] **F4388** (medium, bug, easy/small) `adaptive_refinement.rs:145` — 3D split children are not centred on the parent, so the centre of mass shifts
- [ ] **F4391** (medium, bug, easy/trivial) `adaptive_sph.rs:106` — Public VariableHKernel::dw_dr returns a wrong derivative
- [ ] **F4392** (medium, bug, easy/trivial) `adaptive_sph.rs:209` — compute_omega grad-h factor has the wrong dimensions
- [ ] **F4396** (medium, stub, medium/small) `astrophysics_sph.rs:252` — compute_all_tree is direct O(N²); theta unused
- [ ] **F4398** (medium, bug, medium/small) `astrophysics_sph.rs:479` — apply_pair never heats particle j and the Balsara switch does nothing
- [ ] **F4401** (medium, perf, medium/medium) `astrophysics_sph.rs:642` — AdiabticSph::step runs two all-pairs O(N²) passes per step
- [ ] **F4402** (medium, panic, easy/small) `astrophysics_sph.rs:990` — GalaxyDiskModel::sample_particles panics or loops forever on edge inputs
- [ ] **F4403** (medium, panic, easy/trivial) `astrophysics_sph.rs:998` — sample_particles panics when scale_height=0 and loops forever when sigma0=0
- [ ] **F4404** (medium, bug, easy/trivial) `boundary_sph.rs:188` — compute_boundary_pressure double-counts the hydrostatic term
- [ ] **F4407** (medium, bug, easy/trivial) `climate_sph.rs:88` — OceanLayerSph pressure is 10x too small and the constructor panics on mismatched lengths
- [ ] **F4408** (medium, panic, easy/trivial) `climate_sph.rs:386` — ElNino::step panics when delay rounds to zero; off-by-one delay
- [ ] **F4417** (medium, perf, easy/medium) `coupling/types.rs:320` — Coupling transfers scan every SPH particle for each DEM particle or FEM node
- [ ] **F4418** (medium, bug, medium/small) `coupling/types.rs:332` — hydro_force_on_dem counts buoyancy twice and its drag has the wrong units
- [ ] **F4419** (medium, bug, easy/trivial) `coupling/types.rs:914` — RigidBodyCoupling::sph_coupling_force leaves out the fluid mass m_j
- [ ] **F4420** (medium, bug, medium/small) `coupling/types.rs:1344` — RigidBody6Dof divides a world-frame torque by body-frame inertia
- [ ] **F4421** (medium, bug, easy/small) `cryogenic_sph/types.rs:127` — Thermal velocities are uniform, not Maxwellian, and panic at T=0
- [ ] **F4422** (medium, perf, easy/small) `cryogenic_sph/types.rs:141` — CryogenicSphSimulation::run uses all-pairs density and force loops
- [ ] **F4423** (medium, fabrication, medium/medium) `cryogenic_sph/types.rs:837` — GPE solver claims split-step Fourier but uses explicit FTCS
- [ ] **F4424** (medium, stub, medium/small) `dem_sph.rs:185` — DEM friction coefficient never used; torques always zero
- [ ] **F4425** (medium, fabrication, medium/medium) `dem_sph.rs:233` — DemSimulation ignores friction, torques are always zero, and it panics when particles is pushed directly
- [ ] **F4426** (medium, bug, easy/trivial) `dem_sph.rs:353` — dem_cfl_timestep multiplies by the radius, so dt has units of m·s
- [ ] **F4429** (medium, fabrication, medium/small) `dfsph.rs:712` — correct_density_error_warm_start is identical to the cold solve
- [ ] **F4440** (medium, fabrication, easy/small) `dfsph_solver.rs:1026` — apply_no_slip negates velocity anywhere; FreeSlip is a no-op
- [ ] **F4453** (medium, bug, medium/small) `elastic_sph/functions.rs:845` — hourglass_correction subtracts a velocity from a displacement
- [ ] **F4454** (medium, bug, easy/trivial) `elastic_sph/types.rs:263` — Johnson-Cook radial return gives 1.5x too much plastic strain
- [ ] **F4457** (medium, missing-feature, hard/medium) `electromagnetic_sph.rs:283` — Coulomb treecode advertised but theta_bh unused; force evaluated twice per pair
- [ ] **F4458** (medium, bug, easy/trivial) `electromagnetic_sph.rs:573` — Dipole–dipole force uses r⁵ with unit vectors (extra 1/r)
- [ ] **F4459** (medium, stub, medium/small) `electromagnetic_sph.rs:707` — ElectrokineticSph never computes EDL potential; step is effectively a no-op
- [ ] **F4460** (medium, bug, easy/small) `electromagnetic_sph.rs:970` — PlasmaSheath field interpolation returns −φ instead of −∇φ
- [ ] **F4461** (medium, api, medium/medium) `error.rs:28` — Unified crate Error is dead: no API returns it
- [ ] **F4466** (medium, missing-feature, medium/medium) `explosion_sph.rs:1153` — Blast driver has no time integrator; 1-D mass with 3-D kernel
- [ ] **F4467** (medium, bug, easy/small) `explosion_sph.rs:1256` — Probit-to-probability uses tanh instead of the normal CDF
- [ ] **F4470** (medium, bug, medium/small) `fracture_sph/types.rs:120` — Damage increment independent of dt (fabricated strain rate)
- [ ] **F4471** (medium, perf, medium/medium) `fracture_sph/types.rs:629` — Fracture update_stress and fragment labelling are O(n²); no momentum update
- [ ] **F4474** (medium, fabrication, medium/medium) `free_surface.rs:346` — WaterSurface dam-break integrates gravity only, no pressure
- [ ] **F4475** (medium, perf, easy/medium) `free_surface.rs:430` — Free-surface routines are O(n²) and droplet IDs are nondeterministic
- [ ] **F4477** (medium, bug, easy/small) `free_surface.rs:792` — Surface reconstruction silently truncates grid to 100 cells/axis
- [ ] **F4478** (medium, bug, easy/small) `free_surface.rs:907` — Paraboloid curvature fit is not least squares (2–4× error)
- [ ] **F4480** (medium, bug, easy/trivial) `free_surface_sph.rs:165` — Free-surface 'outward' normal actually points inward
- [ ] **F4481** (medium, fabrication, hard/large) `geothermal_sph.rs:230` — Steam quality and supercritical density are invented formulas
- [ ] **F4484** (medium, perf, easy/small) `geothermal_sph.rs:1182` — Geothermal step: O(n²), per-particle Vec allocation, per-pair property recomputation
- [ ] **F4485** (medium, bug, medium/small) `granular/types.rs:64` — Rolling resistance uses total particle force; friction mu unused
- [ ] **F4487** (medium, bug, medium/small) `granular/types.rs:552` — GranularCollision damping anti-dissipative; k_n lacks sqrt(R*)
- [ ] **F4488** (medium, perf, medium/medium) `granular/types.rs:685` — DEM contact detection O(n²) in three paths; SpatialHashDem unused
- [ ] **F4490** (medium, bug, easy/trivial) `granular/types.rs:746` — Friction torque on particle j has wrong sign
- [ ] **F4491** (medium, bug, easy/small) `ice_sph.rs:644` — 'Isotropic' fabric and fragment directions not uniform on sphere
- [ ] **F4492** (medium, bug, easy/trivial) `ice_sph.rs:699` — Ice fabric eigenvalues taken from the tensor diagonal
- [ ] **F4493** (medium, fabrication, hard/large) `ice_sph.rs:1014` — ice_sph_step has no SPH, ignores gravity, dimensionally wrong damping
- [ ] **F4497** (medium, stub, easy/small) `iisph.rs:573` — iisph_flip_blend ignores velocities and blend factor (no-op)
- [ ] **F4498** (medium, perf, easy/small) `immersed_boundary_sph.rs:365` — IBM interpolation and spreading scan all fluid particles per marker
- [ ] **F4501** (medium, perf, medium/medium) `incompressible_sph.rs:327` — ISPH step O(iter·n²): all-pairs loops in solve, predict, divergence
- [ ] **F4504** (medium, fabrication, medium/small) `incompressible_sph.rs:701` — ISPH CSF surface tension uses |∇c| as curvature
- [ ] **F4505** (medium, bug, easy/small) `isph.rs:470` — ISPH config validation passes NaN, skips lengths, wrong message
- [ ] **F4508** (medium, bug, easy/trivial) `kernel.rs:1318` — sph_laplacian_1d returns the negated Laplacian
- [ ] **F4512** (medium, fabrication, medium/small) `lagrangian_particles.rs:283` — Deformation gradient fabricated from velocities; no reference configuration
- [ ] **F4516** (medium, bug, easy/trivial) `lubrication_sph.rs:380` — Van der Waals uses centre distance and plate-pressure formula
- [ ] **F4517** (medium, fabrication, easy/small) `lubrication_sph.rs:485` — EHL film thickness hardcodes materials parameter G=4000 and k=1
- [ ] **F4520** (medium, perf, easy/medium) `microfluidics_sph.rs:879` — micro_step runs five all-pairs O(n²) passes
- [ ] **F4521** (medium, bug, medium/medium) `microfluidics_sph.rs:901` — micro_step: wettability BC applied after CSF; phi_e never solved
- [ ] **F4523** (medium, panic, easy/small) `microfluidics_sph.rs:1065` — solve_poisson_boltzmann_1d panics at n=0; slow Gauss-Seidel for a linear system
- [ ] **F4526** (medium, bug, medium/small) `multiphase/types.rs:256` — Interface sharpening is anti-diffusive (wrong signs)
- [ ] **F4527** (medium, bug, easy/trivial) `multiphase/types.rs:688` — add_phase discards previously set surface-tension coefficients
- [ ] **F4528** (medium, bug, easy/small) `multiphase/types.rs:707` — MultiphaseSystem masses hardcoded to 1.0; viscous_force fakes masses
- [ ] **F4530** (medium, perf, easy/medium) `multiphase/types.rs:868` — MultiphaseSystem::step and InterfaceSharpening run all-pairs passes
- [ ] **F4531** (medium, bug, easy/trivial) `multiphase_sph/functions.rs:27` — Colour-gradient loops cut the 2h-support kernel at h
- [ ] **F4534** (medium, fabrication, medium/small) `multiphase_sph/functions.rs:813` — surface_tension_csf fakes neighbour normals from coordinate signs
- [ ] **F4536** (medium, bug, easy/trivial) `multiphase_sph/types.rs:220` — VolumeFractionField upwind update has the wrong sign for v < 0
- [ ] **F4539** (medium, perf, medium/medium) `multiphase_sph/types.rs:825` — Multiphase and reactive per-particle operators scan every particle (O(N²))
- [ ] **F4540** (medium, bug, easy/trivial) `multiphase_sph/types.rs:931` — Cahn-Hilliard initial profile divides tanh by epsilon
- [ ] **F4541** (medium, bug, easy/trivial) `multiphase_sph_ext.rs:169` — Wendland C2 3-D normalisation constant is 8x too large
- [ ] **F4542** (medium, bug, easy/trivial) `multiphase_sph_ext.rs:698` — SPH Fick diffusion has the wrong sign (anti-diffusion)
- [ ] **F4543** (medium, fabrication, easy/small) `multiphase_sph_ext.rs:726` — Colour and level-set advection rates compute an invented product of differences
- [ ] **F4544** (medium, bug, medium/small) `nanoscale_sph.rs:429` — Thermal noise amplitude has units of N·m, RNG is unseeded, momentum not conserved
- [ ] **F4547** (medium, panic, medium/medium) `neighbor/types.rs:266` — LinkedCellList sizes a dense grid from the bounding box, so outliers overflow or OOM
- [ ] **F4550** (medium, perf, medium/medium) `neighbor/types.rs:1112` — Per-step neighbor search uses SipHash HashMap and per-particle Vec allocs
- [ ] **F4552** (medium, perf, easy/small) `neighbor/types.rs:1271` — PeriodicSpatialHash is O(N²) brute force although its docs describe 3×3×3 cells
- [ ] **F4553** (medium, bug, easy/trivial) `neighbor/types.rs:1341` — ZOrderCurve::spread_bits uses 2-D masks; Morton codes collide
- [ ] **F4554** (medium, perf, easy/small) `neighbor/types.rs:1422` — Neighbour builds are serial; sph crate does not use rayon
- [ ] **F4555** (medium, fabrication, easy/small) `neighbor/types.rs:1485` — HashGridStats::from_linked_cell ignores the linked-cell grid it is given
- [ ] **F4556** (medium, bug, medium/small) `neural_sph.rs:88` — NeuralKernel::train_from_data uses fixed-lr gradient descent that diverges or stalls
- [ ] **F4559** (medium, fabrication, medium/small) `neural_sph.rs:382` — PhysicsInformedSph::pde_residual uses vx+vy+vz as 'divergence'
- [ ] **F4561** (medium, bug, easy/trivial) `nuclear_sph/functions.rs:316` — skyrme_force scaled by NUCLEON_MASS, so the mean field is 939x too strong
- [ ] **F4562** (medium, fabrication, easy/small) `nuclear_sph/functions.rs:355` — sph_sound_speed is hard-coded and ignores the Skyrme EOS
- [ ] **F4563** (medium, perf, medium/small) `nuclear_sph/types.rs:429` — update_forces clones every other particle for each particle (N² clones per step)
- [ ] **F4564** (medium, bug, easy/trivial) `ocean_sph.rs:104` — JONSWAP spectrum returns S(f) as S(omega), so Hs is sqrt(2pi) too large
- [ ] **F4565** (medium, bug, easy/trivial) `ocean_sph.rs:467` — M2 tidal acceleration has the wrong sign
- [ ] **F4567** (medium, bug, easy/small) `ocean_sph.rs:953` — seed_from_jonswap sets particle velocity to wave slope a·k (dimensionless)
- [ ] **F4569** (medium, bug, easy/small) `ocean_sph.rs:997` — OceanSimulation::step gates wind on the Stokes flag and ignores viscosity
- [ ] **F4570** (medium, stub, medium/medium) `open_boundary/functions.rs:187` — Turbulence injection is uncorrelated hash noise; length_scale is ignored
- [ ] **F4571** (medium, bug, easy/small) `open_boundary/functions.rs:207` — nscbc_inlet_pressure subtracts rate times pressure with no dt, so result clamps to 0
- [ ] **F4572** (medium, bug, easy/small) `open_boundary/types.rs:511` — CharacteristicBc::apply_inlet returns 2x reference density at rest
- [ ] **F4573** (medium, bug, medium/small) `open_boundary/types.rs:628` — AbsorbingOutlet uses a time difference, not du/dn; normal and vel_prev2 unused
- [ ] **F4577** (medium, bug, easy/trivial) `particle/types.rs:599` — SphParticleSet::apply_pbc hangs on infinite or very large coordinates
- [ ] **F4584** (medium, bug, easy/trivial) `particle/types.rs:1772` — SphParticleSet::divergence_at returns the negated divergence
- [ ] **F4585** (medium, bug, medium/small) `pcisph.rs:98` — compute_delta hard-codes an h/2 lattice and ignores the real particle spacing
- [ ] **F4588** (medium, bug, easy/trivial) `pcisph.rs:774` — pcisph_correction_vector moves compressed particles toward each other
- [ ] **F4589** (medium, perf, medium/medium) `plasma_sph.rs:294` — MhdSph::euler_step makes three all-pairs passes per step
- [ ] **F4593** (medium, bug, easy/small) `plasma_sph.rs:373` — Resistive Laplacian uses W instead of dW/dr, and eta instead of eta/mu0
- [ ] **F4595** (medium, bug, easy/trivial) `plasma_sph.rs:818` — line_cooling is negative (heating) above 1e6 K and discontinuous at log T = 6
- [ ] **F4596** (medium, fabrication, easy/small) `pressure_solvers/types_advanced.rs:148` — PCISPH delta fabricated as h^4/(2n), missing dt and mass
- [ ] **F4600** (medium, fabrication, medium/medium) `pressure_solvers/types_core.rs:211` — PcisphSolverSimple/PcisphIterativeLoop iterate on frozen predicted densities
- [ ] **F4602** (medium, bug, easy/small) `pressure_solvers/types_core.rs:651` — pressure_solvers uses fake ∇W; PcisphCorrector also divides by mass
- [ ] **F4604** (medium, perf, easy/small) `quantum_sph.rs:210` — QuantumSPH loops are all-pairs over an untruncated Gaussian kernel
- [ ] **F4606** (medium, fabrication, medium/small) `quantum_sph.rs:313` — wigner_function returns a Gaussian guess, not the Wigner transform
- [ ] **F4608** (medium, fabrication, medium/small) `quantum_sph.rs:527` — FeynmanPathSPH ignores the external potential when beads move
- [ ] **F4609** (medium, bug, easy/trivial) `reactive_sph.rs:243` — Chapman-Jouguet detonation velocity formula is wrong
- [ ] **F4611** (medium, bug, easy/small) `reactive_sph.rs:277` — reactive_sph_step computes heat after the fuel is depleted, so energy is not conserved
- [ ] **F4612** (medium, fabrication, easy/small) `sediment/functions.rs:902` — bedform_height uses an invented formula that scales dunes with grain size
- [ ] **F4613** (medium, fabrication, medium/small) `sediment/types.rs:187` — melville_pier_scour discards Fr/d50 work; wrong K_yb
- [ ] **F4616** (medium, bug, easy/small) `sediment/types_advanced.rs:525` — van Rijn dune height has d50 and depth swapped (≈20x low)
- [ ] **F4620** (medium, bug, medium/small) `sediment_sph.rs:1130` — ExnerGrid::compute_bed_change does not conserve volume
- [ ] **F4621** (medium, fabrication, easy/trivial) `sediment_sph.rs:1331` — advect_concentration ignores the settling term it documents
- [ ] **F4624** (medium, api, easy/trivial) `sediment_transport.rs:827` — sph_diffusion_rate is anti-diffusive when dw is passed as documented
- [ ] **F4625** (medium, bug, medium/small) `sediment_transport.rs:920` — SedimentParticle::advect_concentration deletes sediment mass
- [ ] **F4630** (medium, bug, easy/trivial) `simulation/types.rs:194` — SphSim pair forces use m_j² on reaction side (Newton III broken)
- [ ] **F4631** (medium, api, easy/medium) `simulation/types.rs:748` — No Params/Config struct validates; SphSimulation::new accepts garbage
- [ ] **F4634** (medium, perf, medium/small) `simulation/types_sim.rs:180` — WcSphSim CSF surface tension does two all-pairs passes
- [ ] **F4635** (medium, bug, easy/small) `simulation/types_sim.rs:288` — Neighbor radius hardcoded to 2h regardless of kernel support
- [ ] **F4640** (medium, fabrication, hard/large) `simulation/types_sim.rs:656` — TwoPhaseSimState never couples phases; sigma_ab unused
- [ ] **F4643** (medium, bug, medium/small) `soil_sph.rs:407` — effective_stress uses the compression-positive Terzaghi form on tension-positive stress
- [ ] **F4644** (medium, bug, easy/trivial) `soil_sph.rs:461` — coefficient_of_consolidation ignores Cc and inverts m_v
- [ ] **F4649** (medium, fabrication, medium/small) `soil_sph.rs:1139` — update_liquefaction: hydrostatic=0, current σ'v; liquefied soil reports r_u=0
- [ ] **F4650** (medium, perf, easy/small) `sph_acoustics.rs:191` — AcousticSPH::step is O(N²) all-pairs with sqrt before the cutoff
- [ ] **F4651** (medium, fabrication, medium/small) `sph_acoustics.rs:463` — Dipole/quadrupole sources: doubled frequency, no directivity
- [ ] **F4661** (medium, bug, easy/trivial) `sph_multiscale.rs:161` — coarsen_region invents density from a sphere of radius h
- [ ] **F4662** (medium, bug, easy/trivial) `sph_multiscale.rs:193` — refine_region children lattice off-centre; COM not conserved
- [ ] **F4663** (medium, bug, easy/trivial) `sph_multiscale.rs:325` — interpolate_velocity kernel lacks 1/h³, so mixed-h weights are wrong
- [ ] **F4664** (medium, perf, easy/small) `sph_thermal.rs:172` — compute_heat_flux is O(N²), with a sqrt per pair and 7 Vec gathers per call
- [ ] **F4668** (medium, stub, medium/medium) `surface.rs:752` — build_anisotropic_kernels is an isotropic stub (G = k_r I)
- [ ] **F4669** (medium, bug, easy/trivial) `surface.rs:786` — implicit_surface_tension_force uses unit normal, dropping the |∇c| delta
- [ ] **F4670** (medium, fabrication, easy/small) `surface.rs:972` — MC triangle count uses a popcount heuristic, not the real table
- [ ] **F4675** (medium, bug, easy/small) `thermal_sph/functions.rs:109` — Enthalpy H(T) jumps at liquidus; inverse non-monotonic
- [ ] **F4676** (medium, bug, easy/trivial) `thermal_sph/types.rs:133` — Bubble departure_frequency ignores bubble diameter and has wrong units
- [ ] **F4677** (medium, bug, easy/trivial) `thermal_sph/types.rs:165` — Zuber critical heat flux formula wrong (~6.7x too low)
- [ ] **F4678** (medium, fabrication, easy/small) `thermal_sph/types.rs:172` — Leidenfrost temperature hardcoded as 1.15·T_sat
- [ ] **F4681** (medium, bug, easy/small) `thermal_sph/types.rs:676` — AnisotropicHeatConduction::dtemp_dt returns W/kg, not K/s
- [ ] **F4683** (medium, fabrication, medium/small) `thermal_sph/types.rs:1204` — NucleateBoilingSph vapor_fraction = Ja/100; q_onb hardcoded 1000
- [ ] **F4686** (medium, fabrication, easy/small) `timestepping.rs:165` — rk2_step/rk4_step_single are fake RK with constant acceleration
- [ ] **F4687** (medium, panic, easy/small) `timestepping.rs:366` — MultiRateTimeStepper panics on n_levels=0 / length mismatch
- [ ] **F4688** (medium, fabrication, easy/small) `timestepping.rs:572` — RungeKuttaSph rk4_step_single/rk2_step use constant acceleration
- [ ] **F4689** (medium, fabrication, easy/small) `timestepping.rs:1135` — AdamsBashforth2 position update collapses to forward Euler
- [ ] **F4694** (medium, perf, easy/small) `transport_sph/types.rs:142` — SoluteTransport runs n_species² full O(N²) passes with deep clones
- [ ] **F4696** (medium, fabrication, medium/medium) `transport_sph/types.rs:551` — Transport 'Laplacians' use W instead of grad W: wrong units and magnitude
- [ ] **F4697** (medium, perf, medium/medium) `transport_sph/types.rs:840` — Transport operators are O(N²) all-pairs with no cutoff or neighbour list
- [ ] **F4702** (medium, bug, easy/small) `turbulence/functions.rs:1132` — lambda2_criterion cubic q has wrong sign; eigenvalues mirrored
- [ ] **F4708** (medium, fabrication, easy/small) `turbulence_sph.rs:121` — velocity_gradient_sph uses h^3 as particle volume (fabricated mass)
- [ ] **F4709** (medium, fabrication, medium/medium) `turbulence_sph.rs:336` — energy_spectrum_estimate is an O(N^2) structure function, not E(k)
- [ ] **F4710** (medium, fabrication, easy/small) `turbulence_sph.rs:393` — SmagorinskyModel::update_particles uses |v|/h as strain rate
- [ ] **F4711** (medium, bug, easy/trivial) `turbulence_sph.rs:513` — Vreman eddy viscosity scales as h*|grad u| (m/s), not h^2*|grad u|
- [ ] **F4713** (medium, fabrication, hard/medium) `turbulence_sph.rs:800` — Dynamic Smagorinsky is a ratio heuristic; Leonard stress uses the unfiltered product
- [ ] **F4717** (medium, panic, easy/small) `variable_resolution.rs:206` — coarsen_particles asserts on empty input, gives NaN at zero mass, panics on bad index
- [ ] **F4718** (medium, bug, easy/small) `variable_resolution.rs:269` — Consistency check uses a truncated Gaussian, hardcoded rho=1000 and O(N^2) pairs
- [ ] **F4721** (medium, bug, medium/small) `viscoelastic_sph.rs:311` — FENE-P Peterlin factor compares stress (Pa) directly with extensibility b
- [ ] **F4722** (medium, bug, easy/small) `viscosity.rs:126` — Morris and Cleary viscosity accelerations are half the correct value
- [ ] **F4723** (medium, bug, easy/trivial) `viscosity.rs:216` — sps_stress_tensor SGS kinetic energy k is about 700x too large
- [ ] **F4729** (medium, bug, easy/small) `volcanic_sph.rs:170` — Lava SPH uses the 2D cubic-spline normalisation on 3D positions
- [ ] **F4731** (medium, bug, easy/small) `volcanic_sph.rs:746` — Magma chamber compressibility combined harmonically; Mogi uplift is dimensionless
- [ ] **F4732** (medium, bug, medium/medium) `volcanic_sph.rs:1012` — Plume momentum double-counts gravity, so buoyant plumes always decelerate
- [ ] **F4733** (medium, perf, medium/medium) `wcsph/functions.rs:15` — WCSPH hot loops use dyn kernel and SipHash neighbour grid
- [ ] **F4735** (medium, bug, easy/trivial) `wcsph/functions.rs:385` — speed_of_sound_tait is too large by a factor of sqrt(gamma)
- [ ] **F4737** (medium, fabrication, medium/medium) `wcsph/functions.rs:991` — 'MLS first-order' correction is Shepard and omits the self term (~20% density inflation)
- [ ] **F4738** (medium, bug, easy/trivial) `wcsph/functions.rs:1006` — MLS/Shepard density corrections omit the self-kernel term
- [ ] **F4743** (medium, bug, easy/small) `wcsph/types.rs:50` — WcsphSolver::compute_adaptive_dt ignores the WCSPH speed of sound
- [ ] **F4744** (medium, missing-feature, medium/medium) `wcsph_solver.rs:192` — WcsphSolverAdvanced has no step(); kernel, compressibility and limiter fields are never read
- [ ] **F4745** (medium, bug, easy/trivial) `wcsph_solver.rs:589` — Predictor-corrector correct() drops a_old, integrating only half the acceleration

### Low (37) — tracked in the external ledger
F4395, F4399, F4416, F4432, F4433, F4437, F4447, F4495, F4513, F4514, F4515, F4535, F4551, F4575, F4576, F4579, F4586, F4587, F4603, F4610, F4626, F4632, F4633, F4637, F4653, F4657, F4659, F4667, F4673, F4680, F4682, F4684, F4685, F4725, F4728, F4736, F4742

<!-- AUDIT-2026-10-12 END -->
