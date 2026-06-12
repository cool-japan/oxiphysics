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
- [ ] ISPH pressure-Poisson via PCG/AMG — **Goal:** ‖∇·v‖ reduced 100× vs
      50 Jacobi sweeps at equal cost; 3D dam break stable. **Design:** replace
      Jacobi PPE with matrix-free PCG using core AMG. Cummins–Rudman 1999.
      *(Shipped: Jacobi-only `IsphPressure`; open delta: the Krylov/AMG path,
      reusing `oxiphysics-core` AMG — no new external dependency.)*

### Viscosity & wall treatment
- [ ] Implicit/semi-implicit viscosity solver — **Goal:** stable lid-driven
      cavity at Re=1 with dt 100× the explicit viscous limit. **Design:**
      matrix-free CG on implicit velocity-Laplacian operator. Takahashi 2015.
      *(Code-verified missing 2026-06-11.)*
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
