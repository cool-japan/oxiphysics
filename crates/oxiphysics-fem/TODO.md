# oxiphysics-fem TODO

Last updated: 2026-05-17 / v0.1.1

## Phase 1: Foundation
- [x] Define core types and traits (CsrMatrix, TetrahedralMesh, etc.)
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Linear static analysis (LinearTetrahedron, LinearElasticMaterial, LinearStaticAnalysis)
- [x] Nonlinear solvers: Newton-Raphson, BFGS, Riks arc-length
- [x] Adaptive refinement: h-, p-, hp-refinement
- [x] Beam, shell, truss element types
- [x] Constitutive models: NeoHookean, MooneyRivlin, J2Plasticity, Perzyna, ThermoElastic, Orthotropic
- [x] Dynamic FEM: explicit time integration, modal analysis, eigenvalue analysis, wave propagation
- [x] Coupled physics: thermo-mechanical, FEM-LBM, electrochemical FEM
- [x] Extended formulations: XFEM fracture, discontinuous Galerkin, isogeometric analysis, spectral FEM
- [x] Failure & damage: cohesive zone, damage mechanics, fatigue, buckling
- [x] Multi-scale & ROM: homogenization, multiscale FEM, reduced order modeling
- [x] Meshfree methods, level-set interface tracking
- [x] Geomechanics FEM, soil FEM, probabilistic/stochastic FEM, reliability FEM
- [x] Biomechanics FEM, fluid-structure interaction (FSI)
- [x] Topology optimization FEM
- [x] Crystal plasticity FEM
- [x] Sparse linear algebra: CsrMatrix (CSR), PcgSolver (PCG)
- [x] Integration tests (4,561 tests, 0 stubs, 124 source files)
- [x] Performance benchmarks (basic)

## Phase 3: Polish
- [x] Documentation (rustdoc on all public items)
- [x] Extended examples (doc-tests in `parallel_solver` module)
- [x] Benchmark suite expansion (`perf_bench` — `BenchHarness`, SpMV/PCG/GMRES/assembly/Ke timing harness)
- [x] Parallel sparse solver integration (`parallel_solver` — Rayon-parallel SpMV, PCG, GMRES, assembler)
- [x] Further SIMD / parallel assembly optimization (parallel Rayon assembler + PCG dot products)

## Phase 4: Algebraic Multigrid & Large-Problem Solvers (v0.2.0)

> **Goal:** Add a real algebraic-multigrid (AMG) solver alongside the existing PCG+Jacobi, enabling 10⁶-DOF linear elasticity / Poisson problems to converge in O(N) work. Today `parallel_solver.rs` ships PCG with a diagonal preconditioner — effective at ≤10⁵ DOF but asymptotically poor. `error_estimation.rs` has `MultiGridAdaptive` for *error analysis*; it is not a solver.

### 4.1 Ruge-Stüben classical AMG
- [x] `solvers/amg/classical.rs` — strong-connection graph, C/F splitting, direct interpolation `P` (planned 2026-04-24; bundles 4.1–4.1e)
  - **Goal:** A self-contained `AmgClassical` hierarchy builder + V/W-cycle driver capable of solving the 3D Poisson equation on 128³ tets (≈ 2 M DOF) with residual reduction ≥ 0.1 per V-cycle.
  - **Design:** Strong-connection graph (θ=0.25), Ruge-Stüben C/F splitting (two-pass), direct interpolation `P`, Galerkin `A_c = P^T A P` via SpMM, GS/SGS smoothers (2 pre/post sweeps), V-cycle + W-cycle driver, coarse-level PCG solve when DOF < 500.
  - **Files:** `src/solvers/amg/mod.rs`, `src/solvers/amg/classical.rs`, `src/solvers/amg/smoothers.rs`, `src/solvers/amg/cycle.rs`, `src/solvers/amg/graph.rs`, `src/solvers/amg/galerkin.rs`, `src/solvers/mod.rs` (re-export)
  - **Tests:** `unit::strong_connection_1d_poisson`, `unit::cf_splitting_partitions`, `unit::galerkin_triple_preserves_symmetry`, `unit::symmetric_gs_reduces_residual`, `integration::vcycle_poisson_2d_32x32` (reduction ≤ 0.1/cycle), `integration::wcycle_converges_harder_problem`
- [x] Galerkin coarse-grid operator `A_c = P^T A P` (planned 2026-04-24; part of 4.1 classical.rs)
- [x] Gauss-Seidel and symmetric GS smoothers — forward, backward, symmetric sweeps (planned 2026-04-24; part of 4.1 classical.rs)
- [x] V-cycle and W-cycle drivers — `AmgSolver::v_cycle(level, b, x)` recursion down to coarse direct solve (planned 2026-04-24; part of 4.1 classical.rs)
- [x] Coarse-level solve: fall back to existing `PcgSolver` when `n_dofs < n_coarse_threshold` (typ. ≤ 500) (planned 2026-04-24; part of 4.1 classical.rs)

### 4.2 Smoothed-aggregation AMG (alternative path)
- [x] `solvers/amg/smoothed_aggregation.rs` — aggregate-based coarsening via strength-of-connection + greedy aggregation (planned 2026-04-24; bundles 4.2–4.2c)
  - **Goal:** Alternative coarsening path via greedy aggregation + rigid-body near-null-space, using the same V/W-cycle driver. SA typically outperforms classical on systems-of-PDEs (3D elasticity).
  - **Design:** SA strength metric `|A[i,j]|² ≥ θ²|A[i,i]A[j,j]|` (θ=0.08), two-pass greedy aggregation, tentative prolongator `P̃` from 6 rigid-body modes (3 translations + 3 rotations) per aggregate orthonormalized via QR, Jacobi-smoothed `P = (I − ω D⁻¹ A) P̃` (ω = 4/(3ρ), ρ via 20 power iterations). Cycle driver reuses `cycle.rs` and `smoothers.rs`.
  - **Files:** `src/solvers/amg/smoothed_aggregation.rs`, `src/solvers/amg/aggregation.rs`, `src/solvers/amg/near_null_space.rs`
  - **Tests:** `unit::aggregation_covers_all_nodes`, `unit::rbm_zero_residual`, `integration::sa_vcycle_elasticity_32cube` (reduction ≤ 0.15/cycle)
- [x] Tentative prolongator from near-null-space — rigid body modes for elasticity (planned 2026-04-24; part of 4.2)
- [x] Jacobi-smoothed prolongator `P = (I - ω D⁻¹ A) P̃` (planned 2026-04-24; part of 4.2)
- [x] Plug into the same V/W-cycle driver as classical AMG (planned 2026-04-24; part of 4.2)

### 4.3 Preconditioned Krylov as outer
- [x] `PcgWithAmg` — use the AMG V-cycle as the preconditioner inside PCG (planned 2026-04-24; bundles 4.3a,b)
  - **Goal:** Wrap existing PCG + GMRES with an AMG V-cycle preconditioner; typically cuts iteration counts 10–50× vs Jacobi-PCG on 10⁶-DOF problems.
  - **Design:** `Preconditioner` trait (`fn apply(&self, r: &[f64], z: &mut [f64])`), `AmgPreconditioner` running one V-cycle, `PcgWithAmg` wrapping `ParallelPcgSolver` (line 323), `GmresWithAmg` wrapping `ParallelGmresSolver` (line 457 — confirmed present).
  - **Files:** `src/solvers/amg/preconditioner.rs` (NEW), `src/parallel_solver.rs` (add Preconditioner trait, ~100 LoC diff)
  - **Tests:** `integration::pcg_amg_poisson_64cube` (≤ 15 outer PCG iters), `integration::gmres_amg_advection_diffusion`
- [x] `GmresWithAmg` — for non-symmetric systems (advection-dominated, some coupled physics) (planned 2026-04-24; part of 4.3)

### 4.4 Parallel assembly scale-up
- [x] Element-coloring assembly so `CsrMatrix::assemble` scales on 32+ cores (planned 2026-04-24; bundles 4.4a,b)
  - **Goal:** `CsrMatrix::assemble_colored` scales linearly to 32+ cores; `spmv_chunked` tiles SpMV to fit L3 caches.
  - **Design:** Greedy vertex coloring of element dual graph (elements as nodes, edge if shared DOF), Rayon-parallel assembly within each color (no locks), ≤ 8 colors typical. Chunked SpMV partitions rows into L3-sized chunks (default 256 rows).
  - **Files:** `src/parallel_solver.rs` (add `assemble_colored`, `spmv_chunked`), `src/solvers/assembly_coloring.rs` (NEW, ~400 LoC)
  - **Tests:** `unit::coloring_valid`, `integration::colored_assembly_matches_serial`, `integration::spmv_chunked_bit_exact`, `bench::colored_assembly_scaling`
- [x] Chunk-interleaved SpMV to cut L3 contention on Zen4 / Sapphire Rapids (planned 2026-04-24; part of 4.4)

### 4.5 Benchmark & validation
- [x] `perf_bench::amg` — Poisson 3D unit cube, 128³ tets (≈ 2M DOF), CG-Jacobi vs CG-AMG iteration count + wall time (planned 2026-04-24; bundles 4.5a,b,c)
  - **Goal:** Criterion benchmark suite + in-repo convergence regression.
  - **Design:** `perf_bench::amg_poisson` (128³, PcgAmg ≤ 15 iters, PcgJacobi ≥ 200 at tol 1e-8), `perf_bench::elasticity_ibeam` (10⁶ DOF, records setup+solve breakdown), `integration::amg_convergence_regression` (32³ Poisson, ratio ≤ 0.1/cycle averaged over iters 2–10).
  - **Files:** `benches/amg_poisson.rs` (NEW), `benches/amg_elasticity.rs` (NEW), `../../oxiphysics/tests/validation_fem_amg.rs` (NEW)
- [x] Linear elasticity I-beam, 10⁶ DOF, record setup + solve breakdown (planned 2026-04-24; part of 4.5)
- [x] Convergence test: AMG residual reduction ≥ 0.1 per V-cycle on 3D Poisson (planned 2026-04-24; part of 4.5)
- [x] Compare against published AMG reference benchmarks (pure-Rust validation; no external C libraries) — completed 2026-05-11
  - Validates mesh-independence (iters bounded as h→0), 2× tolerance vs Stuben 2001 published counts, AMG≥5× speedup over PCG+Jacobi on 32³ Poisson.

### 4.6 Out of scope for Phase 4 (parked)
- Domain decomposition (Schwarz, Schur, BDDC) — separate phase
- GPU AMG — depends on `oxiphysics-gpu` Phase 4 landing first
- Adaptive AMG (Bootstrap AMG, αSA) — research-grade, not needed for v0.2.0
