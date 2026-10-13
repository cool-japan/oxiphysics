# oxiphysics-fem TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 110,420 SLoC | 4,620 tests

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Define core types and traits (CsrMatrix, TetrahedralMesh, etc.)
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation
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
- [x] Integration tests (4,620 tests, 0 stubs, 124 source files)
- [x] Performance benchmarks (basic)

### Phase 3: Polish
- [x] Documentation (rustdoc on all public items)
- [x] Extended examples (doc-tests in `parallel_solver` module)
- [x] Benchmark suite expansion (`perf_bench` — `BenchHarness`, SpMV/PCG/GMRES/assembly/Ke timing harness)
- [x] Parallel sparse solver integration (`parallel_solver` — Rayon-parallel SpMV, PCG, GMRES, assembler)
- [x] Further SIMD / parallel assembly optimization (parallel Rayon assembler + PCG dot products)

### Phase 4: Algebraic Multigrid & Large-Problem Solvers

> **Goal:** Add a real algebraic-multigrid (AMG) solver alongside the existing PCG+Jacobi, enabling 10⁶-DOF linear elasticity / Poisson problems to converge in O(N) work. Today `parallel_solver.rs` ships PCG with a diagonal preconditioner — effective at ≤10⁵ DOF but asymptotically poor. `error_estimation.rs` has `MultiGridAdaptive` for *error analysis*; it is not a solver.

#### 4.1 Ruge-Stüben classical AMG
- [x] `solvers/amg/classical.rs` — strong-connection graph, C/F splitting, direct interpolation `P` (planned 2026-04-24; bundles 4.1–4.1e)
  - **Goal:** A self-contained `AmgClassical` hierarchy builder + V/W-cycle driver capable of solving the 3D Poisson equation on 128³ tets (≈ 2 M DOF) with residual reduction ≥ 0.1 per V-cycle.
  - **Design:** Strong-connection graph (θ=0.25), Ruge-Stüben C/F splitting (two-pass), direct interpolation `P`, Galerkin `A_c = P^T A P` via SpMM, GS/SGS smoothers (2 pre/post sweeps), V-cycle + W-cycle driver, coarse-level PCG solve when DOF < 500.
  - **Files:** `src/solvers/amg/mod.rs`, `src/solvers/amg/classical.rs`, `src/solvers/amg/smoothers.rs`, `src/solvers/amg/cycle.rs`, `src/solvers/amg/graph.rs`, `src/solvers/amg/galerkin.rs`, `src/solvers/mod.rs` (re-export)
  - **Tests:** `unit::strong_connection_1d_poisson`, `unit::cf_splitting_partitions`, `unit::galerkin_triple_preserves_symmetry`, `unit::symmetric_gs_reduces_residual`, `integration::vcycle_poisson_2d_32x32` (reduction ≤ 0.1/cycle), `integration::wcycle_converges_harder_problem`
- [x] Galerkin coarse-grid operator `A_c = P^T A P` (planned 2026-04-24; part of 4.1 classical.rs)
- [x] Gauss-Seidel and symmetric GS smoothers — forward, backward, symmetric sweeps (planned 2026-04-24; part of 4.1 classical.rs)
- [x] V-cycle and W-cycle drivers — `AmgSolver::v_cycle(level, b, x)` recursion down to coarse direct solve (planned 2026-04-24; part of 4.1 classical.rs)
- [x] Coarse-level solve: fall back to existing `PcgSolver` when `n_dofs < n_coarse_threshold` (typ. ≤ 500) (planned 2026-04-24; part of 4.1 classical.rs)

#### 4.2 Smoothed-aggregation AMG (alternative path)
- [x] `solvers/amg/smoothed_aggregation.rs` — aggregate-based coarsening via strength-of-connection + greedy aggregation (planned 2026-04-24; bundles 4.2–4.2c)
  - **Goal:** Alternative coarsening path via greedy aggregation + rigid-body near-null-space, using the same V/W-cycle driver. SA typically outperforms classical on systems-of-PDEs (3D elasticity).
  - **Design:** SA strength metric `|A[i,j]|² ≥ θ²|A[i,i]A[j,j]|` (θ=0.08), two-pass greedy aggregation, tentative prolongator `P̃` from 6 rigid-body modes (3 translations + 3 rotations) per aggregate orthonormalized via QR, Jacobi-smoothed `P = (I − ω D⁻¹ A) P̃` (ω = 4/(3ρ), ρ via 20 power iterations). Cycle driver reuses `cycle.rs` and `smoothers.rs`.
  - **Files:** `src/solvers/amg/smoothed_aggregation.rs`, `src/solvers/amg/aggregation.rs`, `src/solvers/amg/near_null_space.rs`
  - **Tests:** `unit::aggregation_covers_all_nodes`, `unit::rbm_zero_residual`, `integration::sa_vcycle_elasticity_32cube` (reduction ≤ 0.15/cycle)
- [x] Tentative prolongator from near-null-space — rigid body modes for elasticity (planned 2026-04-24; part of 4.2)
- [x] Jacobi-smoothed prolongator `P = (I - ω D⁻¹ A) P̃` (planned 2026-04-24; part of 4.2)
- [x] Plug into the same V/W-cycle driver as classical AMG (planned 2026-04-24; part of 4.2)

#### 4.3 Preconditioned Krylov as outer
- [x] `PcgWithAmg` — use the AMG V-cycle as the preconditioner inside PCG (planned 2026-04-24; bundles 4.3a,b)
  - **Goal:** Wrap existing PCG + GMRES with an AMG V-cycle preconditioner; typically cuts iteration counts 10–50× vs Jacobi-PCG on 10⁶-DOF problems.
  - **Design:** `Preconditioner` trait (`fn apply(&self, r: &[f64], z: &mut [f64])`), `AmgPreconditioner` running one V-cycle, `PcgWithAmg` wrapping `ParallelPcgSolver` (line 323), `GmresWithAmg` wrapping `ParallelGmresSolver` (line 457 — confirmed present).
  - **Files:** `src/solvers/amg/preconditioner.rs` (NEW), `src/parallel_solver.rs` (add Preconditioner trait, ~100 LoC diff)
  - **Tests:** `integration::pcg_amg_poisson_64cube` (≤ 15 outer PCG iters), `integration::gmres_amg_advection_diffusion`
- [x] `GmresWithAmg` — for non-symmetric systems (advection-dominated, some coupled physics) (planned 2026-04-24; part of 4.3)

#### 4.4 Parallel assembly scale-up
- [x] Element-coloring assembly so `CsrMatrix::assemble` scales on 32+ cores (planned 2026-04-24; bundles 4.4a,b)
  - **Goal:** `CsrMatrix::assemble_colored` scales linearly to 32+ cores; `spmv_chunked` tiles SpMV to fit L3 caches.
  - **Design:** Greedy vertex coloring of element dual graph (elements as nodes, edge if shared DOF), Rayon-parallel assembly within each color (no locks), ≤ 8 colors typical. Chunked SpMV partitions rows into L3-sized chunks (default 256 rows).
  - **Files:** `src/parallel_solver.rs` (add `assemble_colored`, `spmv_chunked`), `src/solvers/assembly_coloring.rs` (NEW, ~400 LoC)
  - **Tests:** `unit::coloring_valid`, `integration::colored_assembly_matches_serial`, `integration::spmv_chunked_bit_exact`, `bench::colored_assembly_scaling`
- [x] Chunk-interleaved SpMV to cut L3 contention on Zen4 / Sapphire Rapids (planned 2026-04-24; part of 4.4)

#### 4.5 Benchmark & validation
- [x] `perf_bench::amg` — Poisson 3D unit cube, 128³ tets (≈ 2M DOF), CG-Jacobi vs CG-AMG iteration count + wall time (planned 2026-04-24; bundles 4.5a,b,c)
  - **Goal:** Criterion benchmark suite + in-repo convergence regression.
  - **Design:** `perf_bench::amg_poisson` (128³, PcgAmg ≤ 15 iters, PcgJacobi ≥ 200 at tol 1e-8), `perf_bench::elasticity_ibeam` (10⁶ DOF, records setup+solve breakdown), `integration::amg_convergence_regression` (32³ Poisson, ratio ≤ 0.1/cycle averaged over iters 2–10).
  - **Files:** `benches/amg_poisson.rs` (NEW), `benches/amg_elasticity.rs` (NEW), `../../oxiphysics/tests/validation_fem_amg.rs` (NEW)
- [x] Linear elasticity I-beam, 10⁶ DOF, record setup + solve breakdown (planned 2026-04-24; part of 4.5)
- [x] Convergence test: AMG residual reduction ≥ 0.1 per V-cycle on 3D Poisson (planned 2026-04-24; part of 4.5)
- [x] Compare against published AMG reference benchmarks (pure-Rust validation; no external C libraries) — completed 2026-05-11
  - Validates mesh-independence (iters bounded as h→0), 2× tolerance vs Stuben 2001 published counts, AMG≥5× speedup over PCG+Jacobi on 32³ Poisson.

### Verified-existing infrastructure (code-verified 2026-06-11 — basis for the roadmap below; not re-opened)
- [x] Ruge-Stüben classical AMG + smoothed-aggregation AMG (Phase 4 above)
- [x] PCG / GMRES Krylov solvers, AMG-preconditioned (`parallel_solver.rs`, `solvers/amg/`)
- [x] SUPG / PSPG / VMS fluid stabilization (`fluid_fem_stabilized.rs`)
- [x] Taylor-Hood mixed elements + `inf_sup_check` (`mixed_elements`)
- [x] `MortarContact` / `MortarContactElement` contact coupling
- [x] Zienkiewicz-Zhu patch-recovery error estimator (`error_estimator.rs`)
- [x] `GoalOrientedEstimator` — combines pre-supplied primal/adjoint vectors (the adjoint *solve* is the v0.3.0 delta)
- [x] `LanczosMethod` eigensolver (`eigenvalue.rs`)
- [x] Isogeometric analysis module, standalone (`isogeometric/` — geometry coupling is the v0.3.0 delta)
- [x] XFEM fracture + cohesive-zone + damage mechanics (Phase 2)
- [x] Homogenization module (`homogenization/` — nested FE² validation loop is the v1.0 delta)

## v0.2.0 — Solvers & saddle-point

### Matrix-free high-order operators
- [x] Matrix-free sum-factorized high-order operators (done 2026-06-12)
  - **Goal:** p=4 Poisson matvec at <40 bytes/DOF, ≥2× faster than assembled SpMV at equal DOF.
  - **Design:** tensor-product (sum-factorization) element-operator evaluation, no global matrix. Kronbichler-Kormann 2012.
  - **Delta:** MISSING — no matrix-free path today; complements the assembled CSR + AMG baseline rather than replacing it.
  - **Validation:** memory ≤ 40 B/DOF measured via allocation counter; matvec wall-time ≤ 0.5× assembled SpMV at p=4, equal DOF; result bit-comparable to assembled within 1e-12.
  - **Files (proposed):** `src/solvers/matrix_free/mod.rs`, `src/solvers/matrix_free/sum_factorization.rs`, `src/solvers/matrix_free/operator.rs`
  - **Tests (proposed):** `unit::sumfac_matches_dense_p2`, `integration::matfree_poisson_p4_matches_assembled`, `bench::matfree_vs_spmv_p4`

### Saddle-point preconditioning
- [x] Block-diagonal / Schur preconditioner + Chebyshev-accelerated smoother (shipped 2026-06-13)
    - **Goal:** Stokes GMRES iteration count ≤10% variation across 3 mesh refinements with block-Schur PC
    - **Design:** Chebyshev smoother (3-term recurrence, spectral radius via power_iteration_spectral_radius) + block-diagonal Silvester-Wathen PC (AMG velocity block + pressure mass matrix Schur)
    - **Files:** NEW src/solvers/amg/chebyshev_smoother.rs, NEW src/solvers/block_schur_pc.rs; MODIFY src/solvers/amg/mod.rs, src/solvers/mod.rs
    - **Tests:** Chebyshev smoothing factor < GS; block-Schur GMRES mesh-independence ±10%; parity vs dense-direct on small Stokes
    - **Risk:** Saddle-point GMRES coupling + pressure-Schur spectral equivalence — validate ±10% gate directly
  - **Goal:** Taylor-Hood Stokes GMRES iterations mesh-independent (±10% over 3 refinements).
  - **Design:** block-diagonal/Schur preconditioner + Chebyshev-accelerated smoother on the existing AMG.
  - **Delta:** Ruge-Stüben + SA AMG, PCG/GMRES, SUPG/PSPG/VMS, Taylor-Hood `inf_sup_check`, and `MortarContact` all ship — only the Chebyshev smoother + block-Schur wrapper are new.
  - **Validation:** GMRES iteration count for the Stokes lid-driven cavity varies ≤ ±10% across 3 uniform refinements; pressure mass-matrix Schur approximation converges.
  - **Files (proposed):** `src/solvers/amg/chebyshev_smoother.rs`, `src/solvers/block_schur_pc.rs`
  - **Tests (proposed):** `unit::chebyshev_reduces_high_frequency_error`, `integration::stokes_gmres_mesh_independent`

### Eigensolvers
- [x] LOBPCG shift-invert eigensolver
  - **Goal:** first 20 plate modes within 0.5% of analytic; ≥3× fewer factorizations than Lanczos for 20 modes.
  - **Design:** LOBPCG with shift-invert spectral transform. Knyazev 2001.
  - **Delta:** `LanczosMethod` exists; LOBPCG + shift-invert are new.
  - **Validation:** simply-supported square-plate modes vs analytic within 0.5%; factorization count ≤ ⅓ of Lanczos for 20 modes.
  - **Files (proposed):** `src/eigenvalue/lobpcg.rs`, `src/eigenvalue/shift_invert.rs`
  - **Tests (proposed):** `unit::lobpcg_rayleigh_ritz_orthonormal`, `integration::plate_modes_within_half_percent`

## v0.3.0 — Formulations

### Adaptivity
- [ ] Goal-oriented adjoint (DWR) driving hp-adaptivity
  - **Goal:** point-stress QoI error reduced 10× at ≤2× DOF vs uniform refinement.
  - **Design:** solve the dual problem + DWR weighting feeding the existing refinement machinery. Becker-Rannacher 2001.
  - **Delta:** ZZ patch recovery exists and `GoalOrientedEstimator` combines pre-supplied primal/adjoint vectors — the actual adjoint *solve* is new.
  - **Validation:** L-shaped-domain point-stress QoI error ≤ 0.1× the uniform-refinement error at ≤ 2× DOF.
  - **Files (proposed):** `src/error_estimator/dwr.rs`, `src/error_estimator/adjoint_solve.rs`
  - **Tests (proposed):** `integration::dwr_lshape_stress_qoi_10x`

### Geometry coupling
- [ ] Isogeometric analysis upgrade (geometry-coupled NURBS)
  - **Goal:** IGA plate-with-hole stress concentration within 2% of analytic at coarse mesh.
  - **Design:** trimmed/multi-patch C¹ NURBS elements consuming oxiphysics-geometry NURBS. Hughes 2005.
  - **Delta:** `isogeometric/` exists standalone; geometry coupling + multi-patch are new.
  - **Validation:** stress concentration factor at the hole within 2% of the Kirsch analytic solution on a coarse multi-patch mesh.
  - **Files (proposed):** `src/isogeometric/multi_patch.rs`, `src/isogeometric/geometry_bridge.rs`
  - **Tests (proposed):** `integration::iga_plate_with_hole_scf_within_2pct`

### Incompressible robustness
- [ ] Mixed/incompressible robustness (B-bar/F-bar + Raviart-Thomas + Nitsche)
  - **Goal:** Cook's membrane at ν=0.4999 volumetric-locking-free, tip within 2% of reference.
  - **Design:** B-bar/F-bar projection, RT elements, Nitsche interface coupling.
  - **Delta:** `mixed_elements` inf-sup + `MortarContact`/`MortarContactElement` exist; B-bar/F-bar, RT, and Nitsche are new.
  - **Validation:** Cook's membrane tip displacement at ν=0.4999 within 2% of the near-incompressible reference; no checkerboard pressure.
  - **Files (proposed):** `src/mixed_elements/bbar_fbar.rs`, `src/mixed_elements/raviart_thomas.rs`, `src/mixed_elements/nitsche.rs`
  - **Tests (proposed):** `integration::cook_membrane_locking_free_nu_0p4999`

### Fracture
- [ ] Phase-field fracture (AT1/AT2)
  - **Goal:** single-edge-notch tension load-displacement peak within 10% of reference; mesh-objective crack path.
  - **Design:** Miehe staggered scheme with irreversible history field. Miehe 2010.
  - **Delta:** XFEM + cohesive zone + damage exist; phase-field fracture is MISSING.
  - **Validation:** SENT load-displacement peak within 10% of reference; crack-path Hausdorff error insensitive to mesh refinement (objectivity).
  - **Files (proposed):** `src/fracture/phase_field/mod.rs`, `src/fracture/phase_field/staggered.rs`, `src/fracture/phase_field/history_field.rs`
  - **Tests (proposed):** `integration::sent_peak_load_within_10pct`, `integration::phase_field_crack_path_mesh_objective`

## v1.0 — Production & validation
- [ ] NAFEMS benchmark set (LE1, LE10, LE11)
  - **Goal:** target stresses/displacements within NAFEMS tolerance (2–5%).
  - **Design:** standard linear-elastic verification cases (LE1 elliptic membrane, LE10 thick plate, LE11 solid cylinder/taper/sphere).
  - **Validation:** LE1 σ_yy at point D within tolerance; LE10 σ_yy at point D; LE11 σ_zz at point A — all ≤ 5% error.
  - **Files (proposed):** `../../oxiphysics/tests/validation_fem_nafems.rs`
- [ ] Cook's membrane + Taylor bar
  - **Goal:** Cook tip converges to 23.9 ±1%; Taylor-bar final length within 3% of experiment.
  - **Design:** tapered-panel convergence study + impact elastoplasticity.
  - **Validation:** Cook tip vertical displacement → 23.9 ±1% under refinement; Taylor-bar mushroomed length within 3% of the Wilkins experiment.
  - **Files (proposed):** `../../oxiphysics/tests/validation_fem_cook_taylor.rs`
- [ ] FE² computational homogenization validation
  - **Goal:** effective modulus lands within Hashin-Shtrikman bounds for a 2-phase RVE.
  - **Design:** nested macro-micro driver over the existing `homogenization/` module.
  - **Validation:** effective bulk + shear moduli of a 2-phase RVE bracketed by the HS bounds across volume fractions 0.1–0.9.
  - **Files (proposed):** `../../oxiphysics/tests/validation_fem_fe2.rs`

## Cross-cutting (all milestones)
- [ ] Validation harness as local cargo targets — all v1.0 cases run via `cargo test -p oxiphysics --test validation_fem_*` and a `scripts/` runner; no new CI workflow YAML (pure-Rust, local-only per policy).
- [ ] Convergence-regression baselines checked into the repo (iteration counts / error tables) so AMG, saddle-point, and adaptivity changes are guarded against drift.
- [ ] Criterion benchmark coverage for every new solver path (matrix-free matvec, Chebyshev smoother, LOBPCG) alongside the existing `perf_bench` harness.
- [ ] No-unwrap + no-warnings maintained across new solver/formulation modules; deterministic-mode reductions where solver results feed regression assertions.

## Implementation order & dependencies
1. **Matrix-free sum-factorized operators (v0.2.0)** — independent of AMG; complements the assembled CSR path for high-order elements.
2. **Block-Schur PC + Chebyshev smoother (v0.2.0)** — layers on the existing AMG; enables mesh-independent Stokes/saddle-point solves.
3. **LOBPCG shift-invert (v0.2.0)** — built atop the existing `LanczosMethod`; reuse the AMG V-cycle as its preconditioner.
4. **Goal-oriented adjoint / DWR (v0.3.0)** — extends ZZ recovery + `GoalOrientedEstimator` with a real dual solve.
5. **Mixed/incompressible robustness (v0.3.0)** — B-bar/F-bar/RT/Nitsche on the existing `mixed_elements`.
6. **Phase-field fracture (v0.3.0)** — staggered solver; independent of the existing XFEM/cohesive paths.
7. **IGA geometry coupling (v0.3.0)** — consumes oxiphysics-geometry NURBS (multi-patch).
8. **NAFEMS LE1/LE10/LE11 (v1.0)** — linear-elastic anchors; need only baseline elements.
9. **Cook's membrane + Taylor bar (v1.0)** — Cook depends on the v0.3.0 B-bar/F-bar work for ν → 0.5.
10. **FE² homogenization validation (v1.0)** — nested driver over the existing `homogenization/` module, reusing the AMG-preconditioned solver at the micro scale.

## References
- Matrix-free / sum-factorization — Kronbichler & Kormann, *Comput. Fluids* 63 (2012) 135–147.
- Chebyshev polynomial smoothing — Adams, Brezina, Hu & Tuminaro, *J. Comput. Phys.* 188 (2003) 593–610.
- Block / Schur preconditioning — Elman, Silvester & Wathen, *Finite Elements and Fast Iterative Solvers* (2014).
- LOBPCG — Knyazev, *SIAM J. Sci. Comput.* 23 (2001) 517–541.
- DWR adjoint adaptivity — Becker & Rannacher, *Acta Numerica* 10 (2001) 1–102.
- Isogeometric analysis — Hughes, Cottrell & Bazilevs, *Comput. Methods Appl. Mech. Eng.* 194 (2005) 4135–4195.
- F-bar method — de Souza Neto et al., *Int. J. Solids Struct.* 33 (1996) 3277–3296.
- Raviart-Thomas elements — Raviart & Thomas, *Lect. Notes Math.* 606 (1977) 292–315.
- Nitsche's method — Nitsche, *Abh. Math. Sem. Univ. Hamburg* 36 (1971) 9–15.
- Phase-field fracture — Miehe, Welschinger & Hofacker, *Int. J. Numer. Methods Eng.* 83 (2010) 1273–1311.
- NAFEMS benchmarks — NAFEMS linear-elastic verification suite (LE1, LE10, LE11).
- Cook's membrane — Cook, Malkus & Plesha, *Concepts and Applications of Finite Element Analysis*.
- Taylor bar — Wilkins & Guinan, *J. Appl. Phys.* 44 (1973) 1200–1206.
- FE² homogenization — Feyel & Chaboche, *Comput. Methods Appl. Mech. Eng.* 183 (2000) 309–330.

## Risks & open questions
- Matrix-free operators need carefully ordered tensor contractions / quadrature to beat optimized SpMV — bench early against the assembled baseline.
- Chebyshev smoothing needs robust spectral-radius bounds across element types — auto-estimate the eigenvalue interval per level.
- LOBPCG can stagnate without a good preconditioner — reuse the AMG V-cycle as the LOBPCG preconditioner.
- DWR requires a well-posed dual problem for non-self-adjoint operators — document the supported QoI classes.
- B-bar/F-bar volumetric projection interacts with the existing `mixed_elements` inf-sup — verify there is no double-stabilization.
- Phase-field fracture length-scale ℓ must resolve the mesh (h < ℓ/2) — add a mesh-adequacy check before solving.
- FE² nested solves are expensive — reuse the AMG-preconditioned solver at the micro scale and cache micro tangents.

## Deferred / research track
- [~] Domain decomposition (Schwarz, Schur complement, BDDC) — Ready when: distributed/multi-node execution is in scope (single-node AMG covers current 10⁶-DOF targets).
- [~] GPU AMG — Ready when: `oxiphysics-gpu` zero-readback compute path lands (shared ecosystem work).
- [~] Adaptive AMG (Bootstrap AMG, αSA) — Ready when: a benchmark shows classical/SA AMG iteration counts degrading on a target problem class; research-grade, not needed for v0.2.0.

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 363 (critical 8, high 122, medium 177, low 56); 57 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (8)
- [x] **F1037** (critical, bug, easy/small) `analysis.rs:724` — Newmark-beta solves displacement-form K_eff but treats result as acceleration
- [x] **F1038** (critical, bug, medium/small) `analysis.rs:788` — Wilson-theta solves displacement system, uses result as acceleration
- [ ] **F1140** (critical, bug, medium/small) `dynamics.rs:151` — HHT-α parameters sign-flipped (γ<½, energy growth), not HHT
- [ ] **F1142** (critical, bug, medium/medium) `dynamics.rs:193` — NewmarkIntegrator::step is explicit yet claims unconditional stability
- [ ] **F1152** (critical, bug, medium/medium) `eigenvalue.rs:382` — solve_smallest/ModalAnalysis return k-step Lanczos Ritz values, not smallest modes
- [ ] **F1174** (critical, bug, medium/small) `explicit_fem/types.rs:135` — 3D node-to-segment gap is unsigned, so contact force is always zero
- [x] **F1346** (critical, bug, medium/medium) `solvers/functions.rs:378` — MINRES recurrence wrong; returns garbage flagged converged=true
- [ ] **F1378** (critical, bug, easy/small) `thermal_fem.rs:163` — Tet shape gradients use J^-T rows; fails linear patch test

### High (122)
- [x] **F1024** (high, bug, medium/medium) `adaptive_mesh/functions.rs:565` — coarsen_mesh deletes leaf elements, punching holes in mesh
- [x] **F1028** (high, bug, easy/small) `adaptive_mesh/types.rs:481` — refine_element duplicates shared edge-midpoint nodes between neighbours
- [x] **F1032** (high, api, medium/medium) `analysis.rs:94` — LinearStaticAnalysis::solve infallible; hides CG failure and NaN
- [ ] **F1034** (high, bug, medium/medium) `analysis.rs:429` — Modal RQI misses lowest modes; mode shape fabricated on singular shift
- [ ] **F1036** (high, bug, medium/medium) `analysis.rs:690` — Newmark/Wilson solve wrong effective system; acceleration scaled by βΔt²
- [x] **F1042** (high, perf, easy/small) `assembly.rs:145` — apply_dirichlet_bc sweeps whole matrix per fixed DOF
- [x] **F1048** (high, bug, easy/trivial) `beam_fem.rs:89` — Rectangular torsion constant overestimated: 0.21 term divided by 3
- [x] **F1054** (high, perf, easy/small) `boundary.rs:885` — Dirichlet application scans whole CSR per BC: O(n_bc·nnz)
- [ ] **F1055** (high, bug, easy/small) `boundary_element.rs:99` — Laplace double-layer kernel sign flipped; BEM flux and interior u negated
- [ ] **F1058** (high, bug, medium/small) `boundary_element.rs:650` — Elastostatic BEM self-G 3.6–4.6× too small; rigid-body diag incomplete
- [x] **F1063** (high, bug, easy/small) `boundary_element.rs:1026` — BEM-FEM condense: O(n^4) inverse and wrong H·G⁻¹ order
- [ ] **F1069** (high, bug, medium/small) `buckling/types.rs:88` — Higher buckling modes re-find mode 1 (shift = prev*1.01)
- [ ] **F1071** (high, fabrication, medium/small) `cohesive/types.rs:145` — PPR cohesive law ad hoc: fracture energy ≠ Γn (off by ~1e5)
- [ ] **F1072** (high, bug, easy/small) `cohesive/types.rs:366` — Xu-Needleman normal traction nonzero at zero opening
- [ ] **F1073** (high, bug, easy/small) `composite_fem.rs:330` — Ply strain rotated with stress transform; shear term doubled
- [ ] **F1074** (high, fabrication, medium/small) `composite_fem.rs:454` — first_ply_failure fabricates ply strengths from elastic moduli
- [ ] **F1077** (high, bug, medium/small) `constitutive.rs:317` — Mooney-Rivlin stress is nonzero at zero strain (no pressure/bulk term)
- [ ] **F1080** (high, bug, medium/medium) `contact/functions.rs:443` — Node-face contact reports nodes against their own faces
- [ ] **F1100** (high, bug, easy/small) `contact_mechanics.rs:231` — Mindlin G* convention mixed: k_t 4× too stiff, δt 2× too small
- [ ] **F1102** (high, bug, easy/trivial) `contact_mechanics.rs:368` — JKR contact radius missing 3/4 factor; elastic_energy delta wrong
- [x] **F1105** (high, fabrication, hard/large) `coupled_fem/types.rs:911` — ThermoMechanicalCoupling mechanics is fabricated displacement accumulation
- [x] **F1109** (high, bug, easy/small) `coupled_fem/types.rs:1793` — AcousticCavity1D 'Lax-Wendroff' step is unconditionally unstable FTCS
- [ ] **F1112** (high, bug, medium/small) `coupled_physics_fem.rs:693` — CoupledSolver: sign mismatch between schemes; fake converged/zero residuals
- [ ] **F1114** (high, bug, medium/small) `coupled_physics_fem.rs:873` — CoupledSolver: CG on nonsymmetric block, sign mismatch, hard-coded convergence
- [ ] **F1115** (high, bug, medium/medium) `crystal_plasticity.rs:367` — HCP: zero pyramidal systems created; prismatic slip uses c-axis
- [ ] **F1117** (high, fabrication, easy/trivial) `damage/types.rs:212` — GursonModel::yield_function always returns 0.0
- [ ] **F1126** (high, bug, medium/medium) `discontinuous_galerkin.rs:823` — DG diffusion volume term anti-diffusive; interior-penalty face terms discarded
- [ ] **F1129** (high, bug, easy/trivial) `dynamic/functions.rs:826` — modal_newmark_1dof a4=γ/(β²Δt) makes damped response diverge
- [ ] **F1132** (high, fabrication, medium/small) `dynamic/types.rs:80` — Bathe step_diagonal sub-step 2 is not 3-point backward Euler
- [ ] **F1134** (high, bug, easy/small) `dynamic/types.rs:397` — Wilson-θ uses 3·a1 instead of a2=2·a1 on velocity terms
- [ ] **F1137** (high, fabrication, easy/small) `dynamic/types.rs:806` — DynamicSolver::generalized_alpha: diagonal-only, double-counted fake solve
- [ ] **F1138** (high, bug, easy/small) `dynamic/types.rs:826` — DynamicSolver::generalized_alpha diagonal-only and double-counts predictor terms
- [ ] **F1141** (high, fabrication, medium/medium) `dynamics.rs:161` — Newmark 'implicit' step is explicit; critical_dt claims unconditional stability
- [ ] **F1143** (high, fabrication, medium/medium) `dynamics.rs:402` — BatheIntegrator::step uses an invented 0.9/0.1 blend instead of Bathe–Noh
- [ ] **F1144** (high, bug, easy/trivial) `dynamics.rs:452` — lump_mass_hrz divides total by n_dim; lumped mass n_dim× too small
- [ ] **F1146** (high, bug, easy/trivial) `dynamics.rs:701` — estimate_omega_max uses K_ii/M_ii, so the 'safe' explicit dt is unstable
- [ ] **F1148** (high, bug, medium/small) `eigen/functions.rs:265` — subspace_iteration finds largest eigenpairs and never re-orthonormalizes
- [ ] **F1149** (high, bug, medium/small) `eigen/functions.rs:455` — shift_invert_power returns num_modes copies of one eigenpair
- [ ] **F1155** (high, fabrication, medium/small) `eigenvalue.rs:780` — FrequencyResponse ignores damping_ratio; LU has no pivoting
- [ ] **F1158** (high, bug, medium/small) `electrochemical_fem.rs:148` — Nernst-Planck migration sign wrong, upwind weight inverted, velocity ignored
- [ ] **F1159** (high, bug, easy/small) `electrochemical_fem.rs:173` — Nernst-Planck LHS/RHS mass matrices mismatch; concentration grows ~50%/step
- [ ] **F1160** (high, fabrication, hard/large) `electrochemical_fem.rs:797` — BatteryFem 'DFN' step only does solid diffusion; surface BC wrong
- [ ] **F1162** (high, bug, easy/small) `electrochemical_fem.rs:1191` — PEM O2 'zero-flux' membrane BC forces c=0; ORR sink mis-scaled
- [ ] **F1166** (high, bug, easy/small) `electromechanics/types.rs:1056` — GlobalPiezoAssembler::apply_dirichlet ignores prescribed value
- [ ] **F1172** (high, bug, medium/medium) `error_estimation.rs:578` — refine_marked connectivity disagrees with bisected coords; non-conforming mesh
- [ ] **F1176** (high, bug, easy/small) `explicit_fem/types.rs:222` — Contact Lagrange multipliers clamped to zero forever; Lagrange gives zero force
- [ ] **F1177** (high, fabrication, medium/medium) `explicit_fem/types.rs:829` — ConWep overpressure fit is non-physical and does not decay
- [ ] **F1179** (high, bug, easy/small) `fatigue_fem.rs:441` — Rainflow ignores ASTM start-point half-cycle rule and drops plateau peaks
- [ ] **F1183** (high, bug, easy/trivial) `fluid_fem/fluid_fem_stokes.rs:442` — 1-D NS viscous term uses FD scaling with FEM-integrated mass
- [ ] **F1191** (high, bug, easy/trivial) `fracture/types.rs:528` — SENB K_I uses W² instead of W^{3/2}
- [ ] **F1194** (high, bug, easy/trivial) `homogenization/advanced.rs:759` — Levin thermal expansion wrongly multiplied by volume fraction
- [ ] **F1195** (high, bug, easy/trivial) `homogenization/bounds.rs:300` — Mori-Tanaka Eshelby shear entries missing Voigt factor 2
- [ ] **F1196** (high, fabrication, medium/small) `homogenization/bounds.rs:333` — eshelby_tensor_prolate returns sphere tensor for any aspect ratio
- [ ] **F1205** (high, bug, easy/small) `hyperelastic/types.rs:231` — Yeoh/Gent/Arruda-Boyce carry nonzero stress in undeformed state
- [ ] **F1212** (high, fabrication, hard/medium) `isogeometric/types.rs:541` — DegreeElevation uses curve samples as control points
- [ ] **F1215** (high, bug, easy/trivial) `isogeometric/types.rs:1814` — p_refine_knot makes every interior knot C^-1
- [ ] **F1216** (high, fabrication, hard/medium) `isogeometric/types.rs:1901` — Mortar coupling matrix ignores spans, sizes, and interface location
- [ ] **F1218** (high, bug, easy/trivial) `isogeometric_fem.rs:797` — IGA B-matrix uses J⁻¹ instead of J⁻ᵀ for physical derivatives
- [ ] **F1219** (high, perf, medium/large) `isogeometric_fem.rs:825` — IGA element matrices span all control points; O(N³) assembly
- [ ] **F1220** (high, perf, medium/medium) `isogeometric_fem.rs:947` — IGA elements span all control points: O(N³) dense assembly
- [ ] **F1223** (high, bug, easy/small) `lbm_coupling.rs:243` — transfer_velocity_to_lbm double-weights delta kernel; imposes ω·v
- [ ] **F1225** (high, fabrication, medium/medium) `lbm_coupling.rs:354` — IBM pressure force direction = sign of slip velocity, not normal
- [ ] **F1226** (high, bug, easy/small) `level_set.rs:60` — Reinitialization Godunov gradient omits grid spacing
- [ ] **F1227** (high, bug, easy/trivial) `level_set.rs:302` — Level-set curvature cross-derivative term halved
- [ ] **F1231** (high, bug, easy/trivial) `mesh.rs:306` — Tet dihedral angle returns π−φ instead of φ
- [ ] **F1237** (high, fabrication, hard/medium) `meshfree_fem/functions.rs:689` — apply_nitsche_bc is identical to the penalty method
- [ ] **F1239** (high, fabrication, medium/medium) `meshfree_fem/types.rs:154` — Bridging-domain weights no-op; coupling matrix always zero
- [ ] **F1241** (high, bug, medium/small) `meshless.rs:35` — MLS/RKPM silently return zero shape functions on singular moment matrix
- [ ] **F1247** (high, bug, medium/medium) `modal.rs:449` — harmonic_analysis: 5 fixed block-GS sweeps with CG on indefinite K−ω²M
- [ ] **F1248** (high, bug, medium/medium) `modal.rs:456` — harmonic_analysis sign error and CG on indefinite K−ω²M
- [ ] **F1252** (high, bug, medium/small) `modal.rs:1460` — Response spectrum double-steps SDOF and cancels ground forcing
- [ ] **F1253** (high, fabrication, medium/medium) `multiphysics_fem.rs:224` — Thermal load vector is nodal D·ε_th, not ∫BᵀDε_th dV
- [ ] **F1254** (high, bug, easy/trivial) `multiphysics_fem.rs:417` — PZT-5A piezo coupling entries in wrong matrix positions
- [x] **F1260** (high, fabrication, hard/large) `multiscale_fem.rs:305` — FE² MicroSolver solves no RVE; returns Voigt rule-of-mixtures
- [x] **F1262** (high, fabrication, medium/medium) `multiscale_fem.rs:928` — Craig-Bampton SubstructureMethod is entirely synthetic data
- [ ] **F1266** (high, fabrication, medium/small) `nonlinear/functions.rs:991` — Pseudo-arclength corrector and Riks predictor are heuristic placeholders
- [ ] **F1269** (high, fabrication, medium/small) `nonlinear_dynamics.rs:164` — compute_mle is mean log trajectory speed, not a Lyapunov exponent
- [ ] **F1274** (high, bug, medium/small) `nonlinear_fem.rs:373` — J2 return mapping double-counts plastic strain and breaks kinematic hardening
- [ ] **F1275** (high, bug, easy/trivial) `nonlinear_fem.rs:398` — J2 return mapping accumulates equivalent plastic strain twice
- [ ] **F1276** (high, bug, easy/small) `nonlinear_fem.rs:572` — Newton linear solve divides by 1e-30 for negative pivots
- [ ] **F1277** (high, bug, medium/small) `nonlinear_fem.rs:775` — RiksSolver uses frozen-k residual and clamps negative stiffness
- [ ] **F1280** (high, bug, medium/small) `nonlinear_fem.rs:1075` — B-bar hex shape-derivative sign table mismatches node ordering
- [ ] **F1284** (high, fabrication, medium/medium) `nonlinear_solver.rs:450` — LoadStepper scales whole residual; arc-length reference load degenerate
- [ ] **F1285** (high, fabrication, medium/medium) `nonlinear_solver.rs:520` — ArcLengthMethod corrector degenerate: du_r = −du_t, no load vector
- [x] **F1303** (high, fabrication, hard/large) `poromechanics.rs:194` — ConsolidationSolver allocates unused dense K_uu; no u-p coupling
- [ ] **F1304** (high, fabrication, medium/medium) `poromechanics.rs:214` — ConsolidationSolver::step ignores K_uu/K_up; displacements never solved
- [x] **F1305** (high, bug, easy/trivial) `poromechanics.rs:286` — Terzaghi time factor uses ln instead of log10
- [ ] **F1308** (high, bug, easy/trivial) `probabilistic_fem/types.rs:162` — MonteCarloFem LCG uniform covers only [0, 0.5)
- [ ] **F1309** (high, bug, medium/small) `probabilistic_fem/types.rs:340` — PCE 'tensor quadrature' evaluates only the 3 diagonal points
- [ ] **F1313** (high, bug, easy/small) `reduced_order.rs:546` — ROM Newmark step: zero acceleration, wrong effective RHS
- [ ] **F1314** (high, fabrication, hard/medium) `reduced_order.rs:978` — Balanced truncation uses fabricated Gramians and wrong projection
- [x] **F1316** (high, fabrication, hard/medium) `reduced_order_fem.rs:755` — balanced_truncation never uses A; Gramians are BBᵀ and CᵀC
- [ ] **F1323** (high, bug, medium/small) `reliability_fem/types.rs:478` — SORM curvatures mix x/u spaces and skip the rotation
- [x] **F1327** (high, fabrication, medium/medium) `shell/types.rs:347` — RectPlateElement stiffness is uncoupled diagonal, not a plate element
- [ ] **F1328** (high, fabrication, medium/small) `shell/types.rs:353` — RectPlateElement stiffness is diagonal-only, no rigid-body modes
- [ ] **F1329** (high, fabrication, hard/medium) `shell/types.rs:465` — 'DKT' flat-shell bending gives zero stiffness to w DOFs
- [ ] **F1330** (high, bug, easy/trivial) `shell_fem.rs:604` — Cylinder external-pressure buckling formulas overestimate (unconservative)
- [x] **F1332** (high, bug, easy/small) `soil_fem.rs:160` — principal_stresses drops shear; Mohr-Coulomb never sees shear
- [ ] **F1333** (high, bug, easy/small) `soil_fem.rs:176` — Mohr-Coulomb yield uses wrong friction sign for tension-positive stress
- [ ] **F1339** (high, panic, hard/medium) `solvers/amg/smoothed_aggregation.rs:193` — SA-AMG elasticity path indexes node coords by DOF (panic)
- [x] **F1347** (high, fabrication, medium/small) `solvers/functions.rs:675` — Dense V-cycle overwrites pre-smoothing and coarse correction
- [x] **F1348** (high, perf, medium/small) `solvers/functions.rs:1176` — amg_pcg rebuilds AMG hierarchy every iteration; non-SPD preconditioner
- [ ] **F1356** (high, perf, medium/medium) `sparse/types.rs:61` — ICC(0) factorization is O(n²) and breakdown clamp explodes factor
- [x] **F1361** (high, bug, easy/trivial) `spectral_fem/types.rs:456` — SEM load_rhs divides out Jacobian; solution scaled by 1/J
- [x] **F1371** (high, bug, easy/small) `stochastic_fem/types.rs:808` — PCE Gauss-Hermite projection mis-scaled; 5-point weights wrong too
- [ ] **F1372** (high, fabrication, easy/small) `thermal/functions.rs:108` — stefan_neumann_lambda ignores ΔT; Stefan number is cp/L
- [ ] **F1375** (high, bug, easy/trivial) `thermal/functions.rs:521` — enthalpy_update divides heat input by density twice
- [ ] **F1376** (high, bug, easy/trivial) `thermal/types.rs:358` — 1-D explicit thermal step double-counts convective loss
- [ ] **F1380** (high, fabrication, easy/small) `thermal_fem.rs:1108` — θ-method implicit/CN step drops off-diagonal K (diagonal approximation)
- [ ] **F1381** (high, bug, easy/trivial) `thermal_fem.rs:1187` — Linearised radiation BC RHS has wrong sign/terms
- [ ] **F1382** (high, fabrication, easy/small) `thermal_fem.rs:1261` — Tet thermal nodal load split equally; nonzero resultant, wrong units
- [ ] **F1384** (high, fabrication, easy/small) `thermal_stress.rs:517` — Thermal stresses ignore displacements, reference temp; TwoWay unused
- [ ] **F1397** (high, perf, medium/medium) `topology_optimization_ext.rs:91` — All-pairs O(n^2) topology filters; dense n×n filter matrix
- [ ] **F1398** (high, bug, medium/medium) `topology_optimization_ext.rs:311` — MMA subproblem swaps L/U, returns NaN for negative sensitivities
- [ ] **F1407** (high, bug, easy/trivial) `viscoelastic_fem/types.rs:366` — SLS relaxation time formula wrong (τ_r longer than τ_c)
- [ ] **F1408** (high, bug, medium/small) `viscoelastic_fem/types.rs:519` — Viscoelastic bar uses total displacement as strain increment
- [ ] **F1409** (high, bug, easy/small) `viscoelastic_fem/types.rs:865` — PronyFemElement stress uses strain increment, not total strain
- [ ] **F1411** (high, fabrication, medium/small) `wave_propagation.rs:366` — lamb_wave_symmetric residual identically zero; also antisymmetric form
- [ ] **F1412** (high, bug, medium/small) `xfem/functions.rs:16` — Heaviside and LevelSet phi use only last segment's infinite line
- [ ] **F1418** (high, bug, medium/small) `xfem/functions.rs:856` — sub_triangulate_triangle: crossed polygons and ignores interface through vertex
- [ ] **F1427** (high, fabrication, easy/small) `xfem/types.rs:107` — ki_from_displacement drops elastic factor, wrong K_I units
- [ ] **F1429** (high, fabrication, medium/medium) `xfem/types.rs:153` — needs_heaviside enriches every node behind tip, ignores phi

### Medium (177)
- [ ] **F1019** (medium, bug, medium/medium) `adaptive_fem.rs:444` — Longest-edge bisection leaves hanging nodes (no conformity closure)
- [ ] **F1022** (medium, panic, easy/trivial) `adaptive_fem.rs:1096` — estimated_remaining_iters divides by log2(1)=0 then overflows
- [ ] **F1026** (medium, bug, medium/medium) `adaptive_mesh/types.rs:51` — Tet refine duplicates midpoint nodes, never splits neighbours, O(n²)
- [ ] **F1027** (medium, bug, easy/small) `adaptive_mesh/types.rs:97` — AdaptiveMesh::smooth moves boundary nodes, shrinking the domain
- [ ] **F1030** (medium, fabrication, medium/medium) `additive_manufacturing_fem.rs:321` — ThermalHistoryFem::step has no heat conduction; conductivity unused
- [ ] **F1031** (medium, bug, easy/trivial) `additive_manufacturing_fem.rs:627` — Rosenthal melt-pool depth is dimensionally wrong (uses k, not ρc_p)
- [ ] **F1035** (medium, bug, medium/medium) `analysis.rs:651` — Duplicated dense solvers mishandle singularity: silent NaN/zero, public panic
- [ ] **F1039** (medium, bug, easy/small) `analysis.rs:985` — Basquin life: hidden σ_ref and n_ref prefactor give N_ref²
- [ ] **F1043** (medium, perf, easy/small) `assembly.rs:165` — apply_dirichlet_bc scans the whole CSR per fixed DOF
- [ ] **F1044** (medium, panic, easy/small) `assembly.rs:292` — cholesky_solve panics on non-SPD/unconstrained stiffness
- [ ] **F1045** (medium, bug, easy/trivial) `assembly.rs:686` — assemble_surface_traction doc promises area weighting it lacks
- [ ] **F1047** (medium, test, easy/small) `assembly.rs:1262` — No live multi-element FEM patch test; existing one checks symmetry only
- [ ] **F1049** (medium, bug, easy/trivial) `beam_fem.rs:396` — 3-D beam xz-bending coupling signs violate right-hand θy
- [ ] **F1050** (medium, fabrication, easy/small) `biomechanics_fem.rs:187` — Static optimisation ignores n_exp; always equal activation
- [ ] **F1051** (medium, bug, easy/small) `biomechanics_fem.rs:440` — Ligament toe stress jumps at ε=0; artery inner displacement formula wrong
- [ ] **F1056** (medium, fabrication, easy/small) `boundary_element.rs:439` — TellesTransform is not Telles: γ discarded, Jacobian inconsistent
- [ ] **F1057** (medium, stub, medium/small) `boundary_element.rs:603` — BEM assembly builds a Gauss rule then ignores it
- [ ] **F1059** (medium, bug, medium/small) `boundary_element.rs:673` — Elastostatic BEM rigid-body fill and self-G block are wrong
- [ ] **F1060** (medium, bug, medium/small) `boundary_element.rs:714` — Helmholtz BEM reuses Laplace rigid-body trick; self-G imaginary part missing
- [ ] **F1061** (medium, bug, easy/small) `boundary_element.rs:830` — BEM gauss_solve skips singular pivots and returns zeros silently
- [ ] **F1062** (medium, bug, easy/trivial) `boundary_element.rs:948` — DualBem::compute_sif is 4x too large and uses crack length as r
- [ ] **F1064** (medium, bug, medium/small) `boundary_element.rs:1037` — BEM-FEM condensation computes H·G⁻¹ instead of G⁻¹·H
- [ ] **F1065** (medium, panic, medium/small) `boundary_element.rs:1188` — Octree subdivision recurses forever on coincident points; leaves use monopole
- [ ] **F1066** (medium, bug, easy/small) `boundary_element.rs:1246` — Barnes-Hut uses a monopole even for near leaves; error unbounded
- [ ] **F1068** (medium, bug, easy/trivial) `buckling/functions.rs:619` — constrained_beam_thermal_buckling hardcodes area = 1 m²
- [ ] **F1070** (medium, bug, easy/small) `cohesive/functions.rs:23` — Cohesive stiffness halves interface area and hardwires normal to global x
- [ ] **F1075** (medium, api, easy/small) `constitutive.rs:55` — LinearElasticMaterial::new accepts ν=0.5, E≤0, NaN
- [ ] **F1076** (medium, bug, hard/medium) `constitutive.rs:251` — NeoHookean/J2 tangents inconsistent with stress(); J2 has no history
- [ ] **F1078** (medium, panic, easy/small) `constitutive.rs:488` — NeoHookean stress panics on collapsed/degenerate deformation
- [ ] **F1081** (medium, bug, medium/small) `contact/functions.rs:467` — detect_node_face_contact pairs every face with its own vertices
- [ ] **F1082** (medium, fabrication, medium/medium) `contact/functions.rs:691` — hertz_subsurface_stress is an invented, ν-free formula
- [ ] **F1085** (medium, bug, medium/small) `contact/types.rs:499` — EllipticalHertz approach factor 0.5 wrong; circular limit fails
- [ ] **F1088** (medium, perf, easy/small) `contact_fem/functions.rs:325` — Contact stiffness assembled into dense n_dofs×n_dofs Vec
- [ ] **F1089** (medium, fabrication, medium/small) `contact_fem/types.rs:53` — Mortar gap integral lacks shape-function weights; residual is a scalar
- [ ] **F1091** (medium, bug, easy/trivial) `contact_fem/types.rs:154` — Augmented-Lagrangian multiplier never released when gap opens
- [ ] **F1094** (medium, fabrication, medium/small) `contact_mech/types.rs:187` — Maugis-Dugdale λ, pull-off and zero-load radii wrong or heuristic
- [ ] **F1095** (medium, bug, easy/trivial) `contact_mech/types.rs:359` — Hertz subsurface max shear formula wrong: 0.785p0 at the surface
- [ ] **F1096** (medium, bug, easy/trivial) `contact_mech/types.rs:387` — Hertz surface von Mises: wrong σr at r=0 and wrong σθ
- [ ] **F1097** (medium, bug, easy/small) `contact_mech/types.rs:408` — max_subsurface_von_mises returns max shear 0.31p0, not 0.62p0
- [ ] **F1099** (medium, bug, easy/small) `contact_mech/types.rs:776` — NodeToSegment gap unsigned; segment-segment distance wrong
- [ ] **F1101** (medium, bug, medium/small) `contact_mechanics.rs:257` — Mindlin k_t, compliance and fretting energy are mutually inconsistent
- [ ] **F1104** (medium, fabrication, hard/large) `coupled_fem/types.rs:567` — AleFsiSolver: hardcoded hydrostatic pressure, no fluid solve or coupling iteration
- [ ] **F1106** (medium, panic, easy/small) `coupled_fem/types.rs:1027` — ConformingTransfer panics for n_src<2; dense weights; non-conservative
- [ ] **F1107** (medium, panic, easy/small) `coupled_fem/types.rs:1316` — Time-step functions assert!-panic on user dt
- [ ] **F1108** (medium, fabrication, medium/small) `coupled_fem/types.rs:1435` — PodRom::from_snapshots is Gram-Schmidt of first n snapshots, not POD
- [ ] **F1111** (medium, bug, easy/trivial) `coupled_physics_fem.rs:134` — Thermoelastic heat source has wrong sign: expansion heats instead of cools
- [ ] **F1116** (medium, fabrication, hard/medium) `crystal_plasticity.rs:612` — PolycrystalFem::taylor_factor computes |l_x·l_y|, not a Taylor factor
- [ ] **F1118** (medium, bug, easy/trivial) `damage/types.rs:365` — LemaitreCDM energy release rate uses S instead of E
- [ ] **F1119** (medium, bug, medium/small) `damage/types.rs:633` — Damage-plasticity return map ignores back stress and kinematic modulus
- [ ] **F1121** (medium, bug, easy/small) `damage_mechanics.rs:470` — Voce isotropic hardening restarts every increment, never saturates
- [ ] **F1123** (medium, bug, easy/trivial) `discontinuous_galerkin.rs:56` — legendre_poly_deriv returns the wrong sign at x = −1
- [ ] **F1124** (medium, bug, easy/trivial) `discontinuous_galerkin.rs:576` — WENO3 candidate stencils wrong; linear data not reconstructed
- [ ] **F1127** (medium, fabrication, medium/small) `discontinuous_galerkin.rs:978` — euler_roe_flux is a Rusanov flux, not a Roe flux
- [ ] **F1128** (medium, bug, medium/small) `dynamic/functions.rs:340` — Tet critical time step ~2.2x unconservative (bar speed, edge length)
- [ ] **F1130** (medium, bug, easy/small) `dynamic/functions.rs:1538` — modal_superposition projects initial conditions without mass matrix
- [ ] **F1131** (medium, bug, easy/trivial) `dynamic/functions.rs:1731` — cqc_rho numerator uses (ζi+ζj) instead of (ζi + r·ζj)
- [ ] **F1135** (medium, perf, medium/small) `dynamic/types.rs:664` — Newmark step rebuilds K_eff and BC-scans every step
- [ ] **F1136** (medium, perf, medium/medium) `dynamic/types.rs:678` — Implicit integrators rebuild K_eff from HashMap triplets every step
- [ ] **F1145** (medium, bug, easy/small) `dynamics.rs:604` — EnergyBalance damping work halves v_curr; no initial-energy baseline
- [ ] **F1150** (medium, bug, easy/small) `eigen/functions.rs:618` — Unpivoted LU with 1e-12 nudge used on indefinite shifted matrices
- [ ] **F1154** (medium, fabrication, medium/small) `eigenvalue.rs:763` — FrequencyResponse ignores damping_ratio; fails at resonance
- [ ] **F1156** (medium, bug, easy/small) `eigenvalue.rs:787` — FrequencyResponse stores damping_ratio but never uses it
- [ ] **F1164** (medium, fabrication, easy/small) `electromagnetic.rs:181` — EddyCurrent::quality_factor is identically 1 for all inputs
- [ ] **F1167** (medium, bug, easy/trivial) `element.rs:607` — tet_radius_ratio uses max_edge/2 as circumradius and clamps to 1
- [ ] **F1168** (medium, bug, easy/small) `element.rs:779` — spr_recovery returns exactly zero stress for rank-deficient patches
- [ ] **F1169** (medium, api, easy/medium) `error.rs:10` — Crate Error is a String bag; four disjoint error enums
- [ ] **F1171** (medium, perf, easy/small) `error_estimation.rs:271` — ZZ element indicator scans all node patches per element: O(N²)
- [ ] **F1178** (medium, bug, easy/trivial) `explicit_fem/types.rs:905` — Blast reflected-pressure coefficient grows without bound
- [ ] **F1181** (medium, bug, easy/trivial) `fluid_fem/fluid_fem_stabilized.rs:283` — VMS tau_C divides by density, giving wrong units
- [ ] **F1182** (medium, fabrication, medium/medium) `fluid_fem/fluid_fem_stokes.rs:134` — Stokes module has no saddle-point solve; cavity matrix is dense
- [ ] **F1186** (medium, perf, medium/medium) `fluid_structure.rs:710` — MonolithicFsi::assemble_tangent supports dense blocks only
- [ ] **F1187** (medium, fabrication, medium/small) `fracture/functions.rs:127` — Strain-energy-density criterion is a copy of the energy-release criterion
- [ ] **F1188** (medium, fabrication, medium/small) `fracture/functions.rs:176` — NASGRO crack-opening function collapses to R, removing closure
- [ ] **F1189** (medium, bug, easy/small) `fracture/functions.rs:1076` — xfem_domain_integral strain energy built from wrong quantities
- [ ] **F1190** (medium, bug, medium/small) `fracture/types.rs:270` — Williams higher-order terms violate traction-free crack faces
- [ ] **F1193** (medium, bug, easy/trivial) `homogenization/advanced.rs:235` — HS Young's modulus 'bounds' built from Hill-average shear modulus
- [ ] **F1197** (medium, fabrication, medium/medium) `homogenization/extended.rs:78` — 2-D asymptotic homogenization ignores microstructure geometry
- [ ] **F1198** (medium, fabrication, medium/small) `homogenization/extended.rs:260` — Fiber composite nu23 uses an invented formula hidden by a clamp
- [ ] **F1200** (medium, bug, medium/small) `homogenization/extended.rs:563` — Differential-scheme ODE factor is identically 1 (uses K, not K_i)
- [ ] **F1201** (medium, bug, easy/small) `homogenization/extended.rs:669` — Puck mode-A index lacks (1−p·R/S)² factor; no compressive branch
- [ ] **F1202** (medium, fabrication, hard/large) `homogenization/schemes.rs:445` — Computational homogenization driver just returns the Hill average
- [ ] **F1206** (medium, fabrication, medium/small) `hyperelastic/types.rs:430` — Varga model energy is a √(3Ī1) proxy, not Σλ_i
- [ ] **F1207** (medium, bug, easy/small) `inverse_problems.rs:297` — Gauss-Newton uses absolute FD step and fixed damping; silent stagnation
- [ ] **F1211** (medium, fabrication, hard/large) `isogeometric/types.rs:229` — IgaAdaptive 'ZZ estimator' is a geometric second difference
- [ ] **F1214** (medium, missing-feature, medium/small) `isogeometric/types.rs:850` — IGA contact gap uses 11×11 grid search, no closest-point projection
- [ ] **F1222** (medium, stub, medium/medium) `lbm_coupling.rs:28` — Dead config knobs: sub_iterations, interpolation_order, lumping, bulk_viscosity
- [ ] **F1232** (medium, perf, medium/medium) `mesh.rs:614` — coarsen_by_quality clones all nodes per element and punches holes
- [ ] **F1234** (medium, perf, medium/medium) `meshfree_fem/functions.rs:415` — EFG assembles dense ndof² matrix, solved by dense elimination
- [ ] **F1235** (medium, fabrication, medium/medium) `meshfree_fem/functions.rs:602` — generate_voronoi_cells is neither Voronoi nor a valid tiling
- [ ] **F1236** (medium, bug, medium/small) `meshfree_fem/functions.rs:639` — EFG BCs constrain MLS coefficients; 'Nitsche' is just penalty
- [ ] **F1242** (medium, perf, medium/medium) `meshless.rs:377` — efg_stiffness: O(n²) per Gauss point, dense matrix, FD gradients
- [ ] **F1243** (medium, fabrication, easy/small) `mixed_elements.rs:119` — inf_sup_check returns Gershgorin bound, not inf-sup constant
- [ ] **F1244** (medium, panic, easy/small) `modal.rs:76` — Modal eigen solvers panic or fabricate modes on degenerate input
- [ ] **F1249** (medium, bug, easy/small) `modal.rs:562` — Modal superposition critically/overdamped branch is wrong
- [ ] **F1251** (medium, bug, easy/small) `modal.rs:1315` — Modal truncation uses extracted-mode mass as total; static_correction ignored
- [ ] **F1255** (medium, bug, easy/trivial) `multiphysics_fem.rs:535` — Skempton B formula wrong when Biot coefficient α≠1
- [ ] **F1256** (medium, bug, easy/small) `multiphysics_fem.rs:640` — ALE mesh smoothing uses structural node 0 for every fluid node
- [ ] **F1257** (medium, fabrication, easy/small) `multiphysics_fem.rs:1081` — Coffin-Manson damage hardcodes material constants, uses range not amplitude
- [ ] **F1258** (medium, missing-feature, medium/medium) `multiphysics_fem.rs:1166` — CoupledSolverFem advertises three coupling schemes but solves nothing
- [ ] **F1261** (medium, fabrication, medium/small) `multiscale_fem.rs:388` — SelfConsistent scheme returns Voigt–Reuss mean without iterating
- [ ] **F1263** (medium, perf, medium/small) `nonlinear/functions.rs:46` — Newton solvers use unpreconditioned CG with an absolute 1e-12 tolerance
- [ ] **F1264** (medium, perf, easy/small) `nonlinear/functions.rs:503` — modified_newton reassembles the full tangent each iteration, then discards it
- [ ] **F1265** (medium, fabrication, easy/trivial) `nonlinear/functions.rs:579` — lbfgs_1d is fixed-step gradient descent, not L-BFGS
- [ ] **F1267** (medium, fabrication, medium/small) `nonlinear/functions_2.rs:107` — Riks predictor uses fabricated uniform reference load
- [ ] **F1270** (medium, bug, medium/medium) `nonlinear_dynamics.rs:419` — NonlinearFem is linear CST with singular K and silent zero pivots
- [ ] **F1271** (medium, fabrication, medium/medium) `nonlinear_dynamics.rs:437` — NonlinearFem 'Newton' is linear; BCs not applied to K
- [ ] **F1272** (medium, panic, easy/small) `nonlinear_fem.rs:54` — Public mat3_inv asserts; reachable via kinematics/stress APIs
- [ ] **F1273** (medium, bug, easy/trivial) `nonlinear_fem.rs:211` — principal_stretches returns sqrt of C's diagonal, not eigenvalues
- [ ] **F1279** (medium, missing-feature, hard/medium) `nonlinear_fem.rs:944` — EasQ4Element has no condensed stiffness; α never solved, non-standard modes
- [ ] **F1281** (medium, perf, medium/medium) `nonlinear_solver.rs:21` — ResidualFn forces a dense tangent; 'sparse' NewtonKrylov is O(n³)
- [ ] **F1283** (medium, perf, easy/small) `nonlinear_solver.rs:249` — ModifiedNewtonRaphson re-eliminates the cached tangent every iteration
- [ ] **F1288** (medium, bug, easy/small) `nonlinear_solver.rs:1042` — gmres_solve returns Ok(x) when not converged, contrary to its doc
- [ ] **F1293** (medium, fabrication, medium/medium) `optimization_fem.rs:1039` — mbb_beam_topology_opt has no FE solve, loads, or supports
- [ ] **F1294** (medium, duplicate, medium/large) `parallel_solver.rs:53` — Five independent CsrMatrix types across crates
- [ ] **F1295** (medium, perf, medium/medium) `parallel_solver.rs:275` — ParallelAssembler locks a Mutex per nonzero; nondeterministic summation
- [ ] **F1296** (medium, perf, medium/medium) `parallel_solver.rs:284` — Parallel assembly takes a Mutex per nonzero: slow and nondeterministic
- [ ] **F1298** (medium, bug, easy/small) `parallel_solver.rs:492` — Parallel CG/GMRES dot product is a nondeterministic rayon float sum
- [ ] **F1301** (medium, fabrication, easy/small) `piezo/types.rs:300` — Piezo free strain ignores material; d33 computed as e33/c33
- [ ] **F1302** (medium, fabrication, easy/small) `poromechanics.rs:138` — Porous coupling matrix is a constant 1/n, independent of geometry
- [ ] **F1306** (medium, bug, easy/medium) `probabilistic_fem/functions.rs:157` — Dense Cholesky solvers fail silently with zeros or fake pivots
- [ ] **F1307** (medium, perf, easy/small) `probabilistic_fem/types.rs:88` — Tridiagonal bar solved by dense Cholesky in serial MC loop
- [ ] **F1312** (medium, pure-rust, easy/medium) `reduced_order.rs:343` — Hand-rolled deflated power iteration for POD/KL; use symmetric EVD
- [ ] **F1315** (medium, bug, easy/small) `reduced_order_fem.rs:660` — CG used on nonsymmetric Petrov-Galerkin and EIM matrices
- [ ] **F1318** (medium, api, easy/medium) `reliability.rs:224` — Monte Carlo samplers use unseeded thread RNG
- [ ] **F1319** (medium, bug, easy/trivial) `reliability_fem/types.rs:110` — RandomField::sample variance off by factor dx (missing 1/√dx)
- [ ] **F1320** (medium, fabrication, easy/trivial) `reliability_fem/types.rs:119` — RandomField::energy_fraction always returns 1.0
- [ ] **F1321** (medium, bug, medium/small) `reliability_fem/types.rs:195` — FORM β is always non-negative and has no convergence flag
- [ ] **F1322** (medium, fabrication, medium/small) `reliability_fem/types.rs:292` — Ditlevsen bounds use an invented joint failure probability
- [ ] **F1324** (medium, fabrication, medium/small) `reliability_fem/types.rs:511` — tvedt_pf is the Breitung formula under Tvedt's name
- [ ] **F1325** (medium, bug, easy/trivial) `reliability_fem/types.rs:796` — Response-surface quadratic coefficient missing factor 1/2
- [ ] **F1326** (medium, bug, easy/trivial) `shell/functions.rs:22` — Plate natural frequency is π times too high
- [ ] **F1331** (medium, fabrication, easy/small) `shell_fem.rs:912` — tsai_wu_failure_index actually computes a Tsai-Hill-type quadratic
- [ ] **F1334** (medium, bug, easy/small) `soil_fem.rs:772` — Consolidation step silently no-ops when explicit dt is unstable
- [ ] **F1335** (medium, perf, easy/small) `soil_fem.rs:852` — Seepage solver is plain Gauss-Seidel with a weak stop test; unmasked edges pinned to h=0
- [ ] **F1336** (medium, bug, easy/small) `soil_fem.rs:983` — Pile head deflection off by factor β; profile never computed
- [ ] **F1341** (medium, perf, medium/small) `solvers/assembly_coloring.rs:167` — Coloured parallel assembly uses Mutex and HashMap per entry
- [ ] **F1343** (medium, bug, easy/trivial) `solvers/block_schur_pc.rs:125` — Stokes Schur block applies (1/ν)Mp⁻¹ instead of ν·Mp⁻¹
- [ ] **F1345** (medium, api, easy/medium) `solvers/functions.rs:104` — conjugate_gradient: absolute tol, no stats; internal callers blind
- [ ] **F1350** (medium, bug, easy/small) `solvers/functions.rs:1346` — lu_solve skips singular pivots; Schur silently uses S=D
- [ ] **F1351** (medium, perf, medium/medium) `solvers/matrix_free/sum_factorization.rs:568` — Matrix-free apply allocates n_dofs per fold split and 7 Vecs/element
- [ ] **F1353** (medium, bug, medium/small) `sparse/functions.rs:213` — AMG prolongation maps orphan F-points to coarse node 0
- [ ] **F1354** (medium, perf, medium/medium) `sparse/functions.rs:236` — SpGEMM uses SipHash HashMap/HashSet accumulators per row
- [ ] **F1355** (medium, bug, medium/small) `sparse/functions.rs:468` — RCM ignores components; isolated min-degree start makes it a no-op
- [ ] **F1358** (medium, perf, easy/small) `sparse/types.rs:634` — CSR/CSC from_triplets via HashMap; COO to_csr clones
- [ ] **F1362** (medium, fabrication, medium/medium) `spectral_fem/types.rs:656` — BEM H/G matrices are raw kernels: no quadrature, zero diagonals
- [ ] **F1366** (medium, bug, easy/trivial) `stochastic_fem/functions.rs:187` — PCE variance omits n! normalisation of He_n basis
- [ ] **F1367** (medium, bug, easy/small) `stochastic_fem/types.rs:221` — KL from_covariance power iteration from fixed e₀ misses modes
- [ ] **F1370** (medium, bug, easy/trivial) `stochastic_fem/types.rs:622` — Morris elementary effects use stale base output y0
- [ ] **F1374** (medium, bug, medium/small) `thermal/functions.rs:462` — crank_nicolson_transient silently truncates dense K to tridiagonal
- [ ] **F1377** (medium, bug, medium/small) `thermal/types.rs:387` — Implicit/steady thermal solves drop non-adjacent contact resistances
- [ ] **F1383** (medium, api, easy/trivial) `thermal_stress.rs:496` — No reference temperature; absolute T used as thermal strain driver
- [ ] **F1385** (medium, bug, easy/small) `topology_opt/functions.rs:207` — Helmholtz PDE filter is one wrong sweep; constants not preserved
- [ ] **F1386** (medium, fabrication, medium/medium) `topology_opt/types.rs:388` — update_densities_mma is boxed OC, not MMA
- [ ] **F1388** (medium, bug, easy/small) `topology_opt/types.rs:867` — OcOptimizer::has_converged compares to initial density, not last step
- [ ] **F1389** (medium, missing-feature, hard/large) `topology_opt_fem.rs:325` — No FE-coupled topology driver; mbb_beam_topology_opt uses fake oracle
- [ ] **F1390** (medium, bug, easy/trivial) `topology_opt_fem.rs:415` — Q4 element stiffness quartered: Jacobian inverse hard-coded to 1
- [ ] **F1391** (medium, fabrication, medium/medium) `topology_opt_fem.rs:483` — ThermalTopOpt heat_compliance is a series-resistance stand-in
- [ ] **F1393** (medium, fabrication, easy/small) `topology_optimization.rs:609` — run_oc_iterations records fabricated Lagrange multiplier and compliance
- [ ] **F1394** (medium, bug, easy/trivial) `topology_optimization.rs:733` — p_norm_sensitivity has spurious extra 1/σ_ref factor
- [ ] **F1399** (medium, bug, medium/small) `topology_optimization_ext.rs:619` — Level-set reinitialisation uses fixed Δτ=0.5 and central differences
- [ ] **F1402** (medium, bug, easy/small) `truss_frame.rs:543` — support_reactions returns ≈0 residual after penalty BCs
- [ ] **F1403** (medium, bug, easy/small) `truss_frame.rs:552` — support_reactions returns ~0 because it uses penalised K and f
- [ ] **F1404** (medium, bug, easy/trivial) `truss_frame.rs:785` — Portal column moment 2x too large; buckling load has stiffness units
- [ ] **F1405** (medium, bug, easy/small) `truss_frame.rs:797` — Portal sway_buckling_load returns stiffness; column_moments doubled
- [ ] **F1406** (medium, fabrication, medium/medium) `viscoelastic_fem/functions.rs:60` — creep_compliance_integral averages 1/E; retardation times inverted
- [ ] **F1410** (medium, perf, easy/small) `wave_propagation.rs:195` — Explicit wave step: dense O(n²) matvec and unbounded history
- [ ] **F1416** (medium, bug, easy/small) `xfem/functions.rs:638` — Bifurcation base angle has wrong sign; magic 0.4 and ±30° constants
- [ ] **F1419** (medium, bug, easy/small) `xfem/functions.rs:908` — build_enriched_dof_map overwrites Heaviside offset for tip nodes, orphaning DOFs
- [ ] **F1421** (medium, fabrication, easy/small) `xfem/functions.rs:936` — locate_crack_tip_from_level_set ignores phi, fakes interpolation, can panic
- [ ] **F1422** (medium, bug, easy/trivial) `xfem/functions.rs:985` — ctod_irwin returns single-face displacement, half the opening
- [ ] **F1423** (medium, missing-feature, hard/large) `xfem/functions.rs:1006` — XFEM module lacks enriched element stiffness/assembly to feed SIF extraction
- [ ] **F1424** (medium, stub, medium/medium) `xfem/functions_2.rs:16` — reinitialise_level_set is a one-sweep clamp, not signed-distance reinitialisation
- [ ] **F1425** (medium, stub, medium/medium) `xfem/functions_2.rs:27` — reinitialise_level_set is a one-sweep clamp that corrupts valid SDFs
- [ ] **F1430** (medium, bug, medium/small) `xfem/types.rs:164` — update_after_propagation never updates phi; unchecked psi index panics
- [ ] **F1431** (medium, bug, medium/small) `xfem/types.rs:170` — update_after_propagation never updates phi; stale after kink
- [ ] **F1432** (medium, bug, easy/small) `xfem/types.rs:270` — CrackFront3D hardcodes advance direction +X regardless of front geometry
- [ ] **F1433** (medium, bug, easy/trivial) `xfem/types.rs:496` — MixedModeSif::effective computes K_III²/(1-ν)² instead of /(1-ν)

### Low (56) — tracked in the external ledger
F1018, F1023, F1025, F1040, F1041, F1053, F1084, F1086, F1087, F1090, F1093, F1113, F1120, F1125, F1139, F1151, F1157, F1170, F1175, F1180, F1184, F1192, F1203, F1209, F1224, F1228, F1229, F1230, F1240, F1245, F1246, F1250, F1259, F1268, F1278, F1282, F1287, F1289, F1299, F1300, F1310, F1317, F1337, F1338, F1340, F1344, F1352, F1360, F1368, F1387, F1400, F1401, F1413, F1415, F1428, F1435

<!-- AUDIT-2026-10-12 END -->
