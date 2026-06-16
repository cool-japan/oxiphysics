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
