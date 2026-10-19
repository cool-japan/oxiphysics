# oxiphysics-core TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 130,941 SLoC | 5,378 tests

Foundation crate of the OxiPhysics workspace: math types, ODE/PDE solvers,
statistics and optimization, spatial and parallel infrastructure. The forward
roadmap below is code-verified (2026-06-11) and concentrates on numerical
robustness, cross-platform determinism, structure-preserving/stiff integration,
and zero-allocation performance. Everything stays pure Rust: SIMD is
runtime-dispatched Rust (no FFI intrinsics), extended precision is software
double-double, and ecosystem dependencies are marked (oxifft) / (add dep).

**Cross-crate contract** — downstream consumers of this roadmap:
- Adaptive exact predicates (v0.2.0) → `oxiphysics-geometry` exact-arithmetic robust booleans/CSG and `oxiphysics-collision` narrow-phase robustness.
- Deterministic parallel reductions (v0.2.0) and fixed-point mode (v1.0) → workspace determinism program (root TODO Phase 24).
- `Dd` double-double type (v0.2.0) → conceptual basis for GPU double-single (two-float) emulation in `oxiphysics-gpu`.

## Completed (v0.1.0 – v0.1.3)

### Milestone 1: Foundation
- [x] Define core types and traits (`Vec3`, `Quat`, `Mat3`, `Transform`, `AABB`, `types`)
- [x] Implement basic error handling
- [x] Add unit tests (5,378 tests passing)

### Milestone 2: Mathematics & Solvers
- [x] Linear algebra (`linalg`, `tensor`, `sparse`, `simd_math`)
- [x] ODE solvers (`ode`, `numerical_ode`, `adaptive_timestepping`)
- [x] PDE solvers (`pde_solvers`, `finite_difference`, `spectral_methods`)
- [x] Neural ODE integration (`neural_ode`)
- [x] Fractional calculus (`fractional_calculus`)
- [x] Stochastic calculus (`stochastic`)
- [x] Quadrature (`quadrature`)
- [x] Interpolation (`interpolation`)

### Milestone 3: Statistics & Optimization
- [x] Statistics (`statistics`, `probabilistic_models`)
- [x] Monte Carlo methods (`monte_carlo`)
- [x] Bayesian inference (`bayesian_inference`)
- [x] Bayesian optimization (`bayesian_opt`)
- [x] General optimization (`optimization`)
- [x] Game theory (`game_theory`)
- [x] Functional analysis (`functional_analysis`)
- [x] Information geometry (`information_geometry`)
- [x] Machine learning primitives (`machine_learning`)
- [x] Causal inference (`causal_inference`)

### Milestone 4: Advanced Mathematics
- [x] Differential geometry (`differential_geometry`)
- [x] Topology (`topology`)
- [x] Persistent homology (`persistent_homology`)
- [x] Category theory (`category_theory`)
- [x] Chaos theory (`chaos_theory`)
- [x] Dynamic systems (`dynamic_systems`)
- [x] Stability analysis (`stability`)
- [x] Symbolic algebra (`symbolic_algebra`)
- [x] Dual quaternions (`dual_quaternion`)
- [x] Autodiff (`autodiff`)

### Milestone 5: Physics World & Infrastructure
- [x] Physics world (`world`: `Body`, `PhysicsWorld`, `Island`)
- [x] Signal processing (`signal`)
- [x] Spatial structures (`spatial`, `graph`)
- [x] Parallel computing (`parallel`, `parallel_orchestrator`)
- [x] Multi-scale methods (`multiscale_methods`)
- [x] Cache-friendly layouts (`cache_layout`)
- [x] Random number generation (`random`)
- [x] SIMD math (`simd_math`)

### Milestone 6: Polish
- [x] Documentation (rustdoc, 0 warnings)
- [x] Examples
- [x] Benchmarks (criterion)
- [x] 0 stubs (todo!/unimplemented! free in production paths)

### Verified baseline for the forward roadmap (code-checked 2026-06-11)
- [x] Symplectic baseline: `LeapFrog`, velocity-Verlet, Störmer-Verlet in `ode/` — the Yoshida/splitting work below extends these; they are not re-opened
- [x] SoA/cache infrastructure: `ParticleSoA`, `AlignedVec`, Morton ordering in `cache_layout.rs` — the v1.0 arena allocator feeds these existing layouts
- [x] Work-stealing thread pool (`WorkStealingPool` in `parallel/`) — base layer for the ordered deterministic folds below
- [x] Runtime unit handling (`dimensionless.rs`, `unit_conversion.rs`) — runtime-only today; type-level dimensions remain open (v0.3.0)
- [x] Scalar batch math (`Vec3Batch` in `simd_math.rs`) — superseded by runtime-dispatch SIMD at v1.0

## v0.2.0 — Numerical foundations & determinism

### Compensated & deterministic arithmetic
- [x] Compensated + deterministic parallel reductions (planned 2026-06-11)
  - **Files:** NEW src/compensated.rs (neumaier_sum, pairwise_sum, DeterministicReducer w/ fixed 4096 chunking + index-ordered pairwise combine); lib.rs `pub mod compensated;`
  - **Tests:** bit-identical across 1/2/4/8 threads via WorkStealingPool + sequential reference; O(ε) error vs naive on 10⁷ ill-conditioned terms
  - **Risk:** lib.rs concurrently edited by a sibling agent — Edit-only with retry
  - **Goal:** parallel sum of 10⁷ terms bit-identical across 1/2/4/8 threads; rounding error O(ε) not O(nε).
  - **Design:** Kahan/Neumaier accumulator + pairwise tree reduction; fixed-chunk ordered fold layered on existing `WorkStealingPool`. Higham 2002; Ogita-Rump-Oishi 2005.
  - Verified gap: `parallel/` has no ordered reduction.
  - Feeds: workspace determinism program (root TODO Phase 24); deterministic mode of the v1.0 SIMD layer.
  - [ ] Error-free `two_sum` building block + Neumaier serial fold
  - [ ] Pairwise tree reduction with thread-count-independent chunk boundaries
  - [ ] Ordered parallel fold API on `WorkStealingPool` (fixed chunking, stable combine order)
  - [ ] Determinism matrix test: 1/2/4/8 threads × permuted inputs, bit-identical; O(ε) error bound checked against a `Dd` reference sum

- [x] Double-double extended precision type (`Dd`) (planned 2026-06-11)
  - **Files:** NEW src/extended_precision.rs (two_sum/fast_two_sum/split/two_prod via Dekker — NO mul_add/FMA per Phase 24 determinism; Dd add/sub/mul/div/sqrt, PartialOrd)
  - **Tests:** Hilbert-matrix dot residual <1e-28; algebraic identities; round-trip f64
  - **Risk:** none beyond lib.rs shared edit
  - **Goal:** `Dd` reproduces ≈31 decimal digits; Hilbert-matrix dot residual < 1e-28.
  - **Design:** Dekker `two_sum`/`two_prod`, Knuth TwoSum, dd add/mul/div/sqrt. Hida-Li-Bailey QD.
  - Verified gap: no extended-precision type exists in core.
  - Feeds: expansion arithmetic for the exact predicates below; conceptual basis for GPU double-single (two-float) emulation in `oxiphysics-gpu`.
  - [ ] Error-free transforms: Knuth `two_sum`, Dekker `fast_two_sum` / `two_prod` (split-based, FMA-optional)
  - [ ] `Dd` add/sub/mul; div and sqrt via Newton refinement
  - [ ] Operator traits, comparisons, f64 conversions, display/parse
  - [ ] Hilbert-matrix dot-product residual test (< 1e-28) + digit-accuracy property tests against known constants

- [x] Adaptive exact geometric predicates (planned 2026-06-11)
  - **Files:** NEW src/exact_predicates.rs (expansion arithmetic: grow_expansion/expansion_sum/scale_expansion_zeroelim/estimate; orient2d/orient3d w/ stage-A filter + stage-B adaptive + exact fallback; incircle/insphere w/ stage-A filter + exact fallback — C/D adaptive stages = perf-only follow-up)
  - **Tests:** 10⁵ near-degenerate torture (10⁶ behind env flag) — zero sign errors vs Dd oracle; permutation-consistency; naive-f64 failure grid
  - **Risk:** error-bound constants — derive per Shewchuk 1997 §4; consumers: geometry booleans, collision
  - **Goal:** correct sign on 10⁶ near-degenerate coplanar/cocircular inputs where naive f64 misclassifies, 0 errors.
  - **Design:** Shewchuk adaptive-precision `orient2d/3d`, `incircle`, `insphere` via dd expansions + static filter. Shewchuk 1997.
  - Verified gap: missing in core; robustness need is shared with geometry boolean robustness.
  - Consumers: `oxiphysics-geometry` exact-arithmetic robust booleans (geometry v0.2.0); `oxiphysics-collision` narrow phase.
  - [ ] Static forward-error filter (fast path: pure f64 with certified error bound)
  - [ ] Expansion arithmetic (grow/scale/compress) on top of the `Dd` error-free transforms
  - [ ] Adaptive `orient2d` / `orient3d` (staged precision escalation)
  - [ ] Adaptive `incircle` / `insphere`
  - [ ] Adversarial corpus: 10⁶ near-coplanar/cocircular inputs, zero misclassifications vs exact reference

### Structure-preserving integration
- [x] Symplectic family completion (Yoshida 4/6/8 + generic splitting) (planned 2026-06-11)
  - **Files:** NEW src/ode/symplectic.rs (SymplecticComposition: yoshida4/6/8 + strang; reuses LeapFrog::kick/drift); ode module registration
  - **Tests:** Kepler energy drift (runtime-capped ≤10 s debug), measured order ≥3.8/5.5/7.0 via dt-halving; harmonic-oscillator symplecticity
  - **Risk:** Yoshida-6/8 coefficient transcription — cite Yoshida 1990 solution A, cross-check sum(w)=1
  - **Goal:** Yoshida-6 on Kepler conserves energy to 1e-10 over 10⁴ periods; measured order 6.
  - **Design:** triple-jump composition + Strang/operator-splitting wrapper over existing `LeapFrog`/velocity-Verlet. Yoshida 1990.
  - Verified: `LeapFrog`/velocity-Verlet/Störmer-Verlet already ship in `ode/` (see Completed) — only the Yoshida compositions and the generic splitting wrapper are open.
  - [ ] Generic Strang/operator-splitting combinator over separable Hamiltonians
  - [ ] Triple-jump coefficient generation for orders 4/6/8
  - [ ] Kepler-orbit energy-drift test (1e-10 over 10⁴ periods) + measured-order verification (slope ≈ 6)

## v0.3.0 — Geometric & stiff integrators

- [ ] Lie-group integrators on SO(3)/SE(3)
  - **Goal:** free-rigid-top stays on manifold (‖RᵀR−I‖<1e-12), angular momentum conserved to 1e-10.
  - **Design:** Munthe-Kaas RKMK4 + Magnus (2nd/4th) with exp/log maps. Munthe-Kaas 1999; Iserles et al. 2000.
  - Verified gap: no Lie-group integrators in core.
  - [ ] so(3)/se(3) exp/log maps + `dexpinv` truncated series
  - [ ] RKMK4 integrator over a generic Lie-group trait
  - [ ] Magnus expansions (2nd/4th order) for linear time-dependent systems
  - [ ] Free-rigid-top validation: orthogonality drift ‖RᵀR−I‖ < 1e-12, angular momentum conserved to 1e-10

- [ ] IMEX + exponential integrators (oxifft)
  - **Goal:** 1D advection-diffusion-reaction; ETDRK4 4th-order in time, stable at stiffness ratio 10⁶.
  - **Design:** IMEX-ARK/BDF + ETDRK4 with contour/φ-function evaluation (spectral path via oxifft). Cox-Matthews 2002; Kennedy-Carpenter ARK.
  - Verified gap: no IMEX or exponential integrators in core.
  - [ ] φ-function evaluation via complex contour quadrature (stable near small eigenvalues)
  - [ ] ETDRK4 on diagonal/spectral operators (diffusion path via oxifft)
  - [ ] IMEX-ARK tableaux (Kennedy-Carpenter) + IMEX-BDF variants
  - [ ] Advection-diffusion-reaction convergence study: measured 4th order in time, stable at stiffness ratio 10⁶

- [ ] Compile-time dimensional analysis (units type system)
  - **Goal:** dimensionally-inconsistent expression fails to compile; runtime overhead 0 (bench identical to raw f64).
  - **Design:** const-generic SI exponent vector, zero-cost `Quantity<…>` newtype.
  - Verified: `dimensionless.rs`/`unit_conversion.rs` are runtime functions only — no type-level dimensions.
  - [ ] `Quantity` newtype over a const-generic SI exponent vector (m, kg, s, A, K, mol, cd)
  - [ ] Arithmetic impls with compile-time exponent arithmetic + derived-unit aliases
  - [ ] Compile-fail tests for inconsistent expressions; criterion bench proving zero overhead vs raw f64
  - [ ] Bridge to the existing runtime `unit_conversion` tables

## v1.0 — Production & validation

- [ ] Portable SIMD abstraction with runtime dispatch (pure Rust)
  - **Goal:** Vec3 AXPY/dot ≥3× scalar on AVX2; bit-identical to scalar in deterministic mode.
  - **Design:** `f64x4`/`f32x8` wrapper + scalar fallback, `is_x86_feature_detected!` dispatch (pure Rust).
  - Verified: `simd_math.rs` `Vec3Batch` is scalar loops today — no target_feature/std::simd usage.
  - [ ] `f64x4`/`f32x8` wrapper types with scalar fallback lanes
  - [ ] Runtime dispatch (`is_x86_feature_detected!` on x86_64, NEON detection on aarch64), resolved once at startup
  - [ ] Batched kernels: AXPY, dot, normalize, AABB sweeps over `ParticleSoA`
  - [ ] Deterministic-mode equivalence test (bit-identical to scalar) + ≥3× AVX2 criterion gate

- [ ] Arena/bump allocator + SoA hardening
  - **Goal:** zero per-step heap allocations in a 10⁵-particle hot loop (verified via allocation counter).
  - **Design:** typed bump arena + frame-reset pools feeding existing `ParticleSoA`.
  - Verified: `cache_layout.rs` already has `ParticleSoA`/`AlignedVec`/Morton — only the arena layer is open.
  - [ ] Typed bump arena with frame reset (no `unsafe` leaking into the public API)
  - [ ] Frame-scoped scratch pools for solver temporaries
  - [ ] Counting-allocator test harness; zero-allocation assertion over a 10⁵-particle step loop

- [ ] Fixed-point deterministic mode
  - **Goal:** 10⁴-step rigid sim bit-identical across x86_64 and aarch64.
  - **Design:** `Qm.n` type with saturating/rounding ops for cross-platform lockstep.
  - Verified gap: no fixed-point type in core.
  - Feeds: workspace determinism program (root TODO Phase 24) — the lockstep tier beyond ordered-float determinism.
  - [ ] Const-generic `Qm.n` fixed-point type with saturating add/sub/mul/div and explicit rounding modes
  - [ ] Deterministic elementary kernels (table/polynomial) sufficient for rigid-body stepping
  - [ ] Cross-architecture golden-trace test: 10⁴-step rigid sim, byte-identical state hash on x86_64 and aarch64

## Deferred / research track
- [~] GPU double-single (two-float) arithmetic derived from `Dd` — Ready when: `Dd` (v0.2.0) has landed and `oxiphysics-gpu` schedules its f64-emulation pass on wgpu compute.
- [~] Solver-level autograd hooks via scirs2-autograd (add dep) — Ready when: `scirs2-autograd` is added as a workspace dependency (lead consumers: oxiphysics-md ML potentials, oxiphysics-materials differentiable return-mapping); core's existing `autodiff` module covers dual-number needs until then.

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 432 (critical 2, high 118, medium 242, low 70); 70 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (2)
- [x] **F0775** (critical, bug, medium/small) `numerics/functions_2.rs:446` — boys_fn Taylor series diverges for x>8, corrupting md integrals
- [x] **F0885** (critical, bug, medium/small) `sparse/functions.rs:615` — minres_solve wrong recurrence; reports converged with large residual

### High (118)
- [x] **F0518** (high, docs, easy/medium) `crates/oxiphysics-core/README.md:45` — core/geometry/collision README examples use nonexistent APIs
- [x] **F0520** (high, panic, easy/small) `adaptive_timestepping.rs:272` — compute_global_dt panics on disjoint domain dt ranges
- [x] **F0538** (high, bug, easy/small) `causal_inference.rs:578` — FrontdoorCriterion::check: condition 2 no-op, condition 3 rejects valid graphs
- [x] **F0539** (high, fabrication, medium/medium) `causal_inference.rs:621` — Frontdoor adjusted_effect ignores Σ_x' adjustment; confounded estimate
- [x] **F0540** (high, fabrication, medium/small) `causal_inference.rs:892` — 2SLS uses only first instrument; F-stat k=1; panics on bad input
- [x] **F0549** (high, bug, hard/medium) `collision.rs:552` — gjk() reports simplex upper bound and support points as minimum distance
- [x] **F0557** (high, bug, easy/trivial) `complex_analysis.rs:309` — residue() returns leading Laurent coefficient a_{−m}, not residue, for order>1
- [x] **F0558** (high, bug, easy/trivial) `complex_analysis.rs:561` — elliptic_e AGM series subtracts 2^{n+1}c² instead of 2^{n−1}… weights
- [x] **F0560** (high, bug, medium/small) `complex_analysis.rs:693` — bessel_y0/bessel_y1 series wrong: sign flip and -1/(pi x)
- [x] **F0566** (high, fabrication, medium/small) `compressed_sensing.rs:968` — K-SVD discards coefficient update and forbids overcomplete dictionaries
- [x] **F0568** (high, bug, medium/small) `continuation.rs:356` — stability_index and StabilityAnalysis::label count Gershgorin discs for n>2
- [x] **F0569** (high, fabrication, medium/medium) `continuation.rs:498` — branch_switching uses a Jacobian column as the null vector
- [x] **F0572** (high, bug, easy/small) `control.rs:430` — TransferFunction::poles fabricates real roots via unchecked Newton
- [x] **F0576** (high, fabrication, medium/medium) `control.rs:728` — ModelPredictiveControl uses wrong one-step gradient, fixed step
- [x] **F0577** (high, bug, medium/medium) `control.rs:748` — ModelPredictiveControl::compute uses truncated gradient and fixed 50 steps
- [x] **F0583** (high, bug, easy/small) `control_systems/functions.rs:699` — matrix_rank ties pivot row to column index; underestimates rank
- [x] **F0584** (high, fabrication, medium/small) `control_systems/functions.rs:814` — observer_gain_ackermann returns heuristic diagonal gain for n>1
- [x] **F0585** (high, fabrication, medium/medium) `control_systems/functions.rs:881` — step_response/impulse_response use fake equal-weight modal sums, not simulation
- [x] **F0588** (high, bug, easy/trivial) `control_theory.rs:390` — StateSpaceModel::is_stable via power iteration reports unstable systems stable
- [x] **F0589** (high, bug, easy/small) `control_theory.rs:801` — BodeAnalysis::margins never detects −180° crossing (wrapped phase)
- [x] **F0590** (high, bug, easy/small) `control_theory.rs:812` — BodeAnalysis::margins misses −180° crossing because phase is wrapped
- [x] **F0593** (high, fabrication, medium/small) `control_theory.rs:1321` — PolesZeros::step_response_approx fabricates step response; ignores zeros, residues and complex poles
- [x] **F0595** (high, bug, easy/trivial) `differential_geometry/functions.rs:222` — christoffel_symbols_numerical swaps c112/c121, sign-flips Γ¹₁₂ and Γ²₁₁
- [x] **F0597** (high, bug, medium/small) `differential_geometry/functions.rs:350` — mat3_exp: linear scaling exponent overflows; Padé code discarded
- [x] **F0601** (high, bug, medium/small) `dual_quaternion.rs:410` — ScLERP loses pure translations; no shortest-path flip
- [x] **F0604** (high, bug, medium/small) `dynamic_systems.rs:847` — fundamental_solution_matrix transposes Φ every step (computes A·Φᵀ)
- [x] **F0611** (high, bug, medium/small) `dynamical_systems.rs:1106` — DoublePendulum RHS inconsistent with energy; gravity sign flipped
- [ ] **F0613** (high, bug, medium/small) `ergodic_theory.rs:225` — LyapunovSpectrum::update_qr discards the Q frame between steps
- [x] **F0614** (high, fabrication, medium/medium) `ergodic_theory.rs:290` — lyapunov_exponents sums log coordinate increments, not exponents
- [ ] **F0615** (high, fabrication, medium/medium) `ergodic_theory.rs:306` — lyapunov_exponents(trajectory) averages log step sizes, not divergence
- [ ] **F0626** (high, bug, medium/medium) `formal_verification/functions.rs:525` — CEGAR reports Verified after bounded (depth 10) abstract search
- [x] **F0635** (high, bug, easy/small) `fractional_calculus.rs:352` — riesz_derivative right-sided GL sum truncated at j+1 terms
- [x] **F0636** (high, fabrication, medium/medium) `fractional_calculus.rs:379` — fractional_laplacian_1d returns zeros at s=½ and panics at s=1
- [x] **F0638** (high, bug, medium/medium) `fractional_calculus.rs:601` — fde_predictor_corrector adds y0 to GL memory, double-counts history
- [ ] **F0639** (high, bug, medium/medium) `fractional_calculus.rs:610` — FDE predictor-corrector double-counts memory; wrong weights; vacuous test
- [x] **F0646** (high, bug, easy/small) `functional_analysis.rs:561` — smallest_eigenvalue runs forward power iteration, returns largest
- [ ] **F0647** (high, fabrication, easy/small) `functional_analysis.rs:581` — smallest_eigenvalue runs forward power iteration, returns the largest
- [x] **F0652** (high, bug, easy/trivial) `game_theory/types.rs:788` — mixed_nash_2x2 uses wrong payoff entries in indifference formula
- [x] **F0656** (high, bug, easy/trivial) `geometry_primitives.rs:824` — frustum_from_vp plane offsets have wrong sign and are unnormalised
- [x] **F0661** (high, fabrication, hard/large) `information_geometry.rs:34` — StatisticalManifold Fisher metric is a hardcoded Gaussian proxy at x=0
- [x] **F0667** (high, fabrication, medium/medium) `information_geometry.rs:519` — InformationProjection project/reverse_kl ignore input distribution p
- [ ] **F0668** (high, fabrication, medium/medium) `information_geometry.rs:530` — InformationProjection project/reverse_kl ignore the input distribution (NEW)
- [ ] **F0669** (high, fabrication, medium/medium) `information_geometry.rs:560` — StatisticalManifold Fisher metric hardwired to one Gaussian sample (NEW)
- [x] **F0673** (high, fabrication, medium/small) `information_theory/functions.rs:512` — conditional_transfer_entropy is TE(X→Y) minus TE(Z→Y)
- [x] **F0678** (high, bug, easy/small) `information_theory/functions.rs:930` — fisher_rao_gaussian uses an incorrect closed form
- [ ] **F0690** (high, bug, easy/small) `linalg/functions/functions_impl.rs:108` — symmetric_eigenvalues3 has the wrong Cardano q sign, so the spectrum is mirrored about tr/3
- [ ] **F0701** (high, bug, easy/small) `math/geometry.rs:1206` — quat_integrate_rk4 normalises ω and q̇, forcing a 2 rad/s rate and NaN at ω=0
- [ ] **F0707** (high, bug, easy/trivial) `metric_spaces.rs:491` — MetricBallTree kNN heap evicts the nearest neighbour, giving wrong results
- [ ] **F0713** (high, bug, easy/small) `monte_carlo.rs:643` — MarkovChainMC::mala has no Hastings correction and samples a biased target
- [ ] **F0714** (high, bug, easy/small) `monte_carlo.rs:728` — integrate_stratified divides its integral estimate by the number of strata
- [ ] **F0723** (high, bug, easy/small) `network_dynamics.rs:932` — CommunityDetection::modularity applies the null-model penalty only on edges
- [x] **F0728** (high, bug, hard/medium) `neural_ode.rs:596` — AdjointMethod::run integrates J·a with wrong sign and wrong times
- [x] **F0729** (high, fabrication, medium/medium) `neural_ode.rs:776` — TimeSeriesOde::fit computes no gradient; it only shrinks bias_out toward zero
- [x] **F0733** (high, bug, medium/medium) `numerical_linear_algebra/functions.rs:193` — svd_small: wide inputs never orthogonalised, U/V unsorted, V returned as Vt
- [x] **F0737** (high, fabrication, medium/medium) `numerical_linear_algebra/types.rs:134` — Arnoldi/Lanczos: Krylov dimension = k, converged always true, eigenvectors are basis vectors
- [x] **F0738** (high, fabrication, hard/medium) `numerical_linear_algebra/types.rs:211` — lobpcg is fixed-step 0.1 gradient descent, not LOBPCG
- [x] **F0739** (high, fabrication, medium/small) `numerical_linear_algebra/types.rs:340` — Lanczos tridiag_eigenvalues is a non-converging fake Jacobi sweep
- [x] **F0742** (high, bug, hard/medium) `numerical_linear_algebra/types.rs:638` — MINRES recurrence is wrong; a 1×1 system returns b/(1+b²)
- [x] **F0743** (high, bug, medium/medium) `numerical_linear_algebra/types.rs:881` — TFQMR uses a stale A·u and wrong v update; returns converged=true at 1e19
- [x] **F0747** (high, bug, medium/small) `numerical_linear_algebra/types.rs:1321` — MatrixFunctions::logm uses a wrong Padé form and returns about ½·log(A)
- [x] **F0751** (high, bug, easy/small) `numerical_methods/functions.rs:326` — RK45 error estimate omits Dormand–Prince e7·k7 term, so it is O(h)
- [x] **F0759** (high, bug, easy/small) `numerical_methods/functions.rs:938` — qr_eigenvalues forms R·Qᵀ instead of R·Q; arnoldi_eigenvalues inherits it
- [x] **F0761** (high, bug, medium/small) `numerical_methods_ext.rs:332` — Gauss-Hermite n>5 Newton guesses collapse to one duplicated root
- [x] **F0763** (high, bug, easy/trivial) `numerical_methods_ext.rs:864` — SplineInterpolation::eval uses a wrong cubic correction term
- [x] **F0772** (high, bug, easy/small) `numerics/functions.rs:422` — inc_beta continued fraction swapped/wrong terms: returns values above 1
- [x] **F0779** (high, bug, easy/small) `ode/types.rs:101` — DormandPrince/AdaptiveIntegrator evaluate every RK stage at step-start time
- [x] **F0788** (high, bug, easy/small) `optimization/functions.rs:1425` — bfgs inverse-Hessian update is not BFGS and H becomes indefinite
- [x] **F0798** (high, fabrication, medium/medium) `parallel/types.rs:248` — WorkStealingPool workers exit immediately; all tasks run serially
- [x] **F0802** (high, bug, easy/small) `pde/types/types_2.rs:611` — MethodOfLines::integrate_dopri5 mixes DP stages with Fehlberg weights
- [x] **F0828** (high, bug, medium/small) `probabilistic_models.rs:880` — VI ELBO KL formula wrong; CAVI responsibilities never normalised
- [x] **F0831** (high, bug, medium/small) `quadrature.rs:455` — Gauss-Hermite weights are n times too large and nodes collapse for n≥5
- [x] **F0832** (high, bug, easy/small) `quadrature.rs:497` — Gauss-Laguerre uses wrong derivative; nodes collapse, weights sum ≠1
- [x] **F0834** (high, bug, easy/small) `random/functions.rs:617` — regularized_incomplete_beta drops first CF term, so student_t_cdf is wrong
- [x] **F0835** (high, bug, medium/medium) `random/types.rs:198` — SobolSequence capped at 4 dims; dims 0 and 1 identical
- [x] **F0839** (high, fabrication, medium/medium) `random/types.rs:1404` — Sobol direction tables make dimension 1 identical to dimension 0
- [x] **F0840** (high, fabrication, medium/small) `random/types.rs:1518` — fBm 'Hosking' generator ignores Hurst correlation and produces plain BM
- [x] **F0842** (high, fabrication, easy/small) `random_processes.rs:486` — CTMC mean_first_passage_time ignores `from` and returns recurrence time
- [x] **F0843** (high, bug, easy/trivial) `random_processes.rs:633` — stable_sample α=1 branch evaluates tan(π/2), so Cauchy paths explode
- [x] **F0852** (high, fabrication, medium/medium) `sampling/functions.rs:307` — sobol10 direction numbers fabricated; dims 2–8 not stratified
- [x] **F0853** (high, fabrication, easy/small) `sampling/functions.rs:324` — Sobol 10-D direction numbers are invented, not Joe & Kuo
- [x] **F0869** (high, fabrication, medium/small) `signal/types.rs:467` — ButterworthBandpass is ad-hoc RBJ biquad; cutoffs far below −3 dB
- [x] **F0870** (high, fabrication, medium/medium) `signal/types.rs:474` — ButterworthBandpass is an ad-hoc RBJ fudge, not a Butterworth design
- [x] **F0872** (high, bug, easy/trivial) `signal/types.rs:770` — ButterworthLowPass odd-order section has flipped feedback sign
- [x] **F0876** (high, panic, easy/small) `signal_processing.rs:446` — welch_psd hangs on empty input or overlap==segment_len
- [x] **F0882** (high, bug, easy/small) `signal_processing.rs:1076` — db2_reconstruct uses reversed filters, so it does not invert db2_decompose
- [x] **F0883** (high, bug, medium/medium) `signal_processing.rs:1116` — resample loses amplitude, is under-filtered and misaligned
- [x] **F0886** (high, bug, medium/small) `sparse/functions.rs:864` — saddle_point_solve: Schur RHS sign flipped, diagonal-only pressure solve
- [x] **F0887** (high, fabrication, easy/small) `sparse/functions.rs:960` — ritz_values zeroes off-diagonals in one sweep; returns non-eigenvalues
- [x] **F0888** (high, fabrication, medium/medium) `sparse/types.rs:41` — SparseLu::ilu0 is a dense exact LU, not ILU(0)
- [x] **F0893** (high, bug, medium/medium) `sparse/types.rs:757` — IncompleteCholesky::solve_upper does not solve L^T x=b; factor/apply are dense
- [ ] **F0896** (high, bug, easy/small) `spatial/types.rs:412` — BallTree prunes with pivot centre but centroid radius; misses neighbours
- [ ] **F0905** (high, bug, easy/small) `statistics/functions.rs:419` — PCA power iteration from a fixed e0 returns wrong or duplicate components
- [ ] **F0907** (high, duplicate, easy/small) `statistics/functions.rs:682` — shapiro_wilk_w coefficients unnormalised, W approx n instead of <=1
- [ ] **F0908** (high, bug, easy/trivial) `statistics/functions.rs:724` — log_gamma_stirling off by one; Beta pdf 6×, Poisson pmf k× wrong
- [ ] **F0909** (high, duplicate, medium/small) `statistics/functions.rs:921` — t_test_one_sample_full/welch_t_test claim Student-t but use normal p-values
- [ ] **F0910** (high, fabrication, medium/small) `statistics/functions.rs:930` — t-test 'full'/Welch aliases use normal p-values; claimed df unused
- [ ] **F0911** (high, bug, easy/trivial) `statistics/functions.rs:1100` — Anderson-Darling uses ln F instead of ln(1-F) for the upper term
- [ ] **F0916** (high, bug, easy/trivial) `stochastic/types.rs:223` — LangevinDynamics::step noise scaled by 1/√dt; temperature ∝ 1/dt
- [ ] **F0919** (high, bug, easy/trivial) `stochastic/types.rs:458` — SABR Hagan implied vol multiplies by z²/x(z) instead of z/x(z)
- [ ] **F0920** (high, fabrication, medium/small) `stochastic/types.rs:996` — FractionalBrownianMotion::sample claims exact Hosking but uses a one-term AR
- [ ] **F0925** (high, fabrication, medium/small) `stochastic_processes.rs:671` — Levy stable CMS: alpha=1 skew term is zero, skew sign flipped
- [ ] **F0927** (high, bug, easy/small) `symbolic_algebra/functions.rs:108` — diff of x^n uses the f'/f form, so taylor_expand silently zeroes coefficients
- [ ] **F0934** (high, bug, easy/trivial) `tensor/operations.rs:181` — cp_als Khatri-Rao operand order mismatches mode_n_unfold layout
- [ ] **F0935** (high, bug, easy/small) `tensor/operations.rs:230` — cp_als discards lambdas when max_iter reached without convergence
- [ ] **F0952** (high, fabrication, medium/medium) `topological_data.rs:217` — vietoris_rips_h1_approx fabricates H1 birth and death values
- [ ] **F0954** (high, bug, medium/small) `topological_data.rs:515` — VietorisRips::betti_numbers returns β2 = triangle count; β1 ignores triangles
- [x] **F0964** (high, fabrication, medium/medium) `topology_optimization.rs:517` — MultiLoadTopOpt compliance proxy rises with stiffness; no FE solve
- [x] **F0968** (high, stub, easy/small) `topology_optimization.rs:940` — ManufacturingFilter::casting_filter is an identity no-op
- [x] **F0973** (high, bug, easy/trivial) `types/types_3.rs:95` — MassProperties::combine omits parallel-axis shift; NaN COM
- [x] **F0974** (high, bug, easy/small) `types/types_3.rs:184` — First arena handle (0,0) equals BodyHandle::null() and is invalid
- [x] **F0984** (high, fabrication, medium/small) `uncertainty_quantification.rs:926` — GP predictive variance is a fake 'diagonal approximation'
- [ ] **F0987** (high, bug, easy/small) `uncertainty_quantification.rs:1136` — LCG low bits make bootstrap resamples exact permutations
- [x] **F0989** (high, bug, medium/small) `wavelet_transform/functions.rs:137` — Periodic DWT emits duplicate coefficients; Db2–Db6 roundtrip fails
- [x] **F1000** (high, bug, easy/small) `wavelets.rs:860` — SURE risk cancels to n−2k+Σx²; threshold is always max|x|
- [x] **F1007** (high, bug, easy/trivial) `world/types.rs:250` — Quaternion integration treats world-frame angular velocity as body-frame
- [x] **F1013** (high, fabrication, hard/medium) `world/types.rs:926` — PhysicsWorld accepts constraints but step never solves them
- [x] **F1014** (high, bug, easy/small) `world/types.rs:962` — step_with_accumulator hangs when fixed_dt<=0; no substep cap
- [x] **F1015** (high, fabrication, hard/large) `world/types.rs:970` — PhysicsWorld::step ignores stored constraints; add_constraint has no effect

### Medium (242)
- [ ] **F0517** (medium, missing-feature, medium/medium) `crates/oxiphysics-core/Cargo.toml:20` — No serde feature: world state and config cannot be persisted
- [ ] **F0519** (medium, bug, easy/small) `adaptive_timestepping.rs:244` — UnifiedTimeStepper 'subcycling' never activates: global dt is the minimum
- [ ] **F0524** (medium, fabrication, hard/medium) `bayesian_inference.rs:556` — nuts_step is a single-leapfrog HMC step, not NUTS
- [ ] **F0526** (medium, bug, easy/trivial) `bayesian_inference.rs:742` — BayesianLinearRegression::log_evidence uses diagonal log-determinant
- [ ] **F0528** (medium, bug, easy/small) `bayesian_inference.rs:1159` — GP cholesky silently replaces non-PD pivots with 1e-10
- [ ] **F0529** (medium, bug, easy/trivial) `bayesian_opt.rs:63` — Matérn-5/2 kernel uses √10·d/l instead of √5·d/l
- [ ] **F0532** (medium, panic, easy/small) `bayesian_opt.rs:477` — suggest_next panics when n_candidates=0; unseeded RNG; optimize swallows errors
- [ ] **F0533** (medium, bug, easy/small) `bayesian_opt.rs:510` — BayesianOptimizer swallows fit errors, unseeded LHS, panics on zero candidates
- [ ] **F0535** (medium, bug, easy/small) `cache_layout.rs:552` — Morton sort collapses negative coordinates into cell 0, destroying locality
- [ ] **F0537** (medium, bug, easy/small) `causal_inference.rs:457` — Backdoor adjusted_effect ignores tolerance, panics on empty input, silently falls back to E[Y|Z]
- [ ] **F0541** (medium, bug, medium/small) `causal_inference.rs:1082` — PC orientation clears skeleton one-way, producing wrong orientations and cycles
- [ ] **F0542** (medium, perf, easy/small) `causal_inference.rs:1171` — partial_correlation recursion is exponential with full data rescans
- [ ] **F0543** (medium, bug, medium/small) `chaos_theory.rs:441` — Lyapunov spectrum integrates tangent vectors with explicit Euler
- [ ] **F0546** (medium, bug, medium/small) `chaos_theory.rs:1123` — largest_lyapunov is not Rosenstein; mismatched truncation biases exponent
- [ ] **F0550** (medium, bug, easy/trivial) `collision.rs:628` — GJK line/AC-edge search direction points away from origin
- [ ] **F0551** (medium, bug, easy/trivial) `collision.rs:835` — obb_obb_sat returns a separating axis with arbitrary sign
- [ ] **F0552** (medium, bug, medium/small) `collision.rs:1152` — EPA contact point lies in Minkowski space; normal direction mis-documented
- [ ] **F0554** (medium, bug, medium/small) `collision.rs:1412` — shape_cast_sphere_vs_aabb: None on initial overlap, false hits at corners
- [ ] **F0555** (medium, perf, easy/small) `compensated.rs:166` — par_sum copies the input, spawns a pool per call, and effectively runs serially
- [ ] **F0556** (medium, pure-rust, easy/small) `complex/functions.rs:80` — rfft/PSD return zero-padded spectra, contradicting 'N/2+1' doc
- [ ] **F0559** (medium, bug, medium/small) `complex_analysis.rs:606` — bessel_j uses unstable upward recurrence and divergent power series
- [ ] **F0562** (medium, bug, easy/small) `compressed_sensing.rs:310` — Lipschitz estimate collapses to 1e-10 when A·1=0, so ISTA/FISTA diverge
- [ ] **F0565** (medium, fabrication, easy/small) `compressed_sensing.rs:755` — rip_constant tests only uniform vectors on contiguous supports
- [ ] **F0567** (medium, fabrication, medium/small) `compressed_sensing.rs:1089` — MRI CS ignores the DCT sparsity basis it advertises and drops the sine half of k-space
- [ ] **F0571** (medium, bug, easy/small) `control.rs:295` — StateSpace::discretize_zoh uses Bd = dt·B and a truncated Taylor series for Ad
- [ ] **F0573** (medium, bug, easy/small) `control.rs:443` — poles()/zeros() push non-converged Newton iterates as fake real roots
- [ ] **F0574** (medium, bug, easy/small) `control.rs:539` — BodeAnalysis::stability_margins misses wrapped −180° crossing and returns fake 0 dB / 180°
- [ ] **F0578** (medium, bug, easy/small) `control_systems/functions.rs:323` — stability_margins uses ±5%/±5° sample windows instead of crossover interpolation
- [ ] **F0579** (medium, bug, easy/small) `control_systems/functions.rs:342` — stability_margins uses tolerance windows on wrapped phase; biased margins
- [ ] **F0580** (medium, bug, easy/small) `control_systems/functions.rs:396` — nyquist_encirclements integrates only ω>0; wrong count and sign
- [ ] **F0581** (medium, bug, easy/small) `control_systems/functions.rs:452` — discretize_zoh is 2nd-order Taylor; a stable plant becomes an unstable discrete model
- [ ] **F0582** (medium, bug, easy/small) `control_systems/functions.rs:626` — mat_inv_small returns a zero matrix on singular input; LQR and Kalman silently degrade
- [ ] **F0586** (medium, fabrication, medium/medium) `control_systems/types.rs:27` — MpcController::compute is infinite-horizon LQR; horizon misused as iteration cap
- [ ] **F0591** (medium, bug, easy/small) `control_theory.rs:998` — qr_eigenvalues is a self-admitted simplified QR; RootLocus poles unreliable
- [ ] **F0592** (medium, fabrication, easy/trivial) `control_theory.rs:1159` — ZieglerNichols::simc returns ZN PID gains, not SIMC
- [ ] **F0594** (medium, bug, easy/trivial) `convex_analysis.rs:679` — EllipsoidMethod divides by zero for 1-D problems (NaN iterate)
- [ ] **F0596** (medium, bug, easy/small) `differential_geometry/functions.rs:277` — geodesic_curve: explicit Euler, FD step tied to ds
- [ ] **F0598** (medium, bug, easy/trivial) `differential_geometry/functions.rs:397` — so3_log_matrix returns zero at θ=π, ill-conditioned near π
- [ ] **F0599** (medium, panic, medium/medium) `differential_geometry/types.rs:108` — MetricTensorND::inverse panics on singular metric from public APIs
- [ ] **F0600** (medium, bug, easy/trivial) `differential_geometry/types.rs:344` — christoffel_111 uses wrong F_u and wrong sign on g12 term
- [ ] **F0602** (medium, perf, easy/small) `dynamic_systems.rs:338` — correlation_dimension stores and sorts all n²/2 distances
- [ ] **F0603** (medium, bug, medium/small) `dynamic_systems.rs:392` — Rosenstein Lyapunov estimator has no Theiler window
- [ ] **F0605** (medium, fabrication, easy/small) `dynamic_systems.rs:1028` — instantaneous_phase fakes Hilbert transform with cumulative sum
- [ ] **F0606** (medium, fabrication, easy/small) `dynamic_systems.rs:1041` — instantaneous_phase uses a running sum as a fake Hilbert transform
- [ ] **F0609** (medium, bug, medium/medium) `dynamical_systems.rs:328` — Lyapunov tangent flow is Euler with Jacobian at x_{n+1}; panics
- [ ] **F0612** (medium, panic, medium/small) `dynamical_systems.rs:1256` — maximal_lyapunov_exponent panics on renorm=0; Euler tangent
- [ ] **F0616** (medium, bug, easy/small) `ergodic_theory.rs:467` — RQA ignores min_line for DET and counts the identity line
- [ ] **F0620** (medium, bug, easy/trivial) `error.rs:567` — check_md_timestep returns Ok(NaN); check_lbm_mach passes negatives
- [ ] **F0623** (medium, bug, medium/small) `finite_difference.rs:646` — fractional_step_1d Poisson sign makes projection double the divergence
- [ ] **F0628** (medium, perf, easy/small) `formal_verification/types.rs:530` — Model-checker BFS scans every transition per state: O(V·E)
- [ ] **F0629** (medium, bug, easy/trivial) `formal_verification/types.rs:655` — Bisimulation signatures are multisets, splitting bisimilar states
- [ ] **F0632** (medium, bug, hard/medium) `fractional_calculus.rs:121` — Mittag-Leffler power series cancels catastrophically for large negative z
- [ ] **F0633** (medium, bug, easy/small) `fractional_calculus.rs:191` — rl_integral integrates over one extra interval
- [ ] **F0634** (medium, perf, easy/small) `fractional_calculus.rs:239` — rl_derivative recomputes GL coefficients per term: O(n²), series O(n³)
- [ ] **F0637** (medium, fabrication, hard/medium) `fractional_calculus.rs:550` — subdiffusion_pdf returns plain Gaussian, not the fractional propagator
- [ ] **F0640** (medium, bug, medium/small) `fractional_calculus.rs:723` — bagley_torvik applies GL to raw x (no Caputo IC correction)
- [ ] **F0641** (medium, bug, medium/small) `fractional_calculus.rs:812` — Lévy CMS sampler subtracts bogus ξ shift; α=1 branch wrong
- [ ] **F0642** (medium, fabrication, medium/small) `fractional_calculus.rs:933` — subdiffusion_trajectory uses Weibull draws, yields normal diffusion
- [ ] **F0643** (medium, bug, easy/small) `fractional_calculus.rs:962` — Lévy flight reseeds LCG per sample with consecutive seeds
- [ ] **F0645** (medium, bug, easy/trivial) `functional_analysis.rs:463` — BanachSpace L∞ norm uses the signed max
- [ ] **F0648** (medium, fabrication, easy/small) `functional_analysis.rs:626` — resolvent_norm returns 1/‖A−λI‖_F instead of 1/σ_min
- [ ] **F0649** (medium, bug, medium/small) `functional_analysis.rs:905` — euler_lagrange_residual uses a wrong u′ stencil at x_{i−1}
- [ ] **F0654** (medium, bug, easy/small) `game_theory/types.rs:1043` — kalai_smorodinsky can return a Pareto-dominated point
- [ ] **F0655** (medium, bug, medium/small) `geometry_primitives.rs:596` — OrientedBox::intersects_aabb is AABB-of-OBB, not the claimed SAT
- [ ] **F0657** (medium, bug, medium/small) `graph/functions.rs:469` — Edmonds-Karp loops forever when source==sink; dense n×n matrix
- [ ] **F0658** (medium, bug, medium/small) `graph/functions.rs:829` — betweenness_centrality counts one arbitrary shortest path per pair
- [ ] **F0659** (medium, test, easy/medium) `graph/functions.rs:1561` — Tests assert only inside if-let or if-converged branches
- [ ] **F0662** (medium, bug, medium/small) `information_geometry.rs:86` — StatisticalManifold::geodesic raises the Christoffel index twice
- [ ] **F0663** (medium, bug, easy/trivial) `information_geometry.rs:261` — GaussianManifold::geodesic_distance missing the √2 factor
- [ ] **F0664** (medium, fabrication, medium/small) `information_geometry.rs:275` — Gaussian exponential/logarithmic maps are linear, not geodesic
- [ ] **F0670** (medium, bug, easy/small) `information_geometry.rs:623` — invert_matrix silently returns identity on singular input
- [ ] **F0671** (medium, bug, easy/trivial) `information_theory/functions.rs:130` — interaction_information_3 returns the negation of its documented formula
- [ ] **F0674** (medium, fabrication, medium/medium) `information_theory/functions.rs:523` — conditional_transfer_entropy returns TE(X→Y) minus TE(Z→Y)
- [ ] **F0677** (medium, bug, easy/small) `information_theory/functions.rs:908` — natural_gradient_step_categorical leaves the probability simplex
- [ ] **F0679** (medium, bug, easy/trivial) `information_theory/functions.rs:938` — fisher_rao_gaussian uses an incorrect distance formula
- [ ] **F0680** (medium, fabrication, medium/small) `information_theory/functions.rs:1003` — parallel_transport_simplex is not parallel transport
- [ ] **F0682** (medium, bug, medium/small) `interpolation/functions.rs:459` — rbf_tps_fit omits the affine polynomial block and fails for unit-spaced points
- [ ] **F0684** (medium, bug, medium/small) `interval/functions.rs:158` — interval_newton drops boxes when the budget runs out, ignores f_interval
- [ ] **F0691** (medium, bug, medium/small) `linalg/functions/functions_impl.rs:871` — gmres_n never restarts and returns unconverged x as Some
- [ ] **F0692** (medium, bug, easy/trivial) `linalg/functions/functions_impl.rs:968` — ilu0_n replaces negative pivots with EPS, giving garbage factors
- [ ] **F0693** (medium, fabrication, medium/small) `machine_learning/functions.rs:157` — kmeans claims k-means++ but runs deterministic farthest-first
- [ ] **F0694** (medium, bug, easy/small) `machine_learning/functions.rs:282` — PCA power iteration from all-ones loses the second component
- [ ] **F0695** (medium, perf, medium/medium) `machine_learning/functions.rs:326` — best_split is O(d·n²) with Vec and HashMap allocation per threshold
- [ ] **F0696** (medium, bug, easy/small) `machine_learning/functions.rs:976` — GP Cholesky silently clamps non-PD pivots to 1e-300
- [ ] **F0697** (medium, bug, easy/small) `machine_learning/functions.rs:1190` — t-SNE joint probabilities normalised by n², summing to 2/n
- [ ] **F0698** (medium, bug, easy/trivial) `machine_learning/types.rs:328` — DenseLayer::backward computes the input gradient from updated weights
- [ ] **F0699** (medium, fabrication, medium/small) `machine_learning/types.rs:368` — RandomForest computes n_features but never subsamples features
- [ ] **F0700** (medium, bug, easy/trivial) `math/geometry.rs:20` — mat3_to_euler_zyx gimbal-lock branch negates yaw
- [ ] **F0702** (medium, bug, medium/small) `math/linear_algebra.rs:1718` — mat3_log_rotation breaks down near 180° rotations
- [ ] **F0706** (medium, bug, easy/small) `measure_theory.rs:271` — Monte Carlo measure routines use unseeded OS RNG; degenerate bounds panic
- [ ] **F0708** (medium, perf, medium/small) `metric_spaces.rs:645` — Ball tree splits points by insertion order, defeating pruning
- [ ] **F0710** (medium, bug, easy/small) `monte_carlo.rs:63` — inverse_cdf is mislabeled, accurate to only 4.5e-4, and panics
- [ ] **F0711** (medium, bug, medium/small) `monte_carlo.rs:151` — Knuth Poisson sampler returns ≈745 for large λ
- [ ] **F0715** (medium, bug, easy/small) `monte_carlo.rs:779` — integrate_halton silently truncates dimension above 10
- [ ] **F0716** (medium, bug, easy/trivial) `multiscale_methods/functions.rs:134` — cell_problem_1d solves the opposite-sign PDE from its doc
- [ ] **F0717** (medium, bug, easy/trivial) `multiscale_methods/functions.rs:216` — cauchy_born_stress_1d returns force×F, not a stress
- [ ] **F0719** (medium, fabrication, medium/medium) `multiscale_methods/functions.rs:468` — umbrella_free_energy is an ad-hoc scalar, not WHAM as documented
- [ ] **F0721** (medium, perf, easy/small) `network_dynamics.rs:665` — Kuramoto all-to-all RHS is O(n²) where an O(n) identity exists
- [ ] **F0724** (medium, perf, medium/medium) `network_dynamics.rs:965` — greedy_modularity clones the partition and recomputes full Q per pair
- [ ] **F0725** (medium, bug, medium/medium) `network_dynamics.rs:1086` — percolation_curve normalises by surviving nodes and rebuilds the graph per removal
- [ ] **F0726** (medium, perf, medium/medium) `neural_ode.rs:339` — param_grad_contrib does 2P forward passes and 2P param clones per call
- [ ] **F0732** (medium, bug, easy/trivial) `numerical_integration.rs:243` — Clenshaw–Curtis last-term coefficient is wrong for even n
- [ ] **F0735** (medium, bug, easy/small) `numerical_linear_algebra/types.rs:74` — inverse_iteration uses CG on indefinite A-σI; from_triplets doesn't sum duplicates
- [ ] **F0736** (medium, bug, medium/medium) `numerical_linear_algebra/types.rs:107` — Arnoldi/Lanczos: Krylov dim = k, fake converged flag, no Ritz vectors
- [ ] **F0740** (medium, perf, medium/medium) `numerical_linear_algebra/types.rs:476` — ILUT allocates and scans a dense length-n row for every row
- [ ] **F0744** (medium, fabrication, medium/medium) `numerical_linear_algebra/types.rs:946` — SPAIPreconditioner is Jacobi, not even the diagonal-pattern SPAI minimiser
- [ ] **F0745** (medium, bug, medium/small) `numerical_linear_algebra/types.rs:1169` — Nyström returns factors that reconstruct C·Cᵀ, not C·W⁺·Cᵀ
- [ ] **F0750** (medium, bug, easy/small) `numerical_methods/functions.rs:231` — euler_ode/rk2/rk4/solve_rk4 loop forever when dt <= 0
- [ ] **F0753** (medium, fabrication, easy/small) `numerical_methods/functions.rs:564` — deriv_complex_step is a tiny-h central difference and returns 0 at x=0
- [ ] **F0756** (medium, bug, easy/small) `numerical_methods/functions.rs:649` — gauss_legendre/hermite_nodes silently fall back to 3 points; 10-digit constants
- [ ] **F0757** (medium, fabrication, easy/small) `numerical_methods/functions.rs:863` — svd_singular_values truncates to 5 values via power iteration on AᵀA
- [ ] **F0758** (medium, bug, easy/small) `numerical_methods/functions.rs:899` — power_iteration returns |λ| and stalls on its e₀ start vector
- [ ] **F0760** (medium, bug, easy/trivial) `numerical_methods/functions_2.rs:47` — Pseudo-arclength corrector computes Δλ with the wrong formula
- [ ] **F0762** (medium, bug, medium/small) `numerical_methods_ext.rs:784` — solve_linear_system skips singular pivot then divides by it; NaN RBF
- [ ] **F0764** (medium, bug, medium/small) `numerical_methods_ext.rs:932` — Not-a-knot spline boundary rows drop the m2 / m_{n-3} term
- [ ] **F0765** (medium, bug, easy/small) `numerical_ode.rs:149` — Fixed-step integrate loops never terminate for dt<=0
- [ ] **F0767** (medium, fabrication, medium/medium) `numerical_ode.rs:521` — ImplicitEuler/BDF2 Newton computes full FD Jacobian, uses only diagonal
- [ ] **F0768** (medium, bug, medium/small) `numerical_ode.rs:601` — Trapezoidal advertised A-stable but solved by fixed-point iteration
- [ ] **F0769** (medium, bug, easy/small) `numerical_ode.rs:722` — BDF2::integrate applies constant-step coefficients to shortened last step
- [ ] **F0770** (medium, bug, easy/trivial) `numerics/functions.rs:20` — bisect and regula_falsi never validate bracket; report converged
- [ ] **F0774** (medium, bug, easy/small) `numerics/functions_2.rs:54` — gauss_legendre_nw silently returns 5-point rule for n=7,8,10
- [ ] **F0780** (medium, bug, easy/trivial) `ode/types.rs:321` — richardson_extrapolate overwrites Neville table entries before use
- [ ] **F0786** (medium, fabrication, easy/small) `optimization/functions.rs:793` — PSO/GA/DE always report converged=true; empty populations panic
- [ ] **F0787** (medium, bug, medium/small) `optimization/functions.rs:942` — lbfgsb and constrained_min_box never converge at an active bound
- [ ] **F0789** (medium, perf, easy/small) `optimization/types.rs:158` — DifferentialEvolution re-evaluates parent fitness each generation; dim=0 panics
- [ ] **F0790** (medium, perf, easy/trivial) `optimization/types.rs:472` — CMA-ES step re-sorts samples with fresh objective calls, discards result
- [ ] **F0791** (medium, fabrication, medium/small) `optimization/types.rs:523` — TrustRegion 'dogleg' step is fixed full-radius steepest descent
- [ ] **F0793** (medium, perf, medium/large) `parallel/functions.rs:76` — core parallel helpers spawn OS threads per call; solvers have no parallelism
- [ ] **F0795** (medium, perf, easy/small) `parallel/functions.rs:117` — parallel_* helpers deep-copy input and spawn OS threads per call
- [ ] **F0797** (medium, panic, easy/medium) `parallel/functions.rs:296` — NaN-as-Equal comparators make std sort panic or mis-sort
- [ ] **F0800** (medium, missing-feature, medium/medium) `parallel_orchestrator.rs:256` — ParallelOrchestrator runs every wave sequentially; reschedules each call
- [ ] **F0801** (medium, bug, easy/trivial) `pde/types/types_2.rs:138` — NsBurgers1D uses central convective flux and wrong transonic entropy fix
- [ ] **F0803** (medium, bug, medium/small) `pde/types/types_2.rs:846` — FiniteDifference1D implicit/CN freeze Neumann/Periodic/Absorbing boundaries
- [ ] **F0804** (medium, pure-rust, medium/medium) `pde/types/types_3.rs:23` — Hand-rolled O(N²) DFTs and multiple radix-2 copies in core
- [ ] **F0805** (medium, bug, medium/small) `pde/types/types_impl.rs:39` — Multigrid V-cycle uses mismatched transfers and undamped Jacobi
- [ ] **F0806** (medium, pure-rust, easy/medium) `pde/types/types_impl.rs:128` — 'FFT-based' Poisson/spectral solvers use O(N²) DFT; adopt oxifft
- [ ] **F0807** (medium, bug, easy/trivial) `pde/types/types_impl.rs:220` — laplacian_9pt wrong for dx≠dy; dead code; mislabelled Mehrstellen
- [ ] **F0809** (medium, bug, medium/small) `pde/types/types_impl.rs:337` — Dct1D::solve_poisson_dct says Dirichlet but solves Neumann problem
- [ ] **F0811** (medium, bug, medium/medium) `pde/types/types_impl.rs:703` — HeatEquation3D periodic BC is fake; explicit step freezes non-Dirichlet faces
- [ ] **F0812** (medium, bug, easy/trivial) `pde/types/types_impl.rs:950` — FiniteVolume1D::godunov_flux is a copy of roe_flux
- [ ] **F0813** (medium, bug, easy/trivial) `pde/types/types_impl.rs:963` — muscl_reconstruct leaves interface states u_l[1] and u_r[n-1] zero
- [ ] **F0814** (medium, perf, easy/small) `persistent_homology.rs:219` — boundary_column linear-scans all earlier simplices per face (O(n²))
- [ ] **F0817** (medium, bug, medium/medium) `persistent_homology.rs:675` — bottleneck_distance: greedy, wrong metric coordinates, mixes homology dimensions
- [ ] **F0818** (medium, bug, easy/small) `petri_nets.rs:235` — reachability_bfs unbounded: hangs/OOMs on unbounded Petri nets
- [ ] **F0819** (medium, perf, easy/small) `petri_nets.rs:248` — reachability_bfs clones whole net per edge; never terminates on unbounded nets
- [ ] **F0820** (medium, bug, medium/small) `petri_nets.rs:325` — Coverability tree loses OMEGA via 1_000_000 substitution, then hangs
- [ ] **F0821** (medium, bug, medium/medium) `petri_nets.rs:1322` — t_invariants drops mixed-sign basis rows; misses minimal semiflows
- [ ] **F0824** (medium, panic, easy/trivial) `probabilistic_models.rs:369` — baum_welch panics on empty observation sequence
- [ ] **F0825** (medium, bug, easy/small) `probabilistic_models.rs:595` — GaussianProcess::fit masks non-PD kernel matrix via pivot clamping
- [ ] **F0826** (medium, fabrication, medium/medium) `probabilistic_models.rs:731` — DirichletProcess CRP ignores data; stick-breaking weights are fabricated
- [ ] **F0827** (medium, fabrication, medium/small) `probabilistic_models.rs:739` — DirichletProcess::crp_assign ignores the data point x
- [ ] **F0833** (medium, bug, medium/small) `random/functions.rs:460` — sample_truncated_normal unbounded rejection loop hangs on tails/invalid intervals
- [ ] **F0838** (medium, bug, easy/small) `random/types.rs:1291` — LevyProcess::sample CMS has flipped skew, wrong cosine, wrong α=1
- [ ] **F0841** (medium, bug, medium/small) `random_processes.rs:96` — fBm Davies-Harte misses c[n] embedding entry; O(n²) DFTs
- [ ] **F0844** (medium, fabrication, hard/medium) `random_processes.rs:721` — CGMY simulation uses invented jump intensity and exponential jump sizes
- [ ] **F0845** (medium, bug, medium/small) `random_processes.rs:758` — alpha-stable characteristic exponent mixes imaginary term into real part
- [ ] **F0846** (medium, bug, easy/trivial) `renormalization_group.rs:166` — landau_pole inverted sign condition and spurious factor 2
- [ ] **F0847** (medium, fabrication, medium/small) `renormalization_group.rs:682` — conformal_block_2d has wrong series coefficient and single holomorphic factor
- [ ] **F0851** (medium, panic, easy/small) `sampling/functions.rs:275` — Rejection samplers loop forever when target is zero or NaN
- [ ] **F0855** (medium, bug, easy/trivial) `sampling/functions_2.rs:574` — r2_sequence_nd uses plastic constant for every d; 3-D points coplanar
- [ ] **F0857** (medium, bug, easy/trivial) `sampling/functions_2.rs:836` — sampling::sample_kurtosis tends to 0 for every distribution
- [ ] **F0859** (medium, bug, easy/small) `sampling/types.rs:19` — HammersleySequence underflows at d=0; zero-pads dims beyond 10
- [ ] **F0860** (medium, bug, easy/trivial) `sampling/types.rs:218` — SobolSequence dim-3 direction numbers non-standard; poor 2-D stratification
- [ ] **F0864** (medium, pure-rust, medium/large) `signal/functions.rs:39` — Five hand-rolled radix-2 FFTs in core; oxifft not used
- [ ] **F0867** (medium, bug, easy/small) `signal/functions.rs:545` — periodogram_psd mislabels frequencies for non-power-of-two lengths
- [ ] **F0868** (medium, panic, easy/trivial) `signal/types.rs:187` — WelchPsd::compute never terminates when hop == 0
- [ ] **F0871** (medium, api, easy/small) `signal/types.rs:671` — KalmanFilterND::update silently drops measurement; no dimension checks
- [ ] **F0873** (medium, bug, easy/small) `signal_analysis.rs:272` — PowerSpectralDensity::welch returns an all-zero PSD when window > n/2
- [ ] **F0875** (medium, panic, easy/small) `signal_processing.rs:166` — Pow2-only hand-written FFTs panic via public APIs; oxifft unused in core
- [ ] **F0877** (medium, bug, easy/small) `signal_processing.rs:452` — welch_psd/stft/spectrogram map Rectangular etc. to Hann; step underflows
- [ ] **F0878** (medium, bug, easy/trivial) `signal_processing.rs:593` — FirFilter::highpass inverts the wrong tap for even n_taps
- [ ] **F0879** (medium, bug, easy/trivial) `signal_processing.rs:778` — butterworth_lowpass odd-order section uses c=1 instead of c=2
- [ ] **F0880** (medium, bug, easy/small) `signal_processing.rs:846` — chebyshev1_lowpass has a wrong odd section and missing passband gain
- [ ] **F0889** (medium, bug, easy/small) `sparse/types.rs:210` — BlockJacobi stores the un-inverted block when it is singular
- [ ] **F0892** (medium, perf, medium/medium) `sparse/types.rs:531` — CsrMatrix::incomplete_cholesky is O(n^2) through get() column probes
- [ ] **F0894** (medium, perf, medium/medium) `sparse/types.rs:781` — IncompleteCholesky apply O(n²) via get(); factor dense; solve_upper wrong
- [ ] **F0897** (medium, bug, easy/trivial) `spatial/types.rs:447` — BallTree range_query reports each leaf pivot twice
- [ ] **F0898** (medium, perf, easy/small) `spatial/types.rs:1169` — Voxel/grid range_query loops (2r/cell+1)^3 cells; centroid order random
- [ ] **F0901** (medium, bug, easy/small) `statistical_mechanics.rs:132` — Partition function lacks ground-state shift; underflow gives inf/0/NaN
- [ ] **F0902** (medium, perf, easy/small) `statistical_mechanics.rs:708` — GreenKubo::autocorrelation is an O(n^2) direct sum
- [ ] **F0903** (medium, bug, easy/trivial) `statistics/functions.rs:118` — skewness/kurtosis use Bessel std inside population moments
- [ ] **F0906** (medium, bug, easy/small) `statistics/functions.rs:501` — pacf reuses PACF diagonals as AR coefficients instead of Durbin-Levinson
- [ ] **F0912** (medium, bug, medium/small) `statistics/types.rs:97` — PoissonDistribution::sample saturates near 745 for large lambda
- [ ] **F0914** (medium, bug, easy/small) `stochastic/functions.rs:381` — Kramers MFPT divides the barrier by the diffusion coefficient, not kT
- [ ] **F0917** (medium, bug, easy/small) `stochastic/types.rs:242` — step_baoab final half-kick reuses the stale force
- [ ] **F0918** (medium, bug, easy/trivial) `stochastic/types.rs:270` — KleinmanKramers noise amplitude is off by mass and violates equipartition
- [ ] **F0921** (medium, bug, easy/small) `stochastic/types.rs:1104` — LCGs seeded raw; consecutive-seed paths strongly correlated
- [ ] **F0928** (medium, bug, easy/small) `symbolic_algebra/functions.rs:273` — taylor_expand fabricates zero coefficients; differentiates an unsimplified tree
- [ ] **F0929** (medium, bug, medium/medium) `symbolic_algebra/functions.rs:293` — CSE keys on 6-decimal Display string; merges distinct constants
- [ ] **F0931** (medium, panic, easy/small) `symbolic_algebra/functions.rs:728` — sturm_sequence panics on a constant polynomial
- [ ] **F0932** (medium, perf, medium/medium) `tensor/decomposition.rs:107` — general_einsum enumerates all element pairs with per-pair allocations
- [ ] **F0933** (medium, panic, easy/small) `tensor/operations.rs:56` — einsum_2d panics on unknown notation, missing operand, empty input
- [ ] **F0936** (medium, bug, easy/small) `tensor/operations.rs:301` — truncated_svd_left can return the wrong singular subspace
- [ ] **F0937** (medium, perf, medium/small) `tensor/operations.rs:388` — Tucker core uses O(R³N³) loop with allocating DenseTensor::get
- [ ] **F0939** (medium, bug, easy/trivial) `tensor/operations.rs:681` — compliance_from_stiffness returns Kelvin matrix documented as Voigt
- [ ] **F0940** (medium, stub, medium/small) `tensor/types.rs:993` — TensorTrain::from_dense ignores documented tol truncation tolerance
- [ ] **F0943** (medium, bug, medium/small) `time_series_analysis.rs:349` — ADF test uses the no-constant critical value and has no augmentation lags
- [ ] **F0944** (medium, bug, easy/small) `time_series_analysis.rs:371` — find_peaks prominence uses whole-side minima and is O(n²)
- [ ] **F0945** (medium, perf, easy/small) `time_series_analysis.rs:409` — compute_power_spectrum is O(n²) direct DFT; oxifft unused
- [ ] **F0946** (medium, perf, easy/small) `time_series_analysis.rs:420` — compute_power_spectrum is an O(n²) direct DFT; use oxifft::rfft
- [ ] **F0949** (medium, bug, medium/small) `time_series_analysis.rs:740` — Savitzky-Golay edge samples wrong (assumes symmetric window)
- [ ] **F0950** (medium, bug, easy/trivial) `time_series_analysis.rs:819` — ArModel::forecast ignores the series mean removed by Yule-Walker
- [ ] **F0951** (medium, bug, easy/trivial) `topological_data.rs:168` — vietoris_rips_h0 emits one infinite bar even when max_r leaves several components
- [ ] **F0953** (medium, bug, medium/medium) `topological_data.rs:286` — Greedy diagram matching documented as exact; topology.rs version asymmetric
- [ ] **F0956** (medium, bug, easy/small) `topology.rs:445` — CubicalComplex β1 wrong for forests and filled squares
- [ ] **F0957** (medium, fabrication, medium/small) `topology.rs:488` — MorseComplex::morse_index counts lower-valued points instead of Hessian eigenvalues
- [ ] **F0958** (medium, bug, easy/small) `topology.rs:673` — Čech complex filters triangles by circumradius, not minimum enclosing ball
- [ ] **F0959** (medium, bug, medium/medium) `topology.rs:876` — bottleneck_dist asymmetric: unmatched second-diagram points never counted
- [ ] **F0960** (medium, bug, medium/medium) `topology.rs:891` — bottleneck_dist never charges unmatched points of the second diagram
- [ ] **F0961** (medium, perf, easy/small) `topology_optimization.rs:211` — Filter neighbour precomputation is O(n²) over all element pairs
- [ ] **F0962** (medium, bug, easy/trivial) `topology_optimization.rs:239` — Sensitivity filter divides by ΣHρ_f, contradicting its doc and Sigmund
- [ ] **F0963** (medium, missing-feature, hard/large) `topology_optimization.rs:411` — Core topology optimization has no FE solver; only mock physics
- [ ] **F0965** (medium, bug, medium/medium) `topology_optimization.rs:704` — HJ update labelled upwind uses forward differences; fake reinitialize
- [ ] **F0966** (medium, fabrication, medium/medium) `topology_optimization.rs:881` — overhang_filter: angle unused in dirs 1–3; dir-0 formula meaningless
- [ ] **F0967** (medium, fabrication, hard/medium) `topology_optimization.rs:892` — Overhang filter hard-codes rho*0.7, ignores overhang angle
- [ ] **F0969** (medium, fabrication, medium/small) `traits.rs:168` — RungeKutta4 is constant-acceleration kinematics mislabeled as RK4
- [ ] **F0970** (medium, fabrication, medium/medium) `traits.rs:186` — RK4, ImplicitEuler, AB2 and VelocityVerlet integrators reduce to constant-force updates
- [ ] **F0972** (medium, bug, easy/trivial) `types/types.rs:183` — Ritter bounding-sphere centre update overshoots by +0.5
- [ ] **F0976** (medium, bug, easy/trivial) `types/types_impl.rs:183` — Sphere::bounding_sphere Ritter step overshoots centre (radius up to 2×)
- [ ] **F0978** (medium, bug, easy/trivial) `types/types_impl.rs:468` — Transform::from_axis_angle yields NaN rotation for zero axis
- [ ] **F0979** (medium, api, easy/small) `types/types_impl.rs:780` — PhysicsConfig/TimeStep validation passes NaN; builder never validates
- [ ] **F0980** (medium, panic, easy/trivial) `uncertainty_quantification.rs:101` — percentile panics for p > 100 and extrapolates for p < 0
- [ ] **F0983** (medium, bug, medium/small) `uncertainty_quantification.rs:798` — gauss_hermite_nodes silently falls back to 3 points for other orders
- [ ] **F0986** (medium, bug, easy/trivial) `uncertainty_quantification.rs:1081` — FORM β is always ≥0, so Pf is wrong when the origin is in the failure domain
- [ ] **F0990** (medium, bug, medium/small) `wavelet_transform/functions.rs:337` — Best-basis compares to children's own cost with non-additive entropy
- [ ] **F0991** (medium, bug, medium/medium) `wavelet_transform/functions.rs:385` — Real-only CWT makes scalogram oscillate and coherence ill-defined
- [ ] **F0992** (medium, fabrication, medium/medium) `wavelet_transform/functions.rs:410` — icwt hardcodes admissibility constant c_psi = 1.0
- [ ] **F0993** (medium, bug, easy/trivial) `wavelet_transform/functions.rs:841` — Mexican-hat centre frequency ~3x wrong; dt double-counted
- [ ] **F0994** (medium, fabrication, medium/medium) `wavelets.rs:470` — WaveletPacket best_basis_entropy returns a single minimum leaf, not a basis
- [ ] **F0995** (medium, bug, medium/medium) `wavelets.rs:523` — Morlet with omega0=1 is a low-pass filter; frequency mapping wrong
- [ ] **F0996** (medium, perf, medium/small) `wavelets.rs:563` — ContinuousWavelet CWT is direct O(S·N²) with exp+cos per pair
- [ ] **F0997** (medium, bug, medium/small) `wavelets.rs:636` — CWT scale-to-frequency conversion uses wrong units; Morlet ω0=1 is not admissible
- [ ] **F0999** (medium, bug, easy/small) `wavelets.rs:732` — Haar DWT/IDWT mismatch for non-power-of-two n; denoise thresholds one band
- [ ] **F1001** (medium, bug, easy/trivial) `wavelets.rs:880` — denoise_signal thresholds only the finest detail level
- [ ] **F1002** (medium, bug, easy/trivial) `wavelets.rs:1111` — turbulence_spectrum_wavelet assigns frequencies to levels in reverse order
- [ ] **F1003** (medium, stub, medium/medium) `world/functions.rs:169` — resolve_contact ignores friction, penetration depth and angular contact velocity
- [ ] **F1004** (medium, bug, medium/medium) `world/functions.rs:219` — resolve_contact uses linear-only relative velocity; ignores friction and depth
- [ ] **F1005** (medium, bug, easy/trivial) `world/functions.rs:322` — CCD returns the exit time for overlapping spheres and None for penetrating spheres vs planes
- [ ] **F1006** (medium, fabrication, easy/small) `world/types.rs:160` — Body inertia hardcoded mass*0.1; active and generation fields dead
- [ ] **F1008** (medium, bug, easy/small) `world/types.rs:474` — SpatialHashGrid accepts cell_size 0 (hang); nondeterministic pairs
- [ ] **F1016** (medium, bug, easy/small) `world/types.rs:991` — Sleeping bodies never auto-wake; wake branch is dead code
- [ ] **F1017** (medium, bug, easy/small) `world/types.rs:1020` — compute_islands recursive find overflows stack; island order random

### Low (70) — tracked in the external ledger
F0516, F0527, F0530, F0531, F0544, F0547, F0548, F0561, F0563, F0564, F0575, F0587, F0607, F0608, F0617, F0619, F0621, F0622, F0624, F0627, F0630, F0651, F0665, F0666, F0681, F0683, F0685, F0686, F0687, F0703, F0705, F0720, F0730, F0773, F0777, F0778, F0782, F0784, F0794, F0796, F0799, F0810, F0816, F0823, F0836, F0854, F0861, F0862, F0863, F0874, F0884, F0890, F0895, F0913, F0922, F0923, F0924, F0938, F0941, F0942, F0947, F0955, F0971, F0975, F0977, F0981, F0982, F0985, F1010, F1012

<!-- AUDIT-2026-10-12 END -->
