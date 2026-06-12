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
