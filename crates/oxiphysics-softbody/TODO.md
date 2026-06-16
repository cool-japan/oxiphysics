# oxiphysics-softbody TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 88,223 SLoC | 3,473 tests

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation

- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation

- [x] Implement primary algorithms
- [x] Add integration tests
- [x] Performance benchmarks

### Phase 3: Polish

- [x] Documentation
- [x] Examples
- [x] Optimization

### Implemented modules

- [x] Cloth and advanced cloth simulation
- [x] PBD constraint system (distance, bending, volume)
- [x] XPBD solver and XPBD rope
- [x] FEM soft bodies with co-rotational elements
- [x] Rope and RopeSolver
- [x] Cosserat rod model
- [x] Fracture dynamics and crack propagation
- [x] Inflatable bodies and volume preservation
- [x] Muscle simulation
- [x] Numerical soft body integrators
- [x] Particle system
- [x] Surgery simulation and tissue simulation
- [x] Aerodynamics model
- [x] Shape matching
- [x] Hair strands and hair system (fur)
- [x] Membranes and textile / wrinkling
- [x] Metamaterials and origami
- [x] Active matter and bio/biomechanics
- [x] Morphogenesis simulation
- [x] Food physics
- [x] Haptics integration
- [x] Neural deformation
- [x] Position-Based Fluids (PBF)
- [x] Smart materials
- [x] Topology optimization
- [x] Tendon simulation
- [x] Elastic waves and peridynamics
- [x] 3,473 tests passing, 0 stubs

### Code-verified baseline details (2026-06-11)

- [x] MLS-MPM with APIC transfer — `material_point/` ships `MpmGrid` and `MpmParticle` carrying the APIC affine C matrix ("B-matrix in MLS-MPM", `material_point/types.rs:212`)
- [x] MPM elastoplasticity — `SnowParams` (Stomakhin elastoplastic snow) and `drucker_prager_return_mapping` (sand), plus `g2p_transfer_with_dt`
- [x] Small-steps XPBD substepping — `SubstepConfig` (`xpbd/types.rs:705`)
- [x] Classic Neo-Hookean FEM material — `NeoHookeanMaterial` (`fem_soft/materials.rs:252`)
- [x] Co-rotational rotation extraction — approximate iterative "SVD-like" polar decomposition (`fem_soft/corotational.rs:66`); exact branch-free SVD replacement tracked in v0.2.0
- [x] Hair implemented as strand/constraint mass-spring (`SimHairStrand`, `HairConstraint` in `hair.rs`) alongside `cosserat_rods.rs`; Discrete Elastic Rods v2 tracked in v0.3.0

## v0.2.0 — MPM v2 and robust FEM

The MPM headline ("add MLS-MPM") is already shipped — see baseline details above. This milestone makes the existing core couple with rigid bodies, scale past the dense single-thread grid, and survive element inversion.

### Material point method v2

- [ ] CPIC rigid-body coupling for MPM
  - **Goal:** rigid wheel driving on MPM sand with two-way momentum exchange, total momentum error <1% per second.
  - **Design:** Compatible PIC (Hu et al. 2018 MLS-MPM paper): colored grid nodes by rigid-surface side, ghost-velocity projection, impulse feedback to oxiphysics-rigid bodies.
  - **Cross-crate:** two-way impulse exchange with oxiphysics-rigid; prerequisite for the deferred vehicle terramechanics handoff.

- [ ] Sparse blocked MPM grid + parallel P2G/G2P
  - **Goal:** 250k particles ≥30 steps/s on 8 cores (current dense single-thread grid).
  - **Design:** 4³ blocks with occupancy map; rayon (new dep) scatter with block-level coloring to avoid write conflicts; SoA particle layout.
  - **Files:** `material_point/` (grid and transfer paths)

### FEM and cloth robustness

- [x] Invertible FEM elements via true 3×3 SVD
  - **Goal:** fully inverted tet mesh recovers to rest shape; no NaN at J ≤ 0.
  - **Design:** Irving et al. 2004 / Stomakhin 2012 rotation-variant SVD with reflection convention; `fem_soft/corotational.rs:66` currently uses an approximate iterative "SVD-like" polar decomposition — replace with exact branch-free SVD (McAdams 2011).
  - **Files:** `fem_soft/corotational.rs`

- [x] Stable Neo-Hookean
  - **Goal:** Smith et al. 2018 energy as both FEM element and XPBD constraint; no element locking at ν=0.499; classic vs stable energies agree in small strain to 1%.
  - **Design:** add to `fem_soft/materials.rs` (classic `NeoHookeanMaterial` exists at line 252) and as XPBD two-block constraint (deviatoric+hydrostatic, Macklin & Müller 2021).
  - **Files:** `fem_soft/materials.rs`, XPBD constraint set

- [x] Strain limiting for cloth
  - **Goal:** stretch ≤105% of rest under 100x gravity hang.
  - **Design:** per-triangle SVD clamp (Wang 2010) as post-XPBD pass, slotted in per substep — small-steps XPBD substepping already exists (`SubstepConfig`, `xpbd/types.rs:705`; see baseline details).

## v0.3.0 — Guaranteed contact and implicit solvers

Contact with guarantees: both IPC and continuous self-collision consume the cubic-root VF/EE CCD scheduled for oxiphysics-collision v0.2.0.

- [ ] IPC (Incremental Potential Contact)
  - **Goal:** intersection-free guarantee on the IPC benchmark trio (twisted rods, cloth funnel, stacked mats) — zero interpenetrations at any step.
  - **Design:** Li et al. 2020: log-barrier potential + CCD-filtered line search (consumes the new cubic-root VF/EE CCD from oxiphysics-collision), projected-Newton with sparse Hessian (internal CG or scirs2-sparse optional).
  - **Cross-crate:** hard dependency on oxiphysics-collision v0.2.0 "vertex-face/edge-edge continuous CCD with cubic root solver".

- [ ] Projective Dynamics solver
  - **Goal:** 10k-vertex deformable at 60 FPS single-thread with visually plausible convergence (10 local/global iters).
  - **Design:** Bouaziz et al. 2014: local SVD projections + prefactorized global sparse Cholesky (oxiblas new optional dep, or internal LDL).

- [ ] Discrete Elastic Rods for hair v2
  - **Goal:** curly strand (helical rest config) matches DER reference curvature dynamics; 1k strands real-time.
  - **Design:** Bergou 2008/2010 parallel-transport frames, curvature-binormal bending; complements existing `cosserat_rods.rs`; current `hair.rs` is strand/constraint mass-spring (verified `SimHairStrand`, `HairConstraint`).
  - **Files:** `hair.rs`, `cosserat_rods.rs`, new DER module

- [ ] Continuous self-collision + untangling
  - **Goal:** cloth crumple-and-release test ends intersection-free.
  - **Design:** CCD pass (collision crate) + intersection-contour-minimization untangling (Volino & Magnenat-Thalmann 2006 / Baraff 2003 global analysis).
  - **Cross-crate:** same oxiphysics-collision VF/EE CCD as IPC above.

- [ ] Graph-coloring parallel Gauss-Seidel + GPU cloth kernels
  - **Goal:** 65k-particle cloth XPBD step <4 ms on GPU; CPU colored GS ≥4x on 8 cores with identical convergence.
  - **Design:** greedy constraint coloring; WGSL distance/bending kernels via oxiphysics-gpu (wgpu default, cudarc feature-gated non-default).
  - **Cross-crate:** kernels hosted via oxiphysics-gpu.

- [ ] PolyPIC transfer option
  - **Goal:** ≤50% less angular-momentum dissipation than APIC on a spinning-blob test.
  - **Design:** Fu et al. 2017 polynomial PIC modes atop existing APIC transfer (`material_point/types.rs:212`).

## v1.0 — Production & validation

Quantitative validation against analytic and experimental references, robustness fuzzing, and bit-exact reproducibility.

- [ ] Physical validation battery
  - **Goal:** FEM cantilever tip deflection vs Euler-Bernoulli <2%; cloth drape vs catenary profile <2%; MPM sand-column collapse runout vs Lube 2004 experiments <10%; PBF dam-break wavefront vs SPH reference <5%.
  - **Design:** scenario tests with stored tolerances, results in CHANGELOG-tracked baseline file (oxicode).

- [ ] Robustness/fuzzing
  - **Goal:** 1e5 proptest steps over random meshes/params with zero NaN/panic; inverted/degenerate (zero-volume) elements rejected via `error.rs`, never UB.
  - **Design:** proptest generators for tet/tri meshes; NaN-guard debug assertions promoted to checked errors.

- [ ] Determinism + state snapshot
  - **Goal:** serde/oxicode round-trip mid-simulation resumes bit-identically for XPBD/FEM/MPM.
  - **Design:** serialize full solver state (λ accumulators, grid, RNG seeds).

## Deferred / research track

- [~] Vehicle terramechanics MPM handoff — rigid wheel on MPM sand as the high-fidelity tier behind oxiphysics-vehicle's Bekker-Wong terramechanics (vehicle v0.3.0) — Ready when: CPIC rigid-body coupling and the sparse blocked MPM grid (both v0.2.0 above) have shipped and oxiphysics-vehicle defines its per-wheel contact-patch interface.

- [~] Learned constitutive models beyond the existing self-contained `neural_deformation` module — Ready when: scirs2-neural exposes a stable inference API usable as an optional, feature-gated Pure-Rust dependency.

---

History: the pre-2026-06-11 TODO (last updated 2026-06-06) is preserved in full under "Completed"; the forward roadmap derives from the code-verified dynamics design brief of 2026-06-11.
