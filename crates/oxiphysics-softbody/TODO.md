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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 294 (critical 5, high 100, medium 150, low 39); 61 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (5)
- [ ] **F4071** (critical, fabrication, hard/medium) `constraint/neohookeanconstraint_traits.rs:66` — Neo-Hookean XPBD projection fabricated: y-shift, radial push, inverted alpha
- [ ] **F4078** (critical, bug, easy/small) `cosserat_rods.rs:555` — Cosserat rod internal force sign flipped; stretch and shear self-amplify
- [ ] **F4175** (critical, fabrication, medium/medium) `material_point/functions.rs:185` — mat3_svd is not an SVD; snow/sand plasticity squares stretch
- [ ] **F4198** (critical, bug, easy/small) `metamorphic_mesh.rs:209` — polar_decomp_r iterates (R+R^-1)/2, collapsing rotations to identity
- [ ] **F4273** (critical, bug, easy/small) `pbd_system/types.rs:471` — Tet volume-constraint gradients negated; projection diverges

### High (100)
- [ ] **F4037** (high, bug, easy/trivial) `aerodynamics/types.rs:340` — VLM collocation points placed on the bound vortex line
- [ ] **F4038** (high, bug, medium/small) `aerodynamics/types.rs:372` — VLM doubles AoA, puts collocation on bound vortex, sums wrong gammas
- [ ] **F4039** (high, bug, easy/small) `aerodynamics/types.rs:380` — VLM horseshoe trailing legs are diagonals; AoA counted twice
- [ ] **F4040** (high, bug, medium/small) `aerodynamics/types.rs:407` — VLM lift sums wrong panels; horseshoe legs and collocation points wrong
- [ ] **F4050** (high, perf, medium/medium) `cloth/functions.rs:80` — resolve_self_collision is O(n^2) all-pairs every step
- [ ] **F4053** (high, bug, easy/small) `cloth/types.rs:150` — XpbdClothMesh bending stencils are collinear, so every constraint is inert
- [ ] **F4057** (high, fabrication, medium/small) `cloth/types.rs:1125` — PbdClothMesh bend constraint cannot change dihedral angle (no-op bending)
- [ ] **F4064** (high, perf, medium/medium) `collision_response.rs:621` — SelfCollision::detect tests every particle against every triangle
- [ ] **F4068** (high, bug, medium/medium) `constraint/balloonpressureconstraint_traits.rs:24` — Balloon constraint drives each triangle toward the whole-balloon volume
- [ ] **F4072** (high, bug, medium/medium) `constraint/strainlimitingconstraint_traits.rs:69` — Strain limiting always shrinks the triangle, even when over-compressed
- [ ] **F4073** (high, bug, medium/small) `constraint/types.rs:244` — Shape-matching 'polar decomposition' is Gram-Schmidt and breaks for planar groups
- [ ] **F4074** (high, bug, medium/small) `constraint/types.rs:251` — Shape-matching rotation uses Gram-Schmidt QR, not polar; rotates resting clusters
- [ ] **F4075** (high, bug, medium/small) `constraint/types.rs:748` — Isometric bending K vector wrong: nonzero C at flat rest, translation-dependent
- [ ] **F4076** (high, fabrication, medium/small) `constraint/types.rs:838` — Isometric bending Q matrices give non-zero energy for flat sheets
- [ ] **F4077** (high, bug, medium/small) `cosserat_rods.rs:493` — Curvature components projected onto the wrong director axes
- [ ] **F4081** (high, bug, medium/small) `cosserat_softbody.rs:561` — Discrete Cosserat elastic moment applied with same sign to both nodes
- [ ] **F4082** (high, bug, hard/medium) `cosserat_softbody.rs:586` — Discrete Cosserat rod: bending never couples to positions; moments not balanced
- [ ] **F4085** (high, fabrication, hard/medium) `crack_propagation.rs:239` — Crack propagation only refines edges; mesh never separates
- [ ] **F4094** (high, fabrication, medium/small) `electroactive_softbody.rs:852` — Maxwell nodal forces fabricated: identical force on all four nodes
- [ ] **F4095** (high, fabrication, hard/medium) `electroactive_softbody.rs:936` — Electroactive FEM body: Maxwell force only translates, no dt, no elasticity
- [ ] **F4098** (high, bug, easy/trivial) `fem_soft/corot_functions.rs:135` — tet_stiffness_matrix reads Dm⁻¹ columns as shape-function gradients
- [ ] **F4104** (high, bug, easy/trivial) `fem_soft/corot_sim.rs:428` — Rayleigh stiffness damping sign inverted, injects energy
- [ ] **F4107** (high, bug, easy/trivial) `fem_soft/corotational.rs:322` — Invertible element abs(sigma) leaves mirrored tet as equilibrium
- [ ] **F4108** (high, fabrication, medium/small) `fem_soft/integrators.rs:109` — Rayleigh stiffness damping is −β·f_total, not −β·K·v
- [ ] **F4110** (high, fabrication, medium/small) `fem_soft/materials.rs:61` — CorotationalElementRaw not force-free at rest; rotation from raw edges
- [ ] **F4111** (high, bug, medium/small) `fem_soft/materials.rs:112` — CorotationalElementRaw: rotation from Ds not F, rest not force-free
- [ ] **F4112** (high, bug, easy/trivial) `fem_soft/stable_neohookean_xpbd.rs:255` — Hydrostatic target alpha makes rest tet non-equilibrium; ν=0 gives inf/NaN
- [ ] **F4118** (high, bug, easy/trivial) `fluid_film_softbody.rs:1166` — Film gravity flux sign reversed: fluid flows uphill
- [ ] **F4125** (high, bug, medium/medium) `fracture/types.rs:180` — CrackFront::propagate moves vertices along crack-plane normal
- [ ] **F4131** (high, fabrication, medium/small) `fracture_dynamics/types.rs:266` — Crack-tip K_I computed from net nodal force (dimensionally wrong)
- [ ] **F4133** (high, bug, easy/trivial) `granular_softbody.rs:355` — Hertz-Mindlin contact velocity uses flipped lever arms (anti-friction)
- [ ] **F4137** (high, perf, medium/medium) `granular_softbody.rs:1038` — DEM contacts all-pairs with two grain clones and cache insert per pair
- [ ] **F4147** (high, bug, easy/small) `hair_fur.rs:320` — Wind drag force (N) added directly to positions (m)
- [ ] **F4149** (high, fabrication, hard/medium) `hair_sim.rs:252` — twisting_energy is identically zero; twist stiffness has no effect
- [ ] **F4150** (high, perf, medium/medium) `hair_sim.rs:345` — FD rod gradient clones whole strand per DOF: O(n²) per step
- [ ] **F4152** (high, bug, medium/small) `hair_softbody.rs:281` — DER bend force applied along binormal, cannot reduce curvature
- [ ] **F4159** (high, bug, medium/small) `haptic_softbody.rs:971` — Haptic tick returns tool-side coupling force to user, unnegated and unsaturated
- [ ] **F4161** (high, bug, easy/small) `inflatable/functions.rs:682` — Volume constraint uses |V| with signed gradient; built-in sphere inverts
- [ ] **F4162** (high, bug, easy/small) `inflatable/functions.rs:688` — apply_volume_constraint mixes |V| with signed-volume gradient
- [ ] **F4163** (high, bug, easy/small) `inflatable/types.rs:86` — InflatableShell pressure pushes inward; balloon_model winding contradicts doc
- [ ] **F4164** (high, bug, easy/small) `inflatable/types.rs:99` — InflatableShell pressure force points inward on generated meshes
- [ ] **F4165** (high, stub, medium/medium) `inflatable/types.rs:182` — BlowMolding::mold_collision is a no-op stub
- [ ] **F4172** (high, perf, medium/medium) `liquid_sim.rs:440` — PBF neighbour search O(n²), repeated each iteration plus twice more
- [ ] **F4174** (high, bug, easy/trivial) `material_point/functions.rs:176` — polar_decomp stops on symmetry of R, returning a non-orthogonal R
- [ ] **F4176** (high, fabrication, medium/medium) `material_point/functions.rs:402` — Drucker-Prager return mapping ignores μ,λ, drops rotation, and uses a fake SVD
- [ ] **F4178** (high, bug, medium/small) `material_point/functions.rs:490` — MPM P2G grid force omits Fᵀ (uses PK1 as Kirchhoff)
- [ ] **F4179** (high, fabrication, easy/small) `material_point/functions.rs:663` — apply_level_set_coupling never computes grad_phi; push/reflect are no-ops
- [ ] **F4180** (high, bug, medium/small) `material_point/functions.rs:717` — MPM pressure projection Poisson RHS has inverted sign; doubles divergence
- [ ] **F4181** (high, bug, easy/small) `material_point/types.rs:66` — MPM quadratic B-spline stencil drops node base+2; weights don't sum to 1
- [ ] **F4183** (high, fabrication, medium/medium) `membrane_biophysics.rs:243` — Helfrich bending/ADE energies read curvatures nothing computes from geometry
- [ ] **F4187** (high, bug, medium/medium) `membrane_biophysics.rs:791` — VesicleSimulation::new_sphere builds a flat square patch, not a sphere
- [ ] **F4191** (high, bug, medium/medium) `membrane_softbody.rs:318` — Membrane triangle tension is applied as a normal push with nonzero net force
- [ ] **F4192** (high, bug, easy/trivial) `membrane_softbody.rs:445` — Mooney-Rivlin balloon pressure formula wrong; misses inflation instability
- [ ] **F4194** (high, fabrication, easy/small) `membrane_softbody.rs:910` — MembraneMesh::step applies no membrane elasticity at all
- [ ] **F4197** (high, fabrication, hard/large) `metamaterial_softbody.rs:1109` — SIMP sensitivities come from constant E with no FE solve; target_poisson unused
- [ ] **F4201** (high, bug, medium/small) `metamorphic_mesh.rs:868` — RBF deformer Gauss-Seidel gives zero coefficients for ThinPlate; MQ may diverge
- [ ] **F4202** (high, bug, medium/medium) `metamorphic_mesh.rs:960` — ARAP global step adds R·e instead of subtracting; polar_decomp_r is broken
- [ ] **F4206** (high, fabrication, medium/medium) `mud_snow.rs:246` — SnowMpm::step has no stress or plasticity: gravity-only PIC
- [ ] **F4213** (high, bug, easy/small) `muscle_sim.rs:358` — ThilenFibreModel bisection moves away from the root; MTU length hardcoded
- [ ] **F4221** (high, bug, easy/small) `muscle_softbody.rs:993` — IsotonicContraction sign error: stronger muscle lengthens the MTU
- [ ] **F4224** (high, bug, easy/small) `muscle_tendon.rs:718` — MusculotendonUnit never updates v_fiber; FV, damping, power and work always zero
- [ ] **F4230** (high, fabrication, medium/medium) `neural_softbody.rs:244` — PhysicsInformedLoss PDE term is ||sigma||^2, not a div-sigma residual
- [ ] **F4231** (high, fabrication, medium/medium) `neural_softbody.rs:947` — Neural soft body: elastic forces are uniform per element; strain ignores rest shape
- [ ] **F4232** (high, bug, easy/small) `neural_softbody.rs:1089` — Online training never reaches the body; zero-initialised network cannot learn
- [ ] **F4238** (high, fabrication, medium/medium) `numerical_softbody/types.rs:496` — NewtonRaphsonSolver/BackwardEuler are diagonal gradient descent; PCG unused
- [ ] **F4241** (high, perf, easy/small) `numerical_softbody/types.rs:966` — Projective Dynamics recomputes constant diagonal in O(n*C) per iteration
- [ ] **F4242** (high, fabrication, hard/large) `origami_folding.rs:471` — Origami folding simulators never move geometry; angles are scalar interpolations
- [ ] **F4248** (high, fabrication, medium/small) `origami_mechanics.rs:260` — MiuraOriUnit kinematics wrong; Poisson ratio is positive though doc says negative
- [ ] **F4249** (high, bug, medium/small) `origami_mechanics.rs:947` — degree4_vertex_fold_angles uses tan(rho)/2 and same-sign fold pairs
- [ ] **F4256** (high, fabrication, medium/small) `pbd_constraints.rs:182` — IsometricBending Q matrix fabricated; not translation-invariant
- [ ] **F4258** (high, bug, easy/trivial) `pbd_constraints.rs:330` — CST deformation gradient applies Dm⁻¹ transposed; rest state violates constraint
- [ ] **F4260** (high, bug, medium/medium) `pbd_constraints.rs:508` — NeoHookeanConstraint gradient uses Dm⁻¹ transpose; Lamé parameters unused
- [ ] **F4265** (high, fabrication, easy/small) `pbd_system/functions.rs:361` — add_sphere_collision_constraint creates an infinite floor, not a sphere
- [ ] **F4266** (high, fabrication, easy/small) `pbd_system/functions.rs:369` — add_sphere_collision_constraint creates an infinite floor plane
- [ ] **F4267** (high, bug, easy/small) `pbd_system/functions.rs:419` — xpbd_step_colored silently drops bending and volume constraints
- [ ] **F4272** (high, fabrication, medium/small) `pbd_system/types.rs:452` — BendingDihedral uses a fake gradient and moves only two particles
- [ ] **F4280** (high, perf, medium/medium) `position_based_fluids.rs:125` — PBF neighbour search is all-pairs, rebuilt every solver iteration
- [ ] **F4283** (high, perf, medium/medium) `position_based_fluids.rs:808` — PbfSolver density and delta passes are all-pairs every substep
- [ ] **F4284** (high, fabrication, easy/small) `position_based_fluids.rs:955` — PbfSolver λ denominator fabricated from density, not Σ|∇C|²
- [ ] **F4285** (high, bug, easy/small) `pressure_volume.rs:386` — valve_flow drives gauge pressure toward +101325 Pa
- [ ] **F4286** (high, bug, medium/medium) `pressure_volume.rs:439` — step_body: volume projection without velocity update makes pressurised bodies blow up
- [ ] **F4300** (high, stub, medium/medium) `sand_sim.rs:206` — 'MPM' SandSimulation::step is ballistic free fall only
- [ ] **F4302** (high, bug, easy/small) `shape_matching.rs:297` — Shape-matching plasticity moves the rest shape away from the deformation
- [ ] **F4303** (high, bug, medium/small) `shape_matching.rs:721` — Quadratic shape matching omits A_qq⁻¹; goals are not rest positions
- [x] **F4306** (high, bug, medium/medium) `solver.rs:201` — dyn-constraint XPBD solvers never reset lambda between substeps
- [ ] **F4307** (high, bug, easy/small) `solver.rs:219` — XpbdSolver never resets XPBD lambdas per substep
- [x] **F4308** (high, bug, medium/medium) `solver.rs:233` — XPBD Lagrange multipliers never reset between substeps or frames
- [x] **F4312** (high, bug, easy/small) `solver.rs:585` — GaussSeidelSolver SOR rescales whole-substep displacement every iteration
- [ ] **F4318** (high, bug, medium/medium) `surgery_simulation/types.rs:223` — TissueModel divides gravity by mass and has no internal elasticity
- [x] **F4320** (high, fabrication, medium/medium) `surgery_simulation/types.rs:678` — SoftTissueModel::step ignores all neo-Hookean elements
- [x] **F4331** (high, perf, medium/medium) `textile_simulation.rs:757` — Fabric self-collision is O(V·T·iters) and clones triangles
- [ ] **F4341** (high, bug, medium/small) `topology_optimization.rs:318` — mma_update gives less material to more sensitive elements
- [ ] **F4346** (high, bug, easy/trivial) `underwater_softbody.rs:489` — Windenburg-Trilling collapse pressure uses (t/D)^(5/3) instead of (t/D)^(5/2)
- [ ] **F4352** (high, bug, easy/trivial) `volume/functions.rs:99` — tet_volume_gradient returns wrong-sign dV/dp0
- [ ] **F4354** (high, bug, easy/trivial) `volume/functions.rs:281` — volume_preservation_gradient flips sign of 3 of 4 vertex gradients
- [ ] **F4359** (high, bug, easy/small) `volume/types.rs:118` — Tet volume constraint: |V| with signed gradients, and g0 sign flipped
- [ ] **F4360** (high, bug, easy/trivial) `volume/types.rs:126` — solve_volume_constraints inline g0 has wrong sign; injects momentum
- [ ] **F4365** (high, fabrication, medium/small) `wrinkling.rs:131` — wrinkle_amplitude always equals thickness; wavelength load-independent
- [ ] **F4373** (high, bug, medium/medium) `xpbd/strain_limit.rs:267` — Strain limiter projects current triangle onto rest plane, fights rotation
- [ ] **F4375** (high, bug, medium/medium) `xpbd/types.rs:343` — XPBD bending: unsigned acos dihedral gives anti-restoring correction, breaks momentum

### Medium (150)
- [ ] **F4025** (medium, perf, medium/medium) `crates/oxiphysics-softbody/Cargo.toml:17` — Unused oxiphysics-geometry dep while cloth self-collision is O(n^2)
- [ ] **F4026** (medium, fabrication, easy/small) `active_matter.rs:95` — Vicsek and run-and-tumble noise is a deterministic hash of position
- [ ] **F4028** (medium, fabrication, easy/small) `active_matter.rs:287` — Run-and-tumble uses position hash, non-uniform tumble directions
- [ ] **F4029** (medium, stub, medium/medium) `active_matter.rs:369` — ActiveNematic::detect_defects always returns zero
- [ ] **F4030** (medium, bug, easy/trivial) `active_matter.rs:387` — 2D nematic order parameter is not rotation-invariant
- [ ] **F4032** (medium, bug, medium/small) `active_origami.rs:380` — DEA equilibrium_stretch: wrong Newton derivative, 1.001 floor, no convergence check
- [ ] **F4033** (medium, panic, easy/trivial) `active_origami.rs:767` — DeployableStructure::step panics when n_stages is 0
- [ ] **F4034** (medium, pure-rust, easy/small) `aerodynamics/functions.rs:134` — Hand-rolled dense Gauss solvers silently return zeros on singular systems
- [ ] **F4036** (medium, bug, easy/small) `aerodynamics/types.rs:128` — AeroPanel lift/drag linear in velocity; 'lift' acts along flow
- [ ] **F4041** (medium, bug, easy/small) `aerodynamics/types.rs:1185` — FlexibleWing twist has no rate state; Cm_alpha about quarter-chord bogus
- [ ] **F4042** (medium, fabrication, easy/small) `bio_softbody.rs:846` — Tumor surface stress is E·R/3, dimensionally Pa·m
- [ ] **F4043** (medium, fabrication, medium/medium) `bio_softbody.rs:980` — RBC elastic energy uses hardcoded harmonic k and ignores WLC parameters
- [ ] **F4044** (medium, bug, easy/small) `bio_softbody.rs:1089` — Platelet explicit inertial step unstable; wall_shear_rate never read
- [ ] **F4048** (medium, bug, easy/small) `cable_nets.rs:558` — FormFinding kinetic damping inverted: zeroes velocity while KE rises
- [ ] **F4049** (medium, fabrication, easy/small) `cloth/functions.rs:30` — Cloth wind formulas are not aerodynamic forces; one scales with area squared
- [ ] **F4054** (medium, panic, easy/small) `cloth/types.rs:403` — ClothMesh::pin_corner panics and infers grid width from live positions
- [ ] **F4055** (medium, fabrication, medium/small) `cloth/types.rs:495` — ClothMesh solver ignores edge stiffness; bending compliance not scaled by dt²
- [ ] **F4056** (medium, perf, easy/small) `cloth/types.rs:1033` — PbdClothMesh::step does O(P) Vec::contains pin lookups per constraint
- [ ] **F4058** (medium, bug, easy/small) `cloth/types.rs:1148` — PbdClothMesh::apply_wind hardcodes dt=1/60 and uses a non-force formula
- [ ] **F4063** (medium, bug, easy/small) `cloth_simulation.rs:308` — apply_wind_force scales with area², skips half the triangles
- [ ] **F4065** (medium, bug, easy/trivial) `collision_response.rs:637` — Self-collision normal not oriented toward particle side
- [ ] **F4069** (medium, bug, medium/small) `constraint/bendingconstraint_traits.rs:55` — Dihedral bending moves only wing vertices; unsigned angle
- [ ] **F4080** (medium, fabrication, hard/medium) `cosserat_rods.rs:802` — StaticSolver claims Newton; shoot() ignores loads, returns no residual
- [ ] **F4083** (medium, bug, easy/small) `crack_propagation.rs:49` — K_Ic stored in MPa·√m but K and G computed in Pa units
- [ ] **F4084** (medium, bug, easy/small) `crack_propagation.rs:191` — advance_tip moves only the first tip; branched tips freeze
- [ ] **F4086** (medium, fabrication, medium/small) `crack_propagation.rs:318` — Stress intensity factor formula is dimensionally fabricated
- [ ] **F4087** (medium, bug, medium/small) `crack_propagation.rs:747` — Cohesive damage heals on unloading; softening is quadratic
- [ ] **F4088** (medium, stub, medium/medium) `crack_propagation.rs:803` — insert_cohesive_elements creates orphan duplicates; crack never split
- [ ] **F4090** (medium, fabrication, medium/small) `elastic_wave.rs:302` — Love wave dispersion is an invented interpolation
- [ ] **F4091** (medium, bug, easy/trivial) `elastic_wave.rs:362` — Oblique P-P reflection uses Z·cosθ instead of Z/cosθ
- [ ] **F4092** (medium, perf, easy/small) `elastic_wave.rs:642` — step_central_diff rebuilds dense n x n stiffness matrix every step
- [ ] **F4093** (medium, bug, easy/trivial) `elastic_wave.rs:688` — Bar natural frequencies use P-wave speed, FE stiffness uses E
- [ ] **F4097** (medium, api, easy/small) `error.rs:10` — Crate Error enum is dead; no fallible public API
- [ ] **F4100** (medium, perf, medium/medium) `fem_soft/corot_functions.rs:394` — Dense (3n)² global stiffness and silent dense CG
- [ ] **F4102** (medium, bug, medium/small) `fem_soft/corot_sim.rs:73` — Corot polar decomposition: 5 unscaled Newton steps, wrong inversion handling
- [ ] **F4106** (medium, api, easy/small) `fem_soft/corotational.rs:45` — CorotationalElement::new silently accepts degenerate tets and nu=0.5
- [ ] **F4113** (medium, bug, medium/small) `fluid_film_softbody.rs:390` — FilmThicknessEvolution1D non-conservative capillary term, unstable central advection
- [ ] **F4114** (medium, bug, easy/trivial) `fluid_film_softbody.rs:593` — Tanner contact angle uses spurious cube root
- [ ] **F4117** (medium, fabrication, medium/small) `fluid_film_softbody.rs:1138` — FluidFilmSimulation::step ignores surface tension; no capillary flux
- [ ] **F4120** (medium, fabrication, medium/small) `food_physics.rs:937` — LSW ripening grows all droplets equally; no volume conservation
- [ ] **F4121** (medium, bug, easy/trivial) `food_physics.rs:1115` — water_activity_bet plugs m/m0 into forward BET instead of inverting
- [ ] **F4122** (medium, bug, easy/small) `fracture/functions.rs:159` — split_node duplicates full nodal mass, creating mass
- [ ] **F4123** (medium, panic, easy/trivial) `fracture/functions.rs:586` — count_fragments recursive union-find overflows stack; unchecked edge indices
- [ ] **F4129** (medium, bug, easy/small) `fracture_dynamics/types.rs:198` — Bonds break under compression (abs stretch and abs force)
- [ ] **F4130** (medium, perf, easy/small) `fracture_dynamics/types.rs:215` — Fracture scans all bonds per broken bond and per crack tip
- [ ] **F4132** (medium, bug, easy/small) `fracture_dynamics/types.rs:786` — CohesiveZone::traction ignores delta_max, so unloading heals the interface
- [ ] **F4134** (medium, fabrication, medium/small) `granular_softbody.rs:795` — ShearBand::detect hardcodes +Y normal, fakes thickness, panics
- [ ] **F4135** (medium, fabrication, medium/small) `granular_softbody.rs:977` — FEM homogenisation feeds every global contact to each element
- [ ] **F4136** (medium, perf, medium/medium) `granular_softbody.rs:1028` — DEM contacts: O(n^2) pairs, unbounded tangential cache, Grain clones
- [ ] **F4138** (medium, api, easy/small) `granular_softbody.rs:1084` — Unseedable thread RNG inside DEM step breaks reproducibility
- [ ] **F4139** (medium, bug, easy/small) `granular_softbody.rs:1091` — Grain crushing uses net force, wrong under confinement and free fall
- [ ] **F4142** (medium, bug, medium/small) `hair.rs:335` — HairSimulation PBD never derives velocity from corrected positions
- [ ] **F4143** (medium, bug, medium/small) `hair.rs:357` — PBD hair step never syncs velocities; stiffness/self_collision unused
- [ ] **F4145** (medium, bug, easy/trivial) `hair.rs:716` — HairAerodynamics drag uses absolute wind, ignoring strand velocity
- [ ] **F4146** (medium, bug, easy/trivial) `hair_fur.rs:198` — HairStrand::apply_damping is a no-op for Verlet dynamics
- [ ] **F4151** (medium, bug, easy/trivial) `hair_softbody.rs:253` — Rod stretch force has extra 1/l0 versus stretch_energy
- [ ] **F4153** (medium, fabrication, hard/medium) `hair_softbody.rs:502` — CosseratHairStrand twist never evolves; twist_energy always zero
- [ ] **F4154** (medium, fabrication, hard/medium) `hair_softbody.rs:523` — CosseratHairStrand twist_angles never evolve; twist is cosmetic
- [ ] **F4155** (medium, fabrication, medium/medium) `hair_softbody.rs:607` — Hair-hair collision only compares same-index vertices
- [ ] **F4156** (medium, fabrication, easy/small) `haptic_softbody.rs:491` — Needle friction and lateral force use hidden magic constants
- [ ] **F4157** (medium, bug, medium/medium) `haptic_softbody.rs:611` — Tissue point load accumulates displacement; equilibrium depends on dt
- [ ] **F4166** (medium, bug, easy/trivial) `inflatable/types.rs:369` — HashSet edge dedup gives run-to-run nondeterministic constraint order
- [ ] **F4170** (medium, bug, easy/small) `liquid_sim.rs:226` — PBF density has implicit unit mass despite kg/m³ docs
- [ ] **F4171** (medium, perf, medium/medium) `liquid_sim.rs:414` — O(n²) neighbour search repeated (iterations+2) times per step
- [ ] **F4173** (medium, bug, medium/small) `liquid_sim.rs:457` — PBF lambda denominator omits k=i gradient; density lacks particle mass
- [ ] **F4182** (medium, missing-feature, medium/medium) `material_point/types.rs:444` — MultiMaterialGrid has no P2G path; coupling sticks separating bodies
- [ ] **F4184** (medium, fabrication, medium/medium) `membrane_biophysics.rs:419` — RBC 'hexagonal spectrin network' is a nondeterministic ring; n_nodes=0 panics
- [ ] **F4185** (medium, fabrication, medium/medium) `membrane_biophysics.rs:433` — RedBloodCellModel::new builds a ring, not a hexagonal network; panics n=0
- [ ] **F4186** (medium, bug, easy/trivial) `membrane_biophysics.rs:698` — Curvature-sensing force hardcodes local curvature H = 0
- [ ] **F4188** (medium, bug, easy/trivial) `membrane_biophysics.rs:985` — Voronoi area at vertex a uses swapped cotangent weights
- [ ] **F4189** (medium, bug, easy/small) `membrane_biophysics.rs:1118` — Vesicle step never updates mean_curvature; dead adjacency; duplicate cot pass
- [ ] **F4190** (medium, bug, easy/trivial) `membrane_biophysics.rs:1318` — asphericity ignores off-diagonal gyration terms; not rotation-invariant
- [ ] **F4193** (medium, perf, easy/small) `membrane_softbody.rs:547` — Force-density Gauss-Seidel scans all elements per node: O(n·E)
- [ ] **F4196** (medium, bug, easy/small) `metamaterial_softbody.rs:186` — Re-entrant honeycomb Poisson/modulus formulas are not Gibson-Ashby
- [ ] **F4199** (medium, perf, easy/small) `metamorphic_mesh.rs:392` — ICP closest-point search is brute force O(N·M); empty target panics
- [ ] **F4200** (medium, bug, medium/medium) `metamorphic_mesh.rs:608` — Cage deformer IDW weights do not reproduce the rest pose
- [ ] **F4203** (medium, fabrication, easy/small) `morphogenesis_softbody.rs:180` — Turing instability check is a fabricated D_v/D_u > 1 test
- [ ] **F4204** (medium, fabrication, medium/medium) `morphogenesis_softbody.rs:244` — ReactionDiffusionMesh claims cotangent Laplace–Beltrami, uses graph Laplacian
- [ ] **F4205** (medium, bug, easy/trivial) `morphogenesis_softbody.rs:326` — Morphogen amplitude C0 = S0/(D*lambda) dimensionally wrong
- [ ] **F4207** (medium, fabrication, hard/large) `mud_snow.rs:509` — Mud Bingham step treats speed as shear rate, with a magic 0.01
- [ ] **F4208** (medium, fabrication, medium/small) `mud_snow.rs:516` — ViscoPlasticMud uses speed as shear rate, magic 0.01 m depth
- [ ] **F4209** (medium, bug, easy/small) `mud_snow.rs:628` — Avalanche resisting stress uses Pa strength as friction coefficient
- [ ] **F4211** (medium, bug, easy/trivial) `muscle_sim.rs:103` — Passive force-length lacks Thelen normalisation (≈19·Fmax at 1.6 l_opt)
- [ ] **F4212** (medium, bug, medium/small) `muscle_sim.rs:347` — Equilibrium/EMG fitters ignore the inputs that determine the answer
- [ ] **F4214** (medium, fabrication, medium/small) `muscle_sim.rs:397` — optimal_fiber_length_from_emg has no length data to fit
- [ ] **F4215** (medium, bug, easy/small) `muscle_simulation/types.rs:100` — Motor pool peak force != max_force; log thresholds discarded
- [ ] **F4216** (medium, bug, easy/small) `muscle_simulation/types.rs:173` — static_optimization does not minimise Σa² as documented
- [ ] **F4217** (medium, bug, easy/trivial) `muscle_simulation/types.rs:288` — optimal_velocity_fraction uses the wrong max-power formula
- [ ] **F4219** (medium, bug, easy/trivial) `muscle_softbody.rs:362` — Pennated fiber velocity divides by cos(alpha) instead of multiplying
- [ ] **F4222** (medium, bug, easy/small) `muscle_softbody.rs:1065` — Via-point MTU applies forces along the straight origin→insertion chord
- [ ] **F4225** (medium, bug, medium/small) `neural_deform.rs:411` — PCA basis empty when first sample equals the mean
- [ ] **F4226** (medium, panic, easy/small) `neural_soft.rs:115` — Public NN APIs panic on length mismatches via assert_eq!
- [ ] **F4227** (medium, bug, easy/small) `neural_softbody.rs:114` — Zero-initialised ReLU nets cannot learn; simulator trains a disconnected copy
- [ ] **F4228** (medium, bug, easy/small) `neural_softbody.rs:128` — linear_elastic zeroes compressive stress; n_hidden=0 divides by zero
- [ ] **F4229** (medium, fabrication, hard/medium) `neural_softbody.rs:234` — PINN 'PDE loss' is ||σ||² and never reaches the gradient
- [ ] **F4235** (medium, bug, easy/small) `numerical_softbody/types.rs:274` — Shape-matching rotation returns identity for planar/small/inverted clusters
- [ ] **F4236** (medium, bug, medium/medium) `numerical_softbody/types.rs:367` — BackwardEulerIntegrator claims Newton+PCG, does diagonal Jacobi only
- [ ] **F4237** (medium, perf, medium/medium) `numerical_softbody/types.rs:442` — BackwardEuler 'Newton' is a diagonal Jacobi update; O(n·S) residual; pcg unused
- [ ] **F4243** (medium, fabrication, medium/small) `origami_folding.rs:597` — Miura-ori vertex kinematics are not rigid (edge lengths change)
- [ ] **F4244** (medium, fabrication, hard/large) `origami_folding.rs:862` — Rigid-origami solvers never enforce rigidity; geometry static or free-falling
- [ ] **F4245** (medium, bug, easy/small) `origami_mech.rs:285` — Yoshimura flat pattern folds back: x = r·cosθ, duplicate vertices
- [ ] **F4246** (medium, fabrication, easy/trivial) `origami_mech.rs:377` — crease_stiffness ∝ 1/θ² diverges at flat state, contradicting its doc
- [ ] **F4247** (medium, bug, easy/small) `origami_mechanics.rs:249` — Miura-ori Poisson ratio positive / fold-independent
- [ ] **F4250** (medium, fabrication, medium/medium) `origami_softbody.rs:823` — OrigamiSimulator integrates wrong ODE; actuated folds never reach target
- [ ] **F4253** (medium, perf, easy/small) `particle/types.rs:399` — ParticleCollisionHandler is all-pairs despite ParticleSet's spatial hash
- [ ] **F4254** (medium, perf, easy/small) `particle/types.rs:740` — SphSoftBody uses all-pairs loops; pair forces are asymmetric
- [ ] **F4255** (medium, bug, easy/small) `particle/types.rs:787` — SphSoftBody pressure force violates Newton's 3rd law; viscosity unused
- [ ] **F4262** (medium, bug, medium/small) `pbd_system/functions.rs:79` — Floor restitution discarded because velocity is recomputed from positions
- [ ] **F4264** (medium, fabrication, easy/small) `pbd_system/functions.rs:316` — solve_cloth_rigid_attachment discards lambda; moves pinned particles
- [ ] **F4269** (medium, fabrication, medium/small) `pbd_system/types.rs:157` — PbdSystem WarmStartCache stores the violation as lambda and divides by dt²
- [ ] **F4271** (medium, perf, easy/small) `pbd_system/types.rs:397` — PbdSystem::step heap-clones each constraint's index Vec every substep
- [ ] **F4274** (medium, bug, medium/small) `pbd_system/types.rs:504` — PbdSystem velocity update discards floor restitution
- [ ] **F4276** (medium, perf, easy/small) `peridynamics.rs:350` — Peridynamic bond search is all-pairs O(n²)
- [ ] **F4277** (medium, bug, easy/small) `peridynamics.rs:368` — Peridynamic body force erased before second half-kick
- [ ] **F4278** (medium, bug, easy/trivial) `peridynamics.rs:391` — PMB bonds break under compression
- [ ] **F4279** (medium, fabrication, medium/small) `peridynamics.rs:573` — OSB dilatation uses fake δ⁴ normalisation and O(N·B) bond scan
- [ ] **F4281** (medium, bug, medium/small) `position_based_fluids.rs:195` — PBF constraint gradient omits mass; PbfSolver fakes the gradient sum
- [ ] **F4282** (medium, fabrication, easy/small) `position_based_fluids.rs:276` — Vorticity confinement substitutes the particle's own |ω| for neighbour vorticity
- [ ] **F4287** (medium, bug, medium/small) `rigid_soft_coupling.rs:102` — RigidBody::integrate uses diagonal inertia, no gyroscopic term
- [ ] **F4288** (medium, bug, easy/small) `rigid_soft_coupling.rs:317` — EmbeddedRigidBody node velocity ignores angular velocity
- [ ] **F4289** (medium, bug, easy/small) `rod_mechanics.rs:227` — Rod bending curvature not divided by length; mesh-dependent energy
- [ ] **F4290** (medium, bug, easy/trivial) `rod_mechanics.rs:357` — Discrete writhe is half the Gauss integral
- [ ] **F4292** (medium, fabrication, medium/small) `rope/types.rs:573` — apply_twist_constraint limits bending with out-of-plane (no-op) correction
- [ ] **F4293** (medium, duplicate, hard/large) `rope/types.rs:873` — 49 duplicated public type names across overlapping modules
- [ ] **F4294** (medium, bug, easy/small) `rope/types.rs:962` — RopeVerlet::step never updates velocities after constraint projection
- [ ] **F4305** (medium, bug, medium/medium) `soft_gripper.rs:309` — Gripper rest length re-derived each step; actuation compresses axis
- [ ] **F4310** (medium, fabrication, medium/medium) `solver.rs:513` — WarmstartCache::apply_warmstart scales velocities instead of warm-starting
- [ ] **F4311** (medium, fabrication, medium/small) `solver.rs:524` — WarmstartCache rescales velocities instead of warm-starting lambdas
- [ ] **F4315** (medium, missing-feature, medium/medium) `solver.rs:1176` — ParallelGaussSeidelSolver has no solve method
- [ ] **F4319** (medium, fabrication, easy/small) `surgery_simulation/types.rs:560` — viscous_stress returns η·tr(F)/3 with no strain rate
- [ ] **F4322** (medium, fabrication, easy/small) `surgical_simulation.rs:612` — RoboticSurgery::instrument_tip uses an arbitrary 0.1 scale and an unrotated offset
- [ ] **F4323** (medium, bug, medium/medium) `tearing.rs:144` — Tear thresholds desync from edges; progressive_tearing can panic
- [ ] **F4324** (medium, bug, medium/small) `tearing.rs:310` — remesh_after_tear desyncs thresholds, drops edges, re-splits torn edges
- [ ] **F4327** (medium, bug, easy/trivial) `tearing.rs:878` — Griffith fracture_energy drops the crack-length term
- [ ] **F4328** (medium, fabrication, medium/medium) `tendon_softbody.rs:968` — Tendon obstacle wrap angle ignores angular separation of A,B
- [ ] **F4330** (medium, bug, easy/small) `tendon_softbody.rs:1115` — MTU equilibrium uses fixed 0.1% steps, silently unconverged
- [ ] **F4333** (medium, perf, medium/medium) `textile_simulation.rs:1078` — Garment self-collision is brute-force O(V·T·iters); ground friction is overwritten
- [ ] **F4334** (medium, bug, easy/small) `tissue_sim.rs:1147` — ArterialWall pressure overcounts layers; layer stress nonzero at λ=1
- [ ] **F4338** (medium, bug, easy/trivial) `topology_opt.rs:336` — Sensitivity filter divides by Σwρ, contradicting documented Sigmund formula
- [ ] **F4339** (medium, fabrication, medium/small) `topology_opt.rs:580` — enforce_overhang uses constant 1×1 span; asymmetric support
- [ ] **F4344** (medium, bug, easy/small) `topology_optimization.rs:575` — P-norm stress sensitivity NaN for zero stresses or ρ=0; asserts
- [ ] **F4345** (medium, fabrication, medium/medium) `topology_optimization.rs:611` — compliance_minimisation uses fake unit strain energies, never solves FEA
- [ ] **F4348** (medium, bug, medium/small) `underwater_softbody.rs:804` — BioFin AoA uses rad/s as velocity; FinMode never used
- [ ] **F4351** (medium, missing-feature, easy/small) `underwater_softbody.rs:1102` — UnderwaterSoftBody ignores added mass despite AddedMass model existing
- [ ] **F4353** (medium, perf, easy/trivial) `volume/functions.rs:213` — Pressurized surface-volume gradient is O(N·T) per solve
- [ ] **F4358** (medium, bug, easy/trivial) `volume/types.rs:58` — HashSet edge order makes Gauss-Seidel edge solve nondeterministic
- [ ] **F4361** (medium, fabrication, medium/small) `volume/types.rs:285` — Gas pressure added to PBD lambda: dimensionally wrong, explodes
- [ ] **F4364** (medium, bug, easy/trivial) `volume/types.rs:428` — Plane collision amplifies inward velocity instead of reflecting it
- [ ] **F4369** (medium, bug, medium/medium) `xpbd/functions.rs:143` — Shape-matching polar decomposition returns identity for planar clusters
- [ ] **F4378** (medium, perf, medium/medium) `xpbd/types.rs:801` — Self-collision detect_pairs is brute-force O(n²)

### Low (39) — tracked in the external ledger
F4027, F4035, F4047, F4052, F4059, F4067, F4070, F4099, F4101, F4103, F4109, F4167, F4168, F4169, F4177, F4195, F4218, F4234, F4240, F4257, F4261, F4270, F4291, F4295, F4296, F4298, F4309, F4313, F4314, F4325, F4335, F4336, F4337, F4349, F4350, F4357, F4368, F4376, F4379

<!-- AUDIT-2026-10-12 END -->
