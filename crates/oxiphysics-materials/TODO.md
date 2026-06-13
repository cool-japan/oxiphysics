# oxiphysics-materials TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 83,076 SLoC | 4,506 tests

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation
- [x] Implement primary algorithms (75+ modules, 6,241 public items)
- [x] Add integration tests (4,506 tests total, 0 stubs)
- [x] Elastic, plasticity, hyperelastic, viscoelastic, creep, fatigue, fracture, damage
- [x] Composite / fiber / nanocomposite modules
- [x] Geological, geomechanics, geomaterial_models, porous_media
- [x] Biological, biomaterials, polymer_physics, polymer_mechanics
- [x] Smart materials: shape_memory, shape_memory_alloy, magnetocaloric_materials, metamaterials
- [x] Energy: battery_materials, energy_materials, electrochemistry, thermoelectrics, hydrogen_storage
- [x] Electronic/optical: semiconductor, dielectric, optical, optical_materials, quantum_materials, superconductor
- [x] Thermal, radiation, radiation_shielding, nuclear_materials
- [x] Manufacturing: additive_manufacturing, aerospace, alloy, foam, nano/nanomaterials
- [x] Tribology, tribology_ext, corrosion, multiphysics
- [x] EOS, phase_transform, crystal_plasticity, acoustics
- [x] Contact presets: ContactMaterialPair, FrictionCombineRule, ModulusCombineRule, RestitutionCombineRule
- [x] Performance benchmarks (basic)

### Phase 3: Polish
- [x] Documentation (rustdoc on all public items)
- [x] Extended examples (additional usage examples) — see `materials_bench` module
- [x] Benchmark suite expansion (`materials_bench` — `MatBenchHarness`, elastic/neo-Hookean/viscoplastic/fatigue/thermal batch kernel timing)
- [x] SIMD acceleration paths (`simd_paths` — SoA layout, elastic/neo-Hookean/viscoplastic/fatigue/return-mapping batch kernels)

### Verified-existing infrastructure (code-verified 2026-06-11 — basis for the roadmap below; not re-opened)
- [x] `FailureCriteria` trait (`elastic.rs`) — the only material trait today; the unified material-point interface is the v0.2.0 delta
- [x] Analytic `tangent_modulus` for return-mapping (`plasticity.rs`) — the AD path is the v0.2.0 delta
- [x] `update_allen_cahn` point updates + `LatentHeatModel` (`phase_transform`) — the full field solver is the v0.2.0 delta
- [x] Analytical rule-of-mixtures / Halpin-Tsai homogenization (`composite_materials.rs`) — RVE/FFT homogenization is the v0.2.0/v0.3.0 delta
- [x] Crystal-plasticity point model (`crystal_plasticity.rs`) — the CPFFT spectral solver is the v0.3.0 delta
- [x] Rainflow counting (ASTM E1049), `critical_plane_max_shear`, SWT (`fatigue/`) — the multiaxial extension is the v1.0 delta
- [x] Bayesian inference lives in oxiphysics-core — the materials wiring is the v1.0 delta

## v0.2.0 — Homogenization & constitutive infrastructure

### Spectral homogenization
- [x] (planned 2026-06-12) FFT-based spectral homogenization (Moulinec-Suquet) (oxifft)
  - **Goal:** 2-phase effective stiffness within Hashin-Shtrikman bounds; converges at phase contrast up to 1000.
  - **Design:** periodic Lippmann-Schwinger fixed-point + polarization / Eyre-Milton accelerated scheme via oxifft. Moulinec-Suquet 1998.
  - **Delta:** MISSING.
  - **Validation:** effective bulk + shear moduli of a 2-phase RVE inside the HS bounds; fixed-point + accelerated schemes both converge at phase contrast 1000 (accelerated in far fewer iterations).
  - **Files (proposed):** `homogenization/fft/mod.rs`, `homogenization/fft/lippmann_schwinger.rs`, `homogenization/fft/eyre_milton.rs`
  - **Tests (proposed):** `unit::green_operator_symmetry`, `integration::fft_homog_within_hs_bounds`, `integration::accelerated_converges_high_contrast`

### Constitutive trait
- [ ] Unified ConstitutiveModel / UMAT-style trait
  - **Goal:** the 5 existing models (elastic, J2, neo-Hookean, viscoelastic, crystal-plasticity) implement one stress-update + consistent-tangent trait, drivable directly by FEM.
  - **Design:** material-point trait mirroring the Abaqus UMAT contract.
  - **Delta:** only the `FailureCriteria` trait exists in `elastic.rs`; there is no unified material-point interface yet.
  - **Validation:** all 5 models implement the trait and reproduce their existing standalone stress/tangent outputs to 1e-12; an FEM element drives them through the trait.
  - **Files (proposed):** `constitutive/mod.rs` (trait `ConstitutiveModel`), `constitutive/state.rs`, impls wired into existing model modules
  - **Tests (proposed):** `unit::trait_roundtrip_matches_standalone_all_models`, `integration::fem_drives_constitutive_trait`

### Differentiable constitutivity
- [ ] Differentiable return-mapping (scirs2-autograd — add dep)
  - **Goal:** AD consistent tangent matches analytic to 1e-8; enables gradient-based calibration.
  - **Design:** dual-number / autograd path through the J2 return-map.
  - **Delta:** `plasticity.rs` has analytic `tangent_modulus`; the AD path is new.
  - **Validation:** AD-computed consistent tangent matches the analytic `tangent_modulus` to 1e-8 across the yield surface; gradient w.r.t. (E, σ_y, H) matches finite-difference.
  - **Files (proposed):** `plasticity/autodiff_return_map.rs`
  - **Tests (proposed):** `unit::ad_tangent_matches_analytic_1e8`, `unit::ad_param_gradient_matches_fd`

### Phase-field field solver
- [ ] Phase-field solidification / grain-growth field solver
  - **Goal:** dendrite tip velocity matches the Kobayashi reference; grain-growth ⟨R⟩ ∝ t^0.5.
  - **Design:** 2D/3D Allen-Cahn + thermal coupling field solver. Kobayashi 1993.
  - **Delta:** `phase_transform` has `update_allen_cahn` + `LatentHeatModel` at point level; the full field solver is new.
  - **Validation:** anisotropic dendrite tip velocity within reference tolerance; isotropic grain-growth mean radius scales as t^0.5.
  - **Files (proposed):** `phase_transform/field_solver/mod.rs`, `phase_transform/field_solver/allen_cahn_grid.rs`, `phase_transform/field_solver/thermal_coupling.rs`
  - **Tests (proposed):** `integration::dendrite_tip_velocity_kobayashi`, `integration::grain_growth_r_sqrt_t`

## v0.3.0 — Microstructure & data-driven

### RVE homogenization toolkit
- [ ] RVE generation + Hill-Mandel homogenization toolkit
  - **Goal:** moduli converge with RVE size; KUBC/SUBC bracket the PBC result.
  - **Design:** periodic/random microstructure generator + KUBC/SUBC/PBC macro-homogeneity BCs. Hill 1963; Kanit 2003.
  - **Delta:** only analytical rule-of-mixtures / Halpin-Tsai exist in `composite_materials.rs`.
  - **Validation:** effective moduli converge as RVE size grows; KUBC (upper) and SUBC (lower) bracket the PBC estimate at each size.
  - **Files (proposed):** `homogenization/rve/generator.rs`, `homogenization/rve/hill_mandel.rs`, `homogenization/rve/boundary_conditions.rs`
  - **Tests (proposed):** `unit::kubc_subc_bracket_pbc`, `integration::rve_moduli_converge_with_size`

### Crystal-plasticity FFT
- [ ] Crystal-plasticity FFT spectral solver (oxifft)
  - **Goal:** polycrystal macroscopic stress-strain within 5% of reference CPFEM; (111) pole figure qualitatively correct.
  - **Design:** couple the existing CP point model to a CPFFT RVE solver. Lebensohn 2012.
  - **Delta:** `crystal_plasticity.rs` is a point model; the spectral solver is new.
  - **Validation:** polycrystal macroscopic stress-strain curve within 5% of a reference CPFEM solution; (111) pole figure qualitatively reproduces texture evolution.
  - **Files (proposed):** `crystal_plasticity/fft/mod.rs`, `crystal_plasticity/fft/spectral_solver.rs`
  - **Tests (proposed):** `integration::cpfft_stress_strain_within_5pct`, `integration::cpfft_texture_qualitative`

### Data-driven constitutivity
- [ ] Data-driven (model-free) constitutivity (scirs2-neural — add dep)
  - **Goal:** recover linear-elastic response from a stress-strain database with <1% error, no explicit law.
  - **Design:** Kirchdoerfer-Ortiz distance-minimizing solver + optional NN surrogate. Kirchdoerfer-Ortiz 2016.
  - **Delta:** MISSING.
  - **Validation:** recovered stress-strain response from a synthetic linear-elastic database within 1% of the generating modulus, with no explicit constitutive law supplied.
  - **Files (proposed):** `data_driven/mod.rs`, `data_driven/distance_minimizing.rs`, `data_driven/nn_surrogate.rs`
  - **Tests (proposed):** `unit::distance_solver_recovers_linear_elastic`, `integration::data_driven_error_below_1pct`

## v1.0 — Production & validation
- [ ] Bayesian parameter identification (scirs2-optimize — add dep)
  - **Goal:** recover synthetic J2 (E, σ_y, H) within 2% with credible intervals covering truth.
  - **Design:** wire the oxiphysics-core Bayesian infra to material curves.
  - **Delta:** oxiphysics-core has Bayesian inference; the materials wiring is new.
  - **Validation:** posterior means for (E, σ_y, H) within 2% of the synthetic truth; 95% credible intervals cover the true values.
  - **Files (proposed):** `calibration/bayesian_id.rs`
  - **Tests (proposed):** `integration::j2_bayesian_recovers_params_within_2pct`
- [ ] Eshelby inclusion analytic validation
  - **Goal:** interior Eshelby tensor matches the closed form to 1e-6 (spherical / ellipsoidal).
  - **Design:** ellipsoidal-inclusion stress-field benchmark.
  - **Validation:** interior Eshelby tensor S_ijkl matches the closed-form expression to 1e-6 for spherical and several ellipsoidal aspect ratios.
  - **Files (proposed):** `../../oxiphysics/tests/validation_materials_eshelby.rs`
- [ ] Hashin-Shtrikman bounds + multiaxial fatigue validation
  - **Goal:** all homogenized moduli inside the HS bounds; critical-plane / SWT life within factor-2 of test data.
  - **Design:** HS bounds check + multiaxial fatigue extension.
  - **Delta:** rainflow (ASTM E1049), `critical_plane_max_shear`, and SWT exist in `fatigue/` — extension adds Fatemi-Socie / Findley + frequency-domain.
  - **Validation:** every homogenized modulus from the FFT/RVE toolkit lands inside the HS bounds; Fatemi-Socie / Findley / SWT predicted life within a factor of 2 of multiaxial test data.
  - **Files (proposed):** `fatigue/fatemi_socie.rs`, `fatigue/findley.rs`, `fatigue/frequency_domain.rs`, `../../oxiphysics/tests/validation_materials_hs_bounds.rs`
  - **Tests (proposed):** `integration::homogenized_moduli_within_hs_bounds`, `integration::multiaxial_life_within_factor_2`

## Cross-cutting (all milestones)
- [ ] Validation harness as local cargo targets — Eshelby, HS-bounds, and fatigue cases run via `cargo test -p oxiphysics --test validation_materials_*` plus a `scripts/` runner; no new CI workflow YAML (pure-Rust, local-only per policy).
- [ ] FFT paths (Moulinec-Suquet homogenization, CPFFT) go through oxifft only; no rustfft (ecosystem policy).
- [ ] HS-bounds assertion reused as a guardrail across the FFT homogenization, RVE toolkit, and CPFFT solvers so any effective-modulus regression is caught.
- [ ] No-unwrap + no-warnings maintained across new homogenization/constitutive/data-driven modules; deterministic reductions where homogenized moduli feed regression assertions.

## Implementation order & dependencies
1. **FFT homogenization (v0.2.0, oxifft)** — keystone; its periodic Green-operator solve is reused as the PBC backend for the RVE toolkit and as the kinematic engine for CPFFT. Build first.
2. **Unified ConstitutiveModel trait (v0.2.0)** — gates clean FEM material-point coupling; implemented by all 5 existing models.
3. **Phase-field field solver (v0.2.0)** — reuses the existing `update_allen_cahn` point kernels + `LatentHeatModel`; no new dep.
4. **Differentiable return-mapping (v0.2.0)** — implements the trait from (2); gated on scirs2-autograd (see Deferred).
5. **RVE + Hill-Mandel toolkit (v0.3.0)** — uses the FFT periodic solve from (1) as its PBC backend.
6. **Crystal-plasticity FFT (v0.3.0, oxifft)** — needs the existing CP point model plus the FFT Green operator from (1).
7. **Data-driven constitutivity (v0.3.0)** — the distance-minimizing core lands first; the optional NN surrogate is gated on scirs2-neural (see Deferred).
8. **Eshelby validation (v1.0)** — pure analytic correctness anchor for inclusion mechanics; no dep.
9. **Bayesian parameter ID (v1.0)** — wires the oxiphysics-core Bayesian infra; gated on scirs2-optimize (see Deferred).
10. **HS-bounds + multiaxial fatigue (v1.0)** — bounds reuse the homogenized moduli from (1)/(5); extends the existing rainflow/critical-plane/SWT infra.

## References
- FFT homogenization — Moulinec & Suquet, *Comput. Methods Appl. Mech. Eng.* 157 (1998) 69–94.
- Accelerated FFT scheme — Eyre & Milton, *Eur. Phys. J. AP* 6 (1999) 41–47.
- UMAT contract — Abaqus user-subroutine interface (stress update + consistent Jacobian).
- RVE size effects — Kanit et al., *Int. J. Solids Struct.* 40 (2003) 3647–3679.
- Hill-Mandel condition — Hill, *J. Mech. Phys. Solids* 11 (1963) 357–372.
- Crystal-plasticity FFT — Lebensohn et al., *Int. J. Plast.* 32–33 (2012) 59–69.
- Data-driven mechanics — Kirchdoerfer & Ortiz, *Comput. Methods Appl. Mech. Eng.* 304 (2016) 81–101.
- Phase-field solidification — Kobayashi, *Physica D* 63 (1993) 410–423.
- Eshelby inclusion — Eshelby, *Proc. R. Soc. A* 241 (1957) 376–396.
- Hashin-Shtrikman bounds — Hashin & Shtrikman, *J. Mech. Phys. Solids* 11 (1963) 127–140.
- Critical-plane fatigue — Fatemi & Socie, *Fatigue Fract. Eng. Mater. Struct.* 11 (1988) 149–165.
- Findley criterion — Findley, *J. Eng. Ind.* 81 (1959) 301–306.
- SWT parameter — Smith, Watson & Topper, *J. Mater.* 5 (1970) 767–778.
- Rainflow counting — ASTM E1049-85 (standard practice for cycle counting).
- Halpin-Tsai (existing baseline) — Halpin & Kardos, *Polym. Eng. Sci.* 16 (1976) 344–352.
- Mean-field homogenization (context) — Mori & Tanaka, *Acta Metall.* 21 (1973) 571–574.
- Periodic homogenization theory — Suquet, in *Homogenization Techniques for Composite Media* (1987).
- FFT polarization scheme — Michel, Moulinec & Suquet, *Int. J. Numer. Methods Eng.* 52 (2001) 139–160.

## Risks & open questions
- FFT homogenization convergence degrades at very high phase contrast — the Eyre-Milton / polarization accelerator is essential, not optional.
- Voxelized RVE geometry introduces staircase error at curved interfaces — quantify it against the analytic Eshelby field.
- AD through the return-map must handle the active/inactive yield-branch switch cleanly — verify consistent-tangent continuity across yielding.
- The data-driven solver needs sufficient database coverage — define a coverage metric and an explicit failure mode for extrapolation.
- CPFFT texture accuracy is sensitive to the spectral discretization near grain boundaries — validate against reference pole figures.
- HS bounds only bracket isotropic moduli — anisotropic effective tensors need a separate admissibility check.
- The scirs2-* deps are not yet wired — AD, data-driven, and Bayesian items must remain feature-gated until they land (tracked in Deferred).

## Deferred / research track
- [~] scirs2-neural / scirs2-autograd / scirs2-optimize workspace integration (gates differentiable return-mapping, data-driven constitutivity, and Bayesian parameter ID) — Ready when: added as workspace deps (only oxifft, oxiarc-zstd, scirs2-integrate are wired today).
