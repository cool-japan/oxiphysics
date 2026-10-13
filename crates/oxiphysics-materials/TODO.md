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
  - [x] Fix basic Moulinec-Suquet divergence → HS-bounds NaN (reference-modulus + divergence guard) (planned 2026-06-15)
  - **Goal:** 2-phase effective stiffness within Hashin-Shtrikman bounds; converges at phase contrast up to 1000.
  - **Design:** periodic Lippmann-Schwinger fixed-point + polarization / Eyre-Milton accelerated scheme via oxifft. Moulinec-Suquet 1998.
  - **Delta:** MISSING.
  - **Validation:** effective bulk + shear moduli of a 2-phase RVE inside the HS bounds; fixed-point + accelerated schemes both converge at phase contrast 1000 (accelerated in far fewer iterations).
  - **Files (proposed):** `homogenization/fft/mod.rs`, `homogenization/fft/lippmann_schwinger.rs`, `homogenization/fft/eyre_milton.rs`
  - **Tests (proposed):** `unit::green_operator_symmetry`, `integration::fft_homog_within_hs_bounds`, `integration::accelerated_converges_high_contrast`

### Constitutive trait
- [x] (done 2026-06-14) Unified ConstitutiveModel / UMAT-style trait
  - **Goal:** the 5 existing models (elastic, J2, neo-Hookean, viscoelastic, crystal-plasticity) implement one stress-update + consistent-tangent trait, drivable directly by FEM.
  - **Design:** material-point trait mirroring the Abaqus UMAT contract.
  - **Delta:** only the `FailureCriteria` trait exists in `elastic.rs`; there is no unified material-point interface yet.
  - **Validation:** all 5 models implement the trait and reproduce their existing standalone stress/tangent outputs to 1e-12; an FEM element drives them through the trait.
  - **Files (proposed):** `constitutive/mod.rs` (trait `ConstitutiveModel`), `constitutive/state.rs`, impls wired into existing model modules
  - **Tests (proposed):** `unit::trait_roundtrip_matches_standalone_all_models`, `integration::fem_drives_constitutive_trait`
  - **Status (2026-06-14):** DONE. `ConstitutiveModel` trait + `ConstitutiveResponse<S>` live in `constitutive/mod.rs`; state types `J2State` / `ViscoelasticState` in `constitutive/state.rs`. Implemented for: `LinearElastic`, `IsotropicElastic`, `OrthotropicElastic`, `TransverselyIsotropicElastic` (stateless), `J2ReturnMapping` (isotropic hardening, 7 SDVs, delegates to `return_map` + `J2ConsistentTangent::compute_consistent_tangent`), and `GeneralizedMaxwell` (new concrete `voigt_stress_update` exact-integration step + per-branch history state). 11 integration tests in `tests/constitutive_trait.rs` (parity gate at 1e-12 relative trait-vs-concrete for all wired models, virgin step, `n_state_vars`, generic driver) + 5 unit tests; clippy clean; all touched files < 2000 lines.
  - [ ] Follow-on: wire `crystal_plasticity` into `ConstitutiveModel`. The current `crystal_plasticity.rs` is a kinematics/texture point model (Schmid tensor, Taylor factor, slip-system creep) — it exposes NO Voigt-6 stress-update closure and NO consistent 6×6 tangent, so it cannot meet the 1e-12 parity gate today. Requires either (a) a closed `stress_update(strain)->(stress,tangent)` slip-system integrator with an analytic/algorithmic tangent, or (b) a numerical (finite-difference) tangent over a real stress-update. Tracked as a scope-split from the 2026-06-14 trait work.
  - [ ] Follow-on: wire `NeoHookean` (finite-strain) into `ConstitutiveModel` once a small-strain/Voigt adapter or a deformation-gradient-based trait variant is defined (the current trait takes a Voigt-6 small-strain tensor; the hyperelastic models are large-strain PK1/F-based).

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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 232 (critical 1, high 82, medium 119, low 30); 65 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (1)
- [ ] **F2963** (critical, bug, easy/small) `hyperelastic/functions.rs:86` — principal_stretches uses wrong acos sign, returns reflected eigenvalues

### High (82)
- [x] **F2800** (high, bug, medium/small) `acoustics/functions.rs:48` — Multilayer TMM: c=z/1.2, missing imaginary unit, unchecked slice
- [ ] **F2801** (high, bug, easy/small) `acoustics/functions.rs:87` — ISO 9613-1 absorption: fabricated humidity, f instead of f²
- [x] **F2802** (high, bug, easy/small) `acoustics/functions.rs:93` — ISO 9613-1 relaxation terms use f instead of f², wrong humidity
- [x] **F2808** (high, fabrication, medium/medium) `additive_manufacturing/types.rs:637` — ThermalHistory::from_process melt pool is always (r, 0.6r, 0.4r)
- [ ] **F2809** (high, bug, medium/small) `additive_manufacturing/types.rs:806` — Rosenthal melt_pool_width evaluated at xi=0, ~3x too narrow
- [ ] **F2814** (high, fabrication, medium/medium) `aerospace/types.rs:353` — CLT ply failure stress fabricated as Q̄·N/t (no A⁻¹ solve)
- [x] **F2816** (high, bug, medium/small) `aerospace/types.rs:1051` — CompositeLaminateClt::tsai_wu_index uses Q̄·N/t as ply stress
- [x] **F2819** (high, bug, easy/trivial) `aerospace_materials.rs:241` — Q̄66 mixes two formula variants, dropping 2·Q66·c²s²
- [ ] **F2824** (high, fabrication, medium/medium) `anisotropic/types.rs:229` — LaRC fibre kinking discards sigma22, FI=0.18 at sigma11=-Xc
- [x] **F2828** (high, bug, easy/trivial) `anisotropic/types.rs:963` — Hashin matrix compression: unsquared term and sign removed by abs()
- [ ] **F2829** (high, bug, medium/small) `anisotropic/types.rs:1061` — Puck IFF compression adds p_c|σ22| (wrong sign), no Mode C
- [ ] **F2833** (high, bug, easy/trivial) `biological_materials.rs:401` — Donnan swelling pressure ~2000x too small (units, factor 2)
- [x] **F2838** (high, bug, easy/trivial) `biomechanical_materials/functions.rs:162` — solve_cubic_symmetric eigenvalue amplitude is √3 too large
- [ ] **F2843** (high, bug, medium/small) `biomechanics/types.rs:374` — SoftTissueMaterial Cauchy stress: F=I+E, nonzero stress at rest
- [x] **F2846** (high, bug, easy/small) `biomedical_materials.rs:452` — Flory-Rehner elastic term has wrong sign; swelling solve never brackets
- [x] **F2858** (high, bug, easy/trivial) `composite.rs:507` — Mori-Tanaka shear denominator omits G_m (uses ζ_m, not G_m+ζ_m)
- [x] **F2859** (high, bug, medium/small) `composite.rs:753` — puck_failure transverse compression effort can never reach 1
- [x] **F2861** (high, fabrication, medium/medium) `composite.rs:1082` — interlaminar_shear_estimate returns mean σx proxy, not τxz
- [ ] **F2862** (high, bug, easy/small) `composite_failure.rs:246` — Tsai-Hill interaction term uses σ1σ2/(F1F2) instead of σ1σ2/F1²
- [ ] **F2864** (high, bug, medium/small) `composite_failure.rs:499` — Puck IFF Mode C absent: pure transverse compression never fails
- [ ] **F2866** (high, bug, medium/medium) `composite_failure.rs:735` — progressive_failure_sweep: strain N/(A·t), uncoupled, λ capped at 2
- [ ] **F2870** (high, bug, easy/trivial) `composites_advanced.rs:472` — Tsai-Wu F12 divides by sqrt(F11*F22) instead of multiplying
- [ ] **F2889** (high, bug, easy/small) `crystal_plasticity.rs:147` — FCC and BCC slip tables pair non-orthogonal planes and directions
- [ ] **F2899** (high, bug, easy/trivial) `elastic.rs:299` — Transversely isotropic compliance matrix is non-symmetric (wrong reciprocity)
- [ ] **F2900** (high, bug, easy/trivial) `elastic.rs:573` — inv_transpose3 has a wrong cofactor (2,1), corrupting NeoHookean PK1 stress
- [ ] **F2904** (high, bug, easy/trivial) `electrochemistry/functions.rs:72` — Corrosion penetration rate is 1000x too low (A/m² to µA/cm² conversion)
- [ ] **F2907** (high, bug, easy/small) `electrochemistry/types.rs:368` — Randles z_imag has spurious −X_C term; z_real drops Warburg reactance
- [ ] **F2910** (high, bug, medium/medium) `electromagnetic.rs:831` — Jiles-Atherton descending branch drives M away from anhysteretic
- [ ] **F2918** (high, bug, medium/small) `eos/miegruneiseneos_traits.rs:26` — Mie-Grüneisen sound speed adds a pressure to c² (dimension error)
- [ ] **F2920** (high, bug, easy/small) `eos/stiffenedgaseos_traits.rs:16` — Stiffened-gas pressure hardcodes e and uses −p∞ instead of −γp∞
- [ ] **F2923** (high, bug, easy/small) `eos/types.rs:569` — Mie-Grüneisen Hugoniot uses μ in η-form formula: wrong pressures
- [ ] **F2924** (high, bug, medium/small) `eos/types.rs:578` — Mie-Grüneisen Hugoniot mixes μ with an η-form denominator
- [ ] **F2925** (high, bug, medium/small) `eos/vanderwaalseos_traits.rs:27` — VdW density_from_pressure always returns ~10000 kg/m³; T hardcoded to 300 K
- [ ] **F2926** (high, bug, easy/trivial) `eos/vineteos_traits.rs:22` — VinetEos::sound_speed has reversed sign and always returns 0
- [ ] **F2929** (high, bug, medium/small) `fatigue/functions.rs:38` — rain_flow_count skips the start-point rule and drops residual ranges
- [ ] **F2935** (high, bug, easy/small) `foam_materials.rs:1006` — ImpactAttenuation HIC uses m/s² (≈301× too high); peak_g ignores velocity
- [ ] **F2936** (high, bug, easy/small) `fracture/functions.rs:179` — Weight-function SIF dominated by clamped endpoint singularity
- [ ] **F2937** (high, bug, easy/trivial) `fracture/functions.rs:412` — SENB SIF uses W^0.5 instead of W^1.5 (K low by factor W)
- [ ] **F2939** (high, fabrication, medium/small) `fracture/functions.rs:792` — Crack-from-hole SIF polynomial ignores Kt, ~3.4x non-conservative
- [ ] **F2940** (high, bug, easy/small) `fracture/functions.rs:804` — Surface-crack shape factor uses 4.593 instead of 1.464
- [ ] **F2945** (high, bug, easy/trivial) `geomaterial_models.rs:623` — time_for_consolidation uses ln instead of log10
- [ ] **F2946** (high, fabrication, easy/trivial) `geomaterial_models.rs:796` — CRR clean-sand curve is a fabricated fit, not Youd–Idriss 2001
- [ ] **F2949** (high, bug, easy/small) `geomaterials.rs:766` — SPT fines correction negative and discontinuous; stress unit doc wrong
- [ ] **F2951** (high, bug, easy/trivial) `geomechanics.rs:135` — Drucker-Prager q-p yield uses √J2-form constants
- [ ] **F2952** (high, bug, easy/trivial) `geomechanics.rs:150` — Triaxial volumetric strain has sign error on σ3 term
- [ ] **F2954** (high, bug, easy/trivial) `geomechanics.rs:327` — Hazen permeability off by factor 10^6
- [ ] **F2966** (high, bug, easy/trivial) `hyperelastic/functions.rs:424` — Ogden principal stresses lack deviatoric projection: σ≠0 at identity
- [ ] **F2967** (high, bug, medium/small) `hyperelastic/types.rs:173` — Eigen fallback zeroes isochoric energy for uniaxial F (Gent/Yeoh/AB)
- [ ] **F2969** (high, fabrication, medium/medium) `hyperelastic/types.rs:630` — Plane-stress return mapping claims iteration but leaves σzz≠0
- [ ] **F2978** (high, bug, easy/small) `lib.rs:591` — Hashin-Shtrikman lower bulk bound always equals the upper bound
- [ ] **F2980** (high, fabrication, easy/small) `magnetocaloric_materials.rs:184` — numerical_maxwell_entropy returns ΔM, not ΔS; μ0 missing
- [ ] **F2983** (high, fabrication, hard/large) `metamaterials/functions.rs:56` — "SIMP topology optimisation" is a neighbour-count automaton ignoring target
- [ ] **F2984** (high, bug, easy/trivial) `metamaterials/functions.rs:113` — Masters–Evans Poisson ratio has wrong sign
- [ ] **F2989** (high, bug, easy/small) `metamaterials/types.rs:710` — Re-entrant honeycomb Poisson ratio sign inverted at three sites
- [ ] **F2990** (high, fabrication, medium/small) `metamaterials/types.rs:940` — Acoustic cloak shell matches neither ideal nor reduced Cummer–Schurig parameters
- [ ] **F2992** (high, bug, easy/trivial) `multiphysics.rs:582` — SMA cooling martensite fraction is non-monotone and discontinuous
- [ ] **F2994** (high, bug, easy/small) `nano/types.rs:84` — Laminate D11 omits ply z-offsets (parallel-axis term)
- [ ] **F2996** (high, bug, easy/trivial) `nano/types.rs:126` — VirialStress kinetic term multiplies by atom mass; RVE atom count linear
- [ ] **F2997** (high, fabrication, medium/small) `nano/types.rs:168` — hardy_stress_xx sums force×velocity with magic 1e-9
- [ ] **F3000** (high, bug, easy/trivial) `nano/types.rs:544` — Taylor hardening off by 1e6 (GPa×1e9 labelled MPa)
- [ ] **F3003** (high, fabrication, hard/medium) `nanocomposites.rs:157` — "Mori-Tanaka" nanocomposite moduli use arbitrary fudge factors
- [ ] **F3005** (high, fabrication, easy/small) `nuclear_materials.rs:114` — nrt_dpa is dimensionally meaningless
- [ ] **F3008** (high, fabrication, medium/medium) `nuclear_materials_advanced.rs:569` — RG 1.99 chemistry factor and USE drop are invented fits
- [ ] **F3011** (high, fabrication, easy/small) `optical.rs:379` — CIE 1931 colour-matching fits mis-transcribed; x̄ uses ȳ lobe
- [ ] **F3014** (high, bug, medium/small) `optical_materials.rs:836` — Bragg mirror TMM drops imaginary unit; reflectance wrong off quarter-wave
- [ ] **F3016** (high, bug, easy/trivial) `phase_transform/types.rs:172` — BainiteTransformation::time_to_fraction always returns NaN
- [ ] **F3019** (high, fabrication, hard/medium) `phase_transform/types.rs:1135` — Brinson SMA: ξ_s never evolves, fudge 0.1 term, compounding heating
- [ ] **F3021** (high, bug, easy/trivial) `phase_transform/types.rs:1375` — CalphadModel::gibbs_energy adds T·S_mix instead of subtracting
- [ ] **F3026** (high, bug, medium/small) `phase_transform/types.rs:1895` — Phase-field solidification has no melting point; liquid always solidifies
- [ ] **F3027** (high, bug, easy/trivial) `plasticity.rs:360` — J2 return_map stores tensorial shear plastic strain in engineering-shear state
- [ ] **F3028** (high, bug, medium/small) `plasticity.rs:839` — J2 consistent tangent dimensionally wrong; missing 3/2 and cross terms
- [ ] **F3029** (high, bug, medium/small) `plasticity.rs:908` — Drucker-Prager apex check inverted; cone multiplier denominator wrong
- [ ] **F3039** (high, bug, easy/small) `radiation/functions.rs:423` — Perpendicular-rectangle view factor formula wrong: 0.093 vs exact 0.200
- [ ] **F3043** (high, bug, easy/small) `radiation/types.rs:245` — MC RNG maps to [0,0.5): view factor and transport biased
- [ ] **F3046** (high, fabrication, medium/small) `radiation/types.rs:797` — Slab MC transport ignores direction and weight; scattering only kills photons
- [ ] **F3048** (high, bug, easy/trivial) `radiation_shielding.rs:251` — Saturation activity multiplied by decay constant (units Bq/s)
- [ ] **F3049** (high, bug, easy/small) `semiconductor.rs:204` — fermi_dirac_half up to 161% wrong, discontinuous at eta=5
- [ ] **F3054** (high, bug, easy/small) `shape_memory_alloy.rs:219` — Tanaka stress coefficients wrong sign; a_M doubled; hysteresis window erased
- [ ] **F3057** (high, stub, medium/medium) `shape_memory_alloy.rs:472` — Brinson detwinning never occurs; shape-memory recovery ratio always 0
- [ ] **F3072** (high, bug, medium/small) `superconductor.rs:171` — GL solver mixes xi in metres with unit grid; profile never forms
- [ ] **F3083** (high, bug, easy/trivial) `tribology/types.rs:443` — Gear Hertz contact stress misses zone factor 2/sinφ (~2.4× low)
- [ ] **F3096** (high, bug, medium/medium) `viscoelastic/types.rs:1139` — Generalized Maxwell update relaxes toward E_a·ε: stiffens instead of relaxing

### Medium (119)
- [ ] **F2803** (medium, bug, easy/trivial) `acoustics/functions.rs:265` — Matching-layer reflection sums imaginary terms as reals
- [ ] **F2804** (medium, fabrication, easy/small) `acoustics/functions.rs:443` — weighted_sound_reduction_index is a mean offset, not ISO 717-1
- [ ] **F2807** (medium, bug, medium/small) `acoustics/types.rs:520` — PhononicCrystal1D Bragg frequency uses thickness-averaged speed
- [ ] **F2810** (medium, fabrication, medium/small) `additive_manufacturing/types.rs:937` — DistortionCompensation::iterate just rescales its own output
- [ ] **F2811** (medium, fabrication, medium/medium) `additive_manufacturing_materials.rs:271` — eagar_tsai is not Eagar-Tsai: depth fabricated, width off by sqrt(e)
- [ ] **F2813** (medium, bug, easy/trivial) `aerospace/types.rs:344` — CfrpLaminate::effective_ex divides A12² by A11, not A22
- [ ] **F2815** (medium, bug, easy/small) `aerospace/types.rs:821` — Paris-law life returns 0 or finite when below threshold
- [ ] **F2817** (medium, bug, easy/trivial) `aerospace/types.rs:1412` — CfrpPly Q̄16 first term uses c²s² instead of c³s
- [ ] **F2821** (medium, fabrication, medium/small) `aerospace_materials.rs:839` — required_tps_thickness_mm ignores temperatures; constant 1.5√(αt)
- [ ] **F2822** (medium, bug, easy/small) `anisotropic/functions.rs:49` — pole_figure folds lower hemisphere into last ring, unequal-area bins, panics on 0 bins
- [ ] **F2823** (medium, panic, easy/trivial) `anisotropic/functions.rs:61` — pole_figure panics for zero bins; lower hemisphere piles up
- [ ] **F2825** (medium, fabrication, medium/small) `anisotropic/types.rs:241` — matrix_cracking_index ignores tau23 and compressive sigma22, and returns sqrt despite doc
- [ ] **F2830** (medium, fabrication, medium/small) `anisotropic/types.rs:1071` — Puck fibre failure: invented 1% σ2 factor; compression uses Xt
- [ ] **F2832** (medium, fabrication, medium/small) `battery_materials.rs:656` — predict_cycle_life inconsistent with current_capacity model
- [ ] **F2834** (medium, bug, easy/small) `biomaterials.rs:425` — Flory-Rehner elastic potential uses an invented product form
- [ ] **F2836** (medium, bug, easy/trivial) `biomaterials.rs:509` — CellMechanics effective modulus is half the Hertz reduced modulus
- [ ] **F2837** (medium, bug, easy/trivial) `biomaterials.rs:902` — DentalEnamel Vickers hardness 10⁴× too large; flaw size extra 1/π
- [ ] **F2842** (medium, bug, medium/small) `biomechanical_materials/types.rs:1173` — FungHyperelastic S is 2x dW/dE; scalar tangent fabricated
- [ ] **F2845** (medium, fabrication, easy/trivial) `biomechanics/types.rs:1044` — stroke_work ignores MAP, returns t_max*SV/2
- [ ] **F2847** (medium, fabrication, medium/small) `biomedical_materials.rs:613` — ArterialWall residual_strain_inner: radius cancels, wrong sign
- [ ] **F2848** (medium, bug, easy/small) `biomedical_materials.rs:700` — DegradablePolymer::run hangs for dt=0; overshoots total time
- [ ] **F2849** (medium, fabrication, easy/small) `ceramic_materials.rs:285` — sintering_densification has no rate prefactor; returns green density
- [ ] **F2854** (medium, stub, easy/small) `combination.rs:201` — ContactMaterialPair damping is a hardcoded 0.0 placeholder
- [ ] **F2855** (medium, bug, easy/small) `combination.rs:240` — Thermal contact conductance adds W/(m·K) to W/(m²·K)
- [ ] **F2863** (medium, bug, easy/trivial) `composite_failure.rs:485` — PuckCriterion Mode A omits (1−p⁺Yt/S21) factor; fails below Yt
- [ ] **F2865** (medium, bug, easy/small) `composite_failure.rs:700` — Ply degradation applied twice (Q and thickness) and compounds
- [ ] **F2868** (medium, bug, medium/small) `composites_advanced.rs:396` — from_n_symmetric drops A16/A26; no general ABD solve
- [ ] **F2869** (medium, bug, easy/small) `composites_advanced.rs:405` — from_n_symmetric ignores A16/A26 coupling; docstring index order wrong
- [ ] **F2871** (medium, bug, medium/small) `composites_advanced.rs:634` — PuckCriterion::iff_mode_b formula is dimensionally wrong
- [ ] **F2872** (medium, bug, medium/small) `composites_advanced.rs:640` — PuckCriterion::iff_mode_b is dimensionally wrong and lacks mode C
- [ ] **F2873** (medium, panic, easy/trivial) `composites_advanced.rs:672` — ProgressiveDamage::apply_load_step panics on stress/ply length mismatch
- [ ] **F2875** (medium, bug, easy/small) `construction/types.rs:1843` — Prestressed Mn uses yb as dp and invented width/rho_p
- [ ] **F2877** (medium, bug, easy/small) `corrosion.rs:564` — SCC time_to_fracture counts sub-threshold growth as zero time
- [ ] **F2880** (medium, bug, easy/trivial) `creep/types.rs:523` — NabarroHerringCreep::compliance missing 1/(k_B·T) factor
- [ ] **F2882** (medium, bug, easy/small) `creep/types.rs:581` — TertiaryCreepModel::rupture_time omits 1/(φ+1), overestimating life 5×
- [ ] **F2886** (medium, bug, easy/trivial) `creep/types.rs:1117` — Norton/power-law creep gives NaN for compressive stress
- [ ] **F2887** (medium, fabrication, easy/trivial) `creep/types.rs:1134` — NortonCreep::time_to_rupture returns 1/strain-rate
- [ ] **F2888** (medium, bug, easy/small) `creep/types.rs:1338` — Perzyna step loses accumulated strain; negative strain softens
- [ ] **F2890** (medium, bug, easy/small) `crystal_plasticity.rs:231` — Bunge orientation matrix has φ1 and φ2 swapped
- [ ] **F2891** (medium, fabrication, hard/medium) `crystal_plasticity.rs:280` — taylor_factor sums resolved shears and is not a Taylor factor
- [ ] **F2892** (medium, fabrication, hard/large) `crystal_plasticity.rs:293` — taylor_factor is a sum of |RSS|, not a Taylor factor
- [ ] **F2894** (medium, perf, medium/medium) `damage/types.rs:218` — Non-local averaging is all-pairs O(n²) with no cutoff
- [ ] **F2897** (medium, fabrication, medium/small) `damage/types.rs:991` — Puck fibre-failure criterion ignores moduli, uses invented 0.01 term
- [ ] **F2898** (medium, bug, medium/small) `dielectric_materials.rs:92` — Wave impedance/reflection ignore conductivity and phase
- [ ] **F2903** (medium, perf, medium/small) `elastic.rs:1062` — Orthotropic/TI stress_update re-inverts the 6x6 compliance on every call
- [ ] **F2911** (medium, bug, easy/trivial) `electromagnetic.rs:924` — Shielding reflection loss uses absolute σ in the relative-σ formula
- [ ] **F2912** (medium, fabrication, medium/small) `energy_materials.rs:448` — Solar-cell current_density ignores stored r_series
- [ ] **F2913** (medium, fabrication, medium/medium) `energy_materials.rs:502` — Shockley-Queisser limit is an invented parabola
- [ ] **F2914** (medium, fabrication, medium/small) `energy_materials.rs:906` — pka_lindhard is not the Lindhard/Robinson partition
- [ ] **F2916** (medium, panic, easy/small) `eos/functions.rs:71` — Public EOS fits panic on length mismatch; invented default fit on empty data
- [ ] **F2922** (medium, bug, medium/small) `eos/types.rs:507` — VanDerWaals::compressibility picks arbitrary root below Tc
- [ ] **F2930** (medium, bug, medium/small) `fatigue/functions.rs:105` — Rainflow lacks starting-point rule and skips residual reversals
- [ ] **F2931** (medium, bug, medium/small) `fatigue/functions.rs:331` — critical_plane_max_shear scans only x-y normals, ignoring σz, τyz, τxz
- [ ] **F2932** (medium, bug, medium/small) `fatigue/functions.rs:547` — Dirlik damage truncated at 6 sigma, underestimates damage
- [ ] **F2938** (medium, bug, easy/trivial) `fracture/functions.rs:493` — Max-hoop-stress kink angle wrong root for negative K_II
- [ ] **F2942** (medium, bug, easy/small) `fracture/types.rs:358` — RCurve::instability_load returns None for every sub-critical slope
- [ ] **F2943** (medium, bug, easy/small) `geological.rs:357` — Bolton dilatancy ignores stress argument; e_max/e_min hardcoded
- [ ] **F2944** (medium, bug, easy/trivial) `geological.rs:656` — Hertz-Mindlin stiffness uses wrong effective moduli
- [ ] **F2947** (medium, stub, easy/small) `geomaterial_models.rs:1027` — Frozen-soil/permafrost models hardcode porosity 0.4 and bulk modulus
- [ ] **F2948** (medium, bug, easy/trivial) `geomaterials.rs:574` — average_doc_approx late-time formula wrong and discontinuous
- [ ] **F2950** (medium, bug, easy/trivial) `geomechanics.rs:113` — Hoek-Brown mixes intact m_i with rock-mass s
- [ ] **F2953** (medium, bug, easy/trivial) `geomechanics.rs:203` — Biot modulus ignores grain compressibility; Skempton B exceeds 1
- [ ] **F2957** (medium, bug, medium/small) `homogenization/fft/eyre_milton.rs:340` — Per-voxel 6x6 solve each iteration; incomplete isotropy test
- [ ] **F2958** (medium, perf, medium/small) `homogenization/fft/lippmann_schwinger.rs:180` — FFT plan rebuilt and buffers reallocated on every call
- [ ] **F2964** (medium, bug, easy/trivial) `hyperelastic/functions.rs:253` — tangent_shear_modulus_fd returns shear stress, not modulus
- [ ] **F2965** (medium, fabrication, hard/medium) `hyperelastic/functions.rs:308` — 'Full' Neo-Hookean tangent ignores b, depends only on J
- [ ] **F2970** (medium, api, medium/medium) `hyperelastic/types.rs:721` — Two MieGruneisenEos/JwlEos types with swapped constructor argument order
- [ ] **F2971** (medium, bug, easy/trivial) `hyperelastic/types.rs:816` — Fung uniaxial stress is 2/3 of the incompressible value
- [ ] **F2972** (medium, bug, medium/small) `hyperelastic/types.rs:941` — Blatz-Ko energy has residual pressure at F=I
- [ ] **F2975** (medium, missing-feature, easy/medium) `lib.rs:117` — Serde derives only on Material; material data cannot be persisted
- [ ] **F2979** (medium, panic, easy/trivial) `lib.rs:810` — Material ranking helpers panic on NaN property values
- [ ] **F2985** (medium, fabrication, easy/small) `metamaterials/functions.rs:134` — Honeycomb relative density and modulus are not Gibson–Ashby
- [ ] **F2986** (medium, bug, easy/trivial) `metamaterials/functions.rs:183` — Maxwell conductivity: v1 is actually inclusion fraction
- [ ] **F2987** (medium, fabrication, easy/small) `metamaterials/functions.rs:197` — Mechanical cloak ignores radii; hard-coded steel-to-rubber ramp
- [ ] **F2988** (medium, bug, medium/small) `metamaterials/types.rs:111` — Mechanical SSH dispersion, gap and edge mode have wrong units
- [ ] **F2993** (medium, bug, easy/small) `multiphysics.rs:844` — Newton solver silently zeroes singular directions; no failure flag
- [ ] **F2995** (medium, bug, easy/trivial) `nano/types.rs:98` — Tsai-Wu safety factor returns 1/(failure index), not quadratic root
- [ ] **F2998** (medium, bug, easy/small) `nano/types.rs:294` — Mori-Tanaka bulk modulus applies concentration factor twice
- [ ] **F2999** (medium, bug, easy/trivial) `nano/types.rs:383` — Armchair and chiral CNT diameters use a_cc instead of √3·a_cc
- [ ] **F3001** (medium, bug, easy/trivial) `nano/types.rs:982` — Hall-Petch mixes k in MPa·µm^0.5 with grain size in nm
- [ ] **F3002** (medium, bug, easy/trivial) `nano_materials.rs:120` — Hamaker sphere-sphere force uses sphere-plate formula, ignores a2
- [ ] **F3004** (medium, panic, easy/small) `nanomaterials.rs:988` — sample_cnt_bundle panics for n_max<5, NaN for 0 samples, unseedable
- [ ] **F3006** (medium, bug, medium/small) `nuclear_materials_advanced.rs:268` — Booth FGR uses pellet radius and wrong Xe activation energy
- [ ] **F3007** (medium, fabrication, hard/medium) `nuclear_materials_advanced.rs:453` — RIS 'Cr depletion' model predicts 2x enrichment
- [ ] **F3009** (medium, panic, easy/medium) `optical.rs:312` — Public constructors panic via assert_eq! on mismatched lengths
- [ ] **F3010** (medium, bug, easy/trivial) `optical.rs:334` — ColorFromSpectrum::to_xyz claims trapezoidal rule but ignores spacing
- [ ] **F3012** (medium, bug, easy/small) `optical.rs:521` — total_emissivity is unweighted sample mean, not Planck-weighted integral
- [ ] **F3015** (medium, bug, easy/trivial) `phase_transform/functions.rs:61` — landau_free_energy quartic sign makes potential unbounded
- [ ] **F3018** (medium, bug, medium/small) `phase_transform/types.rs:1122` — Brinson temperature update non-idempotent and direction-unaware
- [ ] **F3022** (medium, bug, easy/small) `phase_transform/types.rs:1433` — Martensite reversion compounds per step; result depends on dt
- [ ] **F3024** (medium, bug, easy/trivial) `phase_transform/types.rs:1616` — FirstOrderTransition latent heat is twice the correct value
- [ ] **F3025** (medium, bug, easy/trivial) `phase_transform/types.rs:1656` — Landau heat-capacity jump omits factor Tc
- [ ] **F3031** (medium, bug, easy/small) `polymer_mechanics.rs:130` — FJC force uses contour length as Kuhn length (N× too small)
- [ ] **F3032** (medium, bug, easy/small) `polymers_materials.rs:592` — Rouse/reptation use 0.1 g/mol monomer; reptation lacks factor 3
- [ ] **F3034** (medium, bug, easy/trivial) `porous_media.rs:74` — Forchheimer solves q² instead of q|q|; wrong for positive gradient
- [ ] **F3036** (medium, duplicate, medium/large) `presets/types.rs:9` — 66 duplicated public type names incl. two ExtendedMaterial types
- [ ] **F3044** (medium, fabrication, medium/medium) `radiation/types.rs:329` — RadiationNetwork documents full radiosity solve; only non-physical simple flux exists
- [ ] **F3047** (medium, panic, easy/small) `radiation_shielding.rs:136` — Domain asserts in public physics functions abort the process
- [ ] **F3053** (medium, fabrication, easy/small) `shape_memory.rs:205` — forward_transformation_stress ignores temperature despite documented Clausius-Clapeyron
- [ ] **F3055** (medium, bug, easy/trivial) `shape_memory_alloy.rs:253` — Tanaka/Liang-Rogers branch windows ignore the model's own C_M/C_A
- [ ] **F3056** (medium, bug, easy/small) `shape_memory_alloy.rs:259` — TanakaModel phase_fraction erases hysteresis; stress-coupling sign inverted
- [ ] **F3058** (medium, bug, easy/small) `shape_memory_alloy.rs:571` — Brinson update_state zeroes the thermoelastic stress term
- [ ] **F3059** (medium, perf, easy/small) `simd_paths.rs:210` — Advertised SIMD batch kernels index unchecked-length pub Vecs
- [ ] **F3061** (medium, bug, easy/small) `simd_paths.rs:337` — neo_hookean_stress_batch: J<0 yields NaN; skipped points keep stale stress
- [ ] **F3064** (medium, bug, medium/small) `smart_materials/types/types_2.rs:82` — SmartStructure modal step uses unstable forward Euler
- [ ] **F3066** (medium, bug, easy/small) `smart_materials/types/types_2.rs:715` — Flory-Rehner elastic potential inverted: stiffer gels swell more
- [ ] **F3067** (medium, bug, easy/small) `smart_materials/types/types_impl.rs:275` — HydrogelSwelling mixing potential swaps polymer and solvent fractions
- [ ] **F3068** (medium, fabrication, easy/trivial) `smart_materials/types/types_impl.rs:573` — SmartComposite::actuation_curvature discards its temperature argument
- [ ] **F3069** (medium, fabrication, medium/small) `smart_materials/types/types_impl.rs:731` — HydrogelActuator 'Flory-Rehner' swelling is an invented exponential
- [ ] **F3071** (medium, bug, medium/small) `smart_materials/types/types_impl.rs:914` — Preisach ferroelectric ignores P_r: remanent equals saturation
- [ ] **F3073** (medium, bug, easy/small) `superconductor.rs:339` — Inter-vortex force off by factor lambda/r; crude K1 approximation
- [ ] **F3075** (medium, bug, easy/small) `thermal/types.rs:205` — HAZ peak temperature and half-width are dimensionally inconsistent
- [ ] **F3079** (medium, panic, easy/trivial) `thermoelectrics.rs:293` — MaterialZT::interpolate(NaN) underflows index and panics
- [ ] **F3082** (medium, fabrication, medium/medium) `tribology/types.rs:339` — Elliptical Hertz ellipticity always 1; reduces to circular contact
- [ ] **F3085** (medium, fabrication, easy/small) `tribology/types.rs:948` — GW normal_load uses constant deformation σ, ignoring separation
- [ ] **F3087** (medium, bug, medium/small) `tribology/types.rs:1524` — FrictionOscillator never re-sticks; stick-slip phenomenon lost
- [ ] **F3089** (medium, bug, easy/small) `tribology_ext.rs:235` — EHL film thickness drops pressure-viscosity parameter G (~90× low)
- [ ] **F3090** (medium, bug, easy/small) `tribology_ext.rs:274` — Flash temperature formulas are dimensionally wrong (not Kelvin)
- [ ] **F3093** (medium, bug, easy/small) `viscoelastic/types.rs:799` — Burgers storage/loss moduli combine real moduli, arbitrary 1 Pa cutoff

### Low (30) — tracked in the external ledger
F2806, F2850, F2852, F2860, F2867, F2879, F2884, F2885, F2917, F2919, F2921, F2927, F2928, F2941, F2960, F2961, F2962, F2968, F2973, F2974, F2976, F2977, F3023, F3037, F3038, F3042, F3051, F3063, F3078, F3091

<!-- AUDIT-2026-10-12 END -->
