# oxiphysics-vehicle TODO

**Version:** 0.1.3 | **Updated:** 2026-06-11 | **Status:** v0.1.x complete — production baseline; forward roadmap active
**Scale:** 53,985 SLoC | 2,687 tests

## Completed (v0.1.0 – v0.1.3)

### Phase 1: Foundation

- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

### Phase 2: Core Implementation

- [x] Implement primary algorithms (Pacejka, Fiala, Ackermann, etc.)
- [x] Add integration tests (2,687 tests passing)
- [x] Performance benchmarks

### Phase 3: Polish

- [x] Documentation
- [x] Examples
- [x] Optimization

### Phase 4: Advanced Features

- [x] Pacejka magic-formula tire model
- [x] Fiala & linear tire models
- [x] Tire wear and thermal modelling
- [x] Active suspension, suspension analysis & optimization
- [x] ABS, traction control, ESC
- [x] Aerodynamics (downforce, drag, wind loading)
- [x] Drivetrain (engine curves, gearbox, differential)
- [x] Electric vehicle (motor, energy recovery, fuel cell, charging)
- [x] Lap simulator, race-line, race simulation
- [x] Autonomous driving (path, sensors, driver model)
- [x] Motorcycle and aircraft dynamics
- [x] NVH / ride quality
- [x] Cooling & thermal management
- [x] Telemetry

### Future / Post-0.1 (all shipped within v0.1.x)

- [x] Real-time hardware-in-the-loop (HiL) interfaces (`hil` — `HilInterface` trait, `SimHilBridge`, `HilSignalLogger`)
- [x] Co-simulation with FEM chassis deformation (`fem_chassis_cosim` — Craig-Bampton modal reduction, Newmark-β integration, coupling bridge)
- [x] GPU-parallel multi-vehicle batch simulation (`gpu_multi_vehicle` — SoA state, `MultiVehicleBatch`, `VehicleParams`, Rayon parallel step kernel)

### Code-verified baseline details (2026-06-11)

- [x] Simplified Pacejka tier — 4-coefficient B/C/D/E `PacejkaCoeffs` (`tire.rs:38-46`); retained as the "arcade" tier beneath the v0.2.0 full Magic Formula
- [x] Tire load/camber/aligning-torque extensions — `LoadSensitivity`, `CamberThrust`, `SelfAligningTorque` (`tire.rs:403-504`)
- [x] Transient brush foundations — `TransientBrushModel` (`tire.rs:774`) and `ContactPatchPressure` (`tire.rs:907`)
- [x] Fiala as a brush model with saturation — `FialaTire` (`tire.rs:190`)
- [x] First-order tire relaxation length — `TireRelaxation` (`tire.rs:320-340`)
- [x] CVT and DCT transmissions — `powertrain_advanced.rs:1008`, `powertrain_advanced.rs:1110`
- [x] Crank-resolved combustion engine and turbo — Wiebe-based `CombustionEngine` (`powertrain_advanced.rs:65`) and `TurboCharger` with compressor map (`powertrain_advanced.rs:436`)
- [x] Analytic suspension kinematics oracles — `MacPhersonStrut`, `DoubleWishbone`, `MultiLink` (`suspension_analysis.rs:127-239`)
- [x] Lane-change path types — `autonomous_driving/types.rs` (path generation only; the validation harness is a v1.0 item)
- [x] HiL timing instrumentation — `HilTimingStats` (real-time budget baseline for the v1.0 performance gates)
- [x] SoA multi-vehicle batch — `MultiVehicleBatch` (`gpu_multi_vehicle.rs:275`) with Rayon step kernel

## v0.2.0 — Tire fidelity and transmission completion

The tire stack moves from the simplified 4-coefficient tier to measurement-grade models; the transmission family gains its missing automatic element; suspension goes multibody.

### Tire stack

- [ ] Full Magic Formula (PAC2002 / MF-Tyre 6.x coefficient set)
  - **Goal:** ~50-coefficient model reproduces published MF reference curves (Fy/Fx/Mz vs slip, 4 loads, camber sweep) to <1%; TIR-file import.
  - **Design:** Pacejka 2012 full pure-slip + combined cosine-weighting formulation with load/camber/pressure scaling; current `PacejkaCoeffs` is simplified 4-coeff B/C/D/E (`tire.rs:38-46`) — keep as the "arcade" tier; subsume existing `LoadSensitivity`, `CamberThrust`, `SelfAligningTorque` (`tire.rs:403-504`); TIR parser in oxiphysics-io.
  - **Files:** `tire.rs`, new full-MF coefficient module
  - **Cross-crate:** TIR-file parser lands in oxiphysics-io.

- [ ] Unified parabolic-pressure brush tire model
  - **Goal:** single brush formulation with parabolic pressure distribution, explicit adhesion/sliding boundary, and pneumatic trail output validated against the Fiala closed form to 1e-6.
  - **Design:** unify the existing pieces — `TransientBrushModel` (`tire.rs:774`), `ContactPatchPressure` (`tire.rs:907`), and `FialaTire` brush-with-saturation (`tire.rs:190`); see baseline details — this item is the remaining unification work, not a new model.
  - **Files:** `tire.rs`

- [ ] MF-Swift-style rigid-ring transient tire
  - **Goal:** valid to ~80 Hz and short-wavelength obstacles; cleat-strike force response matches rigid-ring reference within 10%.
  - **Design:** rigid belt ring on residual stiffness + enveloping via elliptical cams (Schmeitz); feeds the existing first-order relaxation model below 8 Hz (`TireRelaxation`, `tire.rs:320-340` — already shipped, see baseline details; relaxation length itself is not new work).
  - **Files:** `tire.rs`, new rigid-ring module

### Transmission and suspension

- [ ] Torque converter
  - **Goal:** K-factor/torque-ratio vs speed-ratio model with lockup clutch; stall-test and coast-down match characteristic curves to 2%.
  - **Design:** standard capacity-factor lookup + lockup state machine; CVT and DCT already exist (`powertrain_advanced.rs:1008`, `:1110` — see baseline details) — the converter is the missing automatic element; integrate with `Gearbox` shift logic.
  - **Files:** `powertrain_advanced.rs`, gearbox module

- [ ] Multibody suspension from real linkage geometry
  - **Goal:** double-wishbone/MacPherson built from hardpoints reproduce camber/toe vs wheel-travel curves of the analytic models to <0.05°, then add compliance the analytic models cannot.
  - **Design:** constraint-loop linkages via oxiphysics-articulated loop closure (v0.3.0 dependency) or oxiphysics-constraints joint loops; analytic `MacPhersonStrut`/`DoubleWishbone`/`MultiLink` already exist (`suspension_analysis.rs:127-239`) and serve as the validation oracle.
  - **Files:** `suspension_analysis.rs` (oracle), new hardpoint-linkage module
  - **Cross-crate:** depends on oxiphysics-articulated v0.3.0 "Kinematic loop closure" (or the oxiphysics-constraints joint-loop alternative).

## v0.3.0 — Off-road, rigs, and virtual test lab

Off-road soil mechanics, multi-unit rigs, and a virtual kinematics & compliance laboratory.

- [ ] K&C virtual rig
  - **Goal:** automated bounce/roll/lateral-force/aligning-torque sweeps output a K&C report (wheel-rate, roll-steer, compliance-steer) for any suspension model.
  - **Design:** quasi-static solve harness over the multibody suspension; CSV/JSON export via oxiphysics-io.
  - **Cross-crate:** consumes the v0.2.0 multibody suspension; report export via oxiphysics-io.

- [ ] Terramechanics (Bekker-Wong)
  - **Goal:** drawbar pull vs slip on sandy loam matches Wong reference data ±10%; rut depth from pressure-sinkage to ±15%.
  - **Design:** Bekker pressure-sinkage (k_c, k_φ, n) + Janosi-Hanamoto shear; per-wheel contact patch discretization; optional handoff to softbody MPM sand for high fidelity (tracked in Deferred below).
  - **Files:** new terramechanics module; per-wheel contact-patch interface shared with the Deferred MPM tier

- [ ] Trailers and articulated rigs
  - **Goal:** tractor + semi-trailer (fifth wheel) and car+caravan simulate jackknife and trailer-sway onset speeds matching linearized stability analysis ±5%.
  - **Design:** hitch as spherical/revolute constraint into the rigid/constraints crates; multi-unit `Vehicle` composition API.
  - **Cross-crate:** hitch joints provided by oxiphysics-rigid / oxiphysics-constraints.

- [ ] Mean-value engine model unification
  - **Goal:** manifold-filling mean-value mode (intake/exhaust plenum dynamics, turbo lag) switchable against the crank-resolved model; step-throttle boost response within 5% of the crank-resolved reference.
  - **Design:** extends the existing Wiebe-based `CombustionEngine` + `TurboCharger` compressor map (`powertrain_advanced.rs:65`, `:436` — see baseline details); the remaining work is the mean-value mode itself plus the mode switch.
  - **Files:** `powertrain_advanced.rs`

## v1.0 — Production & validation

Regulatory-style maneuver harnesses, golden-scenario validation, and hard real-time gates.

- [ ] ISO 3888-2 double-lane-change + NHTSA fishhook scenario harness
  - **Goal:** scripted maneuvers with pass/fail metrics (cone corridor adherence, two-wheel-lift flag at <2 in, ESC-on vs ESC-off differential documented); runs headless via local script.
  - **Design:** scenario DSL over driver_model + path modules; existing `stability/` ESC/ABS/TCS configs exercised; only lane-change *path types* exist today (`autonomous_driving/types.rs`) — no validation harness yet.
  - **Files:** `autonomous_driving/`, `stability/`, driver model modules

- [ ] Steady-state and transient handling validation
  - **Goal:** understeer gradient on constant-radius (ISO 4138) within ±5% of bicycle-model analytic; step-steer yaw-rate overshoot/settling vs linear 2-DoF reference ±10%; 100-0 km/h braking distance vs μ analytic ±3%.
  - **Design:** golden-scenario suite with stored tolerances.
  - **Tests:** ISO 4138 constant-radius, step-steer, and straight-line braking golden scenarios with stored baselines

- [ ] Real-time and batch performance gates
  - **Goal:** full-fidelity single vehicle ≥10x real-time at 1 kHz (HiL budget per existing `HilTimingStats`); 1000-vehicle SoA batch real-time at 100 Hz (existing `MultiVehicleBatch`, `gpu_multi_vehicle.rs:275`, plus WGSL path).
  - **Design:** criterion gates; document HiL jitter bounds; GPU path stays wgpu-default with cudarc feature-gated non-default; gates run as local scripts (no new CI workflow yamls).

- [ ] Parameter-file stability + fuzz
  - **Goal:** TIR/JSON vehicle parameter round-trip via oxicode bit-exact; proptest over parameter ranges (no NaN, all subsystems reject invalid configs through `error.rs`).
  - **Design:** schema-versioned parameter structs, semver-checked public API.

## Deferred / research track

- [~] High-fidelity MPM terramechanics tier — per-wheel rigid-wheel-on-deformable-sand via oxiphysics-softbody MPM (CPIC coupling) behind the v0.3.0 Bekker-Wong layer — Ready when: oxiphysics-softbody v0.2.0 CPIC rigid-body coupling + sparse blocked MPM grid have shipped and the Bekker-Wong per-wheel contact-patch interface (v0.3.0 above) is defined.

---

History: the pre-2026-06-11 TODO (last updated 2026-06-06) is preserved in full under "Completed"; the forward roadmap derives from the code-verified dynamics design brief of 2026-06-11.
