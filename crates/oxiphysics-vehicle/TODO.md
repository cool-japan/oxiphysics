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

<!-- AUDIT-2026-10-12 BEGIN -->
## Audit 2026-10-12 — verified findings (Phase 32 remediation)

Verified: 178 (critical 2, high 38, medium 107, low 31); 29 further raw findings pending second-round verification. IDs refer to the audit ledger; `[x]` = fixed.

### Critical (2)
- [x] **F4847** (critical, bug, easy/small) `powertrain_advanced.rs:971` — Hybrid battery current sign inverted: motoring charges, regen drains SOC
- [ ] **F4882** (critical, bug, easy/small) `suspension/types.rs:233` — QuarterCarModel::step: wrong unsprung reaction sign, bogus road input, wrong static deflection

### High (38)
- [x] **F4752** (high, bug, easy/small) `aerodynamics/types.rs:1040` — total_lift_coefficient divides by chord (units 1/m), wrong elliptic scale
- [x] **F4754** (high, bug, medium/small) `aircraft_dynamics.rs:950` — SixDofIntegrator ignores p,q,r rate damping and Ixz coupling
- [x] **F4757** (high, bug, easy/small) `autonomous_driving/types.rs:59` — KalmanTracker predict skips F·P·Fᵀ, so velocity is never estimated
- [x] **F4758** (high, fabrication, medium/medium) `autonomous_driving/types.rs:192` — DubinsPath::shortest discards LSL and returns a straight line
- [x] **F4771** (high, fabrication, hard/medium) `autonomous_vehicle/types.rs:322` — ParticleFilterSlam never builds a map; weight update is a no-op
- [x] **F4772** (high, fabrication, medium/medium) `autonomous_vehicle/types.rs:717` — UkfState is a fake UKF: hardcoded +1 m/s predict, no update
- [x] **F4779** (high, bug, easy/small) `brake_system.rs:843` — CornerBrake uses absolute pressure: phantom drag torque and heating
- [x] **F4780** (high, fabrication, medium/small) `brake_system.rs:950` — BrakeSystem ABS modulation is stateless and hold-gated; effectively no-op
- [ ] **F4793** (high, bug, easy/trivial) `electric_vehicle/functions.rs:1059` — step_battery discharge efficiency inverted; r=0 gives NaN current
- [ ] **F4802** (high, fabrication, medium/small) `electric_vehicle/types_ext.rs:623` — ThermalManagement::update drops coolant ODE; heat assumes fixed 1 A
- [x] **F4808** (high, bug, easy/trivial) `fuel_cell.rs:216` — Stack H2 consumption omits n_cells factor; tank drains n_cells× slow
- [ ] **F4809** (high, bug, easy/small) `fuel_cell.rs:397` — FCV range: fabricated tank mass, and the result is seconds, not km
- [x] **F4810** (high, bug, easy/small) `fuel_cell.rs:410` — FcvRange::estimate_range returns seconds/1000 and invents tank mass
- [x] **F4820** (high, fabrication, medium/small) `hvac.rs:376` — OccupantComfort::pmv fabricated: zero radiation, discarded tcl, clothing ignored
- [ ] **F4824** (high, bug, easy/small) `lap_simulator/types.rs:630` — straight_time cruises at average speed; simulate_lap assumes v_lim reached
- [x] **F4828** (high, bug, medium/small) `motorcycle_dynamics.rs:1100` — Motorcycle step double-counts drive/brake force; wheel reaction uses abs()
- [ ] **F4829** (high, bug, medium/medium) `motorcycle_dynamics.rs:1117` — Motorcycle step swaps a/b arms, double-counts drive/brake, abs() wheel torque
- [x] **F4840** (high, bug, easy/small) `path_planning.rs:408` — RRT plan connects new node to goal without collision check
- [x] **F4841** (high, bug, easy/trivial) `path_planning.rs:558` — PathSmoother::smooth_path has no iteration cap and can hang forever
- [x] **F4842** (high, bug, medium/small) `path_planning.rs:694` — DubinsPath::plan: wrong frame angle and wrong p² in all six words
- [x] **F4846** (high, stub, medium/medium) `powertrain_advanced.rs:926` — EmsStrategy::Ecms is a stub; lambda_ecms unused; motor rpm never updated
- [x] **F4853** (high, fabrication, medium/small) `race_line.rs:538` — compute_braking_points always yields zero braking distance
- [x] **F4862** (high, fabrication, easy/small) `race_strategy.rs:939` — strategy_total_time fabricates race length as last pit lap + 30
- [ ] **F4863** (high, fabrication, medium/small) `ride_quality.rs:91` — ISO 2631 Wk/Wb/Wd weightings are invented Gaussian/notch curves
- [ ] **F4867** (high, pure-rust, easy/small) `ride_quality.rs:296` — psd_road_profile is an O(n²) naive DFT; use oxifft::rfft
- [ ] **F4870** (high, bug, easy/small) `sensors/types_advanced.rs:217` — KalmanFusion2D::update_gps covariance update is not (I−KH)P
- [ ] **F4872** (high, bug, easy/small) `sensors/types_advanced.rs:616` — TelemetryBuffer::last_n returns oldest physical slots, not newest samples
- [ ] **F4877** (high, bug, easy/trivial) `steering/types.rs:921` — UndersteerAnalyzer understeer gradient divided by an extra wheelbase
- [ ] **F4884** (high, fabrication, medium/medium) `suspension/types.rs:1078` — DoubleWishboneKinematic camber is rigid translation; camber_gain is identically zero
- [ ] **F4899** (high, bug, medium/small) `tire.rs:174` — PacejkaTire/LinearTire combined forces exceed friction circle (≈1.1–1.4 μFz)
- [ ] **F4908** (high, fabrication, easy/small) `tire_thermal/types.rs:407` — RubberDegradation rate ≈1/s at any temperature; not Arrhenius
- [ ] **F4922** (high, bug, medium/small) `vehicle.rs:348` — Wheel spin clamped >=0 makes reverse gear non-functional
- [ ] **F4927** (high, bug, easy/small) `vehicle_dynamics.rs:137` — Critical/characteristic speed inverted and missing wheelbase
- [ ] **F4929** (high, bug, easy/small) `vehicle_dynamics.rs:313` — NonlinearBicycleModel tire force sign inverted; dynamics diverge
- [ ] **F4930** (high, bug, medium/small) `vehicle_dynamics.rs:439` — DoubleTrackModel::step double-counts TV torque, wrong ay, ignores torques
- [ ] **F4931** (high, fabrication, medium/small) `vehicle_dynamics.rs:537` — RolloverDynamics::j_turn_peak_ay ignores steer/stiffness; returns 0.1·v0²
- [ ] **F4932** (high, bug, easy/trivial) `vehicle_dynamics.rs:637` — steady_state_yaw_rate and UndersteerGradient::lateral_accel formulas wrong
- [ ] **F4937** (high, bug, easy/small) `vehicle_dynamics.rs:945` — WheelSpinControl LSD sends torque to faster wheel, allows negative torque

### Medium (107)
- [ ] **F4746** (medium, api, medium/small) `crates/oxiphysics-vehicle/Cargo.toml:17` — Unused deps oxiphysics-constraints/fem; FEM cosim not wired to FEM
- [ ] **F4748** (medium, panic, easy/medium) `aerodynamics/types.rs:36` — Clamp bounds from public fields panic when min>max or NaN
- [ ] **F4749** (medium, fabrication, easy/small) `aerodynamics/types.rs:235` — 'RBF interpolation' surrogate is Nadaraya-Watson kernel averaging
- [ ] **F4751** (medium, bug, easy/small) `aerodynamics/types.rs:796` — AeroBalance splits downforce by CoG instead of CoP
- [ ] **F4753** (medium, bug, easy/trivial) `aircraft_dynamics.rs:569` — Spiral-mode time constant uses non-physical (Lβ − Nβ) numerator
- [ ] **F4755** (medium, bug, medium/medium) `aircraft_dynamics.rs:1082` — Euler-rate guard flips sign when pitch exceeds 90 degrees
- [ ] **F4756** (medium, missing-feature, easy/small) `autonomous.rs:510` — SpeedPlanner ignores comfort_decel and waypoint speed limits
- [ ] **F4759** (medium, fabrication, medium/medium) `autonomous_driving/types.rs:259` — RRT rrt_star flag never read; NaN distance panics via expect
- [ ] **F4760** (medium, panic, medium/small) `autonomous_driving/types.rs:266` — RRT nearest search panics on NaN; rrt_star flag never read
- [ ] **F4763** (medium, fabrication, medium/medium) `autonomous_driving/types.rs:816` — Parking paths are fabricated geometry that teleports into the slot
- [ ] **F4764** (medium, fabrication, easy/trivial) `autonomous_driving/types.rs:906` — AV step feeds fabricated accelerations to the comfort monitor
- [ ] **F4766** (medium, bug, easy/small) `autonomous_driving/types.rs:1098` — RiskAssessor flags objects behind or beside the ego, triggering AEB
- [ ] **F4767** (medium, panic, easy/trivial) `autonomous_driving/types.rs:1399` — stanley_steer panics on an empty path
- [ ] **F4768** (medium, fabrication, hard/medium) `autonomous_driving/types.rs:1421` — MpcPathTracker ignores horizon and Q/R weights; it is pure pursuit
- [ ] **F4769** (medium, panic, easy/small) `autonomous_vehicle/functions.rs:30` — Zero noise std panics in random_range; uniform noise, not Gaussian
- [ ] **F4775** (medium, api, easy/small) `autonomous_vehicle/types.rs:1093` — Planners/filters use unseeded rand::rng(); RRT panics on empty bounds
- [ ] **F4777** (medium, fabrication, easy/trivial) `brake_system.rs:652` — Regen SOC update uses hardcoded 1 kWh battery
- [ ] **F4778** (medium, fabrication, easy/small) `brake_system.rs:779` — BrakeByWire discards max_pressure; latency and fallback are fake
- [ ] **F4781** (medium, stub, medium/medium) `brake_system.rs:976` — enable_regen/enable_bbw are silent no-ops; step never uses them
- [ ] **F4784** (medium, fabrication, medium/small) `chassis_resonance.rs:82` — natural_frequencies uses only the K diagonal; inertia/damping unused
- [ ] **F4785** (medium, bug, easy/trivial) `chassis_resonance.rs:262` — ISO 8608 PSD uses rad/m Φ(Ω0) values with n0 in cycles/m
- [ ] **F4786** (medium, bug, easy/trivial) `driver_model.rs:188` — trail_braking commands full throttle and 50% brake together
- [ ] **F4787** (medium, bug, easy/trivial) `drivetrain/types.rs:191` — PolynomialTorqueCurve::typical_na peaks at ~5400 N·m for a 2.0 L
- [ ] **F4789** (medium, fabrication, easy/small) `drivetrain/types.rs:1258` — Differential::split_torque ignores lock_factor; limited_slip sign bug
- [ ] **F4790** (medium, bug, easy/trivial) `drivetrain/types.rs:1267` — Differential::limited_slip bias sign wrong for negative or zero speed difference
- [ ] **F4791** (medium, bug, medium/medium) `drivetrain/types.rs:1889` — Wheels beyond index 3 never spin and act as locked brakes
- [ ] **F4792** (medium, panic, easy/trivial) `electric_motor.rs:29` — bilinear_lookup panics on negative peak torque or short maps
- [ ] **F4794** (medium, bug, easy/small) `electric_vehicle/types_core.rs:334` — RangeEstimator hardcodes dt=1 s and ignores DriveCycle.time
- [ ] **F4795** (medium, bug, easy/trivial) `electric_vehicle/types_core.rs:965` — CC/CV charger pack resistance is off by a factor of cells_parallel
- [ ] **F4796** (medium, panic, easy/small) `electric_vehicle/types_core.rs:1375` — EfficiencyMap::efficiency panics on empty, single-point or unsorted tables
- [ ] **F4797** (medium, panic, easy/small) `electric_vehicle/types_core.rs:1400` — EfficiencyMap::efficiency panics on short, descending or ragged maps
- [ ] **F4798** (medium, bug, easy/trivial) `electric_vehicle/types_ext.rs:72` — Flux-weakening id subtracts Wb from A (dimensionally wrong)
- [ ] **F4799** (medium, bug, easy/small) `electric_vehicle/types_ext.rs:402` — V2G charge setpoint has nonsense units; mean_price can panic
- [ ] **F4801** (medium, bug, easy/small) `electric_vehicle/types_ext.rs:541` — MotorEfficiencyMap speed axis normalised by max_power (W)
- [ ] **F4803** (medium, api, easy/medium) `error.rs:9` — Crate error enums unused; no public API returns them
- [ ] **F4804** (medium, bug, easy/small) `fem_chassis_cosim.rs:271` — Cosim: unvalidated attachment panics later; sub_dt silently ignored
- [ ] **F4805** (medium, bug, easy/small) `fem_chassis_cosim.rs:418` — sub_dt field is ignored; sub-stepping contradicts its documentation
- [ ] **F4806** (medium, perf, easy/small) `fem_chassis_cosim.rs:431` — Φη rebuilt twice per step; co-sim log grows without bound
- [ ] **F4812** (medium, fabrication, easy/small) `fuel_systems/types.rs:159` — KERS absorb_energy discards energy; ignores gear ratio and CVT loss
- [ ] **F4813** (medium, fabrication, easy/small) `fuel_systems/types.rs:594` — bsfc_g_per_kwh returns a constant regardless of engine state
- [ ] **F4814** (medium, bug, easy/trivial) `gearbox.rs:194` — Torque converter TR is 1.38 at lock-up; clamp can panic
- [ ] **F4816** (medium, fabrication, easy/small) `gpu_multi_vehicle.rs:418` — step_batch_rayon claims Rayon but calls serial step_batch
- [ ] **F4818** (medium, fabrication, medium/small) `hil.rs:409` — HiL tick records requested dt as timing; overruns never detected
- [ ] **F4819** (medium, fabrication, medium/small) `hil.rs:436` — HiL overrun detection is fake: records nominal dt, never errors
- [ ] **F4821** (medium, fabrication, medium/small) `lap_simulator/functions.rs:116` — simulate_lap reports fabricated tyre degradation and 1-second ERS deployment
- [ ] **F4822** (medium, panic, easy/small) `lap_simulator/types.rs:358` — PitStopStrategy constructors panic for small total_laps
- [ ] **F4823** (medium, bug, easy/trivial) `lap_simulator/types.rs:587` — OptimalLap::v_corner samples lateral grip at fixed 50 m/s
- [ ] **F4825** (medium, fabrication, medium/small) `lap_simulator/types.rs:1330` — RacingLine optimiser grows radii without bound; track width unused
- [ ] **F4827** (medium, api, medium/large) `lib.rs:38` — Root Drivetrain export unusable with RaycastVehicle; duplicate type names
- [ ] **F4830** (medium, stub, hard/medium) `motorcycle_dynamics.rs:1127` — Motorcycle steer, pitch and suspension states never evolve
- [ ] **F4831** (medium, bug, easy/trivial) `noise_vibration.rs:466` — sone_to_phon uses log2 with the log10 coefficient
- [ ] **F4833** (medium, bug, easy/trivial) `path/types.rs:109` — lookahead_point on open paths wraps to segments behind the vehicle
- [ ] **F4834** (medium, panic, easy/small) `path/types.rs:363` — Public constructors assert on caller input and panic
- [ ] **F4835** (medium, api, easy/small) `path/types.rs:529` — Path builders assert on too-few waypoints instead of PathError
- [ ] **F4836** (medium, panic, medium/small) `path/types.rs:738` — MinCurvatureRacingLine::optimize asserts and is O(n²) per iteration
- [ ] **F4837** (medium, perf, medium/small) `path/types.rs:747` — MinCurvatureRacingLine::optimize is O(iters·n²) with 4 clones per point
- [ ] **F4843** (medium, panic, easy/trivial) `path_planning.rs:967` — RoadGraph::add_edge accepts out-of-range target; Dijkstra panics
- [ ] **F4844** (medium, bug, easy/trivial) `powertrain_advanced.rs:141` — CombustionEngine::bmep uses 2π and per-cylinder displacement
- [ ] **F4845** (medium, fabrication, medium/small) `powertrain_advanced.rs:512` — Turbo compressor power fabricated; wastegate has no effect
- [ ] **F4848** (medium, bug, easy/trivial) `powertrain_advanced.rs:1159` — DCT in neutral transmits engine torque with ratio 1.0
- [ ] **F4849** (medium, panic, easy/small) `powertrain_dynamics.rs:57` — bilinear_interp indexes out of bounds for single-point axes
- [ ] **F4851** (medium, bug, medium/small) `race_line.rs:219` — Lap sim single Euler step per segment; braking points always zero
- [ ] **F4852** (medium, bug, easy/small) `race_line.rs:311` — geometric_racing_line apex is the arc centre; arc walked straight
- [ ] **F4855** (medium, fabrication, easy/small) `race_simulation.rs:210` — Point-mass LapSimulator has no braking pass; forward pass off by one
- [ ] **F4856** (medium, bug, medium/small) `race_simulation.rs:224` — Point-mass LapSimulator has no braking pass; sectors are an equal split
- [ ] **F4858** (medium, bug, easy/small) `race_strategy.rs:384` — undercut_possible ignores tyre delta; vsc_pit_bonus returns residual loss
- [ ] **F4859** (medium, bug, easy/small) `race_strategy.rs:431` — predict_lap_time ignores compound pace; softer tyres only wear faster
- [ ] **F4860** (medium, bug, easy/small) `race_strategy.rs:693` — Strategy sim: DRS every lap, compound pace ignored, undercut tautology
- [ ] **F4864** (medium, fabrication, medium/small) `ride_quality.rs:110` — ISO 2631 Wk/Wd/Wb weightings are invented curves
- [ ] **F4866** (medium, bug, medium/small) `ride_quality.rs:224` — IRI quarter-car starts at zero elevation and uses unstable Euler
- [ ] **F4868** (medium, panic, easy/trivial) `sensors/types_advanced.rs:23` — LidarScan zero angular resolution hangs for ~1.8e19 iterations
- [ ] **F4873** (medium, fabrication, medium/small) `sensors/types_advanced.rs:864` — Radar SNR is bare 1/R⁴ and CFAR threshold never rejects
- [ ] **F4881** (medium, fabrication, easy/trivial) `suspension/functions.rs:684` — sinusoidal_road_amplitude ignores position and wavelength
- [ ] **F4883** (medium, fabrication, hard/medium) `suspension/types.rs:459` — DoubleWishbone/McPherson camber, wheel centre and RC are magic formulas
- [ ] **F4885** (medium, fabrication, medium/small) `suspension_analysis.rs:216` — DoubleWishbone RC/IC formulas and kingpin offset are fabricated
- [ ] **F4886** (medium, bug, easy/small) `suspension_analysis.rs:340` — AntiRollBar::roll_gradient uses half-track as roll moment arm
- [ ] **F4888** (medium, fabrication, medium/small) `suspension_analysis.rs:640` — rms_body_acceleration formula is dimensionally fabricated
- [ ] **F4892** (medium, perf, easy/small) `telemetry/types.rs:22` — TelemetryRecorder panics when max_frames=0; O(n) remove(0) ring
- [ ] **F4896** (medium, bug, medium/small) `telemetry/types.rs:720` — LapComparison aligns frames by index, not track distance
- [ ] **F4897** (medium, fabrication, medium/small) `thermal_management.rs:146` — CoolingCircuit has no coolant advection; flow_rate never read
- [ ] **F4900** (medium, bug, medium/small) `tire.rs:254` — FialaTire saturates Fx and Fy independently, resultant up to sqrt2*muFz
- [ ] **F4901** (medium, panic, easy/trivial) `tire.rs:309` — f64::clamp panics on negative or NaN friction/normal force
- [ ] **F4904** (medium, panic, easy/trivial) `tire.rs:842` — Tire force clamp panics on negative or NaN normal force
- [ ] **F4907** (medium, bug, medium/small) `tire_thermal/types.rs:209` — TireThermalModel conductance k*A lacks thickness; rim path copies tread; road_temp fabricated
- [ ] **F4910** (medium, bug, easy/trivial) `tire_thermal/types.rs:945` — TireTemperatureMap::rotate drops fractional rotation; overflows on inf
- [ ] **F4912** (medium, bug, easy/small) `tire_wear.rs:494` — Gay-Lussac law applied to gauge pressure instead of absolute
- [ ] **F4913** (medium, fabrication, easy/small) `tire_wear.rs:681` — WearSimulation distance/remaining_km fabricated; compound durability unused
- [ ] **F4914** (medium, fabrication, medium/small) `track_dynamics.rs:111` — WheelRailContact 'Kalker' creep ignores G·a·b; μ fixed at 0.3
- [ ] **F4916** (medium, fabrication, hard/medium) `track_dynamics.rs:305` — critical_speed_hunting uses hard-coded conicity and an ad-hoc 'Klingel' formula
- [ ] **F4917** (medium, bug, easy/trivial) `track_dynamics.rs:429` — rail_bending_stress ignores documented span; returns N/m³, not Pa
- [ ] **F4918** (medium, bug, easy/trivial) `track_dynamics.rs:460` — hertz_rail_contact_width dimensionally invalid (returns m^1.5)
- [ ] **F4921** (medium, fabrication, hard/medium) `track_dynamics.rs:1109` — BogieDynamics::step: decoupled y/psi, no creep or conicity, yaw moment in N
- [ ] **F4923** (medium, bug, easy/trivial) `vehicle.rs:390` — lateral_g computed from lateral velocity, not forward speed
- [ ] **F4925** (medium, bug, easy/trivial) `vehicle.rs:717` — compute_weight_transfer applies only half the lateral load transfer
- [ ] **F4926** (medium, bug, easy/trivial) `vehicle.rs:818` — vehicle::AntiRollBar sign amplifies roll; contradicts suspension::AntiRollBar
- [ ] **F4928** (medium, bug, easy/trivial) `vehicle_dynamics.rs:178` — steady_state_yaw_rate multiplies understeer term by wheelbase
- [ ] **F4934** (medium, fabrication, medium/small) `vehicle_dynamics.rs:792` — YawMomentDiagram::generate hard-codes yaw rate r=0; vx unused
- [ ] **F4935** (medium, bug, easy/trivial) `vehicle_dynamics.rs:813` — neutral_steer_angle clamps every negative ay to -0.01
- [ ] **F4941** (medium, bug, easy/trivial) `wheel/types.rs:407` — ISO 8608 road class PSDs 16x too small for n0=0.1
- [ ] **F4943** (medium, bug, easy/small) `wheel/types.rs:518` — WheelState::update hard-codes radius/mass; force used as torque
- [ ] **F4944** (medium, bug, easy/small) `wheel/types.rs:527` — WheelState::update spins a resting wheel backward; hard-coded radius/mass
- [ ] **F4946** (medium, panic, easy/trivial) `wheel/types.rs:782` — ContactPatch::from_load panics on zero/NaN width; peak equals mean
- [ ] **F4947** (medium, fabrication, easy/small) `wheel/types.rs:952` — Wheel wobble 'SDOF oscillator' is first-order decay; ωn unused
- [ ] **F4949** (medium, bug, easy/trivial) `wind_loading.rs:69` — relative_wind beta is absolute bearing, heading never subtracted
- [ ] **F4950** (medium, bug, easy/small) `wind_loading.rs:200` — Aero moments lack reference length; side force always zero
- [ ] **F4951** (medium, panic, easy/trivial) `wind_loading.rs:216` — wind_yaw_angle hangs forever on infinite or huge input
- [ ] **F4952** (medium, bug, easy/trivial) `wind_loading.rs:230` — crosswind_stability_factor is dimensionally inconsistent (1/m)

### Low (31) — tracked in the external ledger
F4750, F4762, F4765, F4770, F4783, F4788, F4807, F4826, F4832, F4838, F4839, F4850, F4857, F4861, F4865, F4869, F4871, F4874, F4875, F4876, F4879, F4887, F4889, F4894, F4898, F4905, F4906, F4909, F4911, F4924, F4938

<!-- AUDIT-2026-10-12 END -->
