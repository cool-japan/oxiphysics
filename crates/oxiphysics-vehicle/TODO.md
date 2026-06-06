# oxiphysics-vehicle TODO

Last updated: 2026-06-06 | Version: 0.1.2

## Phase 1: Foundation
- [x] Define core types and traits
- [x] Implement basic error handling
- [x] Add unit tests

## Phase 2: Core Implementation
- [x] Implement primary algorithms (Pacejka, Fiala, Ackermann, etc.)
- [x] Add integration tests (2,687 tests passing)
- [x] Performance benchmarks

## Phase 3: Polish
- [x] Documentation
- [x] Examples
- [x] Optimization

## Phase 4: Advanced Features
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

## Future / Post-0.1
- [x] Real-time hardware-in-the-loop (HiL) interfaces (`hil` — `HilInterface` trait, `SimHilBridge`, `HilSignalLogger`)
- [x] Co-simulation with FEM chassis deformation (`fem_chassis_cosim` — Craig-Bampton modal reduction, Newmark-β integration, coupling bridge)
- [x] GPU-parallel multi-vehicle batch simulation (`gpu_multi_vehicle` — SoA state, `MultiVehicleBatch`, `VehicleParams`, Rayon parallel step kernel)
