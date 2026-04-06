# oxiphysics-vehicle

**Status:** Alpha | **Version:** 0.1.0 | **Tests:** 2,669

Full-featured vehicle dynamics simulation for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.

## Domains

- **Drivetrain:** Engine curves, gearbox, differential, drive layout (FWD/RWD/AWD), powertrain
- **Tires:** Pacejka (magic formula), Fiala, linear tire models; wear, thermal modelling
- **Suspension:** Linear, progressive, active suspension; suspension analysis and optimization
- **Steering:** Ackermann steering geometry
- **Braking & Stability:** ABS, traction control, electronic stability control
- **Aerodynamics:** Downforce, drag, wind loading
- **Electric Vehicles:** Electric motor, energy recovery, battery/fuel cell, charging station, HVAC
- **Motorsport:** Lap simulator, race-line optimization, race simulation, pit strategy
- **Autonomous Driving:** Path planning, sensors, driver model, autonomous vehicle control
- **Special Vehicles:** Motorcycle dynamics, aircraft dynamics
- **NVH / Comfort:** Noise-vibration-harshness, ride quality
- **Thermal:** Cooling system, tire thermal, thermal management
- **Telemetry:** Vehicle state logging and telemetry

## Key Exports

```rust
use oxiphysics_vehicle::{
    // Drivetrain
    Drivetrain, EngineCurve, Gearbox, Differential, DriveLayout,
    // Steering
    AckermannSteering,
    // Suspension
    LinearSuspension, ProgressiveSuspension,
    // Tires
    PacejkaTire, PacejkaCoeffs, FialaTire, LinearTire,
    // Vehicle
    RaycastVehicle, VehicleState,
    // Wheels
    Wheel, WheelConfig,
    // Stability
    AbsController, TractionControl, ElectronicStabilityControl,
};
```

## Modules (35+)

`aerodynamics`, `autonomous_driving`, `autonomous_vehicle`, `brake_system`, `chassis_dynamics`,
`charging_station`, `cooling`, `drivetrain`, `driver_model`, `electric_motor`, `electric_vehicle`,
`energy_recovery`, `fuel_cell`, `fuel_systems`, `gearbox`, `hvac`, `lap_simulator`, `lap_time`,
`motorcycle`, `noise_vibration`, `path`, `powertrain`, `race_line`, `race_simulation`,
`ride_quality`, `sensors`, `stability`, `steering`, `suspension`, `suspension_analysis`,
`suspension_optimization`, `telemetry`, `thermal_management`, `tire`, `tire_thermal`, `tire_wear`,
`track_dynamics`, `wheel`, `wind_loading`

## Statistics

| Metric | Count |
|--------|-------|
| Public items | 3,171 |
| Tests | 2,669 |
| Stubs | 0 |
| Modules | 35+ |

## License

MIT OR Apache-2.0 — Copyright © COOLJAPAN OU (Team Kitasan)
