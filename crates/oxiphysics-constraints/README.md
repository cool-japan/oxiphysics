# oxiphysics-constraints

**Status:** Alpha | **Version:** 0.1.1 | **Tests:** 2,173

Comprehensive constraint solvers and joint systems for the [OxiPhysics](https://github.com/cool-japan/oxiphysics) engine.

## Domains

- **Solvers:** PGS (Projected Gauss-Seidel), TGS (Temporal Gauss-Seidel), island-based solving, warm-starting
- **Joints:** Ball, Fixed, Gear, Motor, Prismatic, Pulley, Rack-and-Pinion, Revolute, Spring, 6-DOF
- **Contact:** Contact constraint theory, CCD (continuous collision detection) constraints, Baumgarte bias
- **Friction:** Coulomb friction, friction constraint solvers
- **Control:** PID motor, servo, motor constraints; robot control, locomotion control
- **Advanced:** Position-Based Dynamics (PBD), variational constraints, holonomic constraints
- **Optimal Control:** Trajectory optimization, optimal control, multi-agent coordination
- **Multibody:** Cable constraints, haptic constraints, game-physics constraint sets, multibody dynamics

## Key Exports

```rust
use oxiphysics_constraints::{
    // Core trait
    Constraint,
    // Contact
    ContactConstraint, CcdConstraint,
    // Solvers
    PgsSolver, TgsSolver,
    // Joints
    BallJoint, FixedJoint, GearJoint, MotorJoint,
    PrismaticJoint, PulleyJoint, RackPinionJoint,
    RevoluteJoint, SpringJoint, SixDofConstraint,
    // Island management
    Island, IslandManager,
    // Motors
    MotorConstraint, MotorPid, MotorSolver,
    // Utilities
    baumgarte_bias,
    // Re-exports
    control_theory::*, friction::*,
};
```

## Modules (18+)

`ccd_constraints`, `contact_theory`, `control_theory`, `friction`, `island_solver`, `islands`,
`joints`, `mechanism`, `motor_constraints`, `optimal_control`, `optimization_constraints`,
`pbd`, `robot_control`, `tgs_solver`, `traits`, `trajectory_optimization`,
`variational_constraints`, `warm_start`

## Statistics

| Metric | Count |
|--------|-------|
| Public items | 3,118 |
| Tests | 2,173 |
| Stubs | 0 |
| Modules | 18+ |

## License

Apache-2.0 — Copyright © COOLJAPAN OU (Team Kitasan)
