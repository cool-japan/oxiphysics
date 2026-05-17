# oxiphysics-articulated

Featherstone articulated-body dynamics (RNEA + ABA) for OxiPhysics.

## Features

- **RNEA** — Recursive Newton-Euler Algorithm: inverse dynamics, computes joint torques from `(q, q̇, q̈)`
- **ABA** — Articulated Body Algorithm: forward dynamics, computes `q̈` from `(q, q̇, τ)`
- **CRBA** — Composite Rigid Body Algorithm: O(n²) joint-space mass matrix
- **OSC** — Operational Space Control: inertia `Λ = (J·M⁻¹·Jᵀ)⁻¹`
- **Centroidal dynamics** — Centroidal Momentum Matrix (CMM)
- **Spatial math** — 6D spatial vectors, spatial inertia, Plücker transforms
- **Joint types** — Revolute, Prismatic, Fixed, FreeFloating, Universal, Spherical, Helical
- **Joint limits** — Soft joint limits via penalty forces (`JointLimit`, `JointLimitSet`)

## Usage

Add to `Cargo.toml`:

```toml
[dependencies]
oxiphysics-articulated = "0.1.1"
```

Construct a model and run inverse or forward dynamics:

```rust
use oxiphysics_articulated::{
    body::RigidBody,
    joint::RevoluteJoint,
    model::ArticulatedModel,
    rnea::rnea,
    aba::aba,
};

// Build a simple two-link arm
let mut model = ArticulatedModel::new();
let link0 = model.add_body(RigidBody::new(/* inertia, parent */));
let link1 = model.add_body(RigidBody::new(/* inertia, parent */));
model.add_joint(RevoluteJoint::new(/* axis, X_T */), link0);
model.add_joint(RevoluteJoint::new(/* axis, X_T */), link1);

let q  = vec![0.0, 0.0];   // joint positions
let qd = vec![0.0, 0.0];   // joint velocities
let qdd = vec![0.0, 0.0];  // joint accelerations

// Inverse dynamics: given motion → torques
let tau = rnea(&model, &q, &qd, &qdd);

// Forward dynamics: given torques → accelerations
let tau_applied = vec![1.0, 0.5];
let qdd_out = aba(&model, &q, &qd, &tau_applied);
```

## Part of OxiPhysics

This crate is part of the [OxiPhysics](https://github.com/cool-japan/oxiphysics) ecosystem — a pure-Rust physics engine built on Featherstone's rigid-body dynamics algorithms.

## License

Apache-2.0
