// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Featherstone articulated-body dynamics for OxiPhysics.
//!
//! This crate implements the two foundational algorithms from
//! Featherstone 2008 "Rigid Body Dynamics Algorithms":
//!
//! - **RNEA** (Recursive Newton-Euler Algorithm) — inverse dynamics:
//!   given motion state `(q, q̇, q̈)`, compute joint torques `τ`.
//!
//! - **ABA** (Articulated Body Algorithm) — forward dynamics:
//!   given state `(q, q̇)` and applied torques `τ`, compute `q̈`.
//!
//! # Crate structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`spatial`] | 6D spatial vectors, spatial inertia, Plücker transforms |
//! | [`joint`] | Joint trait + Revolute, Prismatic, Fixed, FreeFloating, Universal, Spherical, Helical |
//! | [`joint_limits`] | Soft joint limits via penalty forces (`JointLimit`, `JointLimitSet`) |
//! | [`body`] | Rigid body with inertia + parent link |
//! | [`model`] | Kinematic tree (`ArticulatedModel`) |
//! | [`rnea`] | Inverse dynamics — computes torques from motion |
//! | [`aba`] | Forward dynamics — computes accelerations from torques |
//! | [`crba`] | Composite Rigid Body Algorithm — O(n²) mass matrix |
//! | [`centroidal`] | Centroidal Momentum Matrix (CMM) |
//! | [`osc`] | Operational-space inertia Λ = (J·M⁻¹·Jᵀ)⁻¹ |
//!
//! # Examples
//!
//! ## RNEA: gravity torque at horizontal position
//!
//! A unit-mass pendulum of length 1 m rotated to q = π/2 (horizontal)
//! experiences a restoring gravity torque of m·g·l = 9.81 N·m.
//!
//! ```
//! use oxiphysics_articulated::{
//!     rnea::rnea,
//!     model::ArticulatedModel,
//!     body::RigidBody,
//!     joint::RevoluteJoint,
//!     spatial::{SpatialInertia, SpatialTransform},
//! };
//!
//! let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
//! let inertia = SpatialInertia::from_com(1.0, [0.0, -1.0, 0.0], [[0.0; 3]; 3]);
//! let body = RigidBody::new("link", inertia, None, SpatialTransform::IDENTITY);
//! model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
//!
//! // At q = π/2 (horizontal), gravity torque = m·g·l = 9.81 N·m
//! let tau = rnea(&model, &[std::f64::consts::FRAC_PI_2], &[0.0], &[0.0]);
//! assert!((tau[0] - 9.81).abs() < 1e-6, "expected τ≈9.81, got {}", tau[0]);
//! ```
//!
//! ## ABA: zero acceleration at static equilibrium
//!
//! Feeding the holding torque (from RNEA) back into ABA with zero velocity
//! must recover zero joint acceleration — the system is in static balance.
//!
//! ```
//! use oxiphysics_articulated::{
//!     aba::aba,
//!     rnea::rnea,
//!     model::ArticulatedModel,
//!     body::RigidBody,
//!     joint::RevoluteJoint,
//!     spatial::{SpatialInertia, SpatialTransform},
//! };
//!
//! let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
//! let inertia = SpatialInertia::from_com(1.0, [0.0, -1.0, 0.0], [[0.0; 3]; 3]);
//! let body = RigidBody::new("link", inertia, None, SpatialTransform::IDENTITY);
//! model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
//!
//! // Compute holding torque at horizontal via RNEA, then feed into ABA
//! let q = [std::f64::consts::FRAC_PI_2];
//! let tau = rnea(&model, &q, &[0.0], &[0.0]);
//! let q_ddot = aba(&model, &q, &[0.0], &tau);
//! // ABA with holding torque → zero acceleration
//! assert!(q_ddot[0].abs() < 1e-6, "expected q̈≈0, got {}", q_ddot[0]);
//! ```
//!
//! ## CRBA: mass-matrix diagonal is positive for a 1-DOF pendulum
//!
//! For a point mass m at distance l, the inertia about the pivot is m·l².
//! With m = 1 kg and l = 1 m the 1×1 mass matrix element equals 1.0.
//!
//! ```
//! use oxiphysics_articulated::{
//!     crba::compute_mass_matrix_crba,
//!     model::ArticulatedModel,
//!     body::RigidBody,
//!     joint::RevoluteJoint,
//!     spatial::{SpatialInertia, SpatialTransform},
//! };
//!
//! let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
//! // Unit mass at unit length → M[0][0] = m·l² = 1.0
//! let inertia = SpatialInertia::from_com(1.0, [0.0, -1.0, 0.0], [[0.0; 3]; 3]);
//! let body = RigidBody::new("link", inertia, None, SpatialTransform::IDENTITY);
//! model.add_body(body, Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])));
//!
//! let m = compute_mass_matrix_crba(&mut model, &[0.0]);
//! assert_eq!(m.len(), 1, "1×1 mass matrix for 1-DoF pendulum");
//! assert!((m[0][0] - 1.0).abs() < 1e-10, "M[0][0] = m·l² = 1.0, got {}", m[0][0]);
//! ```

#![warn(missing_docs)]
#![deny(unsafe_code)]

pub mod aba;
pub mod body;
pub mod centroidal;
pub mod crba;
pub mod joint;
pub mod joint_limits;
pub mod model;
pub mod osc;
pub mod rnea;
pub mod spatial;

pub use joint::HelicalJoint;
pub use joint::SphericalJoint;
pub use joint::UniversalJoint;
pub use joint_limits::{JointLimit, JointLimitSet, aba_with_limits};
