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
