//! Contact / constraint solver building blocks for the rigid-body world.
//!
//! Currently hosts the [`soft`] sub-module, which provides the
//! frequency/damping-ratio ([`SoftParams`]) parameterization used by the
//! TGS-Soft contact response path of [`crate::world::PhysicsWorld`].

pub mod soft;

pub use soft::SoftParams;
