// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shape trait and raycast types.

use oxiphysics_core::math::{Mat3, Real, Vec3};
use oxiphysics_core::{Aabb, MassProperties};

/// Result of a ray intersection test.
#[derive(Debug, Clone)]
pub struct RayHit {
    /// Point of intersection in world space.
    pub point: Vec3,
    /// Surface normal at the intersection point.
    pub normal: Vec3,
    /// Parameter along the ray (distance = toi * ray_direction.norm()).
    pub toi: Real,
}

/// Trait for geometric shapes used in physics simulation.
///
/// The `'static` bound makes every shape downcastable through [`Shape::as_any`],
/// which the narrow-phase dispatcher relies on for sound, panic-free concrete
/// shape recovery (no raw-pointer transmutes). All concrete shapes are owned
/// value types, so the bound is satisfied automatically.
pub trait Shape: std::fmt::Debug + Send + Sync + 'static {
    /// Compute the axis-aligned bounding box of this shape (in local space).
    fn bounding_box(&self) -> Aabb;

    /// Upcast to [`std::any::Any`] for safe concrete-type recovery.
    ///
    /// Used by the narrow-phase dispatcher to downcast a `&dyn Shape` back to a
    /// concrete shape via [`std::any::Any::downcast_ref`] instead of an `unsafe`
    /// pointer cast. Every implementor returns `self`.
    fn as_any(&self) -> &dyn std::any::Any;

    /// Compute the support point in the given direction (for GJK).
    fn support_point(&self, direction: &Vec3) -> Vec3;

    /// Compute the volume of this shape.
    fn volume(&self) -> Real;

    /// Compute the center of mass in local space.
    fn center_of_mass(&self) -> Vec3;

    /// Compute the inertia tensor for the given mass.
    fn inertia_tensor(&self, mass: Real) -> Mat3;

    /// Compute full mass properties for a given density.
    fn mass_properties(&self, density: Real) -> MassProperties {
        let mass = density * self.volume();
        MassProperties::new(mass, self.center_of_mass(), self.inertia_tensor(mass))
    }

    /// Cast a ray against this shape (in local space).
    /// Returns the first intersection within `max_toi`.
    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit>;
}
