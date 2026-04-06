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
pub trait Shape: std::fmt::Debug + Send + Sync {
    /// Compute the axis-aligned bounding box of this shape (in local space).
    fn bounding_box(&self) -> Aabb;

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
