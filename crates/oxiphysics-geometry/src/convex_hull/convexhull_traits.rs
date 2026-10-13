//! # ConvexHull - Trait Implementations
//!
//! This module contains trait implementations for `ConvexHull`.
//!
//! ## Implemented Traits
//!
//! - `Shape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};

use super::solid::HullSolid;
use super::types::ConvexHull;

impl Shape for ConvexHull {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn bounding_box(&self) -> Aabb {
        if self.vertices.is_empty() {
            return Aabb::new(Vec3::zeros(), Vec3::zeros());
        }
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for v in &self.vertices[1..] {
            min = min.inf(v);
            max = max.sup(v);
        }
        Aabb::new(min, max)
    }
    fn support_point(&self, direction: &Vec3) -> Vec3 {
        self.vertices
            .iter()
            .max_by(|a, b| {
                a.dot(direction)
                    .partial_cmp(&b.dot(direction))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or_else(Vec3::zeros)
    }
    /// Volume of the polyhedron spanned by the vertices (zero if they span no
    /// volume).
    fn volume(&self) -> Real {
        self.with_solid(|solid| solid.map_or(0.0, HullSolid::volume))
    }
    /// Centre of mass of the uniform solid; the vertex mean if the vertices
    /// span no volume.
    fn center_of_mass(&self) -> Vec3 {
        self.with_solid(|solid| solid.map_or_else(|| self.centroid(), HullSolid::center_of_mass))
    }
    /// Inertia tensor of the uniform solid about its centre of mass,
    /// including the products of inertia; zero if the vertices span no
    /// volume.
    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        self.with_solid(|solid| solid.map_or_else(Mat3::zeros, |solid| solid.inertia_tensor(mass)))
    }
    /// Intersect the ray with the hull facets.
    ///
    /// `toi` is in units of `ray_direction`; the normal is the outward facet
    /// normal.  A ray starting inside reports the exit point.
    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        self.with_solid(|solid| solid?.ray_cast(ray_origin, ray_direction, max_toi))
    }
}
