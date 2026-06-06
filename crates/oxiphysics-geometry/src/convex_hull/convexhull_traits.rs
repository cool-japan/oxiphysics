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

use super::types::ConvexHull;

impl Shape for ConvexHull {
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
    fn volume(&self) -> Real {
        self.tetra_decomposition().0
    }
    fn center_of_mass(&self) -> Vec3 {
        self.tetra_decomposition().1
    }
    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let bb = self.bounding_box();
        let size = bb.max - bb.min;
        let k = mass / 12.0;
        Mat3::new(
            k * (size.y * size.y + size.z * size.z),
            0.0,
            0.0,
            0.0,
            k * (size.x * size.x + size.z * size.z),
            0.0,
            0.0,
            0.0,
            k * (size.x * size.x + size.y * size.y),
        )
    }
    fn ray_cast(
        &self,
        _ray_origin: &Vec3,
        _ray_direction: &Vec3,
        _max_toi: Real,
    ) -> Option<RayHit> {
        None
    }
}
