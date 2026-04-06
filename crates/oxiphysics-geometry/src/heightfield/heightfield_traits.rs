//! # HeightField - Trait Implementations
//!
//! This module contains trait implementations for `HeightField`.
//!
//! ## Implemented Traits
//!
//! - `Shape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};

use super::functions::ray_triangle;
#[allow(unused_imports)]
use super::functions::*;
use super::types::HeightField;

impl Shape for HeightField {
    fn bounding_box(&self) -> Aabb {
        let (min_h, max_h) = self.height_bounds();
        Aabb::new(
            Vec3::new(0.0, min_h, 0.0),
            Vec3::new(
                self.scale_x * (self.cols - 1) as Real,
                max_h,
                self.scale_z * (self.rows - 1) as Real,
            ),
        )
    }
    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let mut best_dot = Real::NEG_INFINITY;
        let mut best_point = Vec3::zeros();
        for row in 0..self.rows {
            for col in 0..self.cols {
                let p = Vec3::new(
                    col as Real * self.scale_x,
                    self.height_at(row, col),
                    row as Real * self.scale_z,
                );
                let d = p.dot(direction);
                if d > best_dot {
                    best_dot = d;
                    best_point = p;
                }
            }
        }
        best_point
    }
    fn volume(&self) -> Real {
        if self.rows < 2 || self.cols < 2 {
            return 0.0;
        }
        let cell_area = self.scale_x * self.scale_z;
        let mut total = 0.0;
        for row in 0..(self.rows - 1) {
            for col in 0..(self.cols - 1) {
                let h00 = self.height_at(row, col);
                let h10 = self.height_at(row + 1, col);
                let h01 = self.height_at(row, col + 1);
                let h11 = self.height_at(row + 1, col + 1);
                let avg_h = (h00 + h10 + h01 + h11) / 4.0;
                total += avg_h * cell_area;
            }
        }
        total
    }
    fn center_of_mass(&self) -> Vec3 {
        Vec3::zeros()
    }
    fn inertia_tensor(&self, _mass: Real) -> Mat3 {
        Mat3::zeros()
    }
    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        for row in 0..(self.rows.saturating_sub(1)) {
            for col in 0..(self.cols.saturating_sub(1)) {
                let x0 = col as Real * self.scale_x;
                let x1 = (col + 1) as Real * self.scale_x;
                let z0 = row as Real * self.scale_z;
                let z1 = (row + 1) as Real * self.scale_z;
                let v00 = Vec3::new(x0, self.height_at(row, col), z0);
                let v10 = Vec3::new(x0, self.height_at(row + 1, col), z1);
                let v01 = Vec3::new(x1, self.height_at(row, col + 1), z0);
                let v11 = Vec3::new(x1, self.height_at(row + 1, col + 1), z1);
                for tri in [&[v00, v01, v11], &[v00, v11, v10]] {
                    if let Some(hit) = ray_triangle(
                        ray_origin,
                        ray_direction,
                        max_toi,
                        &tri[0],
                        &tri[1],
                        &tri[2],
                    ) && best.as_ref().is_none_or(|b| hit.toi < b.toi)
                    {
                        best = Some(hit);
                    }
                }
            }
        }
        best
    }
}
