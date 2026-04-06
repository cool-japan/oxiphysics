//! # Compound - Trait Implementations
//!
//! This module contains trait implementations for `Compound`.
//!
//! ## Implemented Traits
//!
//! - `Shape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::shape::{RayHit, Shape};
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Mat3, Real, Vec3};

use super::types::Compound;

impl Shape for Compound {
    fn bounding_box(&self) -> Aabb {
        if self.children.is_empty() {
            return Aabb::new(Vec3::zeros(), Vec3::zeros());
        }
        let mut result: Option<Aabb> = None;
        for (transform, shape) in &self.children {
            let local_aabb = shape.bounding_box();
            let corners = [
                Vec3::new(local_aabb.min.x, local_aabb.min.y, local_aabb.min.z),
                Vec3::new(local_aabb.max.x, local_aabb.min.y, local_aabb.min.z),
                Vec3::new(local_aabb.min.x, local_aabb.max.y, local_aabb.min.z),
                Vec3::new(local_aabb.max.x, local_aabb.max.y, local_aabb.min.z),
                Vec3::new(local_aabb.min.x, local_aabb.min.y, local_aabb.max.z),
                Vec3::new(local_aabb.max.x, local_aabb.min.y, local_aabb.max.z),
                Vec3::new(local_aabb.min.x, local_aabb.max.y, local_aabb.max.z),
                Vec3::new(local_aabb.max.x, local_aabb.max.y, local_aabb.max.z),
            ];
            let mut min = transform.transform_point(&corners[0]);
            let mut max = min;
            for corner in &corners[1..] {
                let p = transform.transform_point(corner);
                min = min.inf(&p);
                max = max.sup(&p);
            }
            let child_aabb = Aabb::new(min, max);
            result = Some(match result {
                Some(r) => r.merge(&child_aabb),
                None => child_aabb,
            });
        }
        result.unwrap_or_else(|| Aabb::new(Vec3::zeros(), Vec3::zeros()))
    }
    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let mut best_dot = Real::NEG_INFINITY;
        let mut best_point = Vec3::zeros();
        for (transform, shape) in &self.children {
            let local_dir = transform.inverse().transform_vector(direction);
            let local_support = shape.support_point(&local_dir);
            let world_support = transform.transform_point(&local_support);
            let d = world_support.dot(direction);
            if d > best_dot {
                best_dot = d;
                best_point = world_support;
            }
        }
        best_point
    }
    fn volume(&self) -> Real {
        self.children.iter().map(|(_, s)| s.volume()).sum()
    }
    fn center_of_mass(&self) -> Vec3 {
        let total_vol: Real = self.children.iter().map(|(_, s)| s.volume()).sum();
        if total_vol < 1e-12 {
            return Vec3::zeros();
        }
        let weighted: Vec3 = self
            .children
            .iter()
            .map(|(t, s)| {
                let local_com = s.center_of_mass();
                let world_com = t.transform_point(&local_com);
                world_com * s.volume()
            })
            .sum();
        weighted / total_vol
    }
    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let total_vol: Real = self.children.iter().map(|(_, s)| s.volume()).sum();
        if total_vol < 1e-12 {
            return Mat3::zeros();
        }
        let com = self.center_of_mass();
        let mut total = Mat3::zeros();
        for (transform, shape) in &self.children {
            let child_mass = mass * shape.volume() / total_vol;
            let child_inertia = shape.inertia_tensor(child_mass);
            let child_com = transform.transform_point(&shape.center_of_mass());
            let r = child_com - com;
            let r2 = r.dot(&r);
            let parallel = Mat3::identity() * r2 - r * r.transpose();
            total += child_inertia + parallel * child_mass;
        }
        total
    }
    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        let mut best: Option<RayHit> = None;
        for (transform, shape) in &self.children {
            let inv = transform.inverse();
            let local_origin = inv.transform_point(ray_origin);
            let local_dir = inv.transform_vector(ray_direction);
            if let Some(hit) = shape.ray_cast(&local_origin, &local_dir, max_toi)
                && best.as_ref().is_none_or(|b| hit.toi < b.toi)
            {
                let world_point = transform.transform_point(&hit.point);
                let world_normal = transform.transform_vector(&hit.normal).normalize();
                best = Some(RayHit {
                    point: world_point,
                    normal: world_normal,
                    toi: hit.toi,
                });
            }
        }
        best
    }
}
