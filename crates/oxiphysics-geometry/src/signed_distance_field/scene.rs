// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SDF scene: composable scene of multiple SDF objects with Boolean operations.

use super::combinators::{sdf_smooth_union, sdf_union};
use super::helpers::{scale3, sub3};
use super::operators::sdf_normal;
use super::primitives::{sdf_box, sdf_capsule, sdf_cylinder, sdf_plane, sdf_sphere, sdf_torus};
use super::ray_marcher::{RayMarchHit, RayMarcher};

// ─────────────────────────────────────────────────────────────────────────────
// SDF Scene
// ─────────────────────────────────────────────────────────────────────────────

/// Composable SDF object with a transform and operation.
#[derive(Debug, Clone)]
pub struct SdfObject {
    /// Object name for debugging.
    pub name: String,
    /// Translation applied before evaluating SDF.
    pub translation: [f64; 3],
    /// Uniform scale factor.
    pub scale: f64,
    /// SDF primitive kind.
    pub kind: SdfKind,
}

/// Discriminated union of SDF primitive types.
#[derive(Debug, Clone)]
pub enum SdfKind {
    /// Sphere with radius.
    Sphere(f64),
    /// Axis-aligned box with half-extents.
    Box([f64; 3]),
    /// Capsule from A to B with radius.
    Capsule([f64; 3], [f64; 3], f64),
    /// Cylinder with radius and half-height.
    Cylinder(f64, f64),
    /// Torus with major and minor radii.
    Torus(f64, f64),
    /// Plane with normal and offset.
    Plane([f64; 3], f64),
}

impl SdfObject {
    /// Construct an SDF sphere object.
    pub fn sphere(name: &str, radius: f64, translation: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            translation,
            scale: 1.0,
            kind: SdfKind::Sphere(radius),
        }
    }

    /// Construct an SDF box object.
    pub fn box_shape(name: &str, half_extents: [f64; 3], translation: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            translation,
            scale: 1.0,
            kind: SdfKind::Box(half_extents),
        }
    }

    /// Evaluate the SDF at world-space point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        // Transform point to local space
        let local = scale3(sub3(p, self.translation), 1.0 / self.scale);
        let raw = match &self.kind {
            SdfKind::Sphere(r) => sdf_sphere(local, *r),
            SdfKind::Box(b) => sdf_box(local, *b),
            SdfKind::Capsule(a, b, r) => sdf_capsule(local, *a, *b, *r),
            SdfKind::Cylinder(r, h) => sdf_cylinder(local, *r, *h),
            SdfKind::Torus(r1, r2) => sdf_torus(local, *r1, *r2),
            SdfKind::Plane(n, d) => sdf_plane(local, *n, *d),
        };
        raw * self.scale
    }
}

/// An SDF scene composed of multiple objects with Boolean operations.
#[derive(Debug, Clone, Default)]
pub struct SdfScene {
    /// Objects in the scene.
    pub objects: Vec<SdfObject>,
    /// Smooth blend radius for scene-level union.
    pub blend_k: f64,
}

impl SdfScene {
    /// Construct an empty scene.
    pub fn new() -> Self {
        Self {
            objects: Vec::new(),
            blend_k: 0.0,
        }
    }

    /// Add an object to the scene.
    pub fn add(&mut self, obj: SdfObject) {
        self.objects.push(obj);
    }

    /// Evaluate the scene SDF (smooth union of all objects) at point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        if self.objects.is_empty() {
            return f64::MAX;
        }
        let mut d = self.objects[0].evaluate(p);
        for obj in &self.objects[1..] {
            let di = obj.evaluate(p);
            d = if self.blend_k > 0.0 {
                sdf_smooth_union(d, di, self.blend_k)
            } else {
                sdf_union(d, di)
            };
        }
        d
    }

    /// Cast a ray through the scene using sphere tracing.
    pub fn ray_cast(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<RayMarchHit> {
        let marcher = RayMarcher::new();
        marcher.march(&|p| self.evaluate(p), origin, dir)
    }

    /// Compute gradient (surface normal) at point `p`.
    pub fn normal(&self, p: [f64; 3]) -> [f64; 3] {
        sdf_normal(&|q| self.evaluate(q), p, 1e-4)
    }
}
