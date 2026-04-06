//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::functions::*;
use crate::shape::Shape;
use oxiphysics_core::Transform;
use std::sync::Arc;

/// Classification of primitive shapes for lightweight compound representations.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum ChildShapeKind {
    /// A sphere with given radius.
    Sphere {
        /// Radius of the sphere.
        radius: f64,
    },
    /// An axis-aligned box with half extents.
    Box {
        /// Half extents \[hx, hy, hz\].
        half_extents: [f64; 3],
    },
    /// A capsule (cylinder with hemispherical caps) along the Y axis.
    Capsule {
        /// Radius of the capsule.
        radius: f64,
        /// Half-height of the cylindrical part.
        half_height: f64,
    },
}
/// A child shape in a lightweight compound, with center position and shape kind.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CompoundChild {
    /// Center position of the child shape.
    pub center: [f64; 3],
    /// The kind of primitive shape.
    pub shape_kind: ChildShapeKind,
}
/// All per-child AABBs plus the merged bounding box.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CompoundAabb {
    /// Individual child AABBs (min, max).
    pub all_aabbs: Vec<([f64; 3], [f64; 3])>,
    /// Merged bounding box over all children.
    pub merged_min: [f64; 3],
    /// Merged bounding box over all children.
    pub merged_max: [f64; 3],
}
/// A 3-D rigid transform: translation + 3×3 rotation matrix (row-major).
///
/// The rotation `rot[i]` is the i-th row of the rotation matrix.
/// The matrix must be orthonormal; no validation is performed.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LocalTransform {
    /// Translation vector.
    pub translation: [f64; 3],
    /// Rotation matrix stored row-major: `rot[row][col]`.
    pub rot: [[f64; 3]; 3],
}
#[allow(dead_code)]
impl LocalTransform {
    /// Identity transform (no rotation, no translation).
    pub fn identity() -> Self {
        Self {
            translation: [0.0; 3],
            rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }
    /// Pure translation (no rotation).
    pub fn from_translation(t: [f64; 3]) -> Self {
        Self {
            translation: t,
            rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }
    /// Transform a point from local space to world space.
    ///
    /// `world_p = rot * local_p + translation`
    pub fn local_to_world(&self, p: [f64; 3]) -> [f64; 3] {
        let mut out = self.translation;
        for i in 0..3 {
            out[i] += self.rot[i][0] * p[0] + self.rot[i][1] * p[1] + self.rot[i][2] * p[2];
        }
        out
    }
    /// Transform a point from world space to local space.
    ///
    /// `local_p = rot^T * (world_p - translation)`
    pub fn world_to_local(&self, p: [f64; 3]) -> [f64; 3] {
        let d = [
            p[0] - self.translation[0],
            p[1] - self.translation[1],
            p[2] - self.translation[2],
        ];
        [
            self.rot[0][0] * d[0] + self.rot[1][0] * d[1] + self.rot[2][0] * d[2],
            self.rot[0][1] * d[0] + self.rot[1][1] * d[1] + self.rot[2][1] * d[2],
            self.rot[0][2] * d[0] + self.rot[1][2] * d[1] + self.rot[2][2] * d[2],
        ]
    }
    /// Transform a direction vector from local space to world space (no translation).
    pub fn local_to_world_dir(&self, v: [f64; 3]) -> [f64; 3] {
        [
            self.rot[0][0] * v[0] + self.rot[0][1] * v[1] + self.rot[0][2] * v[2],
            self.rot[1][0] * v[0] + self.rot[1][1] * v[1] + self.rot[1][2] * v[2],
            self.rot[2][0] * v[0] + self.rot[2][1] * v[1] + self.rot[2][2] * v[2],
        ]
    }
    /// Transform a direction vector from world space to local space (no translation).
    pub fn world_to_local_dir(&self, v: [f64; 3]) -> [f64; 3] {
        [
            self.rot[0][0] * v[0] + self.rot[1][0] * v[1] + self.rot[2][0] * v[2],
            self.rot[0][1] * v[0] + self.rot[1][1] * v[1] + self.rot[2][1] * v[2],
            self.rot[0][2] * v[0] + self.rot[1][2] * v[1] + self.rot[2][2] * v[2],
        ]
    }
}
/// A `CompoundShape` variant that stores each child with an explicit
/// [`LocalTransform`], enabling full 6-DOF placement (translation + rotation).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CompoundShapeEx {
    /// Children: (transform, shape kind).
    pub children: Vec<(LocalTransform, ChildShapeKind)>,
}
impl Default for CompoundShapeEx {
    fn default() -> Self {
        Self::new()
    }
}
#[allow(dead_code)]
impl CompoundShapeEx {
    /// Create an empty compound shape.
    pub fn new() -> Self {
        Self {
            children: Vec::new(),
        }
    }
    /// Add a sphere child with a local transform.
    pub fn add_sphere(&mut self, transform: LocalTransform, radius: f64) {
        self.children
            .push((transform, ChildShapeKind::Sphere { radius }));
    }
    /// Add a box child with a local transform.
    pub fn add_box(&mut self, transform: LocalTransform, half_extents: [f64; 3]) {
        self.children
            .push((transform, ChildShapeKind::Box { half_extents }));
    }
    /// Add a capsule child with a local transform.
    pub fn add_capsule(&mut self, transform: LocalTransform, radius: f64, half_height: f64) {
        self.children.push((
            transform,
            ChildShapeKind::Capsule {
                radius,
                half_height,
            },
        ));
    }
    /// Compute the union AABB of all children in world space.
    ///
    /// Returns `(min, max)` as `[f64; 3]` arrays.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        if self.children.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for (transform, kind) in &self.children {
            let local_he = match kind {
                ChildShapeKind::Sphere { radius } => [*radius, *radius, *radius],
                ChildShapeKind::Box { half_extents } => *half_extents,
                ChildShapeKind::Capsule {
                    radius,
                    half_height,
                } => [*radius, half_height + radius, *radius],
            };
            let corners_local = [
                [-local_he[0], -local_he[1], -local_he[2]],
                [local_he[0], -local_he[1], -local_he[2]],
                [-local_he[0], local_he[1], -local_he[2]],
                [local_he[0], local_he[1], -local_he[2]],
                [-local_he[0], -local_he[1], local_he[2]],
                [local_he[0], -local_he[1], local_he[2]],
                [-local_he[0], local_he[1], local_he[2]],
                [local_he[0], local_he[1], local_he[2]],
            ];
            for corner in &corners_local {
                let w = transform.local_to_world(*corner);
                for i in 0..3 {
                    if w[i] < min[i] {
                        min[i] = w[i];
                    }
                    if w[i] > max[i] {
                        max[i] = w[i];
                    }
                }
            }
        }
        (min, max)
    }
    /// Test if a world-space point is inside any child shape.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for (transform, kind) in &self.children {
            let local_p = transform.world_to_local(p);
            if child_kind_contains(kind, local_p) {
                return true;
            }
        }
        false
    }
    /// Ray cast against all children; returns `(toi, normal)` for the nearest hit.
    pub fn ray_cast(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<(f64, [f64; 3])> {
        let mut best: Option<(f64, [f64; 3])> = None;
        for (transform, kind) in &self.children {
            let local_o = transform.world_to_local(origin);
            let local_d = transform.world_to_local_dir(dir);
            if let Some((t, local_n)) = ray_cast_kind(kind, local_o, local_d, f64::MAX * 0.5) {
                let world_n = transform.local_to_world_dir(local_n);
                if best.as_ref().is_none_or(|(bt, _)| t < *bt) {
                    best = Some((t, world_n));
                }
            }
        }
        best
    }
    /// Total volume of all children.
    pub fn volume(&self) -> f64 {
        self.children
            .iter()
            .map(|(_, k)| CompoundShape::child_volume(k))
            .sum()
    }
    /// Compute the inertia tensor about the world-origin using the parallel axis theorem.
    ///
    /// Assumes uniform density `density`.
    #[allow(clippy::too_many_arguments)]
    pub fn inertia_tensor(&self, density: f64) -> [[f64; 3]; 3] {
        let mut i_xx = 0.0f64;
        let mut i_yy = 0.0f64;
        let mut i_zz = 0.0f64;
        let mut i_xy = 0.0f64;
        let mut i_xz = 0.0f64;
        let mut i_yz = 0.0f64;
        for (transform, kind) in &self.children {
            let vol = CompoundShape::child_volume(kind);
            let m = density * vol;
            let (lxx, lyy, lzz) = CompoundShape::child_local_inertia(kind, m);
            let r = transform.local_to_world([0.0; 3]);
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            i_xx += lxx + m * (r2 - r[0] * r[0]);
            i_yy += lyy + m * (r2 - r[1] * r[1]);
            i_zz += lzz + m * (r2 - r[2] * r[2]);
            i_xy -= m * r[0] * r[1];
            i_xz -= m * r[0] * r[2];
            i_yz -= m * r[1] * r[2];
        }
        [[i_xx, i_xy, i_xz], [i_xy, i_yy, i_yz], [i_xz, i_yz, i_zz]]
    }
}
/// A compound shape made of multiple sub-shapes, each with its own transform.
#[derive(Debug, Clone)]
pub struct Compound {
    /// Sub-shapes with their local transforms.
    pub children: Vec<(Transform, Arc<dyn Shape>)>,
}
impl Compound {
    /// Create a new compound shape.
    pub fn new(children: Vec<(Transform, Arc<dyn Shape>)>) -> Self {
        Self { children }
    }
}
/// A compound shape made of multiple lightweight primitive children.
///
/// Unlike `Compound` which uses `Arc<dyn Shape>`, this uses concrete enum
/// variants for common shapes, avoiding dynamic dispatch and allocations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CompoundShape {
    /// The child shapes.
    pub children: Vec<CompoundChild>,
}
impl Default for CompoundShape {
    fn default() -> Self {
        Self::new()
    }
}
#[allow(dead_code)]
impl CompoundShape {
    /// Create an empty compound shape.
    pub fn new() -> Self {
        Self {
            children: Vec::new(),
        }
    }
    /// Add a sphere child.
    pub fn add_sphere(&mut self, center: [f64; 3], radius: f64) {
        self.children.push(CompoundChild {
            center,
            shape_kind: ChildShapeKind::Sphere { radius },
        });
    }
    /// Add a box child.
    pub fn add_box(&mut self, center: [f64; 3], half_extents: [f64; 3]) {
        self.children.push(CompoundChild {
            center,
            shape_kind: ChildShapeKind::Box { half_extents },
        });
    }
    /// Add a capsule child.
    pub fn add_capsule(&mut self, center: [f64; 3], radius: f64, half_height: f64) {
        self.children.push(CompoundChild {
            center,
            shape_kind: ChildShapeKind::Capsule {
                radius,
                half_height,
            },
        });
    }
    /// Return the number of children.
    pub fn child_count(&self) -> usize {
        self.children.len()
    }
    /// Compute the volume of a single child shape kind.
    fn child_volume(kind: &ChildShapeKind) -> f64 {
        match kind {
            ChildShapeKind::Sphere { radius } => {
                (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius
            }
            ChildShapeKind::Box { half_extents } => {
                8.0 * half_extents[0] * half_extents[1] * half_extents[2]
            }
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => {
                let sphere_vol = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
                let cyl_vol = std::f64::consts::PI * radius * radius * 2.0 * half_height;
                sphere_vol + cyl_vol
            }
        }
    }
    /// Compute total volume of all children (no overlap correction).
    pub fn total_volume(&self) -> f64 {
        self.children
            .iter()
            .map(|c| Self::child_volume(&c.shape_kind))
            .sum()
    }
    /// Compute axis-aligned bounding box of all children.
    ///
    /// Returns `(min, max)` as `[f64; 3]` arrays.
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        if self.children.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for child in &self.children {
            let (cmin, cmax) = Self::child_aabb(child);
            for i in 0..3 {
                if cmin[i] < min[i] {
                    min[i] = cmin[i];
                }
                if cmax[i] > max[i] {
                    max[i] = cmax[i];
                }
            }
        }
        (min, max)
    }
    /// Compute the AABB of a single child.
    fn child_aabb(child: &CompoundChild) -> ([f64; 3], [f64; 3]) {
        let c = child.center;
        match child.shape_kind {
            ChildShapeKind::Sphere { radius } => (
                [c[0] - radius, c[1] - radius, c[2] - radius],
                [c[0] + radius, c[1] + radius, c[2] + radius],
            ),
            ChildShapeKind::Box { half_extents } => (
                [
                    c[0] - half_extents[0],
                    c[1] - half_extents[1],
                    c[2] - half_extents[2],
                ],
                [
                    c[0] + half_extents[0],
                    c[1] + half_extents[1],
                    c[2] + half_extents[2],
                ],
            ),
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => (
                [c[0] - radius, c[1] - half_height - radius, c[2] - radius],
                [c[0] + radius, c[1] + half_height + radius, c[2] + radius],
            ),
        }
    }
    /// Compute volume-weighted center of mass.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let total_vol = self.total_volume();
        if total_vol < 1e-12 {
            return [0.0; 3];
        }
        let mut com = [0.0; 3];
        for child in &self.children {
            let v = Self::child_volume(&child.shape_kind);
            for i in 0..3 {
                com[i] += child.center[i] * v;
            }
        }
        for i in 0..3 {
            com[i] /= total_vol;
        }
        com
    }
    /// Test if a point is inside any child shape.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for child in &self.children {
            if Self::child_contains(child, p) {
                return true;
            }
        }
        false
    }
    /// Check if a single child contains a point.
    fn child_contains(child: &CompoundChild, p: [f64; 3]) -> bool {
        let dx = p[0] - child.center[0];
        let dy = p[1] - child.center[1];
        let dz = p[2] - child.center[2];
        match child.shape_kind {
            ChildShapeKind::Sphere { radius } => dx * dx + dy * dy + dz * dz <= radius * radius,
            ChildShapeKind::Box { half_extents } => {
                dx.abs() <= half_extents[0]
                    && dy.abs() <= half_extents[1]
                    && dz.abs() <= half_extents[2]
            }
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => {
                let clamped_y = dy.clamp(-half_height, half_height);
                let ry = dy - clamped_y;
                dx * dx + ry * ry + dz * dz <= radius * radius
            }
        }
    }
    /// Ray cast against all children, returning `(toi, normal, child_index)` for
    /// the closest hit within `max_toi`.
    pub fn ray_cast(
        &self,
        origin: [f64; 3],
        dir: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3], usize)> {
        let mut best: Option<(f64, [f64; 3], usize)> = None;
        for (idx, child) in self.children.iter().enumerate() {
            if let Some((t, n)) = Self::ray_cast_child(child, origin, dir, max_toi)
                && best.as_ref().is_none_or(|(bt, _, _)| t < *bt)
            {
                best = Some((t, n, idx));
            }
        }
        best
    }
    /// Ray cast against a single child.
    fn ray_cast_child(
        child: &CompoundChild,
        origin: [f64; 3],
        dir: [f64; 3],
        max_toi: f64,
    ) -> Option<(f64, [f64; 3])> {
        let lo = [
            origin[0] - child.center[0],
            origin[1] - child.center[1],
            origin[2] - child.center[2],
        ];
        match child.shape_kind {
            ChildShapeKind::Sphere { radius } => ray_sphere(lo, dir, radius, max_toi),
            ChildShapeKind::Box { half_extents } => ray_box(lo, dir, half_extents, max_toi),
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => ray_capsule(lo, dir, radius, half_height, max_toi),
        }
    }
}
impl CompoundShape {
    /// Compute the inertia tensor of the compound shape about its center of mass.
    ///
    /// Uses the parallel axis theorem. `density` is the uniform mass density.
    pub fn inertia_tensor(&self, density: f64) -> [[f64; 3]; 3] {
        let com = self.center_of_mass();
        let mut i_xx = 0.0f64;
        let mut i_yy = 0.0f64;
        let mut i_zz = 0.0f64;
        let mut i_xy = 0.0f64;
        let mut i_xz = 0.0f64;
        let mut i_yz = 0.0f64;
        for child in &self.children {
            let vol = Self::child_volume(&child.shape_kind);
            let m = density * vol;
            let (lxx, lyy, lzz) = Self::child_local_inertia(&child.shape_kind, m);
            let r = [
                child.center[0] - com[0],
                child.center[1] - com[1],
                child.center[2] - com[2],
            ];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            i_xx += lxx + m * (r2 - r[0] * r[0]);
            i_yy += lyy + m * (r2 - r[1] * r[1]);
            i_zz += lzz + m * (r2 - r[2] * r[2]);
            i_xy -= m * r[0] * r[1];
            i_xz -= m * r[0] * r[2];
            i_yz -= m * r[1] * r[2];
        }
        [[i_xx, i_xy, i_xz], [i_xy, i_yy, i_yz], [i_xz, i_yz, i_zz]]
    }
    /// Local principal inertia components (I_xx, I_yy, I_zz) for a child at its own COM.
    fn child_local_inertia(kind: &ChildShapeKind, mass: f64) -> (f64, f64, f64) {
        match kind {
            ChildShapeKind::Sphere { radius } => {
                let i = 2.0 / 5.0 * mass * radius * radius;
                (i, i, i)
            }
            ChildShapeKind::Box { half_extents } => {
                let [hx, hy, hz] = *half_extents;
                let i_xx = mass / 3.0 * (hy * hy + hz * hz);
                let i_yy = mass / 3.0 * (hx * hx + hz * hz);
                let i_zz = mass / 3.0 * (hx * hx + hy * hy);
                (i_xx, i_yy, i_zz)
            }
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => {
                let r = radius;
                let h = half_height * 2.0;
                let m_cyl = mass * std::f64::consts::PI * r * r * h
                    / (std::f64::consts::PI * r * r * h
                        + (4.0 / 3.0) * std::f64::consts::PI * r * r * r);
                let m_hemi = (mass - m_cyl) / 2.0;
                let i_cyl_xx = m_cyl * (3.0 * r * r + h * h) / 12.0;
                let i_hemi_xx = m_hemi * (2.0 * r * r / 5.0 + (3.0 * half_height / 8.0).powi(2));
                let i_xx = i_cyl_xx + 2.0 * i_hemi_xx;
                let i_yy = m_cyl * r * r / 2.0 + 2.0 * m_hemi * 2.0 * r * r / 5.0;
                (i_xx, i_yy, i_xx)
            }
        }
    }
    /// Compute the bounding sphere (center, radius) of the compound shape.
    pub fn bounding_sphere(&self) -> ([f64; 3], f64) {
        let com = self.center_of_mass();
        let mut max_r = 0.0f64;
        for child in &self.children {
            let child_r = self.child_bounding_radius(child);
            let dist_to_com = {
                let dx = child.center[0] - com[0];
                let dy = child.center[1] - com[1];
                let dz = child.center[2] - com[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            let r = dist_to_com + child_r;
            if r > max_r {
                max_r = r;
            }
        }
        (com, max_r)
    }
    fn child_bounding_radius(&self, child: &CompoundChild) -> f64 {
        match child.shape_kind {
            ChildShapeKind::Sphere { radius } => radius,
            ChildShapeKind::Box { half_extents } => {
                (half_extents[0].powi(2) + half_extents[1].powi(2) + half_extents[2].powi(2)).sqrt()
            }
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => half_height + radius,
        }
    }
    /// Apply a uniform scaling to all child shape positions and sizes.
    pub fn scale(&mut self, factor: f64) {
        for child in &mut self.children {
            child.center[0] *= factor;
            child.center[1] *= factor;
            child.center[2] *= factor;
            match &mut child.shape_kind {
                ChildShapeKind::Sphere { radius } => *radius *= factor,
                ChildShapeKind::Box { half_extents } => {
                    half_extents[0] *= factor;
                    half_extents[1] *= factor;
                    half_extents[2] *= factor;
                }
                ChildShapeKind::Capsule {
                    radius,
                    half_height,
                } => {
                    *radius *= factor;
                    *half_height *= factor;
                }
            }
        }
    }
    /// Translate all child centers by the given offset.
    pub fn translate(&mut self, offset: [f64; 3]) {
        for child in &mut self.children {
            child.center[0] += offset[0];
            child.center[1] += offset[1];
            child.center[2] += offset[2];
        }
    }
    /// Returns a new compound shape merged with another (concatenation of children).
    pub fn merge_with(&self, other: &CompoundShape) -> CompoundShape {
        let mut result = self.clone();
        result.children.extend(other.children.iter().cloned());
        result
    }
    /// Compute the overlap (intersection test) with a sphere at `center` with `radius`.
    ///
    /// Returns true if any child overlaps the query sphere.
    pub fn overlaps_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        for child in &self.children {
            let dx = child.center[0] - center[0];
            let dy = child.center[1] - center[1];
            let dz = child.center[2] - center[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let child_r = self.child_bounding_radius(child);
            if dist < child_r + radius {
                return true;
            }
        }
        false
    }
    /// Recursively ray-cast with early exit once `max_hits` are found.
    pub fn ray_cast_all(
        &self,
        origin: [f64; 3],
        dir: [f64; 3],
        max_toi: f64,
    ) -> Vec<(f64, [f64; 3], usize)> {
        let mut hits: Vec<(f64, [f64; 3], usize)> = Vec::new();
        for (idx, child) in self.children.iter().enumerate() {
            if let Some((t, n)) = Self::ray_cast_child(child, origin, dir, max_toi) {
                hits.push((t, n, idx));
            }
        }
        hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }
}
#[allow(dead_code)]
impl CompoundShape {
    /// Compute merged AABB over all children.
    pub fn merged_aabb(&self) -> ([f64; 3], [f64; 3]) {
        self.aabb()
    }
    /// Ray cast returning `(t, child_index)` for the closest hit within `max_t`.
    ///
    /// Unlike `ray_cast` this omits the normal to match the requested signature.
    pub fn raycast(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_t: f64,
    ) -> Option<(f64, usize)> {
        self.ray_cast(ray_origin, ray_dir, max_t)
            .map(|(t, _n, idx)| (t, idx))
    }
    /// Total volume of all children (sum of child volumes).
    pub fn volume(&self) -> f64 {
        self.total_volume()
    }
    /// Mass-weighted center of mass.
    ///
    /// `masses[i]` is the mass of child `i`.  If `masses` is shorter than
    /// `children`, remaining children have zero mass.
    pub fn center_of_mass_weighted(&self, masses: &[f64]) -> [f64; 3] {
        let total: f64 = masses.iter().copied().take(self.children.len()).sum();
        if total < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0f64; 3];
        for (i, child) in self.children.iter().enumerate() {
            let m = if i < masses.len() { masses[i] } else { 0.0 };
            for k in 0..3 {
                com[k] += m * child.center[k];
            }
        }
        for k in 0..3 {
            com[k] /= total;
        }
        com
    }
    /// Inertia tensor about the compound center of mass using the parallel axis theorem.
    ///
    /// `masses[i]` is the mass of child `i`.
    pub fn inertia_tensor_from_masses(&self, masses: &[f64]) -> [[f64; 3]; 3] {
        let com = self.center_of_mass_weighted(masses);
        let mut i_xx = 0.0f64;
        let mut i_yy = 0.0f64;
        let mut i_zz = 0.0f64;
        let mut i_xy = 0.0f64;
        let mut i_xz = 0.0f64;
        let mut i_yz = 0.0f64;
        for (idx, child) in self.children.iter().enumerate() {
            let m = if idx < masses.len() { masses[idx] } else { 0.0 };
            let (lxx, lyy, lzz) = Self::child_local_inertia(&child.shape_kind, m);
            let r = [
                child.center[0] - com[0],
                child.center[1] - com[1],
                child.center[2] - com[2],
            ];
            let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            i_xx += lxx + m * (r2 - r[0] * r[0]);
            i_yy += lyy + m * (r2 - r[1] * r[1]);
            i_zz += lzz + m * (r2 - r[2] * r[2]);
            i_xy -= m * r[0] * r[1];
            i_xz -= m * r[0] * r[2];
            i_yz -= m * r[1] * r[2];
        }
        [[i_xx, i_xy, i_xz], [i_xy, i_yy, i_yz], [i_xz, i_yz, i_zz]]
    }
    /// Closest surface point to `p` among all children.
    ///
    /// Returns `(closest_point, child_index)`.
    pub fn closest_point(&self, p: [f64; 3]) -> ([f64; 3], usize) {
        let mut best_dist = f64::INFINITY;
        let mut best_pt = p;
        let mut best_idx = 0usize;
        for (i, child) in self.children.iter().enumerate() {
            let cp = Self::child_closest_point(child, p);
            let dx = cp[0] - p[0];
            let dy = cp[1] - p[1];
            let dz = cp[2] - p[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist < best_dist {
                best_dist = dist;
                best_pt = cp;
                best_idx = i;
            }
        }
        (best_pt, best_idx)
    }
    fn child_closest_point(child: &CompoundChild, p: [f64; 3]) -> [f64; 3] {
        let c = child.center;
        match child.shape_kind {
            ChildShapeKind::Sphere { radius } => {
                let dx = p[0] - c[0];
                let dy = p[1] - c[1];
                let dz = p[2] - c[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < 1e-30 {
                    [c[0] + radius, c[1], c[2]]
                } else {
                    let scale = radius / dist;
                    [c[0] + dx * scale, c[1] + dy * scale, c[2] + dz * scale]
                }
            }
            ChildShapeKind::Box { half_extents } => {
                let clamped = [
                    (p[0] - c[0]).clamp(-half_extents[0], half_extents[0]) + c[0],
                    (p[1] - c[1]).clamp(-half_extents[1], half_extents[1]) + c[1],
                    (p[2] - c[2]).clamp(-half_extents[2], half_extents[2]) + c[2],
                ];
                let inside = (0..3).all(|i| {
                    let he = [half_extents[0], half_extents[1], half_extents[2]][i];
                    let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]][i].abs();
                    d <= he
                });
                if inside {
                    let dx_neg = [half_extents[0], half_extents[1], half_extents[2]][0]
                        - (p[0] - c[0]).abs();
                    let dy_neg = [half_extents[0], half_extents[1], half_extents[2]][1]
                        - (p[1] - c[1]).abs();
                    let dz_neg = [half_extents[0], half_extents[1], half_extents[2]][2]
                        - (p[2] - c[2]).abs();
                    if dx_neg <= dy_neg && dx_neg <= dz_neg {
                        let sx = if p[0] >= c[0] { 1.0 } else { -1.0 };
                        [c[0] + half_extents[0] * sx, p[1], p[2]]
                    } else if dy_neg <= dx_neg && dy_neg <= dz_neg {
                        let sy = if p[1] >= c[1] { 1.0 } else { -1.0 };
                        [p[0], c[1] + half_extents[1] * sy, p[2]]
                    } else {
                        let sz = if p[2] >= c[2] { 1.0 } else { -1.0 };
                        [p[0], p[1], c[2] + half_extents[2] * sz]
                    }
                } else {
                    clamped
                }
            }
            ChildShapeKind::Capsule {
                radius,
                half_height,
            } => {
                let dy = p[1] - c[1];
                let clamped_y = dy.clamp(-half_height, half_height);
                let axis_pt = [c[0], c[1] + clamped_y, c[2]];
                let dx = p[0] - axis_pt[0];
                let dpz = p[2] - axis_pt[2];
                let dist_xz = (dx * dx + dpz * dpz).sqrt();
                if dist_xz < 1e-30 {
                    [axis_pt[0] + radius, axis_pt[1], axis_pt[2]]
                } else {
                    let scale = radius / dist_xz;
                    [
                        axis_pt[0] + dx * scale,
                        axis_pt[1],
                        axis_pt[2] + dpz * scale,
                    ]
                }
            }
        }
    }
}
#[allow(dead_code)]
impl CompoundShape {
    /// Remove the child at the given index.  Panics if `index` is out of range.
    pub fn remove_child(&mut self, index: usize) {
        self.children.remove(index);
    }
    /// Remove the child at the given index by swapping with the last element.
    ///
    /// Faster than `remove_child` (O(1) vs O(n)) but changes the order of
    /// remaining children.
    pub fn swap_remove_child(&mut self, index: usize) {
        self.children.swap_remove(index);
    }
    /// Replace the child at `index` with a new sphere.
    pub fn replace_with_sphere(&mut self, index: usize, center: [f64; 3], radius: f64) {
        self.children[index] = CompoundChild {
            center,
            shape_kind: ChildShapeKind::Sphere { radius },
        };
    }
    /// Replace the child at `index` with a new box.
    pub fn replace_with_box(&mut self, index: usize, center: [f64; 3], half_extents: [f64; 3]) {
        self.children[index] = CompoundChild {
            center,
            shape_kind: ChildShapeKind::Box { half_extents },
        };
    }
    /// Return `true` if the compound shape has no children.
    pub fn is_empty(&self) -> bool {
        self.children.is_empty()
    }
    /// Clear all children.
    pub fn clear(&mut self) {
        self.children.clear();
    }
    /// Closest surface point among all children together with the squared
    /// distance and the child index.
    ///
    /// Returns `(closest_point, squared_distance, child_index)`.
    pub fn closest_point_with_dist2(&self, p: [f64; 3]) -> ([f64; 3], f64, usize) {
        let mut best_dist2 = f64::INFINITY;
        let mut best_pt = p;
        let mut best_idx = 0usize;
        for (i, child) in self.children.iter().enumerate() {
            let cp = Self::child_closest_point(child, p);
            let dx = cp[0] - p[0];
            let dy = cp[1] - p[1];
            let dz = cp[2] - p[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_dist2 {
                best_dist2 = d2;
                best_pt = cp;
                best_idx = i;
            }
        }
        (best_pt, best_dist2, best_idx)
    }
    /// Check if two `CompoundShape` instances have any overlapping pair of
    /// children using conservative bounding-sphere overlap tests.
    ///
    /// Returns the index pairs `(i, j)` of all overlapping child pairs.
    pub fn broad_phase_pairs(&self, other: &CompoundShape) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for (i, ci) in self.children.iter().enumerate() {
            let ri = self.child_bounding_radius(ci);
            for (j, cj) in other.children.iter().enumerate() {
                let rj = other.child_bounding_radius(cj);
                let dx = ci.center[0] - cj.center[0];
                let dy = ci.center[1] - cj.center[1];
                let dz = ci.center[2] - cj.center[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < ri + rj {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
    /// Test if two `CompoundShape` instances overlap at all (broad phase).
    pub fn overlaps_compound(&self, other: &CompoundShape) -> bool {
        !self.broad_phase_pairs(other).is_empty()
    }
    /// Compute the centroid of the compound (volume-weighted center) using
    /// per-child densities.
    ///
    /// `densities[i]` is the density for child `i`.  If `densities` is shorter
    /// than `children`, remaining children use density 1.0.
    pub fn centroid_with_densities(&self, densities: &[f64]) -> [f64; 3] {
        let mut total_mass = 0.0f64;
        let mut com = [0.0f64; 3];
        for (i, child) in self.children.iter().enumerate() {
            let rho = if i < densities.len() {
                densities[i]
            } else {
                1.0
            };
            let vol = Self::child_volume(&child.shape_kind);
            let m = rho * vol;
            total_mass += m;
            for k in 0..3 {
                com[k] += m * child.center[k];
            }
        }
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        for k in 0..3 {
            com[k] /= total_mass;
        }
        com
    }
    /// Approximate penetration depth between this compound and a sphere at
    /// `center` with `radius`.
    ///
    /// For each child, computes the signed distance to the child's surface
    /// (negative inside).  Returns the minimum signed distance (most negative
    /// = deepest penetration) together with the child index.
    ///
    /// Returns `None` if there is no penetration.
    pub fn penetration_depth_sphere(&self, center: [f64; 3], radius: f64) -> Option<(f64, usize)> {
        let mut best: Option<(f64, usize)> = None;
        for (i, child) in self.children.iter().enumerate() {
            let signed = self.signed_distance_child(child, center, radius);
            if signed < 0.0 && best.as_ref().is_none_or(|(bd, _)| signed < *bd) {
                best = Some((signed, i));
            }
        }
        best
    }
    /// Signed distance between an external sphere and a single child.
    ///
    /// Negative when the sphere penetrates the child.
    fn signed_distance_child(
        &self,
        child: &CompoundChild,
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> f64 {
        let cp = Self::child_closest_point(child, sphere_center);
        let dx = cp[0] - sphere_center[0];
        let dy = cp[1] - sphere_center[1];
        let dz = cp[2] - sphere_center[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist - sphere_radius
    }
    /// Compute per-child masses from a uniform density.
    pub fn child_masses(&self, density: f64) -> Vec<f64> {
        self.children
            .iter()
            .map(|c| density * Self::child_volume(&c.shape_kind))
            .collect()
    }
    /// Total mass of the compound given uniform density.
    pub fn total_mass(&self, density: f64) -> f64 {
        density * self.total_volume()
    }
    /// Axis-aligned bounding box of a single child (public accessor).
    pub fn child_aabb_public(child: &CompoundChild) -> ([f64; 3], [f64; 3]) {
        Self::child_aabb(child)
    }
    /// Compute the AABB expanded by a margin `margin` on all sides.
    pub fn expanded_aabb(&self, margin: f64) -> ([f64; 3], [f64; 3]) {
        let (mn, mx) = self.aabb();
        (
            [mn[0] - margin, mn[1] - margin, mn[2] - margin],
            [mx[0] + margin, mx[1] + margin, mx[2] + margin],
        )
    }
    /// Test if a sphere (center, radius) overlaps the compound's AABB.
    pub fn sphere_overlaps_aabb(&self, center: [f64; 3], radius: f64) -> bool {
        if self.children.is_empty() {
            return false;
        }
        let (mn, mx) = self.aabb();
        let cx = center[0].clamp(mn[0], mx[0]);
        let cy = center[1].clamp(mn[1], mx[1]);
        let cz = center[2].clamp(mn[2], mx[2]);
        let dx = cx - center[0];
        let dy = cy - center[1];
        let dz = cz - center[2];
        dx * dx + dy * dy + dz * dz <= radius * radius
    }
}
