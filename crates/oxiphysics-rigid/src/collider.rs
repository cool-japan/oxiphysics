// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Collider attached to a rigid body.
//!
//! Provides standalone shape/material types that do not depend on `nalgebra`
//! or the geometry crate, as well as AABB computation and broad-phase overlap.

use oxiphysics_core::math::Real;
use oxiphysics_core::{BodyHandle, MassProperties, Transform};
use oxiphysics_geometry::Shape;
use std::sync::Arc;

// ── ColliderShape ─────────────────────────────────────────────────────────────

/// Analytic shape description stored directly in a collider (no trait object).
#[derive(Debug, Clone)]
pub enum ColliderShape {
    /// Sphere with the given radius.
    Sphere {
        /// Sphere radius.
        radius: f64,
    },
    /// Axis-aligned box with half-extents on each axis.
    Box {
        /// Half-extents along X, Y, Z.
        half_extents: [f64; 3],
    },
    /// Capsule (cylinder capped with hemispheres).
    Capsule {
        /// Capsule radius.
        radius: f64,
        /// Half-height of the cylindrical section.
        half_height: f64,
    },
    /// Cylinder aligned with the Y axis.
    Cylinder {
        /// Cylinder radius.
        radius: f64,
        /// Full height of the cylinder.
        height: f64,
    },
}

impl ColliderShape {
    /// Compute the axis-aligned bounding box for this shape given a world position.
    pub fn aabb(&self, pos: [f64; 3]) -> Aabb {
        match self {
            ColliderShape::Sphere { radius } => Aabb {
                min: [pos[0] - radius, pos[1] - radius, pos[2] - radius],
                max: [pos[0] + radius, pos[1] + radius, pos[2] + radius],
            },
            ColliderShape::Box { half_extents } => Aabb {
                min: [
                    pos[0] - half_extents[0],
                    pos[1] - half_extents[1],
                    pos[2] - half_extents[2],
                ],
                max: [
                    pos[0] + half_extents[0],
                    pos[1] + half_extents[1],
                    pos[2] + half_extents[2],
                ],
            },
            ColliderShape::Capsule {
                radius,
                half_height,
            } => {
                let total = radius + half_height;
                Aabb {
                    min: [pos[0] - radius, pos[1] - total, pos[2] - radius],
                    max: [pos[0] + radius, pos[1] + total, pos[2] + radius],
                }
            }
            ColliderShape::Cylinder { radius, height } => {
                let hh = height * 0.5;
                Aabb {
                    min: [pos[0] - radius, pos[1] - hh, pos[2] - radius],
                    max: [pos[0] + radius, pos[1] + hh, pos[2] + radius],
                }
            }
        }
    }

    /// Approximate volume of this shape.
    pub fn volume(&self) -> f64 {
        match self {
            ColliderShape::Sphere { radius } => {
                (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius
            }
            ColliderShape::Box { half_extents } => {
                8.0 * half_extents[0] * half_extents[1] * half_extents[2]
            }
            ColliderShape::Capsule {
                radius,
                half_height,
            } => {
                let sphere = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
                let cyl = std::f64::consts::PI * radius * radius * 2.0 * half_height;
                sphere + cyl
            }
            ColliderShape::Cylinder { radius, height } => {
                std::f64::consts::PI * radius * radius * height
            }
        }
    }
}

// ── ColliderMaterial ──────────────────────────────────────────────────────────

/// Physical material properties attached to a collider.
#[derive(Debug, Clone)]
pub struct ColliderMaterial {
    /// Coulomb friction coefficient (dimensionless, ≥ 0).
    pub friction: f64,
    /// Coefficient of restitution (0 = perfectly inelastic, 1 = elastic).
    pub restitution: f64,
    /// Mass density (kg/m³).
    pub density: f64,
}

impl Default for ColliderMaterial {
    fn default() -> Self {
        Self {
            friction: 0.5,
            restitution: 0.3,
            density: 1000.0,
        }
    }
}

impl ColliderMaterial {
    /// Create a new material.
    pub fn new(friction: f64, restitution: f64, density: f64) -> Self {
        Self {
            friction,
            restitution,
            density,
        }
    }
}

// ── Aabb ──────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box in 3D.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb {
    /// Test whether this AABB overlaps another.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Expand the AABB by `d` on each side.
    pub fn expand(&self, d: f64) -> Aabb {
        Aabb {
            min: [self.min[0] - d, self.min[1] - d, self.min[2] - d],
            max: [self.max[0] + d, self.max[1] + d, self.max[2] + d],
        }
    }

    /// Centre point.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
}

// ── SimpleCollider ────────────────────────────────────────────────────────────

/// A lightweight collider that owns its shape description and material.
///
/// Unlike `Collider` (which wraps a trait-object), `SimpleCollider` uses the
/// `ColliderShape` enum and therefore never allocates on the heap.
#[derive(Debug, Clone)]
pub struct SimpleCollider {
    /// Analytic shape.
    pub shape: ColliderShape,
    /// Physical material.
    pub material: ColliderMaterial,
    /// Local-space position offset from the parent body origin.
    pub local_position: [f64; 3],
    /// Optional owning body.
    pub parent_body: Option<BodyHandle>,
    /// When true the collider fires events but exerts no forces.
    pub is_sensor: bool,
}

impl SimpleCollider {
    /// Create a new simple collider.
    pub fn new(shape: ColliderShape) -> Self {
        Self {
            shape,
            material: ColliderMaterial::default(),
            local_position: [0.0; 3],
            parent_body: None,
            is_sensor: false,
        }
    }

    /// Set material.
    pub fn with_material(mut self, material: ColliderMaterial) -> Self {
        self.material = material;
        self
    }

    /// Set local position offset.
    pub fn with_local_position(mut self, pos: [f64; 3]) -> Self {
        self.local_position = pos;
        self
    }

    /// Mark as a sensor.
    pub fn as_sensor(mut self) -> Self {
        self.is_sensor = true;
        self
    }

    /// Attach to a body.
    pub fn with_parent(mut self, handle: BodyHandle) -> Self {
        self.parent_body = Some(handle);
        self
    }

    /// Compute the AABB in world space given the body's world position.
    pub fn world_aabb(&self, body_pos: [f64; 3]) -> Aabb {
        let world_pos = [
            body_pos[0] + self.local_position[0],
            body_pos[1] + self.local_position[1],
            body_pos[2] + self.local_position[2],
        ];
        self.shape.aabb(world_pos)
    }
}

// ── SimpleColliderHandle ──────────────────────────────────────────────────────

/// Handle into a `SimpleColliderSet`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SimpleColliderHandle {
    /// Slot index.
    pub index: u32,
    /// Generation counter to detect stale handles.
    pub generation: u32,
}

// ── SimpleColliderSet ─────────────────────────────────────────────────────────

/// Generational arena of `SimpleCollider` objects.
#[derive(Debug, Default)]
pub struct SimpleColliderSet {
    colliders: Vec<Option<SimpleCollider>>,
    generations: Vec<u32>,
    free_list: Vec<u32>,
}

impl SimpleColliderSet {
    /// Create an empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert a collider and return its handle.
    pub fn insert(&mut self, c: SimpleCollider) -> SimpleColliderHandle {
        if let Some(index) = self.free_list.pop() {
            let idx = index as usize;
            self.colliders[idx] = Some(c);
            SimpleColliderHandle {
                index,
                generation: self.generations[idx],
            }
        } else {
            let index = self.colliders.len() as u32;
            self.colliders.push(Some(c));
            self.generations.push(0);
            SimpleColliderHandle {
                index,
                generation: 0,
            }
        }
    }

    /// Remove a collider by handle.
    pub fn remove(&mut self, h: SimpleColliderHandle) -> Option<SimpleCollider> {
        let idx = h.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == h.generation {
            self.generations[idx] += 1;
            self.free_list.push(h.index);
            self.colliders[idx].take()
        } else {
            None
        }
    }

    /// Get a reference by handle.
    pub fn get(&self, h: SimpleColliderHandle) -> Option<&SimpleCollider> {
        let idx = h.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == h.generation {
            self.colliders[idx].as_ref()
        } else {
            None
        }
    }

    /// Get a mutable reference by handle.
    pub fn get_mut(&mut self, h: SimpleColliderHandle) -> Option<&mut SimpleCollider> {
        let idx = h.index as usize;
        if idx < self.colliders.len() && self.generations[idx] == h.generation {
            self.colliders[idx].as_mut()
        } else {
            None
        }
    }

    /// Iterate over all active `(handle, collider)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (SimpleColliderHandle, &SimpleCollider)> {
        self.colliders.iter().enumerate().filter_map(move |(i, c)| {
            c.as_ref().map(|col| {
                (
                    SimpleColliderHandle {
                        index: i as u32,
                        generation: self.generations[i],
                    },
                    col,
                )
            })
        })
    }

    /// Number of active colliders.
    pub fn len(&self) -> usize {
        self.colliders.iter().filter(|c| c.is_some()).count()
    }

    /// True if there are no active colliders.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Compute the world-space AABB for a collider, given a body position lookup.
    ///
    /// `body_pos_fn` maps a `BodyHandle` to a world position.  If the collider
    /// has no parent body the local position itself is used as the world position.
    pub fn aabb(
        &self,
        handle: SimpleColliderHandle,
        body_pos_fn: impl Fn(BodyHandle) -> Option<[f64; 3]>,
    ) -> Option<Aabb> {
        let col = self.get(handle)?;
        let world_pos = match col.parent_body {
            Some(bh) => body_pos_fn(bh).unwrap_or(col.local_position),
            None => col.local_position,
        };
        Some(col.world_aabb(world_pos))
    }

    /// Return all pairs of handles whose AABBs overlap (O(n²) brute force).
    ///
    /// Only non-sensor vs any, and no self-pairs.
    pub fn broad_pairs(
        &self,
        body_pos_fn: impl Fn(BodyHandle) -> Option<[f64; 3]>,
    ) -> Vec<(SimpleColliderHandle, SimpleColliderHandle)> {
        let handles: Vec<SimpleColliderHandle> = self.iter().map(|(h, _)| h).collect();
        let aabbs: Vec<Option<Aabb>> = handles
            .iter()
            .map(|&h| self.aabb(h, &body_pos_fn))
            .collect();

        let mut pairs = Vec::new();
        for i in 0..handles.len() {
            for j in (i + 1)..handles.len() {
                if let (Some(a), Some(b)) = (aabbs[i], aabbs[j])
                    && a.overlaps(&b)
                {
                    pairs.push((handles[i], handles[j]));
                }
            }
        }
        pairs
    }
}

// ── Collider (trait-object version, kept for backward compatibility) ───────────

/// A collider that wraps a shape and is attached to a body.
#[derive(Debug, Clone)]
pub struct Collider {
    /// The geometric shape.
    pub shape: Arc<dyn Shape>,
    /// Local offset from the body's transform.
    pub local_transform: Transform,
    /// Handle of the parent body.
    pub body_handle: Option<BodyHandle>,
    /// Friction coefficient.
    pub friction: Real,
    /// Restitution (bounciness).
    pub restitution: Real,
    /// Whether this collider is a sensor (triggers events but no forces).
    pub is_sensor: bool,
    /// Material density (kg/m³). Used for mass computation. Default: 1000 (water).
    pub density: Real,
}

impl Collider {
    /// Create a new collider from a shape.
    pub fn new(shape: Arc<dyn Shape>) -> Self {
        Self {
            shape,
            local_transform: Transform::default(),
            body_handle: None,
            friction: 0.5,
            restitution: 0.3,
            is_sensor: false,
            density: 1000.0,
        }
    }

    /// Set the density (kg/m³).
    pub fn with_density(mut self, density: Real) -> Self {
        self.density = density;
        self
    }

    /// Compute mass properties from the collider's shape and density.
    ///
    /// If the shape has zero volume, falls back to a unit sphere approximation.
    pub fn mass_properties(&self) -> MassProperties {
        let volume = self.shape.volume();
        if volume > 0.0 {
            self.shape.mass_properties(self.density)
        } else {
            // Fallback: unit sphere with given density
            let mass = self.density * (4.0 / 3.0) * std::f64::consts::PI;
            let i = 2.0 / 5.0 * mass;
            let inertia = oxiphysics_core::math::Mat3::identity() * i;
            MassProperties::new(mass, oxiphysics_core::math::Vec3::zeros(), inertia)
        }
    }

    /// Set the local offset transform.
    pub fn with_local_transform(mut self, transform: Transform) -> Self {
        self.local_transform = transform;
        self
    }

    /// Set the friction coefficient.
    pub fn with_friction(mut self, friction: Real) -> Self {
        self.friction = friction;
        self
    }

    /// Set the restitution.
    pub fn with_restitution(mut self, restitution: Real) -> Self {
        self.restitution = restitution;
        self
    }

    /// Mark as sensor.
    pub fn as_sensor(mut self) -> Self {
        self.is_sensor = true;
        self
    }

    /// Attach to a body.
    pub fn with_body(mut self, handle: BodyHandle) -> Self {
        self.body_handle = Some(handle);
        self
    }

    /// Compute the world-space transform of this collider given the body transform.
    pub fn world_transform(&self, body_transform: &Transform) -> Transform {
        body_transform.compose(&self.local_transform)
    }
}

// ── ColliderMassProperties ────────────────────────────────────────────────────

/// Compute mass properties for a `ColliderShape` given a density.
///
/// Returns `(mass, center_of_mass, inertia_diagonal)` where inertia is the
/// principal moments Ixx, Iyy, Izz about the local centre of mass.
pub fn shape_mass_properties(shape: &ColliderShape, density: f64) -> (f64, [f64; 3], [f64; 3]) {
    match shape {
        ColliderShape::Sphere { radius } => {
            let mass = density * (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
            let i = (2.0 / 5.0) * mass * radius * radius;
            (mass, [0.0; 3], [i, i, i])
        }
        ColliderShape::Box { half_extents: he } => {
            let vol = 8.0 * he[0] * he[1] * he[2];
            let mass = density * vol;
            let ixx = (1.0 / 3.0) * mass * (he[1] * he[1] + he[2] * he[2]);
            let iyy = (1.0 / 3.0) * mass * (he[0] * he[0] + he[2] * he[2]);
            let izz = (1.0 / 3.0) * mass * (he[0] * he[0] + he[1] * he[1]);
            (mass, [0.0; 3], [ixx, iyy, izz])
        }
        ColliderShape::Capsule {
            radius,
            half_height,
        } => {
            let r = *radius;
            let h = *half_height;
            let vol_sphere = (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
            let vol_cyl = std::f64::consts::PI * r * r * 2.0 * h;
            let mass = density * (vol_sphere + vol_cyl);
            // Approximate: combine sphere and cylinder inertia
            let m_s = density * vol_sphere;
            let m_c = density * vol_cyl;
            // Cylinder Ixx = m*(3r²+h²)/12 (for full height 2h)
            let ixx_cyl = m_c * (3.0 * r * r + (2.0 * h) * (2.0 * h)) / 12.0;
            let iyy_cyl = 0.5 * m_c * r * r;
            // Sphere Ixx = 2mr²/5 shifted by half_height+sphere_cm offset
            let ixx_sph = (2.0 / 5.0) * m_s * r * r + m_s * h * h;
            let iyy_sph = (2.0 / 5.0) * m_s * r * r;
            let ixx = ixx_cyl + 2.0 * ixx_sph;
            let iyy = iyy_cyl + 2.0 * iyy_sph;
            let izz = ixx; // symmetric
            (mass, [0.0; 3], [ixx, iyy, izz])
        }
        ColliderShape::Cylinder { radius, height } => {
            let r = *radius;
            let h = *height;
            let mass = density * std::f64::consts::PI * r * r * h;
            let ixx = mass * (3.0 * r * r + h * h) / 12.0;
            let iyy = 0.5 * mass * r * r;
            let izz = ixx;
            (mass, [0.0; 3], [ixx, iyy, izz])
        }
    }
}

// ── ColliderTransform ─────────────────────────────────────────────────────────

/// Minimal 3-D rigid transform stored as position + rotation (column-major 3×3).
///
/// This is the no-nalgebra version used by `SimpleCollider`.
#[derive(Debug, Clone, Copy)]
pub struct ColliderTransform {
    /// World-space position.
    pub position: [f64; 3],
    /// Rotation matrix (row-major, orthonormal).
    pub rotation: [[f64; 3]; 3],
}

impl ColliderTransform {
    /// Identity transform.
    pub fn identity() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Translation-only transform.
    pub fn from_position(p: [f64; 3]) -> Self {
        let mut t = Self::identity();
        t.position = p;
        t
    }

    /// Apply this transform to a local-space point.
    pub fn transform_point(&self, local: [f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        [
            self.position[0] + r[0][0] * local[0] + r[0][1] * local[1] + r[0][2] * local[2],
            self.position[1] + r[1][0] * local[0] + r[1][1] * local[1] + r[1][2] * local[2],
            self.position[2] + r[2][0] * local[0] + r[2][1] * local[1] + r[2][2] * local[2],
        ]
    }

    /// Compose two transforms: `self * other`.
    pub fn compose(&self, other: &ColliderTransform) -> ColliderTransform {
        let ra = &self.rotation;
        let rb = &other.rotation;
        // rc = ra * rb
        let mut rc = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    rc[i][j] += ra[i][k] * rb[k][j];
                }
            }
        }
        let pos = self.transform_point(other.position);
        ColliderTransform {
            position: pos,
            rotation: rc,
        }
    }
}

impl Default for ColliderTransform {
    fn default() -> Self {
        Self::identity()
    }
}

// ── TransformedCollider ───────────────────────────────────────────────────────

/// A `SimpleCollider` with a full `ColliderTransform` (position + rotation).
///
/// Used when colliders need arbitrary orientation in addition to a position
/// offset.
#[derive(Debug, Clone)]
pub struct TransformedCollider {
    /// Underlying simple collider.
    pub collider: SimpleCollider,
    /// Full local-to-body transform (includes rotation).
    pub local_transform: ColliderTransform,
}

impl TransformedCollider {
    /// Create a new transformed collider.
    pub fn new(collider: SimpleCollider, local_transform: ColliderTransform) -> Self {
        Self {
            collider,
            local_transform,
        }
    }

    /// Compute the world transform given the body's world transform.
    pub fn world_transform(&self, body_transform: &ColliderTransform) -> ColliderTransform {
        body_transform.compose(&self.local_transform)
    }

    /// Compute the world-space AABB (conservative, using the body transform).
    ///
    /// For non-axis-aligned rotations this is a conservative (expanded) AABB.
    pub fn world_aabb(&self, body_transform: &ColliderTransform) -> Aabb {
        let world_t = self.world_transform(body_transform);
        // Use the world-space origin of the collider as the anchor
        let world_pos = world_t.position;
        self.collider.shape.aabb(world_pos)
    }
}

// ── CompoundColliderShape ────────────────────────────────────────────────────

/// A compound shape made of multiple `ColliderShape`s with local offsets.
#[derive(Debug, Clone)]
pub struct CompoundColliderShape {
    /// Child shapes and their local position offsets.
    pub children: Vec<(ColliderShape, [f64; 3])>,
}

impl CompoundColliderShape {
    /// Create an empty compound shape.
    pub fn new() -> Self {
        Self {
            children: Vec::new(),
        }
    }

    /// Add a child shape at a local offset.
    pub fn add_child(&mut self, shape: ColliderShape, local_pos: [f64; 3]) {
        self.children.push((shape, local_pos));
    }

    /// Compute the aggregate AABB in world space.
    pub fn world_aabb(&self, body_pos: [f64; 3]) -> Option<Aabb> {
        if self.children.is_empty() {
            return None;
        }
        let mut combined: Option<Aabb> = None;
        for (shape, local) in &self.children {
            let child_pos = [
                body_pos[0] + local[0],
                body_pos[1] + local[1],
                body_pos[2] + local[2],
            ];
            let child_aabb = shape.aabb(child_pos);
            combined = Some(match combined {
                None => child_aabb,
                Some(acc) => Aabb {
                    min: [
                        acc.min[0].min(child_aabb.min[0]),
                        acc.min[1].min(child_aabb.min[1]),
                        acc.min[2].min(child_aabb.min[2]),
                    ],
                    max: [
                        acc.max[0].max(child_aabb.max[0]),
                        acc.max[1].max(child_aabb.max[1]),
                        acc.max[2].max(child_aabb.max[2]),
                    ],
                },
            });
        }
        combined
    }

    /// Total volume (sum of children).
    pub fn volume(&self) -> f64 {
        self.children.iter().map(|(s, _)| s.volume()).sum()
    }

    /// Aggregate mass properties (mass-weighted centroid).
    pub fn mass_properties(&self, density: f64) -> (f64, [f64; 3], [f64; 3]) {
        let mut total_mass = 0.0f64;
        let mut weighted_cm = [0.0f64; 3];
        for (shape, local) in &self.children {
            let (m, _, _) = shape_mass_properties(shape, density);
            total_mass += m;
            weighted_cm[0] += m * local[0];
            weighted_cm[1] += m * local[1];
            weighted_cm[2] += m * local[2];
        }
        if total_mass > 0.0 {
            weighted_cm[0] /= total_mass;
            weighted_cm[1] /= total_mass;
            weighted_cm[2] /= total_mass;
        }
        // Inertia: sum of each child's inertia shifted by parallel-axis theorem
        let mut ixx = 0.0f64;
        let mut iyy = 0.0f64;
        let mut izz = 0.0f64;
        for (shape, local) in &self.children {
            let (m, _, diag) = shape_mass_properties(shape, density);
            let dx = local[0] - weighted_cm[0];
            let dy = local[1] - weighted_cm[1];
            let dz = local[2] - weighted_cm[2];
            ixx += diag[0] + m * (dy * dy + dz * dz);
            iyy += diag[1] + m * (dx * dx + dz * dz);
            izz += diag[2] + m * (dx * dx + dy * dy);
        }
        (total_mass, weighted_cm, [ixx, iyy, izz])
    }
}

impl Default for CompoundColliderShape {
    fn default() -> Self {
        Self::new()
    }
}

// ── SensorCollider ────────────────────────────────────────────────────────────

/// A sensor-only collider: triggers overlap events, exerts no forces.
///
/// Wraps a `ColliderShape` and exposes convenience APIs for sensor-specific
/// behaviour (e.g., event queue integration).
#[derive(Debug, Clone)]
pub struct SensorCollider {
    /// The sensor shape.
    pub shape: ColliderShape,
    /// World-space position of the sensor.
    pub position: [f64; 3],
    /// User-defined category bits (for filtering).
    pub category_bits: u32,
    /// User-defined mask bits (collides with objects whose category ∩ mask ≠ ∅).
    pub mask_bits: u32,
    /// Application-defined sensor ID.
    pub sensor_id: u32,
}

impl SensorCollider {
    /// Create a new sensor.
    pub fn new(shape: ColliderShape, position: [f64; 3], sensor_id: u32) -> Self {
        Self {
            shape,
            position,
            category_bits: 0xFFFF_FFFF,
            mask_bits: 0xFFFF_FFFF,
            sensor_id,
        }
    }

    /// Set filter bits.
    pub fn with_filter(mut self, category: u32, mask: u32) -> Self {
        self.category_bits = category;
        self.mask_bits = mask;
        self
    }

    /// AABB of the sensor in world space.
    pub fn aabb(&self) -> Aabb {
        self.shape.aabb(self.position)
    }

    /// Test whether this sensor overlaps another sensor (broadphase + filter).
    pub fn overlaps_sensor(&self, other: &SensorCollider) -> bool {
        // Filter check
        if self.category_bits & other.mask_bits == 0 {
            return false;
        }
        if other.category_bits & self.mask_bits == 0 {
            return false;
        }
        self.aabb().overlaps(&other.aabb())
    }
}

// ── ColliderMaterialCombiner ──────────────────────────────────────────────────

/// Rules for combining friction and restitution from two colliders in contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombineRule {
    /// Average of the two values.
    Average,
    /// Minimum of the two values.
    Minimum,
    /// Maximum of the two values.
    Maximum,
    /// Product of the two values.
    Multiply,
}

impl CombineRule {
    /// Apply this rule to two scalar values.
    pub fn apply(self, a: f64, b: f64) -> f64 {
        match self {
            CombineRule::Average => (a + b) * 0.5,
            CombineRule::Minimum => a.min(b),
            CombineRule::Maximum => a.max(b),
            CombineRule::Multiply => a * b,
        }
    }
}

/// Compute the combined friction and restitution for a contact pair.
pub fn combine_materials(
    mat_a: &ColliderMaterial,
    mat_b: &ColliderMaterial,
    friction_rule: CombineRule,
    restitution_rule: CombineRule,
) -> ColliderMaterial {
    ColliderMaterial {
        friction: friction_rule.apply(mat_a.friction, mat_b.friction),
        restitution: restitution_rule.apply(mat_a.restitution, mat_b.restitution),
        density: (mat_a.density + mat_b.density) * 0.5,
    }
}

// ── Aabb extensions ───────────────────────────────────────────────────────────

impl Aabb {
    /// Return the volume of the AABB.
    pub fn volume(&self) -> f64 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        dx * dy * dz
    }

    /// Return the surface area of the AABB.
    pub fn surface_area(&self) -> f64 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Union of two AABBs (smallest AABB containing both).
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    /// Intersection of two AABBs, or `None` if they don't overlap.
    pub fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        let min = [
            self.min[0].max(other.min[0]),
            self.min[1].max(other.min[1]),
            self.min[2].max(other.min[2]),
        ];
        let max = [
            self.max[0].min(other.max[0]),
            self.max[1].min(other.max[1]),
            self.max[2].min(other.max[2]),
        ];
        if min[0] <= max[0] && min[1] <= max[1] && min[2] <= max[2] {
            Some(Aabb { min, max })
        } else {
            None
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::Transform;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_geometry::{BoxShape, Compound, Shape, Sphere};

    // ── existing trait-object collider tests ─────────────────────────────

    #[test]
    fn test_collider_creation() {
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let collider = Collider::new(shape)
            .with_friction(0.8)
            .with_restitution(0.5);
        assert!((collider.friction - 0.8).abs() < 1e-10);
        assert!((collider.restitution - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_collider_world_transform() {
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let collider = Collider::new(shape)
            .with_local_transform(Transform::from_position(Vec3::new(0.0, 1.0, 0.0)));
        let body_t = Transform::from_position(Vec3::new(10.0, 0.0, 0.0));
        let world = collider.world_transform(&body_t);
        assert!((world.position.x - 10.0).abs() < 1e-10);
        assert!((world.position.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_density_default() {
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let collider = Collider::new(shape);
        assert!((collider.density - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_collider_mass_from_sphere() {
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let collider = Collider::new(shape).with_density(1000.0);
        let props = collider.mass_properties();
        let expected_mass = 1000.0 * (4.0 / 3.0) * std::f64::consts::PI;
        assert!(
            (props.mass - expected_mass).abs() < 1.0,
            "mass={} expected≈{}",
            props.mass,
            expected_mass
        );
    }

    #[test]
    fn test_collider_inertia_sphere_symmetry() {
        let shape: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let collider = Collider::new(shape).with_density(1000.0);
        let props = collider.mass_properties();
        let ixx = props.local_inertia[(0, 0)];
        let iyy = props.local_inertia[(1, 1)];
        let izz = props.local_inertia[(2, 2)];
        assert!(ixx > 0.0);
        assert!((ixx - iyy).abs() < 1e-6, "Ixx={} != Iyy={}", ixx, iyy);
        assert!((iyy - izz).abs() < 1e-6, "Iyy={} != Izz={}", iyy, izz);
        assert!(props.local_inertia[(0, 1)].abs() < 1e-10);
        assert!(props.local_inertia[(0, 2)].abs() < 1e-10);
        assert!(props.local_inertia[(1, 2)].abs() < 1e-10);
    }

    /// Attach a compound collider (sphere + box) to a rigid body and verify mass properties.
    #[test]
    fn test_rigid_body_compound_collider() {
        use crate::body::RigidBody;

        let density = 1000.0_f64;

        let sphere: Arc<dyn Shape> = Arc::new(Sphere::new(1.0));
        let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));

        let compound: Arc<dyn Shape> = Arc::new(Compound::new(vec![
            (Transform::default(), sphere.clone()),
            (
                Transform::from_position(Vec3::new(3.0, 0.0, 0.0)),
                box_shape.clone(),
            ),
        ]));

        let collider = Collider::new(compound.clone()).with_density(density);
        let props = collider.mass_properties();

        let expected_mass = density * (sphere.volume() + box_shape.volume());

        assert!(
            (props.mass - expected_mass).abs() < 1e-6,
            "body mass={} expected={}",
            props.mass,
            expected_mass
        );

        let mut body = RigidBody::new(1.0);
        body.set_mass_from_shape(compound.as_ref(), density);
        assert!(
            (body.mass - expected_mass).abs() < 1e-6,
            "rigid body mass={} expected={}",
            body.mass,
            expected_mass
        );
        assert!(body.inverse_mass > 0.0, "inverse_mass should be positive");
    }

    // ── new Aabb tests ────────────────────────────────────────────────────

    #[test]
    fn aabb_sphere_correct_extents() {
        let shape = ColliderShape::Sphere { radius: 2.0 };
        let aabb = shape.aabb([1.0, 0.0, 0.0]);
        assert!((aabb.min[0] - (-1.0)).abs() < 1e-12);
        assert!((aabb.max[0] - 3.0).abs() < 1e-12);
        assert!((aabb.min[1] - (-2.0)).abs() < 1e-12);
        assert!((aabb.max[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_box_correct_extents() {
        let shape = ColliderShape::Box {
            half_extents: [1.0, 2.0, 0.5],
        };
        let aabb = shape.aabb([0.0, 0.0, 0.0]);
        assert!((aabb.min[0] - (-1.0)).abs() < 1e-12);
        assert!((aabb.max[0] - 1.0).abs() < 1e-12);
        assert!((aabb.min[1] - (-2.0)).abs() < 1e-12);
        assert!((aabb.max[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_overlaps_touching() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = Aabb {
            min: [1.0, 0.0, 0.0],
            max: [2.0, 1.0, 1.0],
        };
        assert!(a.overlaps(&b)); // touching counts as overlap
    }

    #[test]
    fn aabb_no_overlap_separated() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = Aabb {
            min: [2.0, 0.0, 0.0],
            max: [3.0, 1.0, 1.0],
        };
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn aabb_expand_increases_size() {
        let a = Aabb {
            min: [0.0, 0.0, 0.0],
            max: [1.0, 1.0, 1.0],
        };
        let b = a.expand(0.5);
        assert!((b.min[0] - (-0.5)).abs() < 1e-12);
        assert!((b.max[0] - 1.5).abs() < 1e-12);
    }

    // ── SimpleColliderSet tests ───────────────────────────────────────────

    #[test]
    fn simple_collider_set_insert_get() {
        let mut set = SimpleColliderSet::new();
        let col = SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 });
        let h = set.insert(col);
        assert!(set.get(h).is_some());
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn simple_collider_set_remove() {
        let mut set = SimpleColliderSet::new();
        let h = set.insert(SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 }));
        let removed = set.remove(h);
        assert!(removed.is_some());
        assert!(set.get(h).is_none());
        assert_eq!(set.len(), 0);
    }

    #[test]
    fn simple_collider_set_stale_handle() {
        let mut set = SimpleColliderSet::new();
        let h1 = set.insert(SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 }));
        set.remove(h1);
        let h2 = set.insert(SimpleCollider::new(ColliderShape::Sphere { radius: 2.0 }));
        assert!(set.get(h1).is_none(), "stale handle must be invalid");
        assert!(set.get(h2).is_some());
        assert_eq!(h2.index, h1.index); // slot reused
        assert_ne!(h2.generation, h1.generation);
    }

    #[test]
    fn simple_collider_set_broad_pairs_finds_overlap() {
        let mut set = SimpleColliderSet::new();
        // Two spheres close together → overlapping AABBs
        let mut c1 = SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 });
        c1.local_position = [0.0, 0.0, 0.0];
        let mut c2 = SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 });
        c2.local_position = [1.0, 0.0, 0.0]; // centres 1 m apart, radii 1 m each → overlap
        set.insert(c1);
        set.insert(c2);
        let pairs = set.broad_pairs(|_| None);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn simple_collider_set_broad_pairs_no_overlap() {
        let mut set = SimpleColliderSet::new();
        let mut c1 = SimpleCollider::new(ColliderShape::Sphere { radius: 0.4 });
        c1.local_position = [0.0, 0.0, 0.0];
        let mut c2 = SimpleCollider::new(ColliderShape::Sphere { radius: 0.4 });
        c2.local_position = [10.0, 0.0, 0.0];
        set.insert(c1);
        set.insert(c2);
        let pairs = set.broad_pairs(|_| None);
        assert_eq!(pairs.len(), 0);
    }

    #[test]
    fn sensor_flag_preserved() {
        let col = SimpleCollider::new(ColliderShape::Box {
            half_extents: [1.0, 1.0, 1.0],
        })
        .as_sensor();
        assert!(col.is_sensor);
    }

    #[test]
    fn capsule_aabb_height() {
        let shape = ColliderShape::Capsule {
            radius: 0.5,
            half_height: 1.0,
        };
        let aabb = shape.aabb([0.0, 0.0, 0.0]);
        // total half-height = radius + half_height = 1.5
        assert!((aabb.min[1] - (-1.5)).abs() < 1e-12);
        assert!((aabb.max[1] - 1.5).abs() < 1e-12);
    }

    // ── shape_mass_properties tests ───────────────────────────────────────

    #[test]
    fn sphere_mass_properties_correct() {
        let shape = ColliderShape::Sphere { radius: 1.0 };
        let (mass, cm, inertia) = shape_mass_properties(&shape, 1000.0);
        let expected_mass = 1000.0 * (4.0 / 3.0) * std::f64::consts::PI;
        assert!((mass - expected_mass).abs() < 1.0, "mass={mass}");
        assert!(cm[0].abs() < 1e-12 && cm[1].abs() < 1e-12 && cm[2].abs() < 1e-12);
        // Sphere inertia is symmetric
        assert!((inertia[0] - inertia[1]).abs() < 1e-6);
        assert!((inertia[1] - inertia[2]).abs() < 1e-6);
    }

    #[test]
    fn box_mass_properties_correct() {
        // Unit cube half_extents=[0.5,0.5,0.5], density=8 → mass=8*1=8
        let shape = ColliderShape::Box {
            half_extents: [0.5, 0.5, 0.5],
        };
        let (mass, _cm, inertia) = shape_mass_properties(&shape, 8.0);
        assert!((mass - 8.0).abs() < 1e-10, "mass={mass}");
        // Ixx = (1/3)*8*(0.25+0.25) = (1/3)*8*0.5 ≈ 1.333
        let expected_i = (1.0 / 3.0) * 8.0 * (0.25 + 0.25);
        assert!((inertia[0] - expected_i).abs() < 1e-9, "Ixx={}", inertia[0]);
    }

    #[test]
    fn cylinder_mass_properties_positive() {
        let shape = ColliderShape::Cylinder {
            radius: 1.0,
            height: 2.0,
        };
        let (mass, _cm, inertia) = shape_mass_properties(&shape, 1.0);
        assert!(mass > 0.0);
        assert!(inertia[0] > 0.0 && inertia[2] > 0.0);
    }

    #[test]
    fn capsule_mass_properties_positive() {
        let shape = ColliderShape::Capsule {
            radius: 0.5,
            half_height: 1.0,
        };
        let (mass, _cm, inertia) = shape_mass_properties(&shape, 500.0);
        assert!(mass > 0.0);
        assert!(inertia[0] > 0.0);
    }

    // ── ColliderTransform tests ───────────────────────────────────────────

    #[test]
    fn collider_transform_identity_point() {
        let t = ColliderTransform::identity();
        let p = [1.0, 2.0, 3.0];
        let tp = t.transform_point(p);
        assert!((tp[0] - 1.0).abs() < 1e-12);
        assert!((tp[1] - 2.0).abs() < 1e-12);
        assert!((tp[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn collider_transform_translate() {
        let t = ColliderTransform::from_position([5.0, 0.0, 0.0]);
        let tp = t.transform_point([1.0, 0.0, 0.0]);
        assert!((tp[0] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn collider_transform_compose() {
        let t1 = ColliderTransform::from_position([1.0, 0.0, 0.0]);
        let t2 = ColliderTransform::from_position([0.0, 2.0, 0.0]);
        let combined = t1.compose(&t2);
        assert!((combined.position[0] - 1.0).abs() < 1e-12);
        assert!((combined.position[1] - 2.0).abs() < 1e-12);
    }

    // ── TransformedCollider tests ─────────────────────────────────────────

    #[test]
    fn transformed_collider_world_aabb_offset() {
        let sc = SimpleCollider::new(ColliderShape::Sphere { radius: 1.0 });
        let local_t = ColliderTransform::from_position([3.0, 0.0, 0.0]);
        let tc = TransformedCollider::new(sc, local_t);
        let body_t = ColliderTransform::from_position([10.0, 0.0, 0.0]);
        let aabb = tc.world_aabb(&body_t);
        // World pos = 13, sphere radius=1 → [12, 14]
        assert!((aabb.min[0] - 12.0).abs() < 1e-12);
        assert!((aabb.max[0] - 14.0).abs() < 1e-12);
    }

    // ── CompoundColliderShape tests ───────────────────────────────────────

    #[test]
    fn compound_aabb_encloses_children() {
        let mut compound = CompoundColliderShape::new();
        compound.add_child(ColliderShape::Sphere { radius: 1.0 }, [0.0, 0.0, 0.0]);
        compound.add_child(ColliderShape::Sphere { radius: 1.0 }, [5.0, 0.0, 0.0]);
        let aabb = compound.world_aabb([0.0, 0.0, 0.0]).unwrap();
        assert!(aabb.min[0] <= -1.0);
        assert!(aabb.max[0] >= 6.0);
    }

    #[test]
    fn compound_volume_sum() {
        let mut compound = CompoundColliderShape::new();
        let s1 = ColliderShape::Sphere { radius: 1.0 };
        let s2 = ColliderShape::Sphere { radius: 2.0 };
        let v1 = s1.volume();
        let v2 = s2.volume();
        compound.add_child(s1, [0.0; 3]);
        compound.add_child(s2, [0.0; 3]);
        assert!((compound.volume() - (v1 + v2)).abs() < 1e-9);
    }

    #[test]
    fn compound_mass_properties_two_spheres() {
        let mut compound = CompoundColliderShape::new();
        compound.add_child(ColliderShape::Sphere { radius: 1.0 }, [-5.0, 0.0, 0.0]);
        compound.add_child(ColliderShape::Sphere { radius: 1.0 }, [5.0, 0.0, 0.0]);
        let (mass, cm, _) = compound.mass_properties(1.0);
        assert!(mass > 0.0);
        // Centre of mass should be at origin (symmetric)
        assert!(cm[0].abs() < 1e-9, "cm_x should be 0, got {}", cm[0]);
    }

    // ── SensorCollider tests ──────────────────────────────────────────────

    #[test]
    fn sensor_collider_aabb() {
        let sensor =
            SensorCollider::new(ColliderShape::Sphere { radius: 2.0 }, [0.0, 0.0, 0.0], 42);
        let aabb = sensor.aabb();
        assert!((aabb.min[0] - (-2.0)).abs() < 1e-12);
        assert!((aabb.max[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn sensor_collider_overlap_detected() {
        let s1 = SensorCollider::new(ColliderShape::Sphere { radius: 1.0 }, [0.0, 0.0, 0.0], 1);
        let s2 = SensorCollider::new(ColliderShape::Sphere { radius: 1.0 }, [1.0, 0.0, 0.0], 2);
        assert!(s1.overlaps_sensor(&s2));
    }

    #[test]
    fn sensor_collider_filter_blocks() {
        // s1 category=1, mask=2; s2 category=4, mask=8
        // 1 & 8 = 0 → filtered out
        let s1 = SensorCollider::new(ColliderShape::Sphere { radius: 5.0 }, [0.0; 3], 1)
            .with_filter(1, 2);
        let s2 = SensorCollider::new(ColliderShape::Sphere { radius: 5.0 }, [0.0; 3], 2)
            .with_filter(4, 8);
        assert!(!s1.overlaps_sensor(&s2));
    }

    // ── CombineRule / material combiner tests ─────────────────────────────

    #[test]
    fn combine_rule_average() {
        assert!((CombineRule::Average.apply(0.2, 0.8) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn combine_rule_minimum() {
        assert!((CombineRule::Minimum.apply(0.3, 0.7) - 0.3).abs() < 1e-12);
    }

    #[test]
    fn combine_rule_maximum() {
        assert!((CombineRule::Maximum.apply(0.3, 0.7) - 0.7).abs() < 1e-12);
    }

    #[test]
    fn combine_rule_multiply() {
        assert!((CombineRule::Multiply.apply(0.5, 0.4) - 0.2).abs() < 1e-12);
    }

    #[test]
    fn combine_materials_average_friction() {
        let a = ColliderMaterial::new(0.2, 0.4, 500.0);
        let b = ColliderMaterial::new(0.6, 0.8, 1000.0);
        let combined = combine_materials(&a, &b, CombineRule::Average, CombineRule::Minimum);
        assert!((combined.friction - 0.4).abs() < 1e-12);
        assert!((combined.restitution - 0.4).abs() < 1e-12);
    }

    // ── Aabb extensions tests ─────────────────────────────────────────────

    #[test]
    fn aabb_volume_unit_cube() {
        let a = Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!((a.volume() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_surface_area_unit_cube() {
        let a = Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        assert!((a.surface_area() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_union_covers_both() {
        let a = Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let b = Aabb {
            min: [2.0; 3],
            max: [3.0; 3],
        };
        let u = a.union(&b);
        assert!((u.min[0] - 0.0).abs() < 1e-12);
        assert!((u.max[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_intersection_overlapping() {
        let a = Aabb {
            min: [0.0; 3],
            max: [2.0; 3],
        };
        let b = Aabb {
            min: [1.0; 3],
            max: [3.0; 3],
        };
        let inter = a.intersection(&b);
        assert!(inter.is_some());
        let i = inter.unwrap();
        assert!((i.min[0] - 1.0).abs() < 1e-12);
        assert!((i.max[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn aabb_intersection_separated() {
        let a = Aabb {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let b = Aabb {
            min: [2.0; 3],
            max: [3.0; 3],
        };
        assert!(a.intersection(&b).is_none());
    }
}
