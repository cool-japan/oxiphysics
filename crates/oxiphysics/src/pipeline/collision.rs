// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Broad-phase and narrow-phase collision detection of the pipeline.
//!
//! Every step the pipeline snapshots the collider set into
//! [`ColliderEntry`]s, bounds each collider by its *world-space* AABB (the
//! local bounding box carried through the full collider transform, rotation
//! included), asks the broad phase for candidate pairs and runs the
//! [`NarrowPhaseDispatcher`](oxiphysics_collision::NarrowPhaseDispatcher) on
//! them. Pairs involving a sensor never reach the contact solver; their
//! overlaps are reported separately.

use std::time::{Duration, Instant};

use oxiphysics_collision::narrowphase::dispatch::ShapeType;
use oxiphysics_collision::{BroadPhase, CollisionPair, ContactManifold};
use oxiphysics_core::{Aabb, ColliderHandle, Transform};
use oxiphysics_geometry::{
    BoxShape, Capsule, Compound, Cone, ConvexHull, Cylinder, HeightField, Shape, Sphere,
    TriangleMesh,
};
use oxiphysics_rigid::{BodyState, BodyType, Collider, ColliderSet, RigidBody, RigidBodySet};

use super::PhysicsPipeline;

/// A collider of the current step together with what the pipeline needs to
/// know about it.
pub(super) struct ColliderEntry<'a> {
    /// Handle of the collider in its [`ColliderSet`].
    pub handle: ColliderHandle,
    /// The collider.
    pub collider: &'a Collider,
    /// Shape category used to select the narrow-phase algorithm.
    pub shape_type: ShapeType,
}

impl ColliderEntry<'_> {
    /// The body the collider is attached to, if it has one and it is alive.
    pub(super) fn body<'b>(&self, bodies: &'b RigidBodySet) -> Option<&'b RigidBody> {
        self.collider.body_handle.and_then(|h| bodies.get(h))
    }

    /// Transform of the collider in world space.
    pub(super) fn world_transform(&self, bodies: &RigidBodySet) -> Transform {
        let body_transform = self
            .body(bodies)
            .map(|b| b.transform.clone())
            .unwrap_or_default();
        self.collider.world_transform(&body_transform)
    }
}

/// Snapshot all colliders of `colliders` in iteration order.
///
/// The position of an entry in the returned vector is the collider index used
/// by the broad phase.
pub(super) fn collect_entries(colliders: &ColliderSet) -> Vec<ColliderEntry<'_>> {
    colliders
        .iter()
        .map(|(handle, collider)| ColliderEntry {
            handle,
            collider,
            shape_type: classify_shape(collider.shape.as_ref()),
        })
        .collect()
}

/// Shape category of `shape` for narrow-phase dispatch.
///
/// Shapes without a dedicated category are reported as [`ShapeType::Custom`]
/// and handled by the dispatcher's GJK/EPA fallback.
pub(super) fn classify_shape(shape: &dyn Shape) -> ShapeType {
    let any = shape.as_any();
    if any.is::<Sphere>() {
        ShapeType::Sphere
    } else if any.is::<BoxShape>() {
        ShapeType::Box
    } else if any.is::<Capsule>() {
        ShapeType::Capsule
    } else if any.is::<Cylinder>() {
        ShapeType::Cylinder
    } else if any.is::<Cone>() {
        ShapeType::Cone
    } else if any.is::<ConvexHull>() {
        ShapeType::ConvexHull
    } else if any.is::<TriangleMesh>() {
        ShapeType::TriangleMesh
    } else if any.is::<Compound>() {
        ShapeType::Compound
    } else if any.is::<HeightField>() {
        ShapeType::HeightField
    } else {
        ShapeType::Custom(0)
    }
}

/// `true` for a dynamic body that is not asleep.
fn is_awake_dynamic(body: Option<&RigidBody>) -> bool {
    body.is_some_and(|b| b.is_dynamic() && b.state != BodyState::Sleeping)
}

/// `true` for a dynamic body, awake or asleep.
fn is_dynamic(body: Option<&RigidBody>) -> bool {
    body.is_some_and(RigidBody::is_dynamic)
}

/// `true` for a kinematic body.
fn is_kinematic(body: Option<&RigidBody>) -> bool {
    body.is_some_and(|b| b.body_type == BodyType::Kinematic)
}

/// `true` if the body can move (dynamic or kinematic).
pub(super) fn can_move(body: Option<&RigidBody>) -> bool {
    body.is_some_and(|b| b.body_type != BodyType::Static)
}

/// Whether a solid pair has to be resolved this step.
///
/// A pair needs resolving if an awake dynamic body takes part in it, or if a
/// kinematic body touches a dynamic one (which must be woken up if it sleeps).
/// Pairs of static, kinematic and sleeping bodies can neither move nor be
/// moved, so they are skipped altogether.
pub(super) fn pair_needs_resolution(a: Option<&RigidBody>, b: Option<&RigidBody>) -> bool {
    is_awake_dynamic(a)
        || is_awake_dynamic(b)
        || (is_kinematic(a) && is_dynamic(b))
        || (is_kinematic(b) && is_dynamic(a))
}

/// `true` if both colliders are attached to the same body (they never collide
/// with each other).
pub(super) fn shares_body(a: &Collider, b: &Collider) -> bool {
    a.body_handle.is_some() && a.body_handle == b.body_handle
}

/// World-space AABB of a collider with the given world transform.
pub(super) fn world_aabb(collider: &Collider, world: &Transform) -> Aabb {
    collider.shape.bounding_box().transformed(world)
}

/// A solid contact between two colliders of the current step.
pub(super) struct ContactRecord {
    /// Index of collider A in the entry list.
    pub a: usize,
    /// Index of collider B in the entry list.
    pub b: usize,
    /// Contacts between the two; normals point from B toward A.
    pub manifold: ContactManifold,
}

/// Result of the discrete collision detection pass.
pub(super) struct Detection {
    /// Number of candidate pairs the broad phase reported.
    pub broadphase_pairs: usize,
    /// Solid contacts to feed to the solver.
    pub contacts: Vec<ContactRecord>,
    /// Pairs in which at least one collider is a sensor and the shapes overlap.
    pub sensor_overlaps: Vec<(ColliderHandle, ColliderHandle)>,
    /// Wall-clock time spent in the broad phase.
    pub broadphase_time: Duration,
    /// Wall-clock time spent in the narrow phase.
    pub narrowphase_time: Duration,
}

impl PhysicsPipeline {
    /// Broad phase and narrow phase at the current body positions.
    pub(super) fn detect_collisions(
        &self,
        entries: &[ColliderEntry<'_>],
        bodies: &RigidBodySet,
    ) -> Detection {
        let t_broadphase = Instant::now();
        let transforms: Vec<Transform> =
            entries.iter().map(|e| e.world_transform(bodies)).collect();
        let aabbs: Vec<Aabb> = entries
            .iter()
            .zip(&transforms)
            .map(|(entry, transform)| world_aabb(entry.collider, transform))
            .collect();
        let pairs = self.broadphase.find_pairs(&aabbs);
        let broadphase_time = t_broadphase.elapsed();

        let t_narrowphase = Instant::now();
        let mut contacts = Vec::new();
        let mut sensor_overlaps = Vec::new();
        for pair in &pairs {
            let (Some(entry_a), Some(entry_b)) = (entries.get(pair.a), entries.get(pair.b)) else {
                continue;
            };
            let (Some(transform_a), Some(transform_b)) =
                (transforms.get(pair.a), transforms.get(pair.b))
            else {
                continue;
            };
            if shares_body(entry_a.collider, entry_b.collider) {
                continue;
            }
            let body_a = entry_a.body(bodies);
            let body_b = entry_b.body(bodies);
            let is_sensor_pair = entry_a.collider.is_sensor || entry_b.collider.is_sensor;
            let relevant = if is_sensor_pair {
                can_move(body_a) || can_move(body_b)
            } else {
                pair_needs_resolution(body_a, body_b)
            };
            if !relevant {
                continue;
            }

            let result = self.narrowphase.dispatch(
                entry_a.collider.shape.as_ref(),
                entry_a.shape_type,
                transform_a,
                entry_b.collider.shape.as_ref(),
                entry_b.shape_type,
                transform_b,
                CollisionPair::new(pair.a, pair.b),
            );
            let Some(manifold) = result.manifold else {
                continue;
            };
            if manifold.contacts.is_empty() {
                continue;
            }
            if is_sensor_pair {
                sensor_overlaps.push((entry_a.handle, entry_b.handle));
            } else {
                contacts.push(ContactRecord {
                    a: pair.a,
                    b: pair.b,
                    manifold,
                });
            }
        }

        Detection {
            broadphase_pairs: pairs.len(),
            contacts,
            sensor_overlaps,
            broadphase_time,
            narrowphase_time: t_narrowphase.elapsed(),
        }
    }
}
