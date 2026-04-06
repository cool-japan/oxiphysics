// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core collision detection types.
//!
//! Provides the fundamental value types used throughout the collision pipeline:
//! contact points, contact manifolds, collision pairs, feature IDs, friction /
//! restitution material properties, collision filters, and island/group helpers.

use oxiphysics_core::math::{Real, Vec3};

// ── CollisionPair ─────────────────────────────────────────────────────────────

/// A pair of colliding object indices.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CollisionPair {
    /// Index of the first object.
    pub a: usize,
    /// Index of the second object.
    pub b: usize,
}

impl CollisionPair {
    /// Create a new collision pair.
    #[allow(dead_code)]
    pub fn new(a: usize, b: usize) -> Self {
        Self { a, b }
    }

    /// Return a canonicalized pair with `a <= b`.
    #[allow(dead_code)]
    pub fn canonical(a: usize, b: usize) -> Self {
        if a <= b {
            Self { a, b }
        } else {
            Self { a: b, b: a }
        }
    }

    /// Returns `true` if the given index is one of the pair's members.
    #[allow(dead_code)]
    pub fn contains(&self, idx: usize) -> bool {
        self.a == idx || self.b == idx
    }

    /// Returns the other member of the pair (panics if `idx` is not in the pair).
    #[allow(dead_code)]
    pub fn other(&self, idx: usize) -> usize {
        if self.a == idx {
            self.b
        } else if self.b == idx {
            self.a
        } else {
            panic!(
                "CollisionPair::other: index {idx} is not in pair ({}, {})",
                self.a, self.b
            )
        }
    }
}

// ── FeatureId ─────────────────────────────────────────────────────────────────

/// Identifies a geometric feature (vertex, edge, or face) on a shape.
///
/// Used by the contact manifold generator to determine which features are in
/// contact, enabling warm-starting and persistent manifold tracking.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FeatureId {
    /// A vertex feature with the given local index.
    Vertex(u32),
    /// An edge feature with the given local index.
    Edge(u32),
    /// A face feature with the given local index.
    Face(u32),
    /// Feature identity is unknown or not applicable.
    Unknown,
}

impl FeatureId {
    /// Returns `true` if this is a vertex feature.
    #[allow(dead_code)]
    pub fn is_vertex(self) -> bool {
        matches!(self, FeatureId::Vertex(_))
    }

    /// Returns `true` if this is an edge feature.
    #[allow(dead_code)]
    pub fn is_edge(self) -> bool {
        matches!(self, FeatureId::Edge(_))
    }

    /// Returns `true` if this is a face feature.
    #[allow(dead_code)]
    pub fn is_face(self) -> bool {
        matches!(self, FeatureId::Face(_))
    }

    /// Extract the raw index, or `None` if `Unknown`.
    #[allow(dead_code)]
    pub fn index(self) -> Option<u32> {
        match self {
            FeatureId::Vertex(i) | FeatureId::Edge(i) | FeatureId::Face(i) => Some(i),
            FeatureId::Unknown => None,
        }
    }
}

// ── Contact ───────────────────────────────────────────────────────────────────

/// A single contact point between two objects.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Contact {
    /// Contact point on body A in world space.
    pub point_a: Vec3,
    /// Contact point on body B in world space.
    pub point_b: Vec3,
    /// Contact normal (from B towards A).
    pub normal: Vec3,
    /// Penetration depth (positive means overlapping).
    pub depth: Real,
}

impl Contact {
    /// Create a new contact point.
    #[allow(dead_code)]
    pub fn new(point_a: Vec3, point_b: Vec3, normal: Vec3, depth: Real) -> Self {
        Self {
            point_a,
            point_b,
            normal,
            depth,
        }
    }

    /// The midpoint of the two contact points.
    #[allow(dead_code)]
    pub fn point(&self) -> Vec3 {
        (self.point_a + self.point_b) * 0.5
    }
}

/// An extended contact point that also carries feature IDs and solver impulses.
///
/// This is the richer variant used by the persistent manifold and warm-starting
/// system.  The simpler [`Contact`] is kept for backward compatibility.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RichContact {
    /// Contact point on body A in world space.
    pub point_a: Vec3,
    /// Contact point on body B in world space.
    pub point_b: Vec3,
    /// Contact normal (from B towards A).
    pub normal: Vec3,
    /// Penetration depth (positive means overlapping).
    pub depth: Real,
    /// Feature on shape A that generated this contact.
    pub feature_a: FeatureId,
    /// Feature on shape B that generated this contact.
    pub feature_b: FeatureId,
    /// Accumulated normal impulse (used for warm-starting the constraint solver).
    pub impulse_n: Real,
    /// Accumulated tangent impulse along the first friction direction.
    pub impulse_t1: Real,
    /// Accumulated tangent impulse along the second friction direction.
    pub impulse_t2: Real,
}

impl RichContact {
    /// Create a `RichContact` with feature and impulse fields zeroed.
    #[allow(dead_code)]
    pub fn new(point_a: Vec3, point_b: Vec3, normal: Vec3, depth: Real) -> Self {
        Self {
            point_a,
            point_b,
            normal,
            depth,
            feature_a: FeatureId::Unknown,
            feature_b: FeatureId::Unknown,
            impulse_n: 0.0,
            impulse_t1: 0.0,
            impulse_t2: 0.0,
        }
    }

    /// Create from a plain [`Contact`], zeroing feature IDs and impulses.
    #[allow(dead_code)]
    pub fn from_contact(c: &Contact) -> Self {
        Self::new(c.point_a, c.point_b, c.normal, c.depth)
    }

    /// Downgrade to a plain [`Contact`] (drops feature IDs and impulses).
    #[allow(dead_code)]
    pub fn to_contact(&self) -> Contact {
        Contact {
            point_a: self.point_a,
            point_b: self.point_b,
            normal: self.normal,
            depth: self.depth,
        }
    }

    /// The midpoint of the two contact points.
    #[allow(dead_code)]
    pub fn point(&self) -> Vec3 {
        (self.point_a + self.point_b) * 0.5
    }

    /// Build a pair of orthogonal tangent vectors from the contact normal.
    ///
    /// Returns `(t1, t2)` forming a right-handed frame with `self.normal`.
    #[allow(dead_code)]
    pub fn tangent_basis(&self) -> (Vec3, Vec3) {
        let n = self.normal;
        // Choose a reference vector not parallel to n
        let ref_vec = if n.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let t1 = n.cross(&ref_vec).normalize();
        let t2 = n.cross(&t1);
        (t1, t2)
    }

    /// Returns `true` if the contact is a penetrating contact (depth > 0).
    #[allow(dead_code)]
    pub fn is_penetrating(&self) -> bool {
        self.depth > 0.0
    }
}

// ── ContactManifold ───────────────────────────────────────────────────────────

/// Maximum number of contact points stored per manifold.
pub const MAX_CONTACTS: usize = 4;

/// A manifold holding multiple contact points.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// The collision pair.
    pub pair: CollisionPair,
    /// Contact points in this manifold.
    pub contacts: Vec<Contact>,
    /// Number of frames this manifold has been alive.
    pub age: u32,
}

impl ContactManifold {
    /// Create a new empty contact manifold for the given pair.
    #[allow(dead_code)]
    pub fn new(pair: CollisionPair) -> Self {
        Self {
            pair,
            contacts: Vec::new(),
            age: 0,
        }
    }

    /// Add a contact to the manifold, capping at [`MAX_CONTACTS`].
    ///
    /// When the manifold is full the contact with the shallowest penetration
    /// depth is replaced if the new contact is deeper.
    #[allow(dead_code)]
    pub fn add_contact(&mut self, contact: Contact) {
        if self.contacts.len() < MAX_CONTACTS {
            self.contacts.push(contact);
        } else {
            // Replace shallowest contact if new one is deeper
            if let Some(min_idx) = self
                .contacts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.depth
                        .partial_cmp(&b.depth)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                && contact.depth > self.contacts[min_idx].depth
            {
                self.contacts[min_idx] = contact;
            }
        }
    }

    /// Returns `true` if there are no contacts.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// Number of contact points.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    /// Maximum penetration depth across all contact points.
    #[allow(dead_code)]
    pub fn max_depth(&self) -> Real {
        self.contacts
            .iter()
            .map(|c| c.depth)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Average contact normal. Returns the zero vector if no contacts.
    #[allow(dead_code)]
    pub fn average_normal(&self) -> Vec3 {
        if self.contacts.is_empty() {
            return Vec3::zeros();
        }
        let sum: Vec3 = self
            .contacts
            .iter()
            .map(|c| c.normal)
            .fold(Vec3::zeros(), |a, b| a + b);
        let n = sum / (self.contacts.len() as Real);
        let norm = n.norm();
        if norm > 1e-12 {
            n / norm
        } else {
            Vec3::zeros()
        }
    }

    /// Remove contacts whose depth is below `threshold` (e.g., bodies separated).
    #[allow(dead_code)]
    pub fn prune_shallow(&mut self, threshold: Real) {
        self.contacts.retain(|c| c.depth >= threshold);
    }

    /// Increment the manifold age counter.
    #[allow(dead_code)]
    pub fn tick(&mut self) {
        self.age = self.age.saturating_add(1);
    }
}

// ── RichContactManifold ───────────────────────────────────────────────────────

/// A manifold holding [`RichContact`] points that carry feature IDs and
/// solver impulses.  Used by the warm-starting / persistent manifold system.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RichContactManifold {
    /// The collision pair.
    pub pair: CollisionPair,
    /// Contact points.
    pub contacts: Vec<RichContact>,
    /// Number of frames this manifold has been alive.
    pub age: u32,
}

impl RichContactManifold {
    /// Create a new empty manifold.
    #[allow(dead_code)]
    pub fn new(pair: CollisionPair) -> Self {
        Self {
            pair,
            contacts: Vec::new(),
            age: 0,
        }
    }

    /// Add a contact, capping at [`MAX_CONTACTS`].
    #[allow(dead_code)]
    pub fn add_contact(&mut self, contact: RichContact) {
        if self.contacts.len() < MAX_CONTACTS {
            self.contacts.push(contact);
        } else {
            if let Some(min_idx) = self
                .contacts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.depth
                        .partial_cmp(&b.depth)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                && contact.depth > self.contacts[min_idx].depth
            {
                self.contacts[min_idx] = contact;
            }
        }
    }

    /// Transfer accumulated impulses from `old` matching on feature IDs (warm-starting).
    #[allow(dead_code)]
    pub fn warm_start_from(&mut self, old: &RichContactManifold) {
        for c in &mut self.contacts {
            if let Some(prev) = old
                .contacts
                .iter()
                .find(|p| p.feature_a == c.feature_a && p.feature_b == c.feature_b)
            {
                c.impulse_n = prev.impulse_n;
                c.impulse_t1 = prev.impulse_t1;
                c.impulse_t2 = prev.impulse_t2;
            }
        }
    }

    /// Returns `true` if there are no contacts.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// Number of contact points.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    /// Maximum penetration depth.
    #[allow(dead_code)]
    pub fn max_depth(&self) -> Real {
        self.contacts
            .iter()
            .map(|c| c.depth)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Increment the age counter.
    #[allow(dead_code)]
    pub fn tick(&mut self) {
        self.age = self.age.saturating_add(1);
    }
}

// ── PhysicsMaterial ───────────────────────────────────────────────────────────

/// Physical surface properties used by the constraint solver.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PhysicsMaterial {
    /// Coefficient of restitution in \[0, 1\]; 0 = perfectly inelastic.
    pub restitution: Real,
    /// Static friction coefficient.
    pub friction_static: Real,
    /// Dynamic (kinetic) friction coefficient.
    pub friction_dynamic: Real,
}

impl Default for PhysicsMaterial {
    fn default() -> Self {
        Self {
            restitution: 0.3,
            friction_static: 0.6,
            friction_dynamic: 0.4,
        }
    }
}

impl PhysicsMaterial {
    /// Create a material with explicit parameters.
    #[allow(dead_code)]
    pub fn new(restitution: Real, friction_static: Real, friction_dynamic: Real) -> Self {
        Self {
            restitution,
            friction_static,
            friction_dynamic,
        }
    }

    /// Combine two materials using geometric-mean mixing for friction and
    /// the minimum for restitution (conservative).
    #[allow(dead_code)]
    pub fn combine(a: &Self, b: &Self) -> Self {
        Self {
            restitution: a.restitution.min(b.restitution),
            friction_static: (a.friction_static * b.friction_static).sqrt(),
            friction_dynamic: (a.friction_dynamic * b.friction_dynamic).sqrt(),
        }
    }
}

// ── CollisionFilter ───────────────────────────────────────────────────────────

/// Bitmask-based collision filter.
///
/// Body A collides with body B if and only if:
/// `A.mask & B.category != 0`.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CollisionFilter {
    /// Categories this body belongs to (bit flags).
    pub category: u32,
    /// Which categories this body should collide with (bit flags).
    pub mask: u32,
}

impl Default for CollisionFilter {
    fn default() -> Self {
        Self {
            category: 0xFFFF_FFFF,
            mask: 0xFFFF_FFFF,
        }
    }
}

impl CollisionFilter {
    /// Create a filter with explicit category and mask.
    #[allow(dead_code)]
    pub fn new(category: u32, mask: u32) -> Self {
        Self { category, mask }
    }

    /// Returns `true` if `self` and `other` should collide.
    #[allow(dead_code)]
    pub fn should_collide(&self, other: &CollisionFilter) -> bool {
        (self.mask & other.category) != 0 && (other.mask & self.category) != 0
    }

    /// A filter that collides with everything (default).
    #[allow(dead_code)]
    pub fn all() -> Self {
        Self::default()
    }

    /// A filter that collides with nothing.
    #[allow(dead_code)]
    pub fn none() -> Self {
        Self {
            category: 0,
            mask: 0,
        }
    }
}

// ── CollisionEvent ────────────────────────────────────────────────────────────

/// High-level event emitted by the collision pipeline for game logic.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum CollisionEvent {
    /// Two bodies began touching this frame.
    ContactStarted(CollisionPair),
    /// Two bodies that were touching are now separated.
    ContactEnded(CollisionPair),
    /// A trigger (sensor) was entered.
    TriggerEntered {
        /// Index of the sensor body.
        sensor: usize,
        /// Index of the body that entered the sensor.
        body: usize,
    },
    /// A trigger (sensor) was exited.
    TriggerExited {
        /// Index of the sensor body.
        sensor: usize,
        /// Index of the body that exited the sensor.
        body: usize,
    },
}

impl CollisionEvent {
    /// Return the pair of indices involved in the event.
    #[allow(dead_code)]
    pub fn pair(&self) -> (usize, usize) {
        match self {
            CollisionEvent::ContactStarted(p) | CollisionEvent::ContactEnded(p) => (p.a, p.b),
            CollisionEvent::TriggerEntered { sensor, body }
            | CollisionEvent::TriggerExited { sensor, body } => (*sensor, *body),
        }
    }

    /// Returns `true` if this event starts a contact or trigger.
    #[allow(dead_code)]
    pub fn is_enter(&self) -> bool {
        matches!(
            self,
            CollisionEvent::ContactStarted(_) | CollisionEvent::TriggerEntered { .. }
        )
    }

    /// Returns `true` if this event ends a contact or trigger.
    #[allow(dead_code)]
    pub fn is_exit(&self) -> bool {
        !self.is_enter()
    }
}

// ── IslandId ──────────────────────────────────────────────────────────────────

/// Identifier for a rigid-body island (connected component in the contact graph).
///
/// Bodies in the same island share constraints and are solved together.  The
/// sentinel value [`IslandId::NONE`] marks a body that has not been assigned to
/// any island.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct IslandId(pub u32);

impl IslandId {
    /// Sentinel: body has not been assigned to any island.
    pub const NONE: Self = IslandId(u32::MAX);

    /// Returns `true` if this is the sentinel value.
    #[allow(dead_code)]
    pub fn is_none(self) -> bool {
        self == Self::NONE
    }

    /// Returns `true` if this is a valid island id.
    #[allow(dead_code)]
    pub fn is_some(self) -> bool {
        self != Self::NONE
    }

    /// Raw integer index.
    #[allow(dead_code)]
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl Default for IslandId {
    fn default() -> Self {
        Self::NONE
    }
}

impl std::fmt::Display for IslandId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_none() {
            write!(f, "IslandId::NONE")
        } else {
            write!(f, "IslandId({})", self.0)
        }
    }
}

// ── RigidBodyKind ─────────────────────────────────────────────────────────────

/// Classifies how a rigid body participates in dynamics.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RigidBodyKind {
    /// Fully simulated: responds to forces, collisions, and constraints.
    Dynamic,
    /// Has infinite mass; cannot be moved by forces or impulses.
    Static,
    /// Kinematically driven: moves under scripted velocity but pushes dynamic bodies.
    Kinematic,
    /// Sensor (trigger): generates overlap events but no contact impulses.
    Sensor,
}

impl RigidBodyKind {
    /// Returns `true` if the body can be moved by the solver.
    #[allow(dead_code)]
    pub fn is_movable(self) -> bool {
        matches!(self, RigidBodyKind::Dynamic | RigidBodyKind::Kinematic)
    }

    /// Returns `true` if the body generates contact forces.
    #[allow(dead_code)]
    pub fn generates_contacts(self) -> bool {
        !matches!(self, RigidBodyKind::Sensor)
    }

    /// Returns `true` if the body participates in island building.
    #[allow(dead_code)]
    pub fn is_island_member(self) -> bool {
        matches!(self, RigidBodyKind::Dynamic)
    }
}

// ── BodyGroup ─────────────────────────────────────────────────────────────────

/// A strongly-connected group of rigid bodies linked by persistent contacts.
///
/// Groups are computed from the contact graph each frame to allow the solver to
/// partition work and detect sleeping candidates.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BodyGroup {
    /// Island this group belongs to (may span multiple groups).
    pub island: IslandId,
    /// Indices of bodies in this group.
    pub bodies: Vec<usize>,
    /// Whether every body in the group is at rest (can be put to sleep).
    pub all_sleeping: bool,
}

impl BodyGroup {
    /// Create a new group from a list of body indices.
    #[allow(dead_code)]
    pub fn new(island: IslandId, bodies: Vec<usize>) -> Self {
        Self {
            island,
            bodies,
            all_sleeping: false,
        }
    }

    /// Number of bodies in the group.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    /// Returns `true` if the group has no members.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    /// Returns `true` if `body_idx` is a member of this group.
    #[allow(dead_code)]
    pub fn contains(&self, body_idx: usize) -> bool {
        self.bodies.contains(&body_idx)
    }

    /// Mark all bodies as sleeping/not-sleeping.
    #[allow(dead_code)]
    pub fn set_sleeping(&mut self, sleeping: bool) {
        self.all_sleeping = sleeping;
    }
}

// ── ContactConstraint ─────────────────────────────────────────────────────────

/// A position-level contact constraint ready to be fed to a constraint solver.
///
/// Carries the pair, penetration depth, normal, and warm-start impulse.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactConstraint {
    /// The colliding pair.
    pub pair: CollisionPair,
    /// Contact normal (from B toward A, unit length).
    pub normal: Vec3,
    /// Penetration depth (positive = overlapping).
    pub depth: Real,
    /// World-space contact point (midpoint of the two witness points).
    pub point: Vec3,
    /// Accumulated normal impulse from the previous frame (warm-start).
    pub impulse_n: Real,
    /// Accumulated tangent impulse 1 from the previous frame.
    pub impulse_t1: Real,
    /// Accumulated tangent impulse 2 from the previous frame.
    pub impulse_t2: Real,
    /// Effective restitution coefficient for this contact.
    pub restitution: Real,
    /// Effective friction coefficient for this contact.
    pub friction: Real,
}

impl ContactConstraint {
    /// Create a constraint with zero warm-start impulses and default material.
    #[allow(dead_code)]
    pub fn new(pair: CollisionPair, normal: Vec3, depth: Real, point: Vec3) -> Self {
        Self {
            pair,
            normal,
            depth,
            point,
            impulse_n: 0.0,
            impulse_t1: 0.0,
            impulse_t2: 0.0,
            restitution: 0.3,
            friction: 0.5,
        }
    }

    /// Build a pair of tangent vectors for friction from the contact normal.
    #[allow(dead_code)]
    pub fn tangent_basis(&self) -> (Vec3, Vec3) {
        let n = self.normal;
        let ref_vec = if n.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let t1 = n.cross(&ref_vec).normalize();
        let t2 = n.cross(&t1);
        (t1, t2)
    }

    /// Returns `true` if the contact is currently penetrating.
    #[allow(dead_code)]
    pub fn is_penetrating(&self) -> bool {
        self.depth > 0.0
    }

    /// Zero all accumulated impulses (invalidates warm-start data).
    #[allow(dead_code)]
    pub fn clear_impulses(&mut self) {
        self.impulse_n = 0.0;
        self.impulse_t1 = 0.0;
        self.impulse_t2 = 0.0;
    }
}

// ── ContactVelocityConstraint ─────────────────────────────────────────────────

/// A velocity-level (impulse-based) contact constraint.
///
/// Pre-computed quantities that remain constant during the solver iteration.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactVelocityConstraint {
    /// The colliding pair.
    pub pair: CollisionPair,
    /// Contact normal.
    pub normal: Vec3,
    /// Target relative velocity along the normal (bias + restitution).
    pub velocity_bias: Real,
    /// Effective mass along the normal direction.
    pub effective_mass_n: Real,
    /// Effective mass along tangent 1.
    pub effective_mass_t1: Real,
    /// Effective mass along tangent 2.
    pub effective_mass_t2: Real,
    /// Friction coefficient.
    pub friction: Real,
    /// Accumulated normal impulse (clamped to >= 0 during solve).
    pub impulse_n: Real,
    /// Accumulated tangent impulse 1.
    pub impulse_t1: Real,
    /// Accumulated tangent impulse 2.
    pub impulse_t2: Real,
}

impl ContactVelocityConstraint {
    /// Create a velocity constraint with zeroed impulses.
    #[allow(dead_code)]
    pub fn new(pair: CollisionPair, normal: Vec3) -> Self {
        Self {
            pair,
            normal,
            velocity_bias: 0.0,
            effective_mass_n: 1.0,
            effective_mass_t1: 1.0,
            effective_mass_t2: 1.0,
            friction: 0.5,
            impulse_n: 0.0,
            impulse_t1: 0.0,
            impulse_t2: 0.0,
        }
    }

    /// Returns `true` if there is any accumulated impulse.
    #[allow(dead_code)]
    pub fn has_impulse(&self) -> bool {
        self.impulse_n.abs() > 1e-12
            || self.impulse_t1.abs() > 1e-12
            || self.impulse_t2.abs() > 1e-12
    }

    /// The total impulse magnitude (for debugging/telemetry).
    #[allow(dead_code)]
    pub fn impulse_magnitude(&self) -> Real {
        (self.impulse_n * self.impulse_n
            + self.impulse_t1 * self.impulse_t1
            + self.impulse_t2 * self.impulse_t2)
            .sqrt()
    }

    /// Warm-start: copy impulses from a previous-frame constraint.
    #[allow(dead_code)]
    pub fn warm_start_from(&mut self, prev: &ContactVelocityConstraint) {
        self.impulse_n = prev.impulse_n;
        self.impulse_t1 = prev.impulse_t1;
        self.impulse_t2 = prev.impulse_t2;
    }
}

// ── ContactBatch ──────────────────────────────────────────────────────────────

/// A batch of velocity constraints to be solved together (same island).
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct ContactBatch {
    /// All velocity constraints in the batch.
    pub constraints: Vec<ContactVelocityConstraint>,
    /// The island this batch belongs to.
    pub island: IslandId,
}

impl ContactBatch {
    /// Create an empty batch for the given island.
    #[allow(dead_code)]
    pub fn new(island: IslandId) -> Self {
        Self {
            constraints: Vec::new(),
            island,
        }
    }

    /// Add a constraint to the batch.
    #[allow(dead_code)]
    pub fn push(&mut self, c: ContactVelocityConstraint) {
        self.constraints.push(c);
    }

    /// Number of constraints.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.constraints.len()
    }

    /// Returns `true` if no constraints are present.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }

    /// Total accumulated normal impulse across all constraints.
    #[allow(dead_code)]
    pub fn total_normal_impulse(&self) -> Real {
        self.constraints.iter().map(|c| c.impulse_n).sum()
    }

    /// Clear all impulses (invalidates warm-start data for this batch).
    #[allow(dead_code)]
    pub fn clear_impulses(&mut self) {
        for c in &mut self.constraints {
            c.impulse_n = 0.0;
            c.impulse_t1 = 0.0;
            c.impulse_t2 = 0.0;
        }
    }
}

// ── CollisionMatrix ───────────────────────────────────────────────────────────

/// Dense bitmatrix of which category bits collide with which others.
///
/// Up to 32 categories are supported.  The matrix is symmetric: if `(i, j)` is
/// set then `(j, i)` is also set.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CollisionMatrix {
    pub(super) rows: [u32; 32],
}

impl Default for CollisionMatrix {
    fn default() -> Self {
        Self::all_collide()
    }
}

impl CollisionMatrix {
    /// Create a matrix where every category collides with every other.
    #[allow(dead_code)]
    pub fn all_collide() -> Self {
        Self {
            rows: [u32::MAX; 32],
        }
    }

    /// Create a matrix where no categories collide.
    #[allow(dead_code)]
    pub fn none_collide() -> Self {
        Self { rows: [0u32; 32] }
    }

    /// Enable collision between categories `a` and `b` (symmetric).
    #[allow(dead_code)]
    pub fn enable(&mut self, a: u8, b: u8) {
        let a = (a & 31) as usize;
        let b = (b & 31) as usize;
        self.rows[a] |= 1 << b;
        self.rows[b] |= 1 << a;
    }

    /// Disable collision between categories `a` and `b` (symmetric).
    #[allow(dead_code)]
    pub fn disable(&mut self, a: u8, b: u8) {
        let a = (a & 31) as usize;
        let b = (b & 31) as usize;
        self.rows[a] &= !(1 << b);
        self.rows[b] &= !(1 << a);
    }

    /// Returns `true` if categories `a` and `b` should collide.
    #[allow(dead_code)]
    pub fn should_collide(&self, a: u8, b: u8) -> bool {
        let a = (a & 31) as usize;
        let b = (b & 31) as usize;
        (self.rows[a] >> b) & 1 == 1
    }
}

// ── ContactStats ──────────────────────────────────────────────────────────────

/// Per-frame contact statistics collected by the pipeline.
#[allow(dead_code)]
#[derive(Debug, Default, Clone, Copy)]
pub struct ContactStats {
    /// Total number of broadphase pairs tested.
    pub broadphase_pairs: u32,
    /// Number of pairs that passed to narrowphase.
    pub narrowphase_pairs: u32,
    /// Number of contacts generated by narrowphase.
    pub contacts_generated: u32,
    /// Number of contacts removed due to manifold aging.
    pub contacts_expired: u32,
    /// Number of islands built.
    pub islands: u32,
}

impl ContactStats {
    /// Create a zeroed stats block.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the ratio of narrowphase work to broadphase pairs.
    ///
    /// A value close to 1 means the broadphase is poor at rejection.
    #[allow(dead_code)]
    pub fn narrowphase_ratio(&self) -> f64 {
        if self.broadphase_pairs == 0 {
            return 0.0;
        }
        self.narrowphase_pairs as f64 / self.broadphase_pairs as f64
    }

    /// Reset all counters.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collision_pair_basic() {
        let pair = CollisionPair::new(0, 1);
        assert_eq!(pair.a, 0);
        assert_eq!(pair.b, 1);
    }

    #[test]
    fn test_collision_pair_canonical() {
        let p1 = CollisionPair::canonical(3, 1);
        assert!(p1.a <= p1.b);
        let p2 = CollisionPair::canonical(1, 3);
        assert_eq!(p1, p2);
    }

    #[test]
    fn test_collision_pair_contains() {
        let p = CollisionPair::new(2, 7);
        assert!(p.contains(2));
        assert!(p.contains(7));
        assert!(!p.contains(5));
    }

    #[test]
    fn test_collision_pair_other() {
        let p = CollisionPair::new(2, 7);
        assert_eq!(p.other(2), 7);
        assert_eq!(p.other(7), 2);
    }

    #[test]
    #[should_panic]
    fn test_collision_pair_other_panics() {
        let p = CollisionPair::new(2, 7);
        let _ = p.other(99);
    }

    #[test]
    fn test_feature_id() {
        assert!(FeatureId::Vertex(0).is_vertex());
        assert!(FeatureId::Edge(2).is_edge());
        assert!(FeatureId::Face(5).is_face());
        assert_eq!(FeatureId::Face(5).index(), Some(5));
        assert_eq!(FeatureId::Unknown.index(), None);
    }

    #[test]
    fn test_contact_manifold_basic() {
        let pair = CollisionPair::new(0, 1);
        let mut manifold = ContactManifold::new(pair);
        manifold.add_contact(Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::new(0.0, -0.1, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.1,
        });
        assert_eq!(manifold.contacts.len(), 1);
    }

    fn make_contact(depth: f64) -> Contact {
        Contact {
            point_a: Vec3::zeros(),
            point_b: Vec3::zeros(),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth,
        }
    }

    #[test]
    fn test_manifold_max_contacts_cap() {
        let pair = CollisionPair::new(0, 1);
        let mut manifold = ContactManifold::new(pair);
        for i in 0..8 {
            manifold.add_contact(make_contact(i as f64 * 0.01));
        }
        assert!(manifold.contacts.len() <= MAX_CONTACTS);
    }

    #[test]
    fn test_manifold_max_depth() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(make_contact(0.1));
        m.add_contact(make_contact(0.5));
        m.add_contact(make_contact(0.3));
        assert!((m.max_depth() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_prune_shallow() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(make_contact(0.1));
        m.add_contact(make_contact(-0.05));
        m.prune_shallow(0.0);
        assert_eq!(m.contacts.len(), 1);
        assert!((m.contacts[0].depth - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_average_normal() {
        let pair = CollisionPair::new(0, 1);
        let mut m = ContactManifold::new(pair);
        m.add_contact(make_contact(0.1));
        m.add_contact(make_contact(0.2));
        let avg = m.average_normal();
        assert!((avg.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_manifold_empty_average_normal() {
        let m = ContactManifold::new(CollisionPair::new(0, 1));
        let avg = m.average_normal();
        assert_eq!(avg.norm(), 0.0);
    }

    #[test]
    fn test_rich_contact_tangent_basis_orthogonal() {
        let c = RichContact::new(Vec3::zeros(), Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 0.1);
        let (t1, t2) = c.tangent_basis();
        assert!(c.normal.dot(&t1).abs() < 1e-10, "t1 not perp to normal");
        assert!(c.normal.dot(&t2).abs() < 1e-10, "t2 not perp to normal");
        assert!(t1.dot(&t2).abs() < 1e-10, "t1 not perp to t2");
    }

    #[test]
    fn test_physics_material_combine() {
        let a = PhysicsMaterial::new(0.8, 0.4, 0.3);
        let b = PhysicsMaterial::new(0.5, 0.9, 0.6);
        let c = PhysicsMaterial::combine(&a, &b);
        assert!((c.restitution - 0.5).abs() < 1e-12, "min restitution");
        let expected_fs = (0.4_f64 * 0.9).sqrt();
        assert!((c.friction_static - expected_fs).abs() < 1e-12);
    }

    #[test]
    fn test_collision_filter_default_collides() {
        let f1 = CollisionFilter::default();
        let f2 = CollisionFilter::default();
        assert!(f1.should_collide(&f2));
    }

    #[test]
    fn test_collision_filter_none_no_collide() {
        let f1 = CollisionFilter::none();
        let f2 = CollisionFilter::default();
        assert!(!f1.should_collide(&f2));
    }

    #[test]
    fn test_collision_filter_category_mask() {
        // player (cat 1) collides with wall (cat 2) if its mask includes cat 2
        let player = CollisionFilter::new(0b0001, 0b0010);
        let wall = CollisionFilter::new(0b0010, 0b0001);
        assert!(player.should_collide(&wall));
        // two players: cat 1, mask 0b0010 → they should NOT collide with each other
        let player2 = CollisionFilter::new(0b0001, 0b0010);
        assert!(!player.should_collide(&player2)); // mask 0b0010 & cat 0b0001 == 0
    }

    #[test]
    fn test_collision_event_pair() {
        let ev = CollisionEvent::ContactStarted(CollisionPair::new(3, 7));
        assert_eq!(ev.pair(), (3, 7));
        assert!(ev.is_enter());
        assert!(!ev.is_exit());
    }

    #[test]
    fn test_collision_event_trigger() {
        let ev = CollisionEvent::TriggerEntered {
            sensor: 10,
            body: 2,
        };
        assert_eq!(ev.pair(), (10, 2));
        assert!(ev.is_enter());
        let ev2 = CollisionEvent::TriggerExited {
            sensor: 10,
            body: 2,
        };
        assert!(ev2.is_exit());
    }

    #[test]
    fn test_rich_manifold_warm_start_from() {
        let pair = CollisionPair::new(0, 1);
        let mut old = RichContactManifold::new(pair);
        let mut c = RichContact::new(Vec3::zeros(), Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 0.1);
        c.feature_a = FeatureId::Face(0);
        c.feature_b = FeatureId::Face(1);
        c.impulse_n = 5.0;
        old.add_contact(c);

        let mut new_m = RichContactManifold::new(pair);
        let mut c2 = RichContact::new(Vec3::zeros(), Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 0.12);
        c2.feature_a = FeatureId::Face(0);
        c2.feature_b = FeatureId::Face(1);
        new_m.add_contact(c2);

        new_m.warm_start_from(&old);
        assert!((new_m.contacts[0].impulse_n - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_tick_age() {
        let mut m = ContactManifold::new(CollisionPair::new(0, 1));
        assert_eq!(m.age, 0);
        m.tick();
        m.tick();
        assert_eq!(m.age, 2);
    }

    #[test]
    fn test_rich_contact_from_contact_roundtrip() {
        let c = Contact {
            point_a: Vec3::new(1.0, 0.0, 0.0),
            point_b: Vec3::new(2.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.3,
        };
        let rc = RichContact::from_contact(&c);
        assert_eq!(rc.depth, 0.3);
        assert_eq!(rc.feature_a, FeatureId::Unknown);
        let back = rc.to_contact();
        assert!((back.depth - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_rich_contact_is_penetrating() {
        let c_pen = RichContact::new(Vec3::zeros(), Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), 0.1);
        let c_sep = RichContact::new(Vec3::zeros(), Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0), -0.1);
        assert!(c_pen.is_penetrating());
        assert!(!c_sep.is_penetrating());
    }

    // ── IslandId tests ───────────────────────────────────────────────────────

    #[test]
    fn island_id_none_is_sentinel() {
        let id = IslandId::NONE;
        assert!(id.is_none());
        assert!(!id.is_some());
    }

    #[test]
    fn island_id_valid() {
        let id = IslandId(3);
        assert!(!id.is_none());
        assert!(id.is_some());
        assert_eq!(id.index(), 3);
    }

    #[test]
    fn island_id_default_is_none() {
        let id = IslandId::default();
        assert!(id.is_none());
    }

    #[test]
    fn island_id_ordering() {
        assert!(IslandId(0) < IslandId(1));
        assert!(IslandId(0) < IslandId::NONE);
    }

    #[test]
    fn island_id_display_none() {
        let s = format!("{}", IslandId::NONE);
        assert!(s.contains("NONE"), "display: {s}");
    }

    #[test]
    fn island_id_display_value() {
        let s = format!("{}", IslandId(42));
        assert!(s.contains("42"), "display: {s}");
    }

    // ── RigidBodyKind tests ──────────────────────────────────────────────────

    #[test]
    fn rigid_body_kind_movable() {
        assert!(RigidBodyKind::Dynamic.is_movable());
        assert!(RigidBodyKind::Kinematic.is_movable());
        assert!(!RigidBodyKind::Static.is_movable());
        assert!(!RigidBodyKind::Sensor.is_movable());
    }

    #[test]
    fn rigid_body_kind_generates_contacts() {
        assert!(RigidBodyKind::Dynamic.generates_contacts());
        assert!(RigidBodyKind::Static.generates_contacts());
        assert!(RigidBodyKind::Kinematic.generates_contacts());
        assert!(!RigidBodyKind::Sensor.generates_contacts());
    }

    #[test]
    fn rigid_body_kind_island_member() {
        assert!(RigidBodyKind::Dynamic.is_island_member());
        assert!(!RigidBodyKind::Static.is_island_member());
        assert!(!RigidBodyKind::Kinematic.is_island_member());
        assert!(!RigidBodyKind::Sensor.is_island_member());
    }

    // ── BodyGroup tests ──────────────────────────────────────────────────────

    #[test]
    fn body_group_basic() {
        let g = BodyGroup::new(IslandId(0), vec![1, 2, 3]);
        assert_eq!(g.len(), 3);
        assert!(g.contains(2));
        assert!(!g.contains(99));
        assert!(!g.is_empty());
    }

    #[test]
    fn body_group_empty() {
        let g = BodyGroup::new(IslandId(0), vec![]);
        assert!(g.is_empty());
        assert_eq!(g.len(), 0);
    }

    #[test]
    fn body_group_set_sleeping() {
        let mut g = BodyGroup::new(IslandId(0), vec![1]);
        assert!(!g.all_sleeping);
        g.set_sleeping(true);
        assert!(g.all_sleeping);
        g.set_sleeping(false);
        assert!(!g.all_sleeping);
    }

    // ── ContactConstraint tests ──────────────────────────────────────────────

    #[test]
    fn contact_constraint_new() {
        let pair = CollisionPair::new(0, 1);
        let c = ContactConstraint::new(pair, Vec3::new(0.0, 1.0, 0.0), 0.05, Vec3::zeros());
        assert!(c.is_penetrating());
        assert_eq!(c.impulse_n, 0.0);
        assert_eq!(c.impulse_t1, 0.0);
        assert_eq!(c.impulse_t2, 0.0);
    }

    #[test]
    fn contact_constraint_not_penetrating() {
        let pair = CollisionPair::new(0, 1);
        let c = ContactConstraint::new(pair, Vec3::new(0.0, 1.0, 0.0), -0.01, Vec3::zeros());
        assert!(!c.is_penetrating());
    }

    #[test]
    fn contact_constraint_tangent_basis_orthogonal() {
        let pair = CollisionPair::new(0, 1);
        let c = ContactConstraint::new(pair, Vec3::new(0.0, 1.0, 0.0), 0.1, Vec3::zeros());
        let (t1, t2) = c.tangent_basis();
        assert!(c.normal.dot(&t1).abs() < 1e-10, "t1 not perp to normal");
        assert!(c.normal.dot(&t2).abs() < 1e-10, "t2 not perp to normal");
        assert!(t1.dot(&t2).abs() < 1e-10, "t1 not perp to t2");
    }

    #[test]
    fn contact_constraint_clear_impulses() {
        let pair = CollisionPair::new(0, 1);
        let mut c = ContactConstraint::new(pair, Vec3::new(0.0, 1.0, 0.0), 0.1, Vec3::zeros());
        c.impulse_n = 5.0;
        c.impulse_t1 = 2.0;
        c.impulse_t2 = -1.0;
        c.clear_impulses();
        assert_eq!(c.impulse_n, 0.0);
        assert_eq!(c.impulse_t1, 0.0);
        assert_eq!(c.impulse_t2, 0.0);
    }

    // ── ContactVelocityConstraint tests ──────────────────────────────────────

    #[test]
    fn velocity_constraint_new_zeroed() {
        let c = ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        assert!(!c.has_impulse());
        assert_eq!(c.impulse_magnitude(), 0.0);
    }

    #[test]
    fn velocity_constraint_has_impulse_after_set() {
        let mut c =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        c.impulse_n = 3.0;
        assert!(c.has_impulse());
    }

    #[test]
    fn velocity_constraint_impulse_magnitude() {
        let mut c =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        c.impulse_n = 3.0;
        c.impulse_t1 = 4.0;
        let mag = c.impulse_magnitude();
        assert!((mag - 5.0).abs() < 1e-10, "magnitude={mag}");
    }

    #[test]
    fn velocity_constraint_warm_start_from() {
        let mut prev =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        prev.impulse_n = 7.0;
        prev.impulse_t1 = -2.0;
        prev.impulse_t2 = 1.5;

        let mut curr =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        curr.warm_start_from(&prev);
        assert!((curr.impulse_n - 7.0).abs() < 1e-12);
        assert!((curr.impulse_t1 - (-2.0)).abs() < 1e-12);
        assert!((curr.impulse_t2 - 1.5).abs() < 1e-12);
    }

    // ── ContactBatch tests ────────────────────────────────────────────────────

    #[test]
    fn contact_batch_empty() {
        let batch = ContactBatch::new(IslandId(0));
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
        assert_eq!(batch.total_normal_impulse(), 0.0);
    }

    #[test]
    fn contact_batch_push_and_len() {
        let mut batch = ContactBatch::new(IslandId(1));
        let c = ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        batch.push(c);
        assert_eq!(batch.len(), 1);
        assert!(!batch.is_empty());
    }

    #[test]
    fn contact_batch_total_normal_impulse() {
        let mut batch = ContactBatch::new(IslandId(0));
        let mut c1 =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        c1.impulse_n = 3.0;
        let mut c2 =
            ContactVelocityConstraint::new(CollisionPair::new(1, 2), Vec3::new(0.0, 1.0, 0.0));
        c2.impulse_n = 5.0;
        batch.push(c1);
        batch.push(c2);
        assert!((batch.total_normal_impulse() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn contact_batch_clear_impulses() {
        let mut batch = ContactBatch::new(IslandId(0));
        let mut c =
            ContactVelocityConstraint::new(CollisionPair::new(0, 1), Vec3::new(0.0, 1.0, 0.0));
        c.impulse_n = 9.0;
        batch.push(c);
        batch.clear_impulses();
        assert_eq!(batch.total_normal_impulse(), 0.0);
    }

    // ── CollisionMatrix tests ─────────────────────────────────────────────────

    #[test]
    fn collision_matrix_all_collide() {
        let m = CollisionMatrix::all_collide();
        for a in 0..32u8 {
            for b in 0..32u8 {
                assert!(
                    m.should_collide(a, b),
                    "all_collide: ({a},{b}) should collide"
                );
            }
        }
    }

    #[test]
    fn collision_matrix_none_collide() {
        let m = CollisionMatrix::none_collide();
        for a in 0..32u8 {
            for b in 0..32u8 {
                assert!(
                    !m.should_collide(a, b),
                    "none_collide: ({a},{b}) should not collide"
                );
            }
        }
    }

    #[test]
    fn collision_matrix_enable_symmetric() {
        let mut m = CollisionMatrix::none_collide();
        m.enable(2, 5);
        assert!(m.should_collide(2, 5));
        assert!(m.should_collide(5, 2), "must be symmetric");
        assert!(!m.should_collide(2, 3));
    }

    #[test]
    fn collision_matrix_disable_symmetric() {
        let mut m = CollisionMatrix::all_collide();
        m.disable(3, 7);
        assert!(!m.should_collide(3, 7));
        assert!(!m.should_collide(7, 3), "must be symmetric");
        assert!(m.should_collide(3, 8), "other pairs unaffected");
    }

    #[test]
    fn collision_matrix_self_collision() {
        let mut m = CollisionMatrix::none_collide();
        m.enable(4, 4);
        assert!(m.should_collide(4, 4));
    }

    // ── ContactStats tests ────────────────────────────────────────────────────

    #[test]
    fn contact_stats_default_zeroed() {
        let s = ContactStats::new();
        assert_eq!(s.broadphase_pairs, 0);
        assert_eq!(s.narrowphase_pairs, 0);
        assert_eq!(s.contacts_generated, 0);
    }

    #[test]
    fn contact_stats_narrowphase_ratio_zero_broadphase() {
        let s = ContactStats::new();
        assert_eq!(s.narrowphase_ratio(), 0.0);
    }

    #[test]
    fn contact_stats_narrowphase_ratio() {
        let s = ContactStats {
            broadphase_pairs: 100,
            narrowphase_pairs: 25,
            ..ContactStats::default()
        };
        assert!((s.narrowphase_ratio() - 0.25).abs() < 1e-12);
    }

    #[test]
    fn contact_stats_reset() {
        let mut s = ContactStats {
            broadphase_pairs: 50,
            contacts_generated: 10,
            islands: 3,
            ..ContactStats::default()
        };
        s.reset();
        assert_eq!(s.broadphase_pairs, 0);
        assert_eq!(s.contacts_generated, 0);
        assert_eq!(s.islands, 0);
    }
}
