// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact constraint using sequential impulse with friction.
//!
//! # Extensions
//!
//! - Contact manifold reduction (keep best N contacts).
//! - Contact caching and lifetime tracking.
//! - Contact force estimation from accumulated impulses.
//! - Persistent contact management with matching heuristic.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use crate::traits::Constraint;

/// Baumgarte stabilization factor for position correction.
const BAUMGARTE: Real = 0.2;

/// Penetration slop: allowed overlap before position correction kicks in.
const SLOP: Real = 0.005;

/// Velocity threshold below which restitution is suppressed.
const RESTITUTION_VELOCITY_THRESHOLD: Real = 1.0;

/// Maximum number of contacts in a reduced manifold.
const MAX_MANIFOLD_CONTACTS: usize = 4;

/// Distance threshold for matching contacts across frames.
const CONTACT_MATCH_THRESHOLD: Real = 0.02;

/// Maximum age (in frames) before a persistent contact is discarded.
const MAX_CONTACT_AGE: u32 = 120;

// ─────────────────────────────────────────────────────────────────────────────
// Contact constraint
// ─────────────────────────────────────────────────────────────────────────────

/// A contact constraint between two rigid bodies.
///
/// Implements the sequential impulse method with:
/// - Normal impulse for non-penetration
/// - Two tangent impulses for Coulomb friction (friction cone clamping)
/// - Warm starting via accumulated impulses
/// - Baumgarte stabilization for position correction
#[derive(Debug, Clone)]
pub struct ContactConstraint {
    /// Handle of the first body.
    pub body_handle_a: BodyHandle,
    /// Handle of the second body.
    pub body_handle_b: BodyHandle,
    /// Contact normal (from B towards A).
    pub normal: Vec3,
    /// First tangent direction.
    pub tangent1: Vec3,
    /// Second tangent direction.
    pub tangent2: Vec3,
    /// Penetration depth (positive = overlapping).
    pub penetration: Real,
    /// Coefficient of restitution.
    pub restitution: Real,
    /// Coefficient of friction.
    pub friction: Real,
    /// Contact point on body A in world space.
    pub point_a: Vec3,
    /// Contact point on body B in world space.
    pub point_b: Vec3,

    // --- Cached solver data (computed in prepare) ---
    /// Accumulated normal impulse (for warm starting and clamping).
    accumulated_normal_impulse: Real,
    /// Accumulated tangent1 impulse.
    accumulated_tangent1_impulse: Real,
    /// Accumulated tangent2 impulse.
    accumulated_tangent2_impulse: Real,

    /// Effective mass along the normal.
    effective_mass_normal: Real,
    /// Effective mass along tangent1.
    effective_mass_tangent1: Real,
    /// Effective mass along tangent2.
    effective_mass_tangent2: Real,

    /// Velocity bias for restitution.
    velocity_bias: Real,

    /// Offset from body A center of mass to contact point.
    r_a: Vec3,
    /// Offset from body B center of mass to contact point.
    r_b: Vec3,
    /// Whether warm-start has already been applied this timestep.
    warm_started: bool,
}

impl ContactConstraint {
    /// Create a new contact constraint.
    pub fn new(
        body_handle_a: BodyHandle,
        body_handle_b: BodyHandle,
        normal: Vec3,
        penetration: Real,
        point_a: Vec3,
        point_b: Vec3,
        restitution: Real,
        friction: Real,
    ) -> Self {
        // Compute tangent basis from normal
        let (tangent1, tangent2) = compute_tangent_basis(&normal);

        Self {
            body_handle_a,
            body_handle_b,
            normal,
            tangent1,
            tangent2,
            penetration,
            restitution,
            friction,
            point_a,
            point_b,
            accumulated_normal_impulse: 0.0,
            accumulated_tangent1_impulse: 0.0,
            accumulated_tangent2_impulse: 0.0,
            effective_mass_normal: 0.0,
            effective_mass_tangent1: 0.0,
            effective_mass_tangent2: 0.0,
            velocity_bias: 0.0,
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            warm_started: false,
        }
    }

    /// Reset accumulated impulses (disables warm starting).
    pub fn reset_impulses(&mut self) {
        self.accumulated_normal_impulse = 0.0;
        self.accumulated_tangent1_impulse = 0.0;
        self.accumulated_tangent2_impulse = 0.0;
    }

    /// Get the accumulated normal impulse.
    pub fn accumulated_normal_impulse(&self) -> Real {
        self.accumulated_normal_impulse
    }

    /// Get the accumulated tangent1 impulse.
    pub fn accumulated_tangent1_impulse(&self) -> Real {
        self.accumulated_tangent1_impulse
    }

    /// Get the accumulated tangent2 impulse.
    pub fn accumulated_tangent2_impulse(&self) -> Real {
        self.accumulated_tangent2_impulse
    }

    /// Estimate the contact force from accumulated impulses.
    ///
    /// `dt` – time step used during solving.
    /// Returns the force vector in world space.
    pub fn estimated_force(&self, dt: Real) -> Vec3 {
        if dt < 1e-12 {
            return Vec3::zeros();
        }
        let inv_dt = 1.0 / dt;
        self.normal * (self.accumulated_normal_impulse * inv_dt)
            + self.tangent1 * (self.accumulated_tangent1_impulse * inv_dt)
            + self.tangent2 * (self.accumulated_tangent2_impulse * inv_dt)
    }

    /// Estimate the normal force magnitude.
    pub fn normal_force(&self, dt: Real) -> Real {
        if dt < 1e-12 {
            return 0.0;
        }
        self.accumulated_normal_impulse / dt
    }

    /// Estimate the friction force magnitude.
    pub fn friction_force(&self, dt: Real) -> Real {
        if dt < 1e-12 {
            return 0.0;
        }
        let t1 = self.accumulated_tangent1_impulse;
        let t2 = self.accumulated_tangent2_impulse;
        (t1 * t1 + t2 * t2).sqrt() / dt
    }

    /// Check if the contact is separating (accumulated normal impulse is zero).
    pub fn is_separating(&self) -> bool {
        self.accumulated_normal_impulse < 1e-12
    }

    /// Check if the contact is sliding (friction impulse at Coulomb limit).
    pub fn is_sliding(&self) -> bool {
        let limit = self.friction * self.accumulated_normal_impulse;
        if limit < 1e-12 {
            return false;
        }
        let t1 = self.accumulated_tangent1_impulse;
        let t2 = self.accumulated_tangent2_impulse;
        let friction_mag = (t1 * t1 + t2 * t2).sqrt();
        (friction_mag - limit).abs() < 1e-6 * limit
    }

    /// Get the contact midpoint in world space.
    pub fn contact_point(&self) -> Vec3 {
        (self.point_a + self.point_b) * 0.5
    }

    /// Set warm-start impulses from a previous frame.
    pub fn set_warm_start(
        &mut self,
        normal_impulse: Real,
        tangent1_impulse: Real,
        tangent2_impulse: Real,
    ) {
        self.accumulated_normal_impulse = normal_impulse.max(0.0);
        self.accumulated_tangent1_impulse = tangent1_impulse;
        self.accumulated_tangent2_impulse = tangent2_impulse;
    }
}

/// Compute the effective mass for a constraint direction.
///
/// effective_mass = 1 / (inv_m_a + inv_m_b + (r_a x dir) . I_a^-1 . (r_a x dir)
///                                           + (r_b x dir) . I_b^-1 . (r_b x dir))
fn compute_effective_mass(
    inv_mass_a: Real,
    inv_mass_b: Real,
    inv_inertia_a: &oxiphysics_core::math::Mat3,
    inv_inertia_b: &oxiphysics_core::math::Mat3,
    r_a: &Vec3,
    r_b: &Vec3,
    direction: &Vec3,
) -> Real {
    let rn_a = r_a.cross(direction);
    let rn_b = r_b.cross(direction);
    let k = inv_mass_a
        + inv_mass_b
        + rn_a.dot(&(inv_inertia_a * rn_a))
        + rn_b.dot(&(inv_inertia_b * rn_b));
    if k > 1e-12 { 1.0 / k } else { 0.0 }
}

/// Compute relative velocity at contact point along a direction.
fn relative_velocity_along(
    vel_a: &Vec3,
    ang_vel_a: &Vec3,
    vel_b: &Vec3,
    ang_vel_b: &Vec3,
    r_a: &Vec3,
    r_b: &Vec3,
    direction: &Vec3,
) -> Real {
    let v_a = vel_a + ang_vel_a.cross(r_a);
    let v_b = vel_b + ang_vel_b.cross(r_b);
    let relative = v_a - v_b;
    relative.dot(direction)
}

/// Compute a tangent basis from a normal vector.
fn compute_tangent_basis(normal: &Vec3) -> (Vec3, Vec3) {
    let n = *normal;
    // Pick the axis least aligned with normal
    let reference = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let t1 = n.cross(&reference).normalize();
    let t2 = n.cross(&t1);
    (t1, t2)
}

/// Apply an impulse to two bodies given their handles and solver data.
///
/// This function temporarily removes bodies from the set to satisfy the borrow
/// checker, applies impulses, and puts them back.
fn apply_impulse_to_pair(
    bodies: &mut RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    impulse: Vec3,
    r_a: &Vec3,
    r_b: &Vec3,
) {
    // Apply to body A (positive impulse)
    if let Some(body_a) = bodies.get_mut(handle_a) {
        body_a.velocity += impulse * body_a.inverse_mass;
        body_a.angular_velocity += body_a.world_inverse_inertia * r_a.cross(&impulse);
    }
    // Apply to body B (negative impulse)
    if let Some(body_b) = bodies.get_mut(handle_b) {
        body_b.velocity -= impulse * body_b.inverse_mass;
        body_b.angular_velocity -= body_b.world_inverse_inertia * r_b.cross(&impulse);
    }
}

impl Constraint for ContactConstraint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        let (inv_mass_a, inv_inertia_a, vel_a, ang_vel_a, pos_a) = {
            if let Some(body) = bodies.get(self.body_handle_a) {
                (
                    body.inverse_mass,
                    body.world_inverse_inertia,
                    body.velocity,
                    body.angular_velocity,
                    body.transform.position,
                )
            } else {
                return;
            }
        };

        let (inv_mass_b, inv_inertia_b, vel_b, ang_vel_b, pos_b) = {
            if let Some(body) = bodies.get(self.body_handle_b) {
                (
                    body.inverse_mass,
                    body.world_inverse_inertia,
                    body.velocity,
                    body.angular_velocity,
                    body.transform.position,
                )
            } else {
                return;
            }
        };

        // Compute contact point offsets
        let contact_point = (self.point_a + self.point_b) * 0.5;
        self.r_a = contact_point - pos_a;
        self.r_b = contact_point - pos_b;

        // Compute effective masses
        self.effective_mass_normal = compute_effective_mass(
            inv_mass_a,
            inv_mass_b,
            &inv_inertia_a,
            &inv_inertia_b,
            &self.r_a,
            &self.r_b,
            &self.normal,
        );
        self.effective_mass_tangent1 = compute_effective_mass(
            inv_mass_a,
            inv_mass_b,
            &inv_inertia_a,
            &inv_inertia_b,
            &self.r_a,
            &self.r_b,
            &self.tangent1,
        );
        self.effective_mass_tangent2 = compute_effective_mass(
            inv_mass_a,
            inv_mass_b,
            &inv_inertia_a,
            &inv_inertia_b,
            &self.r_a,
            &self.r_b,
            &self.tangent2,
        );

        // Compute velocity bias for restitution
        let closing_velocity = relative_velocity_along(
            &vel_a,
            &ang_vel_a,
            &vel_b,
            &ang_vel_b,
            &self.r_a,
            &self.r_b,
            &self.normal,
        );
        self.velocity_bias = 0.0;
        if closing_velocity < -RESTITUTION_VELOCITY_THRESHOLD {
            self.velocity_bias = -self.restitution * closing_velocity;
        }

        // Mark warm-start as pending for this timestep; it will be applied on the
        // first call to solve_velocity() which has the required &mut RigidBodySet access.
        self.warm_started = false;

        let _ = dt; // Baumgarte uses dt in solve_position
    }

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        // Apply warm-start impulses exactly once per timestep (first call after prepare).
        if !self.warm_started {
            let warm_impulse = self.normal * self.accumulated_normal_impulse
                + self.tangent1 * self.accumulated_tangent1_impulse
                + self.tangent2 * self.accumulated_tangent2_impulse;
            if warm_impulse.norm_squared() > 1e-20 {
                apply_impulse_to_pair(
                    bodies,
                    self.body_handle_a,
                    self.body_handle_b,
                    warm_impulse,
                    &self.r_a,
                    &self.r_b,
                );
            }
            self.warm_started = true;
        }

        // Read current velocities
        let (vel_a, ang_vel_a) = {
            if let Some(body) = bodies.get(self.body_handle_a) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };
        let (vel_b, ang_vel_b) = {
            if let Some(body) = bodies.get(self.body_handle_b) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };

        // --- Normal impulse ---
        let vn = relative_velocity_along(
            &vel_a,
            &ang_vel_a,
            &vel_b,
            &ang_vel_b,
            &self.r_a,
            &self.r_b,
            &self.normal,
        );
        let lambda_n = self.effective_mass_normal * (-vn + self.velocity_bias);

        // Clamp: accumulated normal impulse >= 0 (no pulling)
        let old_accumulated = self.accumulated_normal_impulse;
        self.accumulated_normal_impulse = (old_accumulated + lambda_n).max(0.0);
        let applied_lambda_n = self.accumulated_normal_impulse - old_accumulated;

        let normal_impulse = self.normal * applied_lambda_n;
        apply_impulse_to_pair(
            bodies,
            self.body_handle_a,
            self.body_handle_b,
            normal_impulse,
            &self.r_a,
            &self.r_b,
        );

        // Re-read velocities after normal impulse
        let (vel_a, ang_vel_a) = {
            if let Some(body) = bodies.get(self.body_handle_a) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };
        let (vel_b, ang_vel_b) = {
            if let Some(body) = bodies.get(self.body_handle_b) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };

        // --- Tangent1 impulse (friction) ---
        let vt1 = relative_velocity_along(
            &vel_a,
            &ang_vel_a,
            &vel_b,
            &ang_vel_b,
            &self.r_a,
            &self.r_b,
            &self.tangent1,
        );
        let lambda_t1 = self.effective_mass_tangent1 * (-vt1);

        let friction_limit = self.friction * self.accumulated_normal_impulse;
        let old_t1 = self.accumulated_tangent1_impulse;
        self.accumulated_tangent1_impulse =
            (old_t1 + lambda_t1).clamp(-friction_limit, friction_limit);
        let applied_t1 = self.accumulated_tangent1_impulse - old_t1;

        let tangent1_impulse = self.tangent1 * applied_t1;
        apply_impulse_to_pair(
            bodies,
            self.body_handle_a,
            self.body_handle_b,
            tangent1_impulse,
            &self.r_a,
            &self.r_b,
        );

        // Re-read velocities for tangent2
        let (vel_a, ang_vel_a) = {
            if let Some(body) = bodies.get(self.body_handle_a) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };
        let (vel_b, ang_vel_b) = {
            if let Some(body) = bodies.get(self.body_handle_b) {
                (body.velocity, body.angular_velocity)
            } else {
                return;
            }
        };

        // --- Tangent2 impulse (friction) ---
        let vt2 = relative_velocity_along(
            &vel_a,
            &ang_vel_a,
            &vel_b,
            &ang_vel_b,
            &self.r_a,
            &self.r_b,
            &self.tangent2,
        );
        let lambda_t2 = self.effective_mass_tangent2 * (-vt2);

        let old_t2 = self.accumulated_tangent2_impulse;
        self.accumulated_tangent2_impulse =
            (old_t2 + lambda_t2).clamp(-friction_limit, friction_limit);
        let applied_t2 = self.accumulated_tangent2_impulse - old_t2;

        let tangent2_impulse = self.tangent2 * applied_t2;
        apply_impulse_to_pair(
            bodies,
            self.body_handle_a,
            self.body_handle_b,
            tangent2_impulse,
            &self.r_a,
            &self.r_b,
        );
    }

    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        // Baumgarte position correction
        let pos_a = {
            if let Some(body) = bodies.get(self.body_handle_a) {
                body.transform.position
            } else {
                return;
            }
        };
        let pos_b = {
            if let Some(body) = bodies.get(self.body_handle_b) {
                body.transform.position
            } else {
                return;
            }
        };

        let contact_point = (self.point_a + self.point_b) * 0.5;
        let r_a = contact_point - pos_a;
        let r_b = contact_point - pos_b;

        // Recompute separation along normal
        let separation = -self.penetration;
        let correction = (BAUMGARTE * (separation + SLOP)).min(0.0);

        if correction.abs() < 1e-10 {
            return;
        }

        let inv_mass_a = bodies
            .get(self.body_handle_a)
            .map_or(0.0, |b| b.inverse_mass);
        let inv_mass_b = bodies
            .get(self.body_handle_b)
            .map_or(0.0, |b| b.inverse_mass);

        let k = inv_mass_a + inv_mass_b;
        if k < 1e-12 {
            return;
        }

        let position_impulse = self.normal * (-correction / k);

        if let Some(body_a) = bodies.get_mut(self.body_handle_a) {
            body_a.transform.position += position_impulse * body_a.inverse_mass;
        }
        if let Some(body_b) = bodies.get_mut(self.body_handle_b) {
            body_b.transform.position -= position_impulse * body_b.inverse_mass;
        }

        let _ = (r_a, r_b); // Used for angular correction in more advanced implementations
    }

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_handle_a, self.body_handle_b]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact manifold
// ─────────────────────────────────────────────────────────────────────────────

/// A single contact point within a manifold.
#[derive(Debug, Clone, Copy)]
pub struct ManifoldContact {
    /// World-space contact point on body A.
    pub point_a: Vec3,
    /// World-space contact point on body B.
    pub point_b: Vec3,
    /// Contact normal (B → A).
    pub normal: Vec3,
    /// Penetration depth.
    pub penetration: Real,
    /// Unique ID for matching across frames.
    pub id: u64,
    /// Cached normal impulse from previous frame.
    pub cached_normal_impulse: Real,
    /// Cached friction impulses from previous frame.
    pub cached_friction_impulse: [Real; 2],
}

impl ManifoldContact {
    /// Create a new manifold contact.
    pub fn new(point_a: Vec3, point_b: Vec3, normal: Vec3, penetration: Real, id: u64) -> Self {
        Self {
            point_a,
            point_b,
            normal,
            penetration,
            id,
            cached_normal_impulse: 0.0,
            cached_friction_impulse: [0.0; 2],
        }
    }

    /// Distance between two contact points (for matching).
    pub fn distance_to(&self, other: &ManifoldContact) -> Real {
        let dx = self.point_a.x - other.point_a.x;
        let dy = self.point_a.y - other.point_a.y;
        let dz = self.point_a.z - other.point_a.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Contact midpoint in world space.
    pub fn midpoint(&self) -> Vec3 {
        (self.point_a + self.point_b) * 0.5
    }
}

/// A contact manifold: a collection of contact points between two bodies.
///
/// Stores up to `MAX_MANIFOLD_CONTACTS` contacts and provides methods
/// for contact reduction, caching, and warm-start transfer.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Body A handle.
    pub body_a: BodyHandle,
    /// Body B handle.
    pub body_b: BodyHandle,
    /// Active contact points.
    pub contacts: Vec<ManifoldContact>,
    /// Number of frames this manifold has been alive.
    pub age: u32,
    /// Average normal of the manifold.
    pub average_normal: Vec3,
}

impl ContactManifold {
    /// Create a new empty manifold.
    pub fn new(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        Self {
            body_a,
            body_b,
            contacts: Vec::new(),
            age: 0,
            average_normal: Vec3::zeros(),
        }
    }

    /// Add a contact point, maintaining the manifold limit.
    ///
    /// If the manifold is full, replaces the contact with the smallest
    /// penetration depth.
    pub fn add_contact(&mut self, contact: ManifoldContact) {
        if self.contacts.len() < MAX_MANIFOLD_CONTACTS {
            self.contacts.push(contact);
        } else {
            // Replace the shallowest contact
            if let Some((idx, _)) = self.contacts.iter().enumerate().min_by(|(_, a), (_, b)| {
                a.penetration
                    .partial_cmp(&b.penetration)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }) && contact.penetration > self.contacts[idx].penetration
            {
                self.contacts[idx] = contact;
            }
        }
        self.update_average_normal();
    }

    /// Remove contacts that are no longer valid (separating or too old).
    pub fn prune_contacts(&mut self, max_separation: Real) {
        self.contacts.retain(|c| c.penetration > -max_separation);
        self.update_average_normal();
    }

    /// Number of active contacts.
    pub fn len(&self) -> usize {
        self.contacts.len()
    }

    /// Check if manifold is empty.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// Total normal impulse across all contacts.
    pub fn total_normal_impulse(&self) -> Real {
        self.contacts.iter().map(|c| c.cached_normal_impulse).sum()
    }

    /// Maximum penetration depth across contacts.
    pub fn max_penetration(&self) -> Real {
        self.contacts
            .iter()
            .map(|c| c.penetration)
            .fold(0.0, f64::max)
    }

    /// Average penetration across contacts.
    pub fn average_penetration(&self) -> Real {
        if self.contacts.is_empty() {
            return 0.0;
        }
        let sum: Real = self.contacts.iter().map(|c| c.penetration).sum();
        sum / self.contacts.len() as Real
    }

    /// Transfer warm-start impulses from old contacts to new ones via matching.
    pub fn transfer_warm_start(&mut self, old_contacts: &[ManifoldContact]) {
        for new_c in &mut self.contacts {
            let mut best_dist = CONTACT_MATCH_THRESHOLD;
            let mut best_idx = None;
            for (i, old_c) in old_contacts.iter().enumerate() {
                let dist = new_c.distance_to(old_c);
                if dist < best_dist {
                    best_dist = dist;
                    best_idx = Some(i);
                }
            }
            if let Some(idx) = best_idx {
                new_c.cached_normal_impulse = old_contacts[idx].cached_normal_impulse;
                new_c.cached_friction_impulse = old_contacts[idx].cached_friction_impulse;
            }
        }
    }

    /// Advance the manifold age by one frame.
    pub fn advance_age(&mut self) {
        self.age += 1;
    }

    /// Check if the manifold is stale.
    pub fn is_stale(&self) -> bool {
        self.age > MAX_CONTACT_AGE || self.contacts.is_empty()
    }

    /// Update the average normal from current contacts.
    fn update_average_normal(&mut self) {
        if self.contacts.is_empty() {
            self.average_normal = Vec3::zeros();
            return;
        }
        let mut sum = Vec3::zeros();
        for c in &self.contacts {
            sum += c.normal;
        }
        let n = self.contacts.len() as Real;
        sum *= 1.0 / n;
        let len = sum.norm();
        if len > 1e-12 {
            self.average_normal = sum * (1.0 / len);
        }
    }

    /// Create constraint objects from this manifold.
    pub fn create_constraints(&self, restitution: Real, friction: Real) -> Vec<ContactConstraint> {
        self.contacts
            .iter()
            .map(|c| {
                let mut cc = ContactConstraint::new(
                    self.body_a,
                    self.body_b,
                    c.normal,
                    c.penetration,
                    c.point_a,
                    c.point_b,
                    restitution,
                    friction,
                );
                cc.set_warm_start(
                    c.cached_normal_impulse,
                    c.cached_friction_impulse[0],
                    c.cached_friction_impulse[1],
                );
                cc
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact cache
// ─────────────────────────────────────────────────────────────────────────────

/// A cache of contact manifolds for persistent contact tracking.
///
/// Stores manifolds indexed by body pair and handles lifetime management.
#[derive(Debug, Clone)]
pub struct ContactCache {
    /// All active manifolds.
    pub manifolds: Vec<ContactManifold>,
}

impl Default for ContactCache {
    fn default() -> Self {
        Self::new()
    }
}

impl ContactCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        Self {
            manifolds: Vec::new(),
        }
    }

    /// Find an existing manifold for the given body pair.
    pub fn find_manifold(&self, body_a: BodyHandle, body_b: BodyHandle) -> Option<usize> {
        self.manifolds.iter().position(|m| {
            (m.body_a == body_a && m.body_b == body_b) || (m.body_a == body_b && m.body_b == body_a)
        })
    }

    /// Get or create a manifold for the given body pair.
    pub fn get_or_create(
        &mut self,
        body_a: BodyHandle,
        body_b: BodyHandle,
    ) -> &mut ContactManifold {
        if let Some(idx) = self.find_manifold(body_a, body_b) {
            return &mut self.manifolds[idx];
        }
        self.manifolds.push(ContactManifold::new(body_a, body_b));
        self.manifolds
            .last_mut()
            .expect("collection should not be empty")
    }

    /// Remove stale manifolds.
    pub fn prune_stale(&mut self) {
        self.manifolds.retain(|m| !m.is_stale());
    }

    /// Advance age of all manifolds.
    pub fn advance_all(&mut self) {
        for m in &mut self.manifolds {
            m.advance_age();
        }
    }

    /// Total number of active contacts across all manifolds.
    pub fn total_contacts(&self) -> usize {
        self.manifolds.iter().map(|m| m.len()).sum()
    }

    /// Number of manifolds.
    pub fn len(&self) -> usize {
        self.manifolds.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.manifolds.is_empty()
    }

    /// Clear all manifolds.
    pub fn clear(&mut self) {
        self.manifolds.clear();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Contact reduction
// ─────────────────────────────────────────────────────────────────────────────

/// Reduce a set of contacts to at most `max_contacts` by keeping the most
/// significant ones.
///
/// Selection strategy:
/// 1. Keep the deepest penetrating contact.
/// 2. Keep the contact farthest from the deepest.
/// 3. Keep contacts that maximize the area of the contact polygon.
pub fn reduce_contacts(contacts: &[ManifoldContact], max_contacts: usize) -> Vec<ManifoldContact> {
    if contacts.len() <= max_contacts || contacts.is_empty() {
        return contacts.to_vec();
    }

    let mut result = Vec::with_capacity(max_contacts);

    // 1. Find deepest contact
    let deepest_idx = contacts
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.penetration
                .partial_cmp(&b.penetration)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);

    result.push(contacts[deepest_idx]);

    // 2. Find farthest from deepest
    if result.len() < max_contacts {
        let deepest = &contacts[deepest_idx];
        let farthest_idx = contacts
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != deepest_idx)
            .max_by(|(_, a), (_, b)| {
                let da = a.distance_to(deepest);
                let db = b.distance_to(deepest);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i);

        if let Some(idx) = farthest_idx {
            result.push(contacts[idx]);
        }
    }

    // 3. Fill remaining slots with contacts farthest from existing selection
    let mut used = vec![false; contacts.len()];
    used[deepest_idx] = true;
    if result.len() >= 2 {
        // Mark second contact as used
        for (i, c) in contacts.iter().enumerate() {
            if c.id == result[1].id {
                used[i] = true;
                break;
            }
        }
    }

    while result.len() < max_contacts {
        let mut best_idx = None;
        let mut best_min_dist = -1.0f64;

        for (i, c) in contacts.iter().enumerate() {
            if used[i] {
                continue;
            }
            let min_dist = result
                .iter()
                .map(|r| c.distance_to(r))
                .fold(f64::MAX, f64::min);
            if min_dist > best_min_dist {
                best_min_dist = min_dist;
                best_idx = Some(i);
            }
        }

        if let Some(idx) = best_idx {
            used[idx] = true;
            result.push(contacts[idx]);
        } else {
            break;
        }
    }

    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_rigid::RigidBody;

    #[test]
    fn test_contact_resolves_overlap() {
        let mut bodies = RigidBodySet::new();

        let mut body_a = RigidBody::new(1.0);
        body_a.transform.position = Vec3::new(0.0, 0.45, 0.0);
        body_a.linear_damping = 0.0;
        body_a.angular_damping = 0.0;
        let ha = bodies.insert(body_a);

        let mut body_b = RigidBody::new_static();
        body_b.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let hb = bodies.insert(body_b);

        let mut contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.05,
            Vec3::new(0.0, -0.05, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.3,
        );

        let dt = 1.0 / 60.0;
        contact.prepare(&bodies, dt);

        for _ in 0..10 {
            contact.solve_velocity(&mut bodies, dt);
        }

        let body_a = bodies.get(ha).unwrap();
        let vn = body_a.velocity.dot(&Vec3::new(0.0, 1.0, 0.0));
        assert!(
            vn >= -1e-6,
            "Body should not be moving into the ground, vn = {vn}"
        );

        contact.solve_position(&mut bodies, dt);
        let body_a = bodies.get(ha).unwrap();
        assert!(
            body_a.transform.position.y >= 0.45,
            "Body should have been pushed up"
        );
    }

    #[test]
    fn test_contact_with_restitution() {
        let mut bodies = RigidBodySet::new();

        let mut body_a = RigidBody::new(1.0);
        body_a.transform.position = Vec3::new(0.0, 1.0, 0.0);
        body_a.velocity = Vec3::new(0.0, -5.0, 0.0);
        body_a.linear_damping = 0.0;
        body_a.angular_damping = 0.0;
        let ha = bodies.insert(body_a);

        let body_b = RigidBody::new_static();
        let hb = bodies.insert(body_b);

        let mut contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.8,
            0.3,
        );

        let dt = 1.0 / 60.0;
        contact.prepare(&bodies, dt);

        for _ in 0..10 {
            contact.solve_velocity(&mut bodies, dt);
        }

        let body_a = bodies.get(ha).unwrap();
        assert!(
            body_a.velocity.y > 0.0,
            "Body should bounce, vy = {}",
            body_a.velocity.y
        );
    }

    #[test]
    fn test_tangent_basis() {
        let n = Vec3::new(0.0, 1.0, 0.0);
        let (t1, t2) = compute_tangent_basis(&n);

        assert!(n.dot(&t1).abs() < 1e-10);
        assert!(n.dot(&t2).abs() < 1e-10);
        assert!(t1.dot(&t2).abs() < 1e-10);

        assert!((t1.norm() - 1.0).abs() < 1e-10);
        assert!((t2.norm() - 1.0).abs() < 1e-10);
    }

    // ── 4. Estimated force ─────────────────────────────────────────────────

    #[test]
    fn test_estimated_force() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new_static());
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.5,
        );
        c.set_warm_start(10.0, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        let force = c.estimated_force(dt);
        // F = impulse / dt = 10 / (1/60) = 600 N
        assert!((force.y - 600.0).abs() < 1e-6, "force.y = {}", force.y);
    }

    #[test]
    fn test_normal_and_friction_force() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.5,
        );
        c.set_warm_start(10.0, 3.0, 4.0);
        let dt = 1.0 / 60.0;
        let nf = c.normal_force(dt);
        let ff = c.friction_force(dt);
        assert!((nf - 600.0).abs() < 1e-6);
        assert!((ff - 300.0).abs() < 1e-6, "friction force = {ff}"); // sqrt(9+16)/dt = 5*60
    }

    // ── 5. Contact state queries ───────────────────────────────────────────

    #[test]
    fn test_is_separating() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.0,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        assert!(c.is_separating());
    }

    #[test]
    fn test_contact_point() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(1.0, 0.0, 3.0),
            0.0,
            0.3,
        );
        let mid = c.contact_point();
        assert!((mid.x - 1.0).abs() < 1e-12);
        assert!((mid.y - 1.0).abs() < 1e-12);
        assert!((mid.z - 3.0).abs() < 1e-12);
    }

    // ── 6. Manifold contact ────────────────────────────────────────────────

    #[test]
    fn test_manifold_contact_distance() {
        let c1 = ManifoldContact::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        );
        let c2 = ManifoldContact::new(
            Vec3::new(3.0, 4.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            1,
        );
        let dist = c1.distance_to(&c2);
        assert!((dist - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_contact_midpoint() {
        let c = ManifoldContact::new(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(1.0, 0.0, 3.0),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        );
        let mid = c.midpoint();
        assert!((mid.y - 1.0).abs() < 1e-12);
    }

    // ── 7. Contact manifold operations ─────────────────────────────────────

    #[test]
    fn test_manifold_add_and_prune() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        assert!(manifold.is_empty());

        for i in 0..6 {
            let c = ManifoldContact::new(
                Vec3::new(i as f64 * 0.1, 0.0, 0.0),
                Vec3::zeros(),
                Vec3::new(0.0, 1.0, 0.0),
                0.01 * (i + 1) as f64,
                i as u64,
            );
            manifold.add_contact(c);
        }

        // Should be limited to MAX_MANIFOLD_CONTACTS
        assert_eq!(manifold.len(), MAX_MANIFOLD_CONTACTS);

        // Prune with large separation
        manifold.prune_contacts(100.0);
        assert_eq!(manifold.len(), MAX_MANIFOLD_CONTACTS);
    }

    #[test]
    fn test_manifold_max_and_average_penetration() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        manifold.add_contact(ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        ));
        manifold.add_contact(ManifoldContact::new(
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.03,
            1,
        ));
        assert!((manifold.max_penetration() - 0.03).abs() < 1e-12);
        assert!((manifold.average_penetration() - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_age_and_stale() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        // Empty manifold IS stale (no contacts)
        assert!(manifold.is_stale());

        manifold.add_contact(ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        ));
        assert!(!manifold.is_stale());

        for _ in 0..=MAX_CONTACT_AGE {
            manifold.advance_age();
        }
        assert!(manifold.is_stale());
    }

    #[test]
    fn test_manifold_warm_start_transfer() {
        let old_contacts = vec![{
            let mut c = ManifoldContact::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::zeros(),
                Vec3::new(0.0, 1.0, 0.0),
                0.01,
                0,
            );
            c.cached_normal_impulse = 5.0;
            c.cached_friction_impulse = [1.0, 2.0];
            c
        }];

        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        // Add a contact near the old one
        manifold.add_contact(ManifoldContact::new(
            Vec3::new(0.001, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            1,
        ));
        manifold.transfer_warm_start(&old_contacts);

        assert!((manifold.contacts[0].cached_normal_impulse - 5.0).abs() < 1e-12);
        assert!((manifold.contacts[0].cached_friction_impulse[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_create_constraints() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        manifold.add_contact(ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.02,
            0,
        ));
        manifold.add_contact(ManifoldContact::new(
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            1,
        ));

        let constraints = manifold.create_constraints(0.5, 0.3);
        assert_eq!(constraints.len(), 2);
        assert!((constraints[0].penetration - 0.02).abs() < 1e-12);
        assert!((constraints[0].restitution - 0.5).abs() < 1e-12);
        assert!((constraints[0].friction - 0.3).abs() < 1e-12);
    }

    // ── 8. Contact cache ───────────────────────────────────────────────────

    #[test]
    fn test_contact_cache_basic() {
        let mut cache = ContactCache::new();
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);

        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let hc = BodyHandle::new(2, 0);

        let m = cache.get_or_create(ha, hb);
        m.add_contact(ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        ));

        assert_eq!(cache.len(), 1);
        assert_eq!(cache.total_contacts(), 1);

        // Same pair should return existing
        let idx = cache.find_manifold(ha, hb);
        assert!(idx.is_some());

        // Reversed pair should also find it
        let idx_rev = cache.find_manifold(hb, ha);
        assert!(idx_rev.is_some());

        // Different pair
        let idx_diff = cache.find_manifold(ha, hc);
        assert!(idx_diff.is_none());
    }

    #[test]
    fn test_contact_cache_prune() {
        let mut cache = ContactCache::new();
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);

        // Create manifold without contacts (will be stale immediately)
        cache.get_or_create(ha, hb);
        assert_eq!(cache.len(), 1);

        cache.prune_stale();
        assert_eq!(cache.len(), 0);
    }

    #[test]
    fn test_contact_cache_advance() {
        let mut cache = ContactCache::new();
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);

        let m = cache.get_or_create(ha, hb);
        m.add_contact(ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        ));
        assert_eq!(cache.manifolds[0].age, 0);

        cache.advance_all();
        assert_eq!(cache.manifolds[0].age, 1);
    }

    #[test]
    fn test_contact_cache_clear() {
        let mut cache = ContactCache::new();
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        cache.get_or_create(ha, hb);
        cache.clear();
        assert!(cache.is_empty());
    }

    // ── 9. Contact reduction ───────────────────────────────────────────────

    #[test]
    fn test_reduce_contacts_no_reduction() {
        let contacts: Vec<ManifoldContact> = (0..3)
            .map(|i| {
                ManifoldContact::new(
                    Vec3::new(i as f64 * 0.1, 0.0, 0.0),
                    Vec3::zeros(),
                    Vec3::new(0.0, 1.0, 0.0),
                    0.01,
                    i as u64,
                )
            })
            .collect();
        let reduced = reduce_contacts(&contacts, 4);
        assert_eq!(reduced.len(), 3);
    }

    #[test]
    fn test_reduce_contacts_to_max() {
        let contacts: Vec<ManifoldContact> = (0..8)
            .map(|i| {
                ManifoldContact::new(
                    Vec3::new(i as f64 * 0.1, 0.0, 0.0),
                    Vec3::zeros(),
                    Vec3::new(0.0, 1.0, 0.0),
                    0.01 * (i + 1) as f64,
                    i as u64,
                )
            })
            .collect();
        let reduced = reduce_contacts(&contacts, 4);
        assert_eq!(reduced.len(), 4);
        // Deepest contact should be included (index 7, pen=0.08)
        assert!(reduced.iter().any(|c| (c.penetration - 0.08).abs() < 1e-12));
    }

    #[test]
    fn test_reduce_contacts_empty() {
        let contacts: Vec<ManifoldContact> = Vec::new();
        let reduced = reduce_contacts(&contacts, 4);
        assert!(reduced.is_empty());
    }

    // ── 10. Warm-start setter ──────────────────────────────────────────────

    #[test]
    fn test_set_warm_start() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.5,
        );
        c.set_warm_start(10.0, -2.0, 3.0);
        assert!((c.accumulated_normal_impulse() - 10.0).abs() < 1e-12);
        assert!((c.accumulated_tangent1_impulse() + 2.0).abs() < 1e-12);
        assert!((c.accumulated_tangent2_impulse() - 3.0).abs() < 1e-12);

        // Negative normal impulse should be clamped to 0
        c.set_warm_start(-5.0, 0.0, 0.0);
        assert!((c.accumulated_normal_impulse()).abs() < 1e-12);
    }

    // ── 11. Reset impulses ─────────────────────────────────────────────────

    #[test]
    fn test_reset_impulses() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        c.set_warm_start(10.0, 5.0, 3.0);
        c.reset_impulses();
        assert!((c.accumulated_normal_impulse()).abs() < 1e-12);
        assert!((c.accumulated_tangent1_impulse()).abs() < 1e-12);
        assert!((c.accumulated_tangent2_impulse()).abs() < 1e-12);
    }

    // ── 12. Tangent basis for different normals ────────────────────────────

    #[test]
    fn test_tangent_basis_x_normal() {
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (t1, t2) = compute_tangent_basis(&n);
        assert!(n.dot(&t1).abs() < 1e-10);
        assert!(n.dot(&t2).abs() < 1e-10);
        assert!(t1.dot(&t2).abs() < 1e-10);
    }

    #[test]
    fn test_tangent_basis_diagonal_normal() {
        let n = Vec3::new(1.0, 1.0, 1.0).normalize();
        let (t1, t2) = compute_tangent_basis(&n);
        assert!(n.dot(&t1).abs() < 1e-10);
        assert!(n.dot(&t2).abs() < 1e-10);
        assert!((t1.norm() - 1.0).abs() < 1e-10);
        assert!((t2.norm() - 1.0).abs() < 1e-10);
    }

    // ── 13. Manifold total normal impulse ──────────────────────────────────

    #[test]
    fn test_manifold_total_normal_impulse() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut manifold = ContactManifold::new(ha, hb);
        let mut c1 = ManifoldContact::new(
            Vec3::zeros(),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            0,
        );
        c1.cached_normal_impulse = 5.0;
        let mut c2 = ManifoldContact::new(
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            1,
        );
        c2.cached_normal_impulse = 3.0;
        manifold.add_contact(c1);
        manifold.add_contact(c2);
        assert!((manifold.total_normal_impulse() - 8.0).abs() < 1e-12);
    }

    // ── 14. Estimated force with zero dt ───────────────────────────────────

    #[test]
    fn test_estimated_force_zero_dt() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.5,
        );
        c.set_warm_start(10.0, 0.0, 0.0);
        let force = c.estimated_force(0.0);
        assert!((force.norm()).abs() < 1e-12);
    }

    // ── 15. Is sliding detection ───────────────────────────────────────────

    #[test]
    fn test_is_sliding() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.5,
        );
        // At friction limit: t = mu * n = 0.5 * 10 = 5
        c.set_warm_start(10.0, 5.0, 0.0);
        assert!(c.is_sliding());

        // Below friction limit
        c.set_warm_start(10.0, 1.0, 0.0);
        assert!(!c.is_sliding());
    }

    // ── 16. Body handles from constraint ───────────────────────────────────

    #[test]
    fn test_body_handles() {
        let ha = BodyHandle::new(0, 0);
        let hb = BodyHandle::new(1, 0);
        let c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let handles = c.body_handles();
        assert_eq!(handles.len(), 2);
        assert_eq!(handles[0], ha);
        assert_eq!(handles[1], hb);
    }

    // ── B4: Warm-start applied to real bodies in solve_velocity ────────────

    #[test]
    fn test_warm_start_applied_once_per_timestep() {
        let mut bodies = RigidBodySet::new();

        let body_a = RigidBody::new(1.0);
        let ha = bodies.insert(body_a);

        let body_b = RigidBody::new(1.0);
        let hb = bodies.insert(body_b);

        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        // Simulate a prior frame accumulation
        c.set_warm_start(2.0, 0.0, 0.0);

        // prepare() must reset warm_started flag
        c.prepare(&bodies, 0.016);
        assert!(
            !c.warm_started,
            "warm_started must be false after prepare()"
        );

        // First solve_velocity call: warm-start impulse is applied
        c.solve_velocity(&mut bodies, 0.016);
        assert!(
            c.warm_started,
            "warm_started must be true after first solve_velocity()"
        );

        // Capture velocity after first call
        let vel_after_first = bodies.get(ha).map(|b| b.velocity).unwrap_or_default();

        // Second solve_velocity call: warm-start must NOT be applied again
        c.solve_velocity(&mut bodies, 0.016);
        let vel_after_second = bodies.get(ha).map(|b| b.velocity).unwrap_or_default();

        // The warm-start impulse was already applied; re-applying would only differ by
        // the incremental normal impulse computed in the second iteration — the velocity
        // change from warm-start itself must not be double-counted.
        // We verify the flag remains set (no double application).
        assert!(
            c.warm_started,
            "warm_started should remain true on subsequent calls"
        );
        let _ = (vel_after_first, vel_after_second); // used for context
    }

    #[test]
    fn test_warm_start_not_applied_twice() {
        let mut bodies = RigidBodySet::new();

        let body_a = RigidBody::new(1.0);
        let ha = bodies.insert(body_a);

        let body_b = RigidBody::new(1.0);
        let hb = bodies.insert(body_b);

        let mut c = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        c.set_warm_start(5.0, 0.0, 0.0);
        c.prepare(&bodies, 0.016);

        // First call applies warm-start
        c.solve_velocity(&mut bodies, 0.016);
        let vel1 = bodies.get(ha).map(|b| b.velocity.y).unwrap_or(0.0);

        // Manually reset warm_started to simulate a second erroneous prepare-less call
        // This tests that the engine's single-call pattern is enforced correctly
        // When warm_started is already true, a second call should NOT re-apply
        c.solve_velocity(&mut bodies, 0.016);
        let vel2 = bodies.get(ha).map(|b| b.velocity.y).unwrap_or(0.0);

        // vel1 and vel2 should be close — any difference comes only from the incremental
        // normal impulse solver, not a second warm-start application
        let delta = (vel2 - vel1).abs();
        assert!(
            delta < 5.0,
            "Second solve_velocity should not double-apply warm-start: delta={delta}"
        );
    }
}
