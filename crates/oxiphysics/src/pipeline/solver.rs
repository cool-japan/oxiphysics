// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sequential-impulse contact solver.
//!
//! The solver turns the contact manifolds of one step into velocity-level
//! constraints and solves them with projected Gauss-Seidel (PGS) iterations:
//!
//! * **Accumulated impulses.** Every contact keeps the *total* normal impulse
//!   `lambda_n >= 0` and friction impulse `lambda_t` it has applied during the
//!   step. An iteration computes an increment, clamps the new *total*
//!   (`lambda_n >= 0`, `|lambda_t| <= mu * lambda_n`) and applies only the
//!   difference. Contacts can therefore never pull bodies together, and
//!   restitution is part of a fixed velocity target instead of being re-applied
//!   by every iteration.
//! * **Warm starting.** The accumulated impulses of the previous step are cached
//!   per contact point ([`WarmStartCache`]) and applied up front, so a resting
//!   contact starts each step (almost) at its solution.
//! * **Split-impulse position correction.** Penetration is removed by a second
//!   pass that works on separate *pseudo* motion which is applied once and then
//!   forgotten. The correction therefore does not depend on the number of
//!   iterations and never adds kinetic energy to the simulation.
//!
//! The contact normal convention is that of `oxiphysics-collision`: normals
//! point from body B toward body A, so a relative normal speed
//! `(v_a - v_b) . n < 0` means the bodies approach each other.

use std::collections::HashMap;

use oxiphysics_collision::ContactManifold;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_core::{BodyHandle, ColliderHandle};
use oxiphysics_rigid::{BodyState, BodyType, RigidBody, RigidBodySet};

/// Fraction of the (slop-reduced) penetration removed per step.
const POSITION_CORRECTION_RATE: Real = 0.2;

/// Penetration \[m\] that is tolerated without correction. Keeping a small
/// overlap alive lets resting contacts persist from step to step instead of
/// flickering between touching and separated.
const PENETRATION_SLOP: Real = 0.005;

/// Largest penetration \[m\] a single step tries to correct, so that deeply
/// interpenetrating bodies are separated gradually rather than violently.
const MAX_POSITION_CORRECTION: Real = 0.2;

/// Approach speed \[m/s\] below which collisions are perfectly inelastic. It
/// keeps resting contacts from bouncing on the speed that gravity adds in one
/// step.
const RESTITUTION_THRESHOLD: Real = 1.0;

/// Smallest inverse effective mass for which a constraint is solved.
const MIN_INVERSE_EFFECTIVE_MASS: Real = 1e-10;

/// Largest distance \[m\] between a contact point and a cached contact point
/// of the previous step for the cached impulse to be reused.
const WARM_START_MATCH_DISTANCE: Real = 0.05;

/// Smallest cosine between the normals of a contact and its cached
/// counterpart for the cached impulse to be reused.
const WARM_START_MIN_NORMAL_COSINE: Real = 0.8;

/// Velocity or displacement change below which a sleeping body is considered
/// undisturbed.
const WAKE_EPSILON: Real = 1e-9;

// ---------------------------------------------------------------------------
// Warm-start cache
// ---------------------------------------------------------------------------

/// Identity of a collider pair, independent of the order in which the pair was
/// reported.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PairKey {
    first: ColliderHandle,
    second: ColliderHandle,
}

impl PairKey {
    /// Canonical key of the pair `(a, b)` and whether the canonical order is
    /// the opposite of `(a, b)`.
    fn new(a: ColliderHandle, b: ColliderHandle) -> (Self, bool) {
        if (a.index, a.generation) <= (b.index, b.generation) {
            (
                Self {
                    first: a,
                    second: b,
                },
                false,
            )
        } else {
            (
                Self {
                    first: b,
                    second: a,
                },
                true,
            )
        }
    }
}

/// Accumulated impulses of one contact point, expressed for the canonical
/// pair order (body A is the pair's `first` collider).
#[derive(Debug, Clone)]
struct CachedContact {
    /// World-space contact point.
    point: Vec3,
    /// Unit contact normal pointing from `second` toward `first`.
    normal: Vec3,
    /// Accumulated normal impulse.
    normal_impulse: Real,
    /// Accumulated friction impulse acting on `first` (world space).
    friction_impulse: Vec3,
}

/// Accumulated contact impulses of the previous step, kept per contact point.
///
/// [`PhysicsPipeline`](super::PhysicsPipeline) feeds the cache back into the
/// solver as an initial guess (warm starting). Contacts are matched by
/// proximity of their contact points, so the cache survives changes in the
/// number or the order of the points of a manifold. It is part of a
/// [`PipelineCheckpoint`](super::PipelineCheckpoint).
#[derive(Debug, Clone, Default)]
pub struct WarmStartCache {
    /// Step size the cached impulses were accumulated over.
    step_dt: Real,
    /// Cached contacts per collider pair.
    contacts: HashMap<PairKey, Vec<CachedContact>>,
}

impl WarmStartCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// `true` if no contact is cached.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }

    /// Number of collider pairs with cached contacts.
    pub fn pair_count(&self) -> usize {
        self.contacts.len()
    }

    /// Total number of cached contact points.
    pub fn contact_count(&self) -> usize {
        self.contacts.values().map(Vec::len).sum()
    }

    /// Discard all cached impulses (the next step starts cold).
    pub fn clear(&mut self) {
        self.step_dt = 0.0;
        self.contacts.clear();
    }

    /// Sum of the cached normal impulses of all contact points \[N s\].
    pub fn normal_impulse_total(&self) -> Real {
        self.contacts
            .values()
            .flat_map(|list| list.iter().map(|c| c.normal_impulse))
            .sum()
    }

    /// Sum of the cached normal impulses of the collider pair `(a, b)` \[N s\],
    /// or `0.0` if the pair has no cached contact.
    pub fn normal_impulse(&self, a: ColliderHandle, b: ColliderHandle) -> Real {
        let (key, _) = PairKey::new(a, b);
        self.contacts
            .get(&key)
            .map_or(0.0, |list| list.iter().map(|c| c.normal_impulse).sum())
    }

    /// Cached normal impulses of the contact points of the pair `(a, b)`.
    pub fn contact_normal_impulses(&self, a: ColliderHandle, b: ColliderHandle) -> Vec<Real> {
        let (key, _) = PairKey::new(a, b);
        self.contacts
            .get(&key)
            .map(|list| list.iter().map(|c| c.normal_impulse).collect())
            .unwrap_or_default()
    }

    /// Index of the cached contact of `key` that continues a contact at
    /// `point` with canonical-order normal `normal`, if any.
    fn find_match(
        &self,
        key: &PairKey,
        point: &Vec3,
        normal: &Vec3,
        taken: &mut [bool],
    ) -> Option<&CachedContact> {
        let list = self.contacts.get(key)?;
        let mut best: Option<(usize, Real)> = None;
        for (idx, cached) in list.iter().enumerate() {
            if taken.get(idx).copied().unwrap_or(true) {
                continue;
            }
            if cached.normal.dot(normal) < WARM_START_MIN_NORMAL_COSINE {
                continue;
            }
            let distance = (cached.point - point).norm();
            if distance > WARM_START_MATCH_DISTANCE {
                continue;
            }
            if best.is_none_or(|(_, best_distance)| distance < best_distance) {
                best = Some((idx, distance));
            }
        }
        let (idx, _) = best?;
        if let Some(flag) = taken.get_mut(idx) {
            *flag = true;
        }
        list.get(idx)
    }
}

// ---------------------------------------------------------------------------
// Solver input
// ---------------------------------------------------------------------------

/// One colliding pair handed to the solver.
pub(super) struct ContactPair<'a> {
    /// Collider A of the pair.
    pub handle_a: ColliderHandle,
    /// Collider B of the pair.
    pub handle_b: ColliderHandle,
    /// Body that carries collider A, if any.
    pub body_a: Option<BodyHandle>,
    /// Body that carries collider B, if any.
    pub body_b: Option<BodyHandle>,
    /// Combined Coulomb friction coefficient of the pair.
    pub friction: Real,
    /// Combined coefficient of restitution of the pair.
    pub restitution: Real,
    /// Contacts of the pair; normals point from B toward A.
    pub manifold: &'a ContactManifold,
}

/// Parameters of one solver invocation.
#[derive(Debug, Clone, Copy)]
pub(super) struct SolverSettings {
    /// Number of velocity iterations (and position iterations).
    pub iterations: u32,
    /// Size of the simulation step the solved velocities will be integrated
    /// over. It is only used to rescale cached impulses when the step size
    /// changes between steps.
    pub step_dt: Real,
}

// ---------------------------------------------------------------------------
// Solver-side bodies
// ---------------------------------------------------------------------------

/// Solver-side copy of a rigid body.
#[derive(Debug, Clone)]
struct SolverBody {
    /// Body this slot mirrors; `None` for the shared static world slot.
    handle: Option<BodyHandle>,
    /// Centre of mass in world space.
    position: Vec3,
    /// Inverse mass (zero for bodies that impulses must not move).
    inv_mass: Real,
    /// World-space inverse inertia tensor (zero for non-dynamic bodies).
    inv_inertia: Mat3,
    /// Linear velocity.
    velocity: Vec3,
    /// Angular velocity.
    angular_velocity: Vec3,
    /// Linear velocity at the start of the solve.
    start_velocity: Vec3,
    /// Angular velocity at the start of the solve.
    start_angular_velocity: Vec3,
    /// Pseudo translation accumulated by the position pass.
    shift: Vec3,
    /// Pseudo rotation (scaled axis) accumulated by the position pass.
    turn: Vec3,
    /// Whether the body is dynamic, i.e. can be changed by the solver.
    dynamic: bool,
    /// Whether the body was asleep when the solve started.
    sleeping: bool,
}

impl SolverBody {
    fn from_body(handle: BodyHandle, body: &RigidBody) -> Self {
        let dynamic = body.is_dynamic();
        let (velocity, angular_velocity) = if body.body_type == BodyType::Static {
            (Vec3::zeros(), Vec3::zeros())
        } else {
            (body.velocity, body.angular_velocity)
        };
        Self {
            handle: Some(handle),
            position: body.transform.position,
            inv_mass: if dynamic { body.inverse_mass } else { 0.0 },
            inv_inertia: if dynamic {
                body.world_inverse_inertia
            } else {
                Mat3::zeros()
            },
            velocity,
            angular_velocity,
            start_velocity: velocity,
            start_angular_velocity: angular_velocity,
            shift: Vec3::zeros(),
            turn: Vec3::zeros(),
            dynamic,
            sleeping: body.state == BodyState::Sleeping,
        }
    }

    /// Immovable stand-in for colliders without a (live) body.
    fn world() -> Self {
        Self {
            handle: None,
            position: Vec3::zeros(),
            inv_mass: 0.0,
            inv_inertia: Mat3::zeros(),
            velocity: Vec3::zeros(),
            angular_velocity: Vec3::zeros(),
            start_velocity: Vec3::zeros(),
            start_angular_velocity: Vec3::zeros(),
            shift: Vec3::zeros(),
            turn: Vec3::zeros(),
            dynamic: false,
            sleeping: false,
        }
    }

    /// Apply `sign * impulse` at the point `arm` away from the centre of mass.
    fn apply_impulse(&mut self, impulse: Vec3, arm: Vec3, sign: Real) {
        self.velocity += impulse * (self.inv_mass * sign);
        self.angular_velocity += self.inv_inertia * arm.cross(&impulse) * sign;
    }

    /// Apply `sign * impulse` to the pseudo motion of the position pass.
    fn apply_pseudo_impulse(&mut self, impulse: Vec3, arm: Vec3, sign: Real) {
        self.shift += impulse * (self.inv_mass * sign);
        self.turn += self.inv_inertia * arm.cross(&impulse) * sign;
    }

    /// Velocity of the material point `arm` away from the centre of mass.
    fn point_velocity(&self, arm: &Vec3) -> Vec3 {
        self.velocity + self.angular_velocity.cross(arm)
    }

    /// Pseudo displacement of the material point `arm` away from the centre.
    fn point_shift(&self, arm: &Vec3) -> Vec3 {
        self.shift + self.turn.cross(arm)
    }

    /// Whether the solve changed this body enough to wake it up.
    fn disturbed(&self) -> bool {
        (self.velocity - self.start_velocity).norm_squared() > WAKE_EPSILON * WAKE_EPSILON
            || (self.angular_velocity - self.start_angular_velocity).norm_squared()
                > WAKE_EPSILON * WAKE_EPSILON
            || self.shift.norm_squared() > WAKE_EPSILON * WAKE_EPSILON
            || self.turn.norm_squared() > WAKE_EPSILON * WAKE_EPSILON
    }
}

/// Solver bodies of one invocation, addressed by dense slot index.
#[derive(Debug, Default)]
struct BodyTable {
    slots: Vec<SolverBody>,
    index: HashMap<BodyHandle, usize>,
    world: Option<usize>,
}

impl BodyTable {
    /// Slot of the body `handle`; colliders without a live body share one
    /// immovable world slot.
    fn slot(&mut self, handle: Option<BodyHandle>, bodies: &RigidBodySet) -> usize {
        if let Some(h) = handle {
            if let Some(&idx) = self.index.get(&h) {
                return idx;
            }
            if let Some(body) = bodies.get(h) {
                let idx = self.slots.len();
                self.slots.push(SolverBody::from_body(h, body));
                self.index.insert(h, idx);
                return idx;
            }
        }
        match self.world {
            Some(idx) => idx,
            None => {
                let idx = self.slots.len();
                self.slots.push(SolverBody::world());
                self.world = Some(idx);
                idx
            }
        }
    }

    /// Copy the solved state back into the bodies.
    ///
    /// Awake bodies always take the solved velocities. A sleeping body is only
    /// touched (and woken up) if the solve actually disturbed it, so resting
    /// contacts do not keep sleeping bodies awake.
    fn write_back(&self, bodies: &mut RigidBodySet) {
        for slot in &self.slots {
            let Some(handle) = slot.handle else {
                continue;
            };
            if !slot.dynamic {
                continue;
            }
            let disturbed = slot.disturbed();
            if slot.sleeping && !disturbed {
                continue;
            }
            let Some(body) = bodies.get_mut(handle) else {
                continue;
            };
            body.velocity = slot.velocity;
            body.angular_velocity = slot.angular_velocity;
            if slot.shift.norm_squared() > 0.0 || slot.turn.norm_squared() > 0.0 {
                body.transform.position += slot.shift;
                body.transform.rotation =
                    Quat::from_scaled_axis(slot.turn) * body.transform.rotation;
                body.update_world_inertia();
            }
            if slot.sleeping {
                body.wake_up();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Constraints
// ---------------------------------------------------------------------------

/// Velocity constraint of one contact point.
#[derive(Debug, Clone)]
struct Constraint {
    /// Solver slot of body A.
    a: usize,
    /// Solver slot of body B.
    b: usize,
    /// Canonical key of the collider pair.
    key: PairKey,
    /// Whether `(a, b)` is the opposite of the canonical pair order.
    flipped: bool,
    /// World-space contact point.
    point: Vec3,
    /// Contact point relative to the centre of mass of A.
    arm_a: Vec3,
    /// Contact point relative to the centre of mass of B.
    arm_b: Vec3,
    /// Unit normal pointing from B toward A.
    normal: Vec3,
    /// Two unit tangents spanning the friction plane.
    tangents: [Vec3; 2],
    /// Effective mass along the normal.
    normal_mass: Real,
    /// Effective masses along the tangents.
    tangent_mass: [Real; 2],
    /// Separation speed the normal velocity is driven to (restitution).
    restitution_target: Real,
    /// Separating displacement the position pass drives the contact to.
    position_target: Real,
    /// Coulomb friction coefficient.
    friction: Real,
    /// Accumulated normal impulse.
    normal_impulse: Real,
    /// Accumulated friction impulses along the tangents.
    tangent_impulse: [Real; 2],
    /// Accumulated pseudo impulse of the position pass.
    position_impulse: Real,
}

/// Orthonormal tangent pair of a unit normal.
fn tangent_basis(normal: &Vec3) -> [Vec3; 2] {
    // At least one component of a unit vector has magnitude <= 1/sqrt(3), so
    // crossing with that axis can never degenerate.
    let reference = if normal.x.abs() <= 0.577_350_269 {
        Vec3::x()
    } else if normal.y.abs() <= 0.577_350_269 {
        Vec3::y()
    } else {
        Vec3::z()
    };
    let t1 = normal.cross(&reference).normalize();
    let t2 = normal.cross(&t1);
    [t1, t2]
}

/// Effective mass `1 / (1/m_a + 1/m_b + ...)` of a constraint along `dir`.
fn effective_mass(a: &SolverBody, b: &SolverBody, arm_a: &Vec3, arm_b: &Vec3, dir: &Vec3) -> Real {
    let torque_arm_a = arm_a.cross(dir);
    let torque_arm_b = arm_b.cross(dir);
    let inverse = a.inv_mass
        + b.inv_mass
        + torque_arm_a.dot(&(a.inv_inertia * torque_arm_a))
        + torque_arm_b.dot(&(b.inv_inertia * torque_arm_b));
    if inverse > MIN_INVERSE_EFFECTIVE_MASS {
        1.0 / inverse
    } else {
        0.0
    }
}

impl Constraint {
    /// Relative velocity of A with respect to B at the contact point.
    fn relative_velocity(&self, a: &SolverBody, b: &SolverBody) -> Vec3 {
        a.point_velocity(&self.arm_a) - b.point_velocity(&self.arm_b)
    }

    /// Apply the accumulated impulses (warm start).
    fn apply_accumulated(&self, slots: &mut [SolverBody]) {
        if self.normal_impulse == 0.0 && self.tangent_impulse == [0.0, 0.0] {
            return;
        }
        let impulse = self.normal * self.normal_impulse
            + self.tangents[0] * self.tangent_impulse[0]
            + self.tangents[1] * self.tangent_impulse[1];
        let Ok([a, b]) = slots.get_disjoint_mut([self.a, self.b]) else {
            return;
        };
        a.apply_impulse(impulse, self.arm_a, 1.0);
        b.apply_impulse(impulse, self.arm_b, -1.0);
    }

    /// One Coulomb friction iteration (friction cone clamp on the accumulated
    /// tangential impulse).
    fn solve_friction(&mut self, slots: &mut [SolverBody]) {
        let Ok([a, b]) = slots.get_disjoint_mut([self.a, self.b]) else {
            return;
        };
        let relative = self.relative_velocity(a, b);
        let mut candidate = [
            self.tangent_impulse[0] - self.tangent_mass[0] * relative.dot(&self.tangents[0]),
            self.tangent_impulse[1] - self.tangent_mass[1] * relative.dot(&self.tangents[1]),
        ];
        let limit = self.friction * self.normal_impulse;
        let magnitude = candidate[0].hypot(candidate[1]);
        if magnitude > limit {
            let scale = if magnitude > 0.0 {
                limit / magnitude
            } else {
                0.0
            };
            candidate[0] *= scale;
            candidate[1] *= scale;
        }
        let delta = [
            candidate[0] - self.tangent_impulse[0],
            candidate[1] - self.tangent_impulse[1],
        ];
        self.tangent_impulse = candidate;
        let impulse = self.tangents[0] * delta[0] + self.tangents[1] * delta[1];
        a.apply_impulse(impulse, self.arm_a, 1.0);
        b.apply_impulse(impulse, self.arm_b, -1.0);
    }

    /// One non-penetration iteration with a clamped accumulated impulse.
    fn solve_normal(&mut self, slots: &mut [SolverBody]) {
        let Ok([a, b]) = slots.get_disjoint_mut([self.a, self.b]) else {
            return;
        };
        let normal_speed = self.relative_velocity(a, b).dot(&self.normal);
        let increment = self.normal_mass * (self.restitution_target - normal_speed);
        let total = (self.normal_impulse + increment).max(0.0);
        let applied = total - self.normal_impulse;
        self.normal_impulse = total;
        let impulse = self.normal * applied;
        a.apply_impulse(impulse, self.arm_a, 1.0);
        b.apply_impulse(impulse, self.arm_b, -1.0);
    }

    /// One iteration of the split-impulse position pass.
    fn solve_position(&mut self, slots: &mut [SolverBody]) {
        let Ok([a, b]) = slots.get_disjoint_mut([self.a, self.b]) else {
            return;
        };
        let separation =
            (a.point_shift(&self.arm_a) - b.point_shift(&self.arm_b)).dot(&self.normal);
        let increment = self.normal_mass * (self.position_target - separation);
        let total = (self.position_impulse + increment).max(0.0);
        let applied = total - self.position_impulse;
        self.position_impulse = total;
        let impulse = self.normal * applied;
        a.apply_pseudo_impulse(impulse, self.arm_a, 1.0);
        b.apply_pseudo_impulse(impulse, self.arm_b, -1.0);
    }
}

/// Factor that converts impulses accumulated over `cached_dt` into impulses
/// for a step of `dt`.
fn warm_start_ratio(cached_dt: Real, dt: Real) -> Real {
    let ratio = dt / cached_dt;
    if cached_dt > 0.0 && ratio.is_finite() && ratio > 0.0 {
        ratio
    } else {
        1.0
    }
}

/// Build the constraints of all contacts, and fill the body table.
fn build_constraints(
    pairs: &[ContactPair<'_>],
    bodies: &RigidBodySet,
    cache: &WarmStartCache,
    ratio: Real,
    table: &mut BodyTable,
) -> Vec<Constraint> {
    let mut constraints = Vec::new();
    let mut taken: HashMap<PairKey, Vec<bool>> = HashMap::new();

    for pair in pairs {
        let a = table.slot(pair.body_a, bodies);
        let b = table.slot(pair.body_b, bodies);
        if a == b {
            continue;
        }
        let (key, flipped) = PairKey::new(pair.handle_a, pair.handle_b);
        let (Some(body_a), Some(body_b)) = (table.slots.get(a), table.slots.get(b)) else {
            continue;
        };

        for contact in &pair.manifold.contacts {
            let Some(normal) = contact.normal.try_normalize(1e-12) else {
                continue;
            };
            let point = contact.point();
            let arm_a = point - body_a.position;
            let arm_b = point - body_b.position;
            let normal_mass = effective_mass(body_a, body_b, &arm_a, &arm_b, &normal);
            if normal_mass == 0.0 {
                continue;
            }
            let tangents = tangent_basis(&normal);
            let tangent_mass = [
                effective_mass(body_a, body_b, &arm_a, &arm_b, &tangents[0]),
                effective_mass(body_a, body_b, &arm_a, &arm_b, &tangents[1]),
            ];

            // Restitution target from the approach speed *before* any warm
            // start impulse changes the velocities.
            let approach =
                (body_a.point_velocity(&arm_a) - body_b.point_velocity(&arm_b)).dot(&normal);
            let restitution_target = if approach < -RESTITUTION_THRESHOLD {
                -pair.restitution.max(0.0) * approach
            } else {
                0.0
            };
            let position_target = POSITION_CORRECTION_RATE
                * (contact.depth - PENETRATION_SLOP).clamp(0.0, MAX_POSITION_CORRECTION);

            let friction = pair.friction.max(0.0);
            let mut constraint = Constraint {
                a,
                b,
                key,
                flipped,
                point,
                arm_a,
                arm_b,
                normal,
                tangents,
                normal_mass,
                tangent_mass,
                restitution_target,
                position_target,
                friction,
                normal_impulse: 0.0,
                tangent_impulse: [0.0, 0.0],
                position_impulse: 0.0,
            };

            let canonical_normal = if flipped { -normal } else { normal };
            let flags = taken
                .entry(key)
                .or_insert_with(|| vec![false; cache.contacts.get(&key).map_or(0, Vec::len)]);
            if let Some(cached) = cache.find_match(&key, &point, &canonical_normal, flags) {
                constraint.normal_impulse = cached.normal_impulse * ratio;
                let friction_impulse = if flipped {
                    -cached.friction_impulse
                } else {
                    cached.friction_impulse
                } * ratio;
                let mut tangent_impulse = [
                    friction_impulse.dot(&constraint.tangents[0]),
                    friction_impulse.dot(&constraint.tangents[1]),
                ];
                let limit = constraint.friction * constraint.normal_impulse;
                let magnitude = tangent_impulse[0].hypot(tangent_impulse[1]);
                if magnitude > limit {
                    let scale = if magnitude > 0.0 {
                        limit / magnitude
                    } else {
                        0.0
                    };
                    tangent_impulse[0] *= scale;
                    tangent_impulse[1] *= scale;
                }
                constraint.tangent_impulse = tangent_impulse;
            }
            constraints.push(constraint);
        }
    }
    constraints
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

/// Solve the contacts of one step.
///
/// Updates the velocities of the bodies in `bodies`, moves overlapping bodies
/// apart by split-impulse position correction, wakes sleeping bodies that were
/// disturbed, and replaces the contents of `cache` by the accumulated impulses
/// of this step.
pub(super) fn solve_contacts(
    cache: &mut WarmStartCache,
    pairs: &[ContactPair<'_>],
    bodies: &mut RigidBodySet,
    settings: SolverSettings,
) {
    let mut table = BodyTable::default();
    let ratio = warm_start_ratio(cache.step_dt, settings.step_dt);
    let mut constraints = build_constraints(pairs, bodies, cache, ratio, &mut table);

    // Warm start every constraint only after all restitution targets have been
    // measured, so the result does not depend on the order of the contacts.
    for constraint in &constraints {
        constraint.apply_accumulated(&mut table.slots);
    }

    let iterations = settings.iterations.max(1);
    for _ in 0..iterations {
        for constraint in constraints.iter_mut() {
            constraint.solve_friction(&mut table.slots);
            constraint.solve_normal(&mut table.slots);
        }
    }
    for _ in 0..iterations {
        for constraint in constraints.iter_mut() {
            constraint.solve_position(&mut table.slots);
        }
    }

    table.write_back(bodies);

    let mut contacts: HashMap<PairKey, Vec<CachedContact>> = HashMap::new();
    for constraint in &constraints {
        let friction_impulse = constraint.tangents[0] * constraint.tangent_impulse[0]
            + constraint.tangents[1] * constraint.tangent_impulse[1];
        let (normal, friction_impulse) = if constraint.flipped {
            (-constraint.normal, -friction_impulse)
        } else {
            (constraint.normal, friction_impulse)
        };
        contacts
            .entry(constraint.key)
            .or_default()
            .push(CachedContact {
                point: constraint.point,
                normal,
                normal_impulse: constraint.normal_impulse,
                friction_impulse,
            });
    }
    cache.contacts = contacts;
    cache.step_dt = settings.step_dt;
}

/// A single impact to be resolved at one contact point.
#[derive(Debug, Clone, Copy)]
pub(super) struct Impact {
    /// Body that carries collider A, if any.
    pub body_a: Option<BodyHandle>,
    /// Body that carries collider B, if any.
    pub body_b: Option<BodyHandle>,
    /// World-space contact point.
    pub point: Vec3,
    /// Unit contact normal pointing from B toward A.
    pub normal: Vec3,
    /// Combined coefficient of restitution.
    pub restitution: Real,
}

/// Resolve a single impact with one normal impulse (used at a time of impact
/// found by continuous collision detection).
///
/// Returns `true` if bodies were approaching and an impulse was applied.
pub(super) fn resolve_impact(bodies: &mut RigidBodySet, impact: &Impact) -> bool {
    let Some(normal) = impact.normal.try_normalize(1e-12) else {
        return false;
    };
    let mut table = BodyTable::default();
    let a = table.slot(impact.body_a, bodies);
    let b = table.slot(impact.body_b, bodies);
    if a == b {
        return false;
    }
    let applied = match table.slots.get_disjoint_mut([a, b]) {
        Ok([body_a, body_b]) => {
            let arm_a = impact.point - body_a.position;
            let arm_b = impact.point - body_b.position;
            let mass = effective_mass(body_a, body_b, &arm_a, &arm_b, &normal);
            let approach =
                (body_a.point_velocity(&arm_a) - body_b.point_velocity(&arm_b)).dot(&normal);
            if mass == 0.0 || approach >= 0.0 {
                false
            } else {
                let magnitude = -(1.0 + impact.restitution.max(0.0)) * approach * mass;
                let impulse = normal * magnitude;
                body_a.apply_impulse(impulse, arm_a, 1.0);
                body_b.apply_impulse(impulse, arm_b, -1.0);
                true
            }
        }
        Err(_) => false,
    };
    if applied {
        table.write_back(bodies);
    }
    applied
}

#[cfg(test)]
mod tests {
    use super::*;

    fn handle(index: u32, generation: u32) -> ColliderHandle {
        ColliderHandle::new(index, generation)
    }

    #[test]
    fn pair_key_is_independent_of_argument_order() {
        let a = handle(3, 0);
        let b = handle(7, 2);
        let (key_ab, flipped_ab) = PairKey::new(a, b);
        let (key_ba, flipped_ba) = PairKey::new(b, a);
        assert_eq!(key_ab, key_ba);
        assert!(!flipped_ab);
        assert!(flipped_ba);
        // Same slot, different generation: the generation breaks the tie.
        let (_, flipped) = PairKey::new(handle(3, 1), handle(3, 0));
        assert!(flipped);
    }

    #[test]
    fn tangent_basis_is_orthonormal_for_any_normal() {
        let normals = [
            Vec3::x(),
            Vec3::y(),
            Vec3::z(),
            -Vec3::y(),
            Vec3::new(1.0, 1.0, 1.0).normalize(),
            Vec3::new(-0.3, 0.9, 0.1).normalize(),
            Vec3::new(0.577, 0.577, 0.5778).normalize(),
        ];
        for normal in normals {
            let [t1, t2] = tangent_basis(&normal);
            assert!((t1.norm() - 1.0).abs() < 1e-12, "{normal:?}");
            assert!((t2.norm() - 1.0).abs() < 1e-12, "{normal:?}");
            assert!(t1.dot(&normal).abs() < 1e-12, "{normal:?}");
            assert!(t2.dot(&normal).abs() < 1e-12, "{normal:?}");
            assert!(t1.dot(&t2).abs() < 1e-12, "{normal:?}");
        }
    }

    #[test]
    fn warm_start_ratio_rescales_and_rejects_nonsense() {
        assert_eq!(warm_start_ratio(1.0 / 60.0, 1.0 / 60.0), 1.0);
        assert!((warm_start_ratio(1.0 / 60.0, 1.0 / 120.0) - 0.5).abs() < 1e-12);
        assert!((warm_start_ratio(1.0 / 120.0, 1.0 / 60.0) - 2.0).abs() < 1e-12);
        // No cached step (empty cache) or a degenerate step: no rescaling.
        assert_eq!(warm_start_ratio(0.0, 1.0 / 60.0), 1.0);
        assert_eq!(warm_start_ratio(1.0 / 60.0, 0.0), 1.0);
        assert_eq!(warm_start_ratio(1.0 / 60.0, f64::NAN), 1.0);
    }

    fn cached(point: Vec3, impulse: Real) -> CachedContact {
        CachedContact {
            point,
            normal: Vec3::y(),
            normal_impulse: impulse,
            friction_impulse: Vec3::zeros(),
        }
    }

    #[test]
    fn find_match_prefers_the_nearest_unused_contact_in_range() {
        let (key, _) = PairKey::new(handle(0, 0), handle(1, 0));
        let mut cache = WarmStartCache::new();
        cache.contacts.insert(
            key,
            vec![
                cached(Vec3::new(0.0, 0.0, 0.0), 1.0),
                cached(Vec3::new(0.02, 0.0, 0.0), 2.0),
                cached(Vec3::new(1.0, 0.0, 0.0), 3.0),
            ],
        );
        let mut taken = vec![false; 3];

        // Nearest to x = 0.019 is the second contact.
        let m = cache
            .find_match(&key, &Vec3::new(0.019, 0.0, 0.0), &Vec3::y(), &mut taken)
            .map(|c| c.normal_impulse);
        assert_eq!(m, Some(2.0));
        // It is now taken: the same query falls back to the first one.
        let m = cache
            .find_match(&key, &Vec3::new(0.019, 0.0, 0.0), &Vec3::y(), &mut taken)
            .map(|c| c.normal_impulse);
        assert_eq!(m, Some(1.0));
        // Nothing within range of x = 0.5.
        assert!(
            cache
                .find_match(&key, &Vec3::new(0.5, 0.0, 0.0), &Vec3::y(), &mut taken)
                .is_none()
        );
        // The far contact matches a nearby query.
        let m = cache
            .find_match(&key, &Vec3::new(1.01, 0.0, 0.0), &Vec3::y(), &mut taken)
            .map(|c| c.normal_impulse);
        assert_eq!(m, Some(3.0));
    }

    #[test]
    fn find_match_rejects_contacts_with_a_different_normal() {
        let (key, _) = PairKey::new(handle(0, 0), handle(1, 0));
        let mut cache = WarmStartCache::new();
        cache.contacts.insert(key, vec![cached(Vec3::zeros(), 1.0)]);
        let mut taken = vec![false];
        // Same place, but the contact now faces sideways.
        assert!(
            cache
                .find_match(&key, &Vec3::zeros(), &Vec3::x(), &mut taken)
                .is_none()
        );
        assert_eq!(taken, vec![false]);
    }

    #[test]
    fn cache_accessors_sum_over_pairs_and_ignore_argument_order() {
        let mut cache = WarmStartCache::new();
        assert!(cache.is_empty());
        let a = handle(0, 0);
        let b = handle(1, 0);
        let c = handle(2, 0);
        let (key_ab, _) = PairKey::new(a, b);
        let (key_ac, _) = PairKey::new(a, c);
        cache.contacts.insert(
            key_ab,
            vec![cached(Vec3::zeros(), 1.0), cached(Vec3::x(), 2.0)],
        );
        cache.contacts.insert(key_ac, vec![cached(Vec3::y(), 4.0)]);

        assert!(!cache.is_empty());
        assert_eq!(cache.pair_count(), 2);
        assert_eq!(cache.contact_count(), 3);
        assert_eq!(cache.normal_impulse(a, b), 3.0);
        assert_eq!(cache.normal_impulse(b, a), 3.0);
        assert_eq!(cache.normal_impulse(b, c), 0.0);
        assert_eq!(cache.contact_normal_impulses(c, a), vec![4.0]);
        assert_eq!(cache.normal_impulse_total(), 7.0);

        cache.clear();
        assert!(cache.is_empty());
        assert_eq!(cache.normal_impulse_total(), 0.0);
    }

    fn floor_and_faller(velocity: Vec3) -> (RigidBodySet, BodyHandle, BodyHandle) {
        let mut bodies = RigidBodySet::new();
        let floor = bodies.insert(RigidBody::new_static());
        let mut faller = RigidBody::new(2.0);
        faller.transform.position = Vec3::new(0.0, 0.5, 0.0);
        faller.velocity = velocity;
        let faller = bodies.insert(faller);
        (bodies, floor, faller)
    }

    #[test]
    fn resolve_impact_reflects_the_approach_speed_by_the_restitution() {
        let (mut bodies, floor, faller) = floor_and_faller(Vec3::new(0.0, -5.0, 0.0));
        let impact = Impact {
            body_a: Some(faller),
            body_b: Some(floor),
            point: Vec3::new(0.0, 0.0, 0.0),
            normal: Vec3::y(),
            restitution: 0.5,
        };
        assert!(resolve_impact(&mut bodies, &impact));
        let v = bodies.get(faller).map(|b| b.velocity).unwrap_or_default();
        assert!((v.y - 2.5).abs() < 1e-12, "rebound speed {}", v.y);
        assert_eq!(bodies.get(floor).map(|b| b.velocity), Some(Vec3::zeros()));
    }

    #[test]
    fn resolve_impact_ignores_separating_bodies_and_degenerate_normals() {
        let (mut bodies, floor, faller) = floor_and_faller(Vec3::new(0.0, 3.0, 0.0));
        let mut impact = Impact {
            body_a: Some(faller),
            body_b: Some(floor),
            point: Vec3::zeros(),
            normal: Vec3::y(),
            restitution: 1.0,
        };
        assert!(!resolve_impact(&mut bodies, &impact), "already separating");
        impact.normal = Vec3::zeros();
        assert!(!resolve_impact(&mut bodies, &impact), "no normal");
        impact.normal = Vec3::y();
        impact.body_b = Some(faller);
        assert!(!resolve_impact(&mut bodies, &impact), "a body with itself");
        let v = bodies.get(faller).map(|b| b.velocity).unwrap_or_default();
        assert_eq!(v, Vec3::new(0.0, 3.0, 0.0));
    }

    #[test]
    fn resolve_impact_between_two_dynamic_bodies_conserves_momentum() {
        let mut bodies = RigidBodySet::new();
        let mut light = RigidBody::new(1.0);
        light.velocity = Vec3::new(-4.0, 0.0, 0.0);
        let light = bodies.insert(light);
        let mut heavy = RigidBody::new(3.0);
        heavy.velocity = Vec3::new(1.0, 0.0, 0.0);
        heavy.transform.position = Vec3::new(-1.0, 0.0, 0.0);
        let heavy = bodies.insert(heavy);
        // A = light on the right moving left, B = heavy on the left moving
        // right; the normal points from B toward A (+x).
        let impact = Impact {
            body_a: Some(light),
            body_b: Some(heavy),
            point: Vec3::new(-0.5, 0.0, 0.0),
            normal: Vec3::x(),
            restitution: 0.0,
        };
        let momentum_before = -4.0 + 3.0 * 1.0;
        assert!(resolve_impact(&mut bodies, &impact));
        let v_light = bodies.get(light).map(|b| b.velocity.x).unwrap_or_default();
        let v_heavy = bodies.get(heavy).map(|b| b.velocity.x).unwrap_or_default();
        assert!((v_light + 3.0 * v_heavy - momentum_before).abs() < 1e-12);
        assert!((v_light - v_heavy).abs() < 1e-12, "perfectly inelastic");
    }
}
