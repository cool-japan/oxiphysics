// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact constraint solvers for [`PhysicsWorld`].
//!
//! Sign convention (shared with [`ContactPair`]): the contact normal `n` points
//! FROM body B TOWARD body A, so the relative normal velocity
//! `v_n = n · (v_A(p) − v_B(p))` is negative while the bodies approach and
//! positive while they separate. A normal impulse `λ n` is applied as `+λ n` to
//! A and `−λ n` to B and is always non-negative (contacts push, never pull).
//!
//! The default [`PhysicsWorld::step`] path is a sequential-impulse solver with
//! accumulated (clamped) normal and Coulomb-friction impulses using the full
//! linear + angular Jacobian, a restitution velocity target computed once from
//! the pre-solve approach speed, and a non-linear Gauss–Seidel position pass
//! (Baumgarte-scaled, energy-free) after integration.

use super::contact_detection::detect_contacts;
use super::contact_lcp::{MAX_BLOCK, solve_lcp};
use super::types::{ContactPair, PhysicsWorld};
use crate::body::{BodyType, RigidBody};
use crate::impulse::SubStepRestitutionCorrector;
use crate::sets::RigidBodySet;
use crate::solver::soft::SoftParams;
use oxiphysics_core::math::{Mat3, Quat, Real, Vec3};
use oxiphysics_core::{BodyHandle, ColliderHandle};
use std::collections::HashMap;

/// Largest positional correction applied to one contact per position
/// iteration (metres). Prevents large overlaps from launching bodies.
const MAX_LINEAR_CORRECTION: Real = 0.2;

/// Effective masses below this are treated as "no response possible".
const MIN_EFFECTIVE_MASS: Real = 1e-30;

/// Two contact points (of the same collider pair) whose body-local anchors on
/// both bodies are closer than this are treated as the same point when
/// warm-starting (metres).
const WARM_START_MATCH_DISTANCE: Real = 0.02;

/// A cached impulse is only reused while the contact normal stays within
/// ~18° of its previous direction.
const WARM_START_MIN_NORMAL_DOT: Real = 0.95;

/// Points of one manifold are block-solved only when their normals agree.
const BLOCK_NORMAL_DOT: Real = 0.999;

/// Orthonormal tangent basis perpendicular to the unit normal `n`.
pub(super) fn tangent_basis(n: &Vec3) -> [Vec3; 2] {
    let helper = if n.x.abs() < 0.577_35 {
        Vec3::x()
    } else if n.y.abs() < 0.577_35 {
        Vec3::y()
    } else {
        Vec3::z()
    };
    let t1 = n.cross(&helper).normalize();
    let t2 = n.cross(&t1);
    [t1, t2]
}

fn vec3(a: [Real; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

/// Mass terms of the two bodies of a contact at a given contact point.
#[derive(Debug, Clone)]
pub(super) struct PairMass {
    inv_ma: Real,
    inv_mb: Real,
    inv_ia: Mat3,
    inv_ib: Mat3,
    r_a: Vec3,
    r_b: Vec3,
}

impl PairMass {
    /// Inverse effective mass `J M⁻¹ Jᵀ` along the unit direction `d`.
    fn inv_effective_mass(&self, d: &Vec3) -> Real {
        let ra_d = self.r_a.cross(d);
        let rb_d = self.r_b.cross(d);
        self.inv_ma
            + self.inv_mb
            + ra_d.dot(&(self.inv_ia * ra_d))
            + rb_d.dot(&(self.inv_ib * rb_d))
    }
}

/// Mass terms of a body (zero inverse mass / inertia unless it is dynamic).
fn inverse_mass_terms(body: &RigidBody) -> (Real, Mat3) {
    if body.body_type == BodyType::Dynamic {
        (body.inverse_mass, body.world_inverse_inertia)
    } else {
        (0.0, Mat3::zeros())
    }
}

/// Build the mass terms for a contact at world point `point`.
fn pair_mass(bodies: &RigidBodySet, a: BodyHandle, b: BodyHandle, point: Vec3) -> Option<PairMass> {
    let body_a = bodies.get(a)?;
    let body_b = bodies.get(b)?;
    let (inv_ma, inv_ia) = inverse_mass_terms(body_a);
    let (inv_mb, inv_ib) = inverse_mass_terms(body_b);
    Some(PairMass {
        inv_ma,
        inv_mb,
        inv_ia,
        inv_ib,
        r_a: point - body_a.transform.position,
        r_b: point - body_b.transform.position,
    })
}

/// Relative velocity `v_A(p) − v_B(p)` at the contact point.
fn relative_velocity(bodies: &RigidBodySet, a: BodyHandle, b: BodyHandle, m: &PairMass) -> Vec3 {
    let point_velocity = |h: BodyHandle, r: &Vec3| -> Vec3 {
        bodies
            .get(h)
            .map(|body| body.velocity + body.angular_velocity.cross(r))
            .unwrap_or_else(Vec3::zeros)
    };
    point_velocity(a, &m.r_a) - point_velocity(b, &m.r_b)
}

/// Apply `+p` to A and `−p` to B at the contact offsets (dynamic bodies only).
fn apply_impulse(bodies: &mut RigidBodySet, a: BodyHandle, b: BodyHandle, m: &PairMass, p: Vec3) {
    if m.inv_ma > 0.0
        && let Some(body) = bodies.get_mut(a)
    {
        body.velocity += p * m.inv_ma;
        body.angular_velocity += m.inv_ia * m.r_a.cross(&p);
    }
    if m.inv_mb > 0.0
        && let Some(body) = bodies.get_mut(b)
    {
        body.velocity -= p * m.inv_mb;
        body.angular_velocity -= m.inv_ib * m.r_b.cross(&p);
    }
}

/// Rotate a body by the small rotation vector `dtheta` (world frame).
fn rotate_body(body: &mut RigidBody, dtheta: Vec3) {
    if dtheta.norm_squared() > 1e-24 {
        let dq = Quat::from_scaled_axis(dtheta);
        body.transform.rotation = dq * body.transform.rotation;
        body.update_world_inertia();
    }
}

/// World position of a body-local anchor.
fn world_anchor(bodies: &RigidBodySet, h: BodyHandle, local: &Vec3) -> Option<Vec3> {
    let body = bodies.get(h)?;
    Some(body.transform.transform_point(local))
}

/// Body-local coordinates of the world point `p`.
fn local_anchor(body: &RigidBody, p: &Vec3) -> Vec3 {
    body.transform.rotation.inverse() * (p - body.transform.position)
}

/// Accumulated Coulomb friction for one contact (two tangent directions, each
/// clamped to `±μ λ_n`), shared by the soft solver paths.
pub(super) fn solve_friction_accumulated(
    bodies: &mut RigidBodySet,
    cp: &ContactPair,
    lambda_n: Real,
    lambda_t: &mut [Real; 2],
) {
    if cp.friction <= 0.0 {
        return;
    }
    let point = vec3(cp.contact_point);
    let Some(m) = pair_mass(bodies, cp.body_a, cp.body_b, point) else {
        return;
    };
    let n = vec3(cp.normal);
    let max_friction = cp.friction * lambda_n.max(0.0);
    for (k, t) in tangent_basis(&n).iter().enumerate() {
        let k_t = m.inv_effective_mass(t);
        if k_t < MIN_EFFECTIVE_MASS {
            continue;
        }
        let vt = relative_velocity(bodies, cp.body_a, cp.body_b, &m).dot(t);
        let old = lambda_t[k];
        let new = (old - vt / k_t).clamp(-max_friction, max_friction);
        lambda_t[k] = new;
        apply_impulse(bodies, cp.body_a, cp.body_b, &m, *t * (new - old));
    }
}

// ── Default (rigid) sequential-impulse solver ─────────────────────────────────

/// One prepared contact constraint of the default solver.
#[derive(Debug, Clone)]
pub(super) struct ContactConstraint {
    body_a: BodyHandle,
    body_b: BodyHandle,
    collider_a: ColliderHandle,
    collider_b: ColliderHandle,
    mass: PairMass,
    normal: Vec3,
    tangents: [Vec3; 2],
    normal_mass: Real,
    tangent_mass: [Real; 2],
    friction: Real,
    /// Target separating velocity (restitution bounce, or the allowed
    /// approach speed `depth / dt` of a speculative contact).
    target: Real,
    /// Detected penetration depth (negative for speculative gaps).
    depth: Real,
    local_anchor_a: Vec3,
    local_anchor_b: Vec3,
    lambda_n: Real,
    lambda_t: [Real; 2],
}

impl ContactConstraint {
    fn prepare(
        cp: &ContactPair,
        bodies: &RigidBodySet,
        restitution_threshold: Real,
        dt: Real,
    ) -> Option<Self> {
        let point = vec3(cp.contact_point);
        let normal = vec3(cp.normal);
        let mass = pair_mass(bodies, cp.body_a, cp.body_b, point)?;
        let k_n = mass.inv_effective_mass(&normal);
        if k_n < MIN_EFFECTIVE_MASS {
            return None;
        }
        let tangents = tangent_basis(&normal);
        let tangent_mass = [
            inv_or_zero(mass.inv_effective_mass(&tangents[0])),
            inv_or_zero(mass.inv_effective_mass(&tangents[1])),
        ];
        let vn = relative_velocity(bodies, cp.body_a, cp.body_b, &mass).dot(&normal);
        let target = if cp.depth < 0.0 {
            // Speculative: allow closing exactly the gap during this step.
            if dt > 0.0 { cp.depth / dt } else { 0.0 }
        } else if vn < -restitution_threshold {
            -cp.restitution * vn
        } else {
            0.0
        };
        let body_a = bodies.get(cp.body_a)?;
        let body_b = bodies.get(cp.body_b)?;
        Some(Self {
            body_a: cp.body_a,
            body_b: cp.body_b,
            collider_a: cp.collider_a,
            collider_b: cp.collider_b,
            mass,
            normal,
            tangents,
            normal_mass: 1.0 / k_n,
            tangent_mass,
            friction: cp.friction,
            target,
            depth: cp.depth,
            local_anchor_a: local_anchor(body_a, &point),
            local_anchor_b: local_anchor(body_b, &point),
            lambda_n: 0.0,
            lambda_t: [0.0; 2],
        })
    }

    /// Seed the accumulated impulses from the previous step and apply them.
    fn warm_start(&mut self, bodies: &mut RigidBodySet, previous: &CachedImpulse) {
        if self.depth < 0.0 || previous.normal.dot(&self.normal) < WARM_START_MIN_NORMAL_DOT {
            return;
        }
        self.lambda_n = previous.lambda_n.max(0.0);
        let max_friction = self.friction * self.lambda_n;
        for k in 0..2 {
            self.lambda_t[k] = previous
                .tangent_impulse
                .dot(&self.tangents[k])
                .clamp(-max_friction, max_friction);
        }
        let p = self.normal * self.lambda_n
            + self.tangents[0] * self.lambda_t[0]
            + self.tangents[1] * self.lambda_t[1];
        apply_impulse(bodies, self.body_a, self.body_b, &self.mass, p);
    }

    /// Snapshot of the converged impulses for warm-starting the next step.
    fn cached_impulse(&self) -> CachedImpulse {
        CachedImpulse {
            local_anchor_a: self.local_anchor_a,
            local_anchor_b: self.local_anchor_b,
            normal: self.normal,
            lambda_n: self.lambda_n,
            tangent_impulse: self.tangents[0] * self.lambda_t[0]
                + self.tangents[1] * self.lambda_t[1],
        }
    }

    /// One Gauss–Seidel friction iteration (two tangents, `|λ_t| ≤ μ λ_n`).
    fn solve_friction(&mut self, bodies: &mut RigidBodySet) {
        let max_friction = self.friction * self.lambda_n;
        for k in 0..2 {
            if self.tangent_mass[k] <= 0.0 {
                continue;
            }
            let t = self.tangents[k];
            let vt = relative_velocity(bodies, self.body_a, self.body_b, &self.mass).dot(&t);
            let old = self.lambda_t[k];
            let new = (old - vt * self.tangent_mass[k]).clamp(-max_friction, max_friction);
            self.lambda_t[k] = new;
            apply_impulse(
                bodies,
                self.body_a,
                self.body_b,
                &self.mass,
                t * (new - old),
            );
        }
    }

    /// One Gauss–Seidel normal iteration (accumulated impulse clamped ≥ 0).
    fn solve_normal(&mut self, bodies: &mut RigidBodySet) {
        let vn = self.normal_velocity(bodies);
        let old = self.lambda_n;
        let new = (old - self.normal_mass * (vn - self.target)).max(0.0);
        self.lambda_n = new;
        apply_impulse(
            bodies,
            self.body_a,
            self.body_b,
            &self.mass,
            self.normal * (new - old),
        );
    }

    /// Relative normal velocity at the contact point.
    fn normal_velocity(&self, bodies: &RigidBodySet) -> Real {
        relative_velocity(bodies, self.body_a, self.body_b, &self.mass).dot(&self.normal)
    }

    /// Whether `other` is a further point of the same manifold (same collider
    /// pair, same normal, both penetrating) and can be block-solved with it.
    fn same_manifold(&self, other: &ContactConstraint) -> bool {
        self.collider_a == other.collider_a
            && self.collider_b == other.collider_b
            && self.body_a == other.body_a
            && self.body_b == other.body_b
            && self.depth >= 0.0
            && other.depth >= 0.0
            && self.normal.dot(&other.normal) > BLOCK_NORMAL_DOT
    }

    /// Current penetration re-projected from the body-local anchors.
    fn current_depth(&self, bodies: &RigidBodySet) -> Option<(Real, Vec3, Vec3)> {
        let pa = world_anchor(bodies, self.body_a, &self.local_anchor_a)?;
        let pb = world_anchor(bodies, self.body_b, &self.local_anchor_b)?;
        Some((self.depth - self.normal.dot(&(pa - pb)), pa, pb))
    }

    /// One non-linear Gauss–Seidel position iteration: pushes the bodies apart
    /// until the penetration is at most `slop`, without touching velocities.
    fn solve_position(&self, bodies: &mut RigidBodySet, beta: Real, slop: Real) {
        if self.depth < 0.0 {
            return;
        }
        let Some((depth, pa, pb)) = self.current_depth(bodies) else {
            return;
        };
        let correction = (beta * (depth - slop)).clamp(0.0, MAX_LINEAR_CORRECTION);
        if correction <= 0.0 {
            return;
        }
        let (Some(body_a), Some(body_b)) = (bodies.get(self.body_a), bodies.get(self.body_b))
        else {
            return;
        };
        let (inv_ma, inv_ia) = inverse_mass_terms(body_a);
        let (inv_mb, inv_ib) = inverse_mass_terms(body_b);
        let mass = PairMass {
            inv_ma,
            inv_mb,
            inv_ia,
            inv_ib,
            r_a: pa - body_a.transform.position,
            r_b: pb - body_b.transform.position,
        };
        let k = mass.inv_effective_mass(&self.normal);
        if k < MIN_EFFECTIVE_MASS {
            return;
        }
        let p = self.normal * (correction / k);
        if mass.inv_ma > 0.0
            && let Some(body) = bodies.get_mut(self.body_a)
        {
            body.transform.position += p * mass.inv_ma;
            rotate_body(body, mass.inv_ia * mass.r_a.cross(&p));
        }
        if mass.inv_mb > 0.0
            && let Some(body) = bodies.get_mut(self.body_b)
        {
            body.transform.position -= p * mass.inv_mb;
            rotate_body(body, -(mass.inv_ib * mass.r_b.cross(&p)));
        }
    }
}

fn inv_or_zero(k: Real) -> Real {
    if k < MIN_EFFECTIVE_MASS { 0.0 } else { 1.0 / k }
}

// ── Manifold blocks ───────────────────────────────────────────────────────────

/// A run of constraints forming one manifold, solved together.
#[derive(Debug, Clone)]
struct ManifoldBlock {
    start: usize,
    len: usize,
    /// Normal effective-mass matrix `K_ij = J_i M⁻¹ J_jᵀ` of the block.
    k: [[Real; MAX_BLOCK]; MAX_BLOCK],
}

/// Group consecutive constraints of the same manifold into blocks of at most
/// [`MAX_BLOCK`] points and precompute their normal effective-mass matrices.
fn build_blocks(constraints: &[ContactConstraint]) -> Vec<ManifoldBlock> {
    let mut blocks = Vec::new();
    let mut start = 0;
    while start < constraints.len() {
        let mut len = 1;
        while start + len < constraints.len()
            && len < MAX_BLOCK
            && constraints[start].same_manifold(&constraints[start + len])
        {
            len += 1;
        }
        let mut k = [[0.0; MAX_BLOCK]; MAX_BLOCK];
        if len > 1 {
            for i in 0..len {
                let ci = &constraints[start + i];
                let ra_i = ci.mass.r_a.cross(&ci.normal);
                let rb_i = ci.mass.r_b.cross(&ci.normal);
                for j in 0..len {
                    let cj = &constraints[start + j];
                    let ra_j = cj.mass.r_a.cross(&cj.normal);
                    let rb_j = cj.mass.r_b.cross(&cj.normal);
                    k[i][j] = ci.normal.dot(&cj.normal) * (ci.mass.inv_ma + ci.mass.inv_mb)
                        + ra_i.dot(&(ci.mass.inv_ia * ra_j))
                        + rb_i.dot(&(ci.mass.inv_ib * rb_j));
                }
            }
        }
        blocks.push(ManifoldBlock { start, len, k });
        start += len;
    }
    blocks
}

/// One iteration over a manifold block: per-point friction, then an exact
/// block solve of the normal impulses (falling back to per-point
/// Gauss–Seidel if no feasible active set is found).
fn solve_block(
    block: &ManifoldBlock,
    constraints: &mut [ContactConstraint],
    bodies: &mut RigidBodySet,
) {
    let points = &mut constraints[block.start..block.start + block.len];
    for c in points.iter_mut() {
        c.solve_friction(bodies);
    }
    if block.len == 1 {
        points[0].solve_normal(bodies);
        return;
    }
    let n = block.len;
    let mut b = [0.0; MAX_BLOCK];
    for i in 0..n {
        let accumulated: Real = (0..n).map(|j| block.k[i][j] * points[j].lambda_n).sum();
        b[i] = points[i].normal_velocity(bodies) - points[i].target - accumulated;
    }
    match solve_lcp(&block.k, &b, n) {
        Some(x) => {
            for (i, c) in points.iter_mut().enumerate() {
                let delta = x[i] - c.lambda_n;
                c.lambda_n = x[i];
                apply_impulse(bodies, c.body_a, c.body_b, &c.mass, c.normal * delta);
            }
        }
        None => {
            for c in points.iter_mut() {
                c.solve_normal(bodies);
            }
        }
    }
}

// ── Warm starting ─────────────────────────────────────────────────────────────

/// Converged impulses of one contact point, kept for the next step.
#[derive(Debug, Clone)]
pub(super) struct CachedImpulse {
    local_anchor_a: Vec3,
    local_anchor_b: Vec3,
    normal: Vec3,
    lambda_n: Real,
    /// Friction impulse as a world vector (re-projected onto the next step's
    /// tangent basis).
    tangent_impulse: Vec3,
}

/// Persistent per-collider-pair contact impulses used to warm-start the
/// default solver. A new contact point inherits the impulse of the previous
/// point of the same collider pair whose body-local anchors (on both bodies)
/// lie within [`WARM_START_MATCH_DISTANCE`].
#[derive(Debug, Default)]
pub(crate) struct ContactImpulseCache {
    entries: HashMap<(ColliderHandle, ColliderHandle), Vec<CachedImpulse>>,
}

impl ContactImpulseCache {
    fn find(&self, c: &ContactConstraint) -> Option<&CachedImpulse> {
        let previous = self.entries.get(&(c.collider_a, c.collider_b))?;
        let max_d2 = WARM_START_MATCH_DISTANCE * WARM_START_MATCH_DISTANCE;
        previous
            .iter()
            .filter(|p| {
                (p.local_anchor_a - c.local_anchor_a).norm_squared() <= max_d2
                    && (p.local_anchor_b - c.local_anchor_b).norm_squared() <= max_d2
            })
            .min_by(|x, y| {
                let dx = (x.local_anchor_a - c.local_anchor_a).norm_squared();
                let dy = (y.local_anchor_a - c.local_anchor_a).norm_squared();
                dx.partial_cmp(&dy).unwrap_or(std::cmp::Ordering::Equal)
            })
    }

    /// Replace the cache with the converged impulses of this step.
    fn store(&mut self, constraints: &[ContactConstraint]) {
        self.entries.clear();
        for c in constraints {
            if c.depth < 0.0 {
                continue;
            }
            self.entries
                .entry((c.collider_a, c.collider_b))
                .or_default()
                .push(c.cached_impulse());
        }
    }
}

// ── Small-steps cached contacts ───────────────────────────────────────────────

/// A contact cached at frame start for [`PhysicsWorld::step_small_steps`].
///
/// The contact point is stored as an anchor in each body's local frame, so
/// the penetration can be RE-PROJECTED from the current body transforms every
/// sub-step (`depth = depth₀ − n·(p_A − p_B)`) along the frozen frame-start
/// normal, without re-running broad/narrowphase.
#[derive(Debug, Clone)]
pub(super) struct CachedContact {
    body_a: BodyHandle,
    body_b: BodyHandle,
    collider_a: ColliderHandle,
    collider_b: ColliderHandle,
    normal: Vec3,
    depth: Real,
    local_anchor_a: Vec3,
    local_anchor_b: Vec3,
    friction: Real,
    restitution: Real,
    /// Per-sub-step restitution derived from `restitution`.
    e_sub: Real,
}

impl CachedContact {
    pub(super) fn from_pair(
        cp: &ContactPair,
        bodies: &RigidBodySet,
        restitution_threshold: Real,
        substeps: usize,
    ) -> Option<Self> {
        let body_a = bodies.get(cp.body_a)?;
        let body_b = bodies.get(cp.body_b)?;
        let point = vec3(cp.contact_point);
        let e_sub =
            SubStepRestitutionCorrector::new(cp.restitution, restitution_threshold, substeps)
                .per_substep_restitution();
        Some(Self {
            body_a: cp.body_a,
            body_b: cp.body_b,
            collider_a: cp.collider_a,
            collider_b: cp.collider_b,
            normal: vec3(cp.normal),
            depth: cp.depth,
            local_anchor_a: local_anchor(body_a, &point),
            local_anchor_b: local_anchor(body_b, &point),
            friction: cp.friction,
            restitution: cp.restitution,
            e_sub,
        })
    }

    /// Re-projected penetration depth and current contact point.
    fn reproject(&self, bodies: &RigidBodySet) -> Option<(Real, Vec3)> {
        let pa = world_anchor(bodies, self.body_a, &self.local_anchor_a)?;
        let pb = world_anchor(bodies, self.body_b, &self.local_anchor_b)?;
        Some((self.depth - self.normal.dot(&(pa - pb)), (pa + pb) * 0.5))
    }

    /// Per-sub-step restitution of this contact.
    pub(super) fn per_substep_restitution(&self) -> Real {
        self.e_sub
    }
}

impl PhysicsWorld {
    /// Run the collider-driven broad/narrowphase on the current state and wake
    /// sleeping bodies that are touched by a moving body.
    pub(super) fn detect_world_contacts(&mut self) -> Vec<ContactPair> {
        let speculative =
            if self.solver_config.use_speculative && self.solver_config.speculative_margin > 0.0 {
                Some(self.solver_config.speculative_margin)
            } else {
                None
            };
        let out = detect_contacts(
            &self.bodies,
            &self.colliders,
            &self.narrowphase,
            &mut self.concave_cache,
            speculative,
        );
        for h in out.wake {
            if let Some(body) = self.bodies.get_mut(h) {
                body.wake_up();
            }
        }
        out.contacts
    }

    /// Velocity phase of the default solver: prepares one constraint per
    /// contact, warm-starts it from the previous step's converged impulses and
    /// runs `velocity_iterations` Gauss–Seidel sweeps.
    pub(super) fn solve_contacts_velocity(
        &mut self,
        contacts: &[ContactPair],
        dt: Real,
    ) -> Vec<ContactConstraint> {
        let threshold = self.solver_config.restitution_threshold;
        let mut constraints: Vec<ContactConstraint> = contacts
            .iter()
            .filter_map(|cp| ContactConstraint::prepare(cp, &self.bodies, threshold, dt))
            .collect();
        // Warm start after every restitution target has been computed from the
        // pre-solve velocities.
        for c in &mut constraints {
            if let Some(previous) = self.impulse_cache.find(c).cloned() {
                c.warm_start(&mut self.bodies, &previous);
            }
        }
        let blocks = build_blocks(&constraints);
        for _ in 0..self.solver_config.velocity_iterations {
            for block in &blocks {
                solve_block(block, &mut constraints, &mut self.bodies);
            }
        }
        self.impulse_cache.store(&constraints);
        constraints
    }

    /// Position phase of the default solver (run after velocity integration).
    pub(super) fn solve_contacts_position(&mut self, constraints: &[ContactConstraint]) {
        let beta = self.solver_config.baumgarte_factor.clamp(0.0, 1.0);
        let slop = self.solver_config.linear_slop.max(0.0);
        for _ in 0..self.solver_config.position_iterations {
            for c in constraints {
                c.solve_position(&mut self.bodies, beta, slop);
            }
        }
    }

    /// One velocity-level solve for a single cached small-steps contact.
    ///
    /// The penetration is re-projected from the current body transforms while
    /// the normal is held at its frame-start value; separated contacts
    /// (negative re-projected depth) are skipped. The normal impulse is
    /// `j = max((target − v_n) / K, 0)` with `target = max(bias, e_sub·(−v_n))`
    /// — the larger of the Baumgarte position bias and the restitution bias,
    /// never their sum — and `K = J M⁻¹ Jᵀ` including the angular terms.
    /// Coulomb friction is then applied with `|j_t| ≤ μ j` per tangent.
    ///
    /// With `with_bias == false` (the relax pass) both the position bias and the
    /// restitution term are disabled, leaving a pure non-penetration projection.
    pub(super) fn solve_cached_contact_velocity(
        &mut self,
        c: &CachedContact,
        e_sub: Real,
        h: Real,
        with_bias: bool,
    ) {
        let Some((depth, point)) = c.reproject(&self.bodies) else {
            return;
        };
        if depth < 0.0 {
            return;
        }
        let Some(m) = pair_mass(&self.bodies, c.body_a, c.body_b, point) else {
            return;
        };
        let n = c.normal;
        let k_n = m.inv_effective_mass(&n);
        if k_n < MIN_EFFECTIVE_MASS {
            return;
        }
        let vn = relative_velocity(&self.bodies, c.body_a, c.body_b, &m).dot(&n);
        let bias = if with_bias && h > 0.0 {
            let beta = self.solver_config.baumgarte_factor;
            (beta / h) * (depth - self.solver_config.linear_slop.max(0.0)).max(0.0)
        } else {
            0.0
        };
        let restitution_term = if vn < 0.0 { -e_sub * vn } else { 0.0 };
        let target = bias.max(restitution_term);
        let j = ((target - vn) / k_n).max(0.0);
        apply_impulse(&mut self.bodies, c.body_a, c.body_b, &m, n * j);
        if c.friction <= 0.0 || j <= 0.0 {
            return;
        }
        let max_friction = c.friction * j;
        for t in tangent_basis(&n) {
            let k_t = m.inv_effective_mass(&t);
            if k_t < MIN_EFFECTIVE_MASS {
                continue;
            }
            let vt = relative_velocity(&self.bodies, c.body_a, c.body_b, &m).dot(&t);
            let jt = (-vt / k_t).clamp(-max_friction, max_friction);
            apply_impulse(&mut self.bodies, c.body_a, c.body_b, &m, t * jt);
        }
    }

    /// Soft (TGS-Soft) velocity solve for one cached small-steps contact.
    ///
    /// Re-projects the penetration from the current body transforms along the
    /// cached frame-start normal and routes it through the shared soft normal
    /// solver (full linear + angular Jacobian), followed by accumulated Coulomb
    /// friction. `lambda` holds the warm-started `[normal, tangent₁, tangent₂]`
    /// accumulated impulses of this contact.
    pub(super) fn solve_cached_contact_velocity_soft(
        &mut self,
        c: &CachedContact,
        soft: &SoftParams,
        e_sub: Real,
        lambda: &mut [Real; 3],
        h: Real,
    ) {
        let Some((depth, point)) = c.reproject(&self.bodies) else {
            return;
        };
        let is_speculative_candidate = self.solver_config.use_speculative
            && depth < 0.0
            && depth >= -self.solver_config.speculative_margin;
        if depth < 0.0 && !is_speculative_candidate {
            return;
        }
        let cp = ContactPair {
            body_a: c.body_a,
            body_b: c.body_b,
            collider_a: c.collider_a,
            collider_b: c.collider_b,
            normal: [c.normal.x, c.normal.y, c.normal.z],
            depth,
            contact_point: [point.x, point.y, point.z],
            friction: c.friction,
            restitution: c.restitution,
        };
        let threshold = self.solver_config.restitution_threshold.min(0.5);
        let [lambda_n, t1, t2] = lambda;
        self.apply_soft_normal_constraint(&cp, soft, e_sub, threshold, lambda_n, h);
        let mut lambda_t = [*t1, *t2];
        solve_friction_accumulated(&mut self.bodies, &cp, *lambda_n, &mut lambda_t);
        *t1 = lambda_t[0];
        *t2 = lambda_t[1];
    }
}
