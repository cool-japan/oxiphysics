// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact dynamics and impulse-based collision resolution.
//!
//! Provides contact manifolds, sequential impulse solvers, friction models,
//! island/sleeping management, warm-starting, Baumgarte stabilization,
//! restitution models, soft contacts, rolling resistance, and material tables.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Internal math helpers (no external deps)
// ---------------------------------------------------------------------------

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n > 1e-12 {
        v3_scale(a, 1.0 / n)
    } else {
        [0.0; 3]
    }
}

#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute an orthonormal tangent frame for a given contact normal.
/// Returns two tangent vectors (t1, t2) perpendicular to `n` and to each other.
pub fn contact_tangent_frame(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let t1 = if n[0].abs() < 0.9 {
        v3_normalise(v3_cross(n, [1.0, 0.0, 0.0]))
    } else {
        v3_normalise(v3_cross(n, [0.0, 1.0, 0.0]))
    };
    let t2 = v3_cross(n, t1);
    (t1, t2)
}

// ---------------------------------------------------------------------------
// ContactManifold
// ---------------------------------------------------------------------------

/// A single point in a contact manifold.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// World-space position of the contact point.
    pub position: [f64; 3],
    /// Contact normal pointing from body B to body A (unit vector).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub penetration: f64,
    /// Accumulated normal impulse (for warm-starting).
    pub accumulated_normal_impulse: f64,
    /// Accumulated tangential impulse in the first tangent direction.
    pub accumulated_tangent1_impulse: f64,
    /// Accumulated tangential impulse in the second tangent direction.
    pub accumulated_tangent2_impulse: f64,
}

impl ContactPoint {
    /// Create a new contact point with zero accumulated impulses.
    pub fn new(position: [f64; 3], normal: [f64; 3], penetration: f64) -> Self {
        Self {
            position,
            normal: v3_normalise(normal),
            penetration,
            accumulated_normal_impulse: 0.0,
            accumulated_tangent1_impulse: 0.0,
            accumulated_tangent2_impulse: 0.0,
        }
    }

    /// Reset all accumulated impulses to zero.
    pub fn reset_impulses(&mut self) {
        self.accumulated_normal_impulse = 0.0;
        self.accumulated_tangent1_impulse = 0.0;
        self.accumulated_tangent2_impulse = 0.0;
    }
}

/// A collection of contact points between two rigid bodies.
///
/// A manifold can hold up to 4 contact points (typical physics engine limit).
/// It stores the combined friction and restitution from the material pair,
/// and the relative velocity at the start of each solver iteration.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Handle/ID of the first body.
    pub body_a: usize,
    /// Handle/ID of the second body.
    pub body_b: usize,
    /// Active contact points (maximum 4).
    pub points: Vec<ContactPoint>,
    /// Combined static friction coefficient.
    pub friction_static: f64,
    /// Combined dynamic (kinetic) friction coefficient.
    pub friction_dynamic: f64,
    /// Combined coefficient of restitution.
    pub restitution: f64,
    /// Frame counter when this manifold was last updated.
    pub last_updated_frame: u64,
}

impl ContactManifold {
    /// Construct an empty manifold between two bodies.
    pub fn new(
        body_a: usize,
        body_b: usize,
        friction_static: f64,
        friction_dynamic: f64,
        restitution: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            points: Vec::with_capacity(4),
            friction_static,
            friction_dynamic,
            restitution,
            last_updated_frame: 0,
        }
    }

    /// Add a contact point, capping at 4.
    pub fn add_point(&mut self, cp: ContactPoint) {
        if self.points.len() < 4 {
            self.points.push(cp);
        } else {
            // Replace point with smallest penetration
            if let Some(idx) = self
                .points
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.penetration
                        .partial_cmp(&b.penetration)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                && cp.penetration > self.points[idx].penetration
            {
                self.points[idx] = cp;
            }
        }
    }

    /// Remove contact points whose penetration has become non-positive (separated).
    pub fn remove_separated(&mut self) {
        self.points.retain(|p| p.penetration > 0.0);
    }

    /// Average contact normal over all contact points.
    pub fn average_normal(&self) -> [f64; 3] {
        if self.points.is_empty() {
            return [0.0, 1.0, 0.0];
        }
        let mut avg = [0.0f64; 3];
        for p in &self.points {
            avg = v3_add(avg, p.normal);
        }
        v3_normalise(v3_scale(avg, 1.0 / self.points.len() as f64))
    }

    /// Maximum penetration depth across all contact points.
    pub fn max_penetration(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.penetration)
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// ImpulseResolver
// ---------------------------------------------------------------------------

/// State for a rigid body used in impulse computation (minimal, solver-side).
#[derive(Debug, Clone)]
pub struct SolverBody {
    /// Linear velocity.
    pub vel: [f64; 3],
    /// Angular velocity.
    pub omega: [f64; 3],
    /// Inverse mass (0 = static/infinite mass).
    pub inv_mass: f64,
    /// Inverse inertia tensor diagonal (local frame, assumed diagonal here).
    pub inv_inertia: [f64; 3],
    /// Center of mass world position.
    pub position: [f64; 3],
}

impl SolverBody {
    /// Create a static (infinite-mass) solver body.
    pub fn static_body(position: [f64; 3]) -> Self {
        Self {
            vel: [0.0; 3],
            omega: [0.0; 3],
            inv_mass: 0.0,
            inv_inertia: [0.0; 3],
            position,
        }
    }

    /// Create a dynamic solver body from physical properties.
    pub fn dynamic(
        position: [f64; 3],
        vel: [f64; 3],
        omega: [f64; 3],
        mass: f64,
        inertia: [f64; 3],
    ) -> Self {
        let inv_mass = if mass > 1e-12 { 1.0 / mass } else { 0.0 };
        let inv_inertia = inertia.map(|i| if i > 1e-12 { 1.0 / i } else { 0.0 });
        Self {
            vel,
            omega,
            inv_mass,
            inv_inertia,
            position,
        }
    }

    /// Apply a linear+angular impulse delta at a contact point offset `r`.
    pub fn apply_impulse(&mut self, impulse: [f64; 3], r: [f64; 3]) {
        // Δv = J * inv_mass
        self.vel[0] += impulse[0] * self.inv_mass;
        self.vel[1] += impulse[1] * self.inv_mass;
        self.vel[2] += impulse[2] * self.inv_mass;
        // Δω = I^{-1} * (r × J)
        let torque = v3_cross(r, impulse);
        self.omega[0] += torque[0] * self.inv_inertia[0];
        self.omega[1] += torque[1] * self.inv_inertia[1];
        self.omega[2] += torque[2] * self.inv_inertia[2];
    }

    /// Velocity of a point at offset `r` from the body's COM.
    pub fn point_velocity(&self, r: [f64; 3]) -> [f64; 3] {
        let omega_cross_r = v3_cross(self.omega, r);
        v3_add(self.vel, omega_cross_r)
    }
}

/// Sequential impulse solver for contact resolution.
///
/// Implements the Gauss-Seidel iterative impulse solver used in most
/// physics engines (e.g. Bullet, Box2D). Supports Hertz-based contact
/// stiffness and a configurable restitution coefficient.
#[derive(Debug, Clone)]
pub struct ImpulseResolver {
    /// Number of solver iterations per frame.
    pub iterations: usize,
    /// Default coefficient of restitution (0 = perfectly inelastic, 1 = elastic).
    pub default_restitution: f64,
    /// Hertz contact stiffness (N/m^(3/2)) — used in soft contact mode.
    pub hertz_stiffness: f64,
    /// Velocity threshold below which restitution is zeroed (avoids jitter).
    pub restitution_velocity_threshold: f64,
}

impl Default for ImpulseResolver {
    fn default() -> Self {
        Self {
            iterations: 10,
            default_restitution: 0.3,
            hertz_stiffness: 1.0e6,
            restitution_velocity_threshold: 0.5,
        }
    }
}

impl ImpulseResolver {
    /// Compute the effective mass along a direction `n` for two bodies at offsets `ra`, `rb`.
    pub fn effective_mass(
        ba: &SolverBody,
        bb: &SolverBody,
        n: [f64; 3],
        ra: [f64; 3],
        rb: [f64; 3],
    ) -> f64 {
        let ra_cross_n = v3_cross(ra, n);
        let rb_cross_n = v3_cross(rb, n);
        let term_a = ba.inv_mass
            + ra_cross_n[0] * ra_cross_n[0] * ba.inv_inertia[0]
            + ra_cross_n[1] * ra_cross_n[1] * ba.inv_inertia[1]
            + ra_cross_n[2] * ra_cross_n[2] * ba.inv_inertia[2];
        let term_b = bb.inv_mass
            + rb_cross_n[0] * rb_cross_n[0] * bb.inv_inertia[0]
            + rb_cross_n[1] * rb_cross_n[1] * bb.inv_inertia[1]
            + rb_cross_n[2] * rb_cross_n[2] * bb.inv_inertia[2];
        let denom = term_a + term_b;
        if denom > 1e-14 { 1.0 / denom } else { 0.0 }
    }

    /// Resolve one contact point with sequential impulse.
    /// Modifies `ba` and `bb` in place.
    pub fn resolve_contact(
        &self,
        ba: &mut SolverBody,
        bb: &mut SolverBody,
        cp: &mut ContactPoint,
        restitution: f64,
        _friction: f64,
    ) {
        let n = cp.normal;
        let ra = v3_sub(cp.position, ba.position);
        let rb = v3_sub(cp.position, bb.position);

        let va = ba.point_velocity(ra);
        let vb = bb.point_velocity(rb);
        let rel_vel = v3_sub(va, vb);
        let vn = v3_dot(rel_vel, n);

        // Only resolve if bodies are approaching
        if vn > 0.0 {
            return;
        }

        let e = if vn.abs() > self.restitution_velocity_threshold {
            restitution
        } else {
            0.0
        };
        let eff_mass = Self::effective_mass(ba, bb, n, ra, rb);
        let j_num = -(1.0 + e) * vn;
        let j = j_num * eff_mass;

        // Clamp accumulated impulse (non-penetration constraint)
        let old = cp.accumulated_normal_impulse;
        cp.accumulated_normal_impulse = (old + j).max(0.0);
        let delta_j = cp.accumulated_normal_impulse - old;

        let impulse = v3_scale(n, delta_j);
        ba.apply_impulse(impulse, ra);
        bb.apply_impulse(v3_scale(impulse, -1.0), rb);
    }

    /// Run all iterations over all manifolds in a slice.
    pub fn solve_manifolds(&self, bodies: &mut [SolverBody], manifolds: &mut [ContactManifold]) {
        for _iter in 0..self.iterations {
            for manifold in manifolds.iter_mut() {
                let ia = manifold.body_a;
                let ib = manifold.body_b;
                if ia >= bodies.len() || ib >= bodies.len() {
                    continue;
                }
                let restitution = manifold.restitution;
                let friction = manifold.friction_dynamic;
                // Split borrow
                let (a_slice, b_slice) = bodies.split_at_mut(ia.max(ib));
                let (ba, bb) = if ia < ib {
                    (&mut a_slice[ia], &mut b_slice[0])
                } else {
                    (&mut b_slice[0], &mut a_slice[ib])
                };
                for cp in manifold.points.iter_mut() {
                    self.resolve_contact(ba, bb, cp, restitution, friction);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// FrictionModel
// ---------------------------------------------------------------------------

/// Friction model type selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FrictionType {
    /// Classic Coulomb cone friction (isotropic).
    CoulombIsotropic,
    /// Anisotropic friction with separate coefficients per tangent direction.
    Anisotropic,
    /// Simplified model: friction force is always dynamic (no static/dynamic distinction).
    KineticOnly,
}

/// Coulomb cone friction model with static/dynamic transition.
///
/// Computes tangential friction impulses subject to `|F_t| ≤ μ |F_n|`.
#[derive(Debug, Clone)]
pub struct FrictionModel {
    /// Friction model variant.
    pub kind: FrictionType,
    /// Static friction coefficient (slip starts when tangential force exceeds this).
    pub mu_static: f64,
    /// Dynamic friction coefficient (force during sliding).
    pub mu_dynamic: f64,
    /// Friction in the second tangent direction (anisotropic mode only).
    pub mu_static_t2: f64,
    /// Dynamic friction in the second tangent direction (anisotropic mode only).
    pub mu_dynamic_t2: f64,
    /// Velocity threshold for static-to-dynamic transition (m/s).
    pub slip_velocity_threshold: f64,
}

impl Default for FrictionModel {
    fn default() -> Self {
        Self {
            kind: FrictionType::CoulombIsotropic,
            mu_static: 0.6,
            mu_dynamic: 0.4,
            mu_static_t2: 0.6,
            mu_dynamic_t2: 0.4,
            slip_velocity_threshold: 0.01,
        }
    }
}

impl FrictionModel {
    /// Check if the contact is sliding given a tangential slip speed.
    pub fn is_sliding(&self, slip_speed: f64) -> bool {
        slip_speed > self.slip_velocity_threshold
    }

    /// Return the active friction coefficient given the slip speed.
    pub fn mu_active(&self, slip_speed: f64) -> f64 {
        if self.is_sliding(slip_speed) {
            self.mu_dynamic
        } else {
            self.mu_static
        }
    }

    /// Apply Coulomb friction impulse clamping to a proposed tangential impulse.
    /// `normal_impulse` is the magnitude of the resolved normal impulse.
    /// `tangential` is the proposed tangential impulse vector.
    /// Returns the clamped tangential impulse vector.
    pub fn clamp_impulse_isotropic(
        &self,
        normal_impulse: f64,
        tangential: [f64; 3],
        slip_speed: f64,
    ) -> [f64; 3] {
        let mu = self.mu_active(slip_speed);
        let limit = mu * normal_impulse.abs();
        let t_norm = v3_norm(tangential);
        if t_norm > limit {
            v3_scale(tangential, limit / t_norm)
        } else {
            tangential
        }
    }

    /// Clamp anisotropic friction: separate limits per tangent axis.
    pub fn clamp_impulse_anisotropic(
        &self,
        normal_impulse: f64,
        jt1: f64,
        jt2: f64,
        slip_speed_t1: f64,
        slip_speed_t2: f64,
    ) -> (f64, f64) {
        let mu1 = if self.is_sliding(slip_speed_t1) {
            self.mu_dynamic
        } else {
            self.mu_static
        };
        let mu2 = if self.is_sliding(slip_speed_t2) {
            self.mu_dynamic_t2
        } else {
            self.mu_static_t2
        };
        let limit1 = mu1 * normal_impulse.abs();
        let limit2 = mu2 * normal_impulse.abs();
        (jt1.clamp(-limit1, limit1), jt2.clamp(-limit2, limit2))
    }
}

// ---------------------------------------------------------------------------
// ContactGraph — island detection and sleeping
// ---------------------------------------------------------------------------

/// A contact island: a set of bodies mutually connected by contacts.
#[derive(Debug, Clone, Default)]
pub struct Island {
    /// Body indices in this island.
    pub bodies: Vec<usize>,
    /// Whether all bodies in the island qualify for sleeping.
    pub all_sleeping: bool,
}

/// A graph of bodies connected by contact edges used for island detection
/// and sleeping management.
///
/// Uses union-find to efficiently merge bodies into islands when contacts form.
#[derive(Debug, Clone, Default)]
pub struct ContactGraph {
    /// Union-find parent array.
    parent: Vec<usize>,
    /// Union-find rank array.
    rank: Vec<usize>,
    /// Sleep timers per body (seconds below threshold).
    sleep_timer: Vec<f64>,
    /// Whether each body is currently sleeping.
    is_sleeping: Vec<bool>,
    /// Number of bodies registered.
    num_bodies: usize,
}

impl ContactGraph {
    /// Create a new contact graph for `n` bodies.
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
            sleep_timer: vec![0.0; n],
            is_sleeping: vec![false; n],
            num_bodies: n,
        }
    }

    /// Reset all islands (called each frame before contact detection).
    pub fn reset(&mut self) {
        for i in 0..self.num_bodies {
            self.parent[i] = i;
            self.rank[i] = 0;
        }
    }

    /// Find root of the island for body `i` (path compression).
    pub fn find(&mut self, i: usize) -> usize {
        if self.parent[i] != i {
            self.parent[i] = self.find(self.parent[i]);
        }
        self.parent[i]
    }

    /// Union the islands of bodies `a` and `b`.
    pub fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }

    /// Build islands from all manifolds.
    pub fn build_islands(&mut self, manifolds: &[ContactManifold]) -> Vec<Island> {
        self.reset();
        for m in manifolds {
            self.union(m.body_a, m.body_b);
        }
        let mut map: HashMap<usize, Island> = HashMap::new();
        for i in 0..self.num_bodies {
            let root = self.find(i);
            let island = map.entry(root).or_default();
            island.bodies.push(i);
            if !self.is_sleeping[i] {
                island.all_sleeping = false;
            }
        }
        map.into_values().collect()
    }

    /// Update sleep timers; body sleeps after `threshold` seconds below `vel_limit`.
    pub fn update_sleeping(
        &mut self,
        body_speeds: &[f64],
        dt: f64,
        vel_limit: f64,
        threshold: f64,
    ) {
        for (i, &speed) in body_speeds.iter().enumerate().take(self.num_bodies) {
            if speed < vel_limit {
                self.sleep_timer[i] += dt;
                if self.sleep_timer[i] >= threshold {
                    self.is_sleeping[i] = true;
                }
            } else {
                self.sleep_timer[i] = 0.0;
                self.is_sleeping[i] = false;
            }
        }
    }

    /// Wake a body explicitly (e.g. due to external force).
    pub fn wake_body(&mut self, i: usize) {
        if i < self.num_bodies {
            self.is_sleeping[i] = false;
            self.sleep_timer[i] = 0.0;
        }
    }

    /// Query whether body `i` is sleeping.
    pub fn is_sleeping(&self, i: usize) -> bool {
        i < self.num_bodies && self.is_sleeping[i]
    }

    /// Count sleeping bodies.
    pub fn num_sleeping(&self) -> usize {
        self.is_sleeping.iter().filter(|&&s| s).count()
    }
}

// ---------------------------------------------------------------------------
// WarmStarting
// ---------------------------------------------------------------------------

/// Warm-starting cache: stores contact impulses from the previous frame
/// and re-applies them at the start of the new frame to accelerate convergence.
#[derive(Debug, Clone, Default)]
pub struct WarmStarting {
    /// Map from (body_a, body_b) pair to cached normal impulse.
    cache: HashMap<(usize, usize), Vec<f64>>,
    /// Aging factor: cached impulses are multiplied by this each frame (0..1].
    pub aging: f64,
}

impl WarmStarting {
    /// Create a warm-starting cache with default aging factor 0.95.
    pub fn new() -> Self {
        Self {
            cache: HashMap::new(),
            aging: 0.95,
        }
    }

    /// Store accumulated impulses from a manifold.
    pub fn store(&mut self, manifold: &ContactManifold) {
        let key = (
            manifold.body_a.min(manifold.body_b),
            manifold.body_a.max(manifold.body_b),
        );
        let impulses: Vec<f64> = manifold
            .points
            .iter()
            .map(|p| p.accumulated_normal_impulse)
            .collect();
        self.cache.insert(key, impulses);
    }

    /// Apply cached impulses to a manifold at the start of a frame.
    /// Scales by `aging` to reduce reliance on stale data.
    pub fn apply(&self, manifold: &mut ContactManifold) {
        let key = (
            manifold.body_a.min(manifold.body_b),
            manifold.body_a.max(manifold.body_b),
        );
        if let Some(cached) = self.cache.get(&key) {
            for (i, cp) in manifold.points.iter_mut().enumerate() {
                if let Some(&ci) = cached.get(i) {
                    cp.accumulated_normal_impulse = (ci * self.aging).max(0.0);
                }
            }
        }
    }

    /// Remove entries not present in the current manifold set (stale contacts).
    pub fn prune(&mut self, active_keys: &[(usize, usize)]) {
        self.cache.retain(|k, _| active_keys.contains(k));
    }

    /// Clear all cached impulses.
    pub fn clear(&mut self) {
        self.cache.clear();
    }
}

// ---------------------------------------------------------------------------
// PositionCorrection
// ---------------------------------------------------------------------------

/// Position correction strategy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PositionCorrectionMethod {
    /// Baumgarte stabilization: adds a bias velocity to push bodies apart.
    Baumgarte,
    /// Non-linear pseudo-velocity projection (more stable for stacked objects).
    PseudoVelocity,
    /// No position correction (rely solely on impulses).
    None,
}

/// Handles penetration correction to prevent drift in the constraint solver.
///
/// Baumgarte stabilization injects a velocity bias `β/Δt * penetration` into
/// the velocity constraint, trading energy accuracy for positional stability.
/// Pseudo-velocity projection solves a separate position-level LCP each frame.
#[derive(Debug, Clone)]
pub struct PositionCorrection {
    /// Selected correction method.
    pub method: PositionCorrectionMethod,
    /// Baumgarte factor β ∈ (0, 1].  Typical: 0.2.
    pub baumgarte_beta: f64,
    /// Allowed slop (ignored penetration depth) to avoid jitter (m).
    pub slop: f64,
    /// Maximum correction per step (m), to prevent over-correction.
    pub max_correction: f64,
}

impl Default for PositionCorrection {
    fn default() -> Self {
        Self {
            method: PositionCorrectionMethod::Baumgarte,
            baumgarte_beta: 0.2,
            slop: 0.005,
            max_correction: 0.2,
        }
    }
}

impl PositionCorrection {
    /// Compute the Baumgarte bias velocity for a contact point.
    /// `penetration` — signed overlap (positive = penetrating), `dt` — time step.
    pub fn baumgarte_bias(&self, penetration: f64, dt: f64) -> f64 {
        let corrected = (penetration - self.slop).max(0.0);
        let bias = (self.baumgarte_beta / dt) * corrected;
        bias.min(self.max_correction / dt)
    }

    /// Compute pseudo-velocity correction magnitude for position projection.
    pub fn pseudo_velocity_correction(&self, penetration: f64, _dt: f64) -> f64 {
        let corrected = (penetration - self.slop).max(0.0);
        corrected.min(self.max_correction)
    }

    /// Apply position correction to a pair of solver bodies at a contact point.
    pub fn correct_position(
        &self,
        ba: &mut SolverBody,
        bb: &mut SolverBody,
        normal: [f64; 3],
        penetration: f64,
        dt: f64,
    ) {
        if penetration <= self.slop {
            return;
        }
        let correction = match self.method {
            PositionCorrectionMethod::Baumgarte => {
                // Baumgarte is applied as a velocity bias inside the solver — position unchanged here
                0.0
            }
            PositionCorrectionMethod::PseudoVelocity => {
                self.pseudo_velocity_correction(penetration, dt)
            }
            PositionCorrectionMethod::None => 0.0,
        };
        if correction > 0.0 {
            let total_inv_mass = ba.inv_mass + bb.inv_mass;
            if total_inv_mass > 1e-12 {
                let scale_a = ba.inv_mass / total_inv_mass;
                let scale_b = bb.inv_mass / total_inv_mass;
                ba.position = v3_add(ba.position, v3_scale(normal, correction * scale_a));
                bb.position = v3_sub(bb.position, v3_scale(normal, correction * scale_b));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// BouncingModel — restitution models
// ---------------------------------------------------------------------------

/// Restitution model variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RestitutionModel {
    /// Newton's hypothesis: `e = v_after / v_before` (kinematic).
    Newton,
    /// Poisson's hypothesis: impulse ratio between compression and restitution phases.
    Poisson,
    /// Energetic (stronge) model: energy-consistent restitution.
    Energetic,
}

/// Computes post-collision velocity ratios under different restitution models.
#[derive(Debug, Clone)]
pub struct BouncingModel {
    /// Active model.
    pub model: RestitutionModel,
    /// Coefficient of restitution `e` ∈ \[0, 1\].
    pub coefficient: f64,
    /// Minimum approach speed (m/s) below which e is zeroed (avoids micro-bouncing).
    pub min_speed: f64,
}

impl Default for BouncingModel {
    fn default() -> Self {
        Self {
            model: RestitutionModel::Newton,
            coefficient: 0.5,
            min_speed: 0.2,
        }
    }
}

impl BouncingModel {
    /// Effective restitution for a given approach speed `v_approach` (positive = approaching).
    pub fn effective_e(&self, v_approach: f64) -> f64 {
        if v_approach < self.min_speed {
            0.0
        } else {
            self.coefficient
        }
    }

    /// Compute post-bounce relative normal velocity.
    /// `v_n` — pre-impact relative normal velocity (negative = approaching).
    pub fn post_impact_velocity(&self, v_n: f64) -> f64 {
        let e = self.effective_e(-v_n);
        match self.model {
            RestitutionModel::Newton => -e * v_n,
            RestitutionModel::Poisson => {
                // Poisson: same result for single-contact case
                -e * v_n
            }
            RestitutionModel::Energetic => {
                // Energetic: v_after^2 = e^2 * v_before^2 (energy-consistent sign)
                if v_n < 0.0 { e * (-v_n) } else { -e * v_n }
            }
        }
    }

    /// Impact impulse magnitude for a given relative approach speed and effective mass.
    pub fn impact_impulse(&self, v_n: f64, eff_mass: f64) -> f64 {
        let e = self.effective_e(-v_n);
        -(1.0 + e) * v_n * eff_mass
    }
}

// ---------------------------------------------------------------------------
// SoftContact — penalty / spring-damper
// ---------------------------------------------------------------------------

/// Soft contact model using a spring-damper penalty force.
///
/// Instead of discrete impulses, this model generates continuous forces
/// proportional to penetration depth and relative closing velocity.
/// Suitable for deformable-body interactions or very soft materials.
#[derive(Debug, Clone)]
pub struct SoftContact {
    /// Contact spring stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Tangential spring stiffness (N/m) — resists sliding during static friction.
    pub tangential_stiffness: f64,
    /// Maximum penalty force (N), to prevent explosion.
    pub max_force: f64,
}

impl Default for SoftContact {
    fn default() -> Self {
        Self {
            stiffness: 1.0e4,
            damping: 1.0e2,
            tangential_stiffness: 5.0e3,
            max_force: 1.0e6,
        }
    }
}

impl SoftContact {
    /// Compute normal penalty force for a penetrating contact.
    /// `penetration` — overlap depth (m), `closing_vel` — relative closing velocity (m/s).
    pub fn normal_force(&self, penetration: f64, closing_vel: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        let f = self.stiffness * penetration - self.damping * closing_vel;
        f.max(0.0).min(self.max_force)
    }

    /// Compute tangential spring force resisting micro-slip.
    /// `tangential_disp` — tangential displacement from rest position (m).
    pub fn tangential_spring_force(&self, tangential_disp: f64) -> f64 {
        (self.tangential_stiffness * tangential_disp).min(self.max_force)
    }

    /// Full contact force vector.
    /// Returns `[Fn * normal + Ft * tangent]` as world-space force.
    pub fn contact_force(
        &self,
        normal: [f64; 3],
        tangent: [f64; 3],
        penetration: f64,
        closing_vel: f64,
        tangential_disp: f64,
    ) -> [f64; 3] {
        let fn_ = self.normal_force(penetration, closing_vel);
        let ft = self.tangential_spring_force(tangential_disp);
        v3_add(v3_scale(normal, fn_), v3_scale(tangent, ft))
    }

    /// Hertz contact model: force proportional to penetration^(3/2).
    /// `e_star` — combined elastic modulus (Pa), `r_star` — combined radius (m).
    pub fn hertz_force(e_star: f64, r_star: f64, penetration: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * e_star * r_star.sqrt() * penetration.powf(1.5)
    }
}

// ---------------------------------------------------------------------------
// RollingResistance
// ---------------------------------------------------------------------------

/// Rolling resistance model for spheres/cylinders rolling on surfaces.
///
/// Computes a rolling friction torque that opposes the rolling angular velocity.
/// Uses the simple linear `T = μ_r * F_n * R` model by default.
#[derive(Debug, Clone)]
pub struct RollingResistance {
    /// Rolling resistance coefficient μ_r (dimensionless). Typical: 0.001–0.01.
    pub mu_rolling: f64,
    /// Spinning resistance coefficient μ_s (for spin about the normal axis).
    pub mu_spinning: f64,
    /// Effective contact radius (m) — ball or wheel radius.
    pub contact_radius: f64,
}

impl Default for RollingResistance {
    fn default() -> Self {
        Self {
            mu_rolling: 0.005,
            mu_spinning: 0.002,
            contact_radius: 0.1,
        }
    }
}

impl RollingResistance {
    /// Rolling friction torque magnitude: T = μ_r * F_n * R.
    pub fn rolling_torque(&self, normal_force: f64) -> f64 {
        self.mu_rolling * normal_force.abs() * self.contact_radius
    }

    /// Spinning resistance torque magnitude: T = μ_s * F_n * R.
    pub fn spinning_torque(&self, normal_force: f64) -> f64 {
        self.mu_spinning * normal_force.abs() * self.contact_radius
    }

    /// Rolling torque vector opposing the rolling angular velocity component.
    /// `omega_roll` — angular velocity projected onto tangent plane (rad/s).
    /// `normal_force` — contact normal force magnitude (N).
    pub fn torque_vector(&self, omega_roll: [f64; 3], normal_force: f64) -> [f64; 3] {
        let omega_norm = v3_norm(omega_roll);
        if omega_norm < 1e-12 {
            return [0.0; 3];
        }
        let omega_dir = v3_scale(omega_roll, 1.0 / omega_norm);
        let t_mag = self.rolling_torque(normal_force);
        v3_scale(omega_dir, -t_mag)
    }

    /// Spinning torque vector opposing spin about the contact normal.
    pub fn spin_torque_vector(
        &self,
        omega_spin: f64,
        normal: [f64; 3],
        normal_force: f64,
    ) -> [f64; 3] {
        let t_mag = self.spinning_torque(normal_force);
        let sign = if omega_spin > 0.0 { -1.0 } else { 1.0 };
        v3_scale(normal, sign * t_mag)
    }
}

// ---------------------------------------------------------------------------
// ContactMaterial — material pairing table
// ---------------------------------------------------------------------------

/// Material ID: a simple integer index into the material table.
pub type MaterialId = usize;

/// Physical contact properties for a single material.
#[derive(Debug, Clone)]
pub struct Material {
    /// Human-readable name.
    pub name: String,
    /// Static friction coefficient.
    pub mu_static: f64,
    /// Dynamic friction coefficient.
    pub mu_dynamic: f64,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Rolling resistance coefficient.
    pub mu_rolling: f64,
    /// Young's modulus (Pa), for Hertz contact.
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            name: "default".into(),
            mu_static: 0.5,
            mu_dynamic: 0.4,
            restitution: 0.3,
            mu_rolling: 0.005,
            youngs_modulus: 2.0e11,
            poisson_ratio: 0.3,
        }
    }
}

/// Combined contact properties for a pair of materials.
#[derive(Debug, Clone)]
pub struct PairedProperties {
    /// Combined static friction (geometric mean).
    pub friction_static: f64,
    /// Combined dynamic friction (geometric mean).
    pub friction_dynamic: f64,
    /// Combined restitution (minimum).
    pub restitution: f64,
    /// Combined rolling resistance (geometric mean).
    pub mu_rolling: f64,
}

/// Table mapping material pairs to combined contact properties.
///
/// Uses the geometric-mean mixing rule for friction and minimum for restitution.
/// Can hold explicit overrides for specific pairs.
#[derive(Debug, Clone, Default)]
pub struct ContactMaterial {
    /// Registered materials.
    materials: Vec<Material>,
    /// Explicit overrides for specific pairs `(a, b)` with `a <= b`.
    overrides: HashMap<(MaterialId, MaterialId), PairedProperties>,
}

impl ContactMaterial {
    /// Create an empty material table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a new material and return its ID.
    pub fn register(&mut self, mat: Material) -> MaterialId {
        let id = self.materials.len();
        self.materials.push(mat);
        id
    }

    /// Get a reference to a material by ID.
    pub fn get(&self, id: MaterialId) -> Option<&Material> {
        self.materials.get(id)
    }

    /// Add an explicit override for a specific pair.
    pub fn set_override(&mut self, a: MaterialId, b: MaterialId, props: PairedProperties) {
        let key = (a.min(b), a.max(b));
        self.overrides.insert(key, props);
    }

    /// Compute combined contact properties for a pair of materials.
    /// Returns explicit override if present, otherwise applies mixing rules.
    pub fn combined(&self, a: MaterialId, b: MaterialId) -> PairedProperties {
        let key = (a.min(b), a.max(b));
        if let Some(ov) = self.overrides.get(&key) {
            return ov.clone();
        }
        let ma = self.materials.get(a).cloned().unwrap_or_default();
        let mb = self.materials.get(b).cloned().unwrap_or_default();
        PairedProperties {
            friction_static: (ma.mu_static * mb.mu_static).sqrt(),
            friction_dynamic: (ma.mu_dynamic * mb.mu_dynamic).sqrt(),
            restitution: ma.restitution.min(mb.restitution),
            mu_rolling: (ma.mu_rolling * mb.mu_rolling).sqrt(),
        }
    }

    /// Number of registered materials.
    pub fn len(&self) -> usize {
        self.materials.len()
    }

    /// Returns true if no materials are registered.
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- ContactPoint ----

    #[test]
    fn contact_point_new_normalises_normal() {
        let cp = ContactPoint::new([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.1);
        let n = cp.normal;
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-10,
            "normal should be unit, got {len}"
        );
    }

    #[test]
    fn contact_point_reset_impulses() {
        let mut cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.05);
        cp.accumulated_normal_impulse = 5.0;
        cp.accumulated_tangent1_impulse = 1.0;
        cp.reset_impulses();
        assert!(cp.accumulated_normal_impulse.abs() < 1e-12);
        assert!(cp.accumulated_tangent1_impulse.abs() < 1e-12);
    }

    // ---- ContactManifold ----

    #[test]
    fn manifold_add_up_to_four_points() {
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        for i in 0..6 {
            m.add_point(ContactPoint::new(
                [i as f64, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                0.01 * (i + 1) as f64,
            ));
        }
        assert_eq!(m.points.len(), 4, "manifold should cap at 4 points");
    }

    #[test]
    fn manifold_remove_separated() {
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01));
        m.add_point(ContactPoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], -0.01));
        m.remove_separated();
        assert_eq!(m.points.len(), 1, "separated points should be removed");
    }

    #[test]
    fn manifold_average_normal_single_point() {
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [1.0, 0.0, 0.0], 0.1));
        let n = m.average_normal();
        assert!((n[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn manifold_max_penetration() {
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.02));
        m.add_point(ContactPoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.05));
        assert!((m.max_penetration() - 0.05).abs() < 1e-12);
    }

    // ---- SolverBody ----

    #[test]
    fn solver_body_apply_impulse_linear() {
        let mut b = SolverBody::dynamic([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        b.apply_impulse([2.0, 0.0, 0.0], [0.0; 3]);
        assert!((b.vel[0] - 2.0).abs() < 1e-12, "vel.x should be 2.0");
    }

    #[test]
    fn solver_body_static_unaffected_by_impulse() {
        let mut b = SolverBody::static_body([0.0; 3]);
        b.apply_impulse([100.0, 0.0, 0.0], [0.0; 3]);
        assert!(b.vel[0].abs() < 1e-12, "static body should not move");
    }

    #[test]
    fn solver_body_point_velocity() {
        let b = SolverBody::dynamic([0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, [1.0; 3]);
        // omega = (0,0,1), r = (1,0,0) → omega×r = (0,1,0)
        // v_point = (1,0,0) + (0,1,0) = (1,1,0)
        let vp = b.point_velocity([1.0, 0.0, 0.0]);
        assert!((vp[0] - 1.0).abs() < 1e-10);
        assert!((vp[1] - 1.0).abs() < 1e-10);
        assert!(vp[2].abs() < 1e-10);
        let _ = b.vel; // suppress unused warning
    }

    // ---- ImpulseResolver ----

    #[test]
    fn impulse_resolver_separating_contact_no_impulse() {
        let resolver = ImpulseResolver::default();
        let mut ba = SolverBody::dynamic([0.0; 3], [0.0, 1.0, 0.0], [0.0; 3], 1.0, [1.0; 3]);
        let mut bb = SolverBody::static_body([0.0, -0.5, 0.0]);
        let mut cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01);
        let old_vy = ba.vel[1];
        resolver.resolve_contact(&mut ba, &mut bb, &mut cp, 0.5, 0.4);
        assert!(
            (ba.vel[1] - old_vy).abs() < 1e-12,
            "separating contact should not change velocity"
        );
    }

    #[test]
    fn impulse_resolver_approaching_contact_generates_impulse() {
        let resolver = ImpulseResolver {
            restitution_velocity_threshold: 0.0,
            ..Default::default()
        };
        let mut ba = SolverBody::dynamic([0.0; 3], [0.0, -1.0, 0.0], [0.0; 3], 1.0, [1.0; 3]);
        let mut bb = SolverBody::static_body([0.0, -0.5, 0.0]);
        let mut cp = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01);
        resolver.resolve_contact(&mut ba, &mut bb, &mut cp, 0.5, 0.4);
        assert!(ba.vel[1] > -1.0, "velocity should increase after bounce");
    }

    #[test]
    fn effective_mass_symmetric() {
        let ba = SolverBody::dynamic([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        let bb = SolverBody::dynamic([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        let n = [0.0, 1.0, 0.0];
        let r = [0.0; 3];
        let em_ab = ImpulseResolver::effective_mass(&ba, &bb, n, r, r);
        let em_ba = ImpulseResolver::effective_mass(&bb, &ba, n, r, r);
        assert!(
            (em_ab - em_ba).abs() < 1e-12,
            "effective mass should be symmetric"
        );
    }

    // ---- FrictionModel ----

    #[test]
    fn friction_model_no_slip_below_threshold() {
        let fm = FrictionModel::default();
        assert!(!fm.is_sliding(0.005), "below threshold should not slide");
        assert!(fm.is_sliding(0.05), "above threshold should slide");
    }

    #[test]
    fn friction_model_mu_active() {
        let fm = FrictionModel::default();
        // Below threshold: static
        assert!((fm.mu_active(0.005) - fm.mu_static).abs() < 1e-12);
        // Above threshold: dynamic
        assert!((fm.mu_active(0.1) - fm.mu_dynamic).abs() < 1e-12);
    }

    #[test]
    fn friction_model_clamp_impulse_isotropic() {
        let fm = FrictionModel::default();
        // Normal impulse = 10 N, mu_dynamic = 0.4 → max tangential = 4 N
        let tangential = [5.0, 0.0, 0.0]; // exceeds limit
        let clamped = fm.clamp_impulse_isotropic(10.0, tangential, 0.1);
        let len =
            (clamped[0] * clamped[0] + clamped[1] * clamped[1] + clamped[2] * clamped[2]).sqrt();
        assert!(
            (len - 4.0).abs() < 1e-10,
            "clamped magnitude should be 4.0, got {len}"
        );
    }

    #[test]
    fn friction_model_clamp_within_limit_unchanged() {
        let fm = FrictionModel::default();
        let tangential = [1.0, 0.0, 0.0]; // within limit of 4.0
        let clamped = fm.clamp_impulse_isotropic(10.0, tangential, 0.1);
        assert!(
            (clamped[0] - 1.0).abs() < 1e-10,
            "within limit should be unchanged"
        );
    }

    // ---- ContactGraph ----

    #[test]
    fn contact_graph_single_body_island() {
        let mut cg = ContactGraph::new(3);
        let islands = cg.build_islands(&[]);
        assert_eq!(islands.len(), 3, "3 isolated bodies → 3 islands");
    }

    #[test]
    fn contact_graph_connected_pair_one_island() {
        let mut cg = ContactGraph::new(2);
        let m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        let islands = cg.build_islands(&[m]);
        assert_eq!(islands.len(), 1, "connected pair → 1 island");
        assert_eq!(islands[0].bodies.len(), 2);
    }

    #[test]
    fn contact_graph_sleeping_bodies() {
        let mut cg = ContactGraph::new(3);
        let speeds = [0.001, 0.001, 5.0];
        cg.update_sleeping(&speeds, 0.016, 0.1, 0.5);
        // After one frame, sleep timers are 0.016 — no body should sleep yet
        assert_eq!(cg.num_sleeping(), 0);
    }

    #[test]
    fn contact_graph_sleep_after_threshold() {
        let mut cg = ContactGraph::new(1);
        let speeds = [0.001];
        // 0.5 s / 0.016 = ~32 frames, run 40 frames
        for _ in 0..40 {
            cg.update_sleeping(&speeds, 0.016, 0.1, 0.5);
        }
        assert!(cg.is_sleeping(0), "slow body should sleep after threshold");
    }

    #[test]
    fn contact_graph_wake_body() {
        let mut cg = ContactGraph::new(1);
        cg.is_sleeping[0] = true;
        cg.sleep_timer[0] = 1.0;
        cg.wake_body(0);
        assert!(!cg.is_sleeping(0), "woken body should not be sleeping");
    }

    // ---- WarmStarting ----

    #[test]
    fn warm_starting_store_and_apply() {
        let mut ws = WarmStarting::new();
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01));
        m.points[0].accumulated_normal_impulse = 3.0;
        ws.store(&m);

        // Reset and apply
        m.points[0].accumulated_normal_impulse = 0.0;
        ws.apply(&mut m);
        // After apply: 3.0 * 0.95 = 2.85
        assert!((m.points[0].accumulated_normal_impulse - 3.0 * ws.aging).abs() < 1e-10);
    }

    #[test]
    fn warm_starting_prune_removes_stale() {
        let mut ws = WarmStarting::new();
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01));
        ws.store(&m);
        // Prune with empty active set
        ws.prune(&[]);
        assert!(ws.cache.is_empty(), "cache should be pruned");
    }

    #[test]
    fn warm_starting_clear() {
        let mut ws = WarmStarting::new();
        let mut m = ContactManifold::new(0, 1, 0.5, 0.4, 0.3);
        m.add_point(ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01));
        ws.store(&m);
        ws.clear();
        assert!(ws.cache.is_empty());
    }

    // ---- PositionCorrection ----

    #[test]
    fn baumgarte_bias_zero_for_slop() {
        let pc = PositionCorrection::default();
        let bias = pc.baumgarte_bias(pc.slop * 0.5, 0.016);
        assert!(bias.abs() < 1e-12, "penetration within slop → zero bias");
    }

    #[test]
    fn baumgarte_bias_positive_for_penetration() {
        let pc = PositionCorrection::default();
        let bias = pc.baumgarte_bias(0.1, 0.016);
        assert!(bias > 0.0, "deep penetration → positive bias");
    }

    #[test]
    fn position_correction_pseudo_velocity() {
        let pc = PositionCorrection {
            method: PositionCorrectionMethod::PseudoVelocity,
            ..Default::default()
        };
        let mut ba = SolverBody::dynamic([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        let mut bb = SolverBody::static_body([0.0, -1.0, 0.0]);
        let normal = [0.0, 1.0, 0.0];
        let pen = 0.05;
        let orig_y = ba.position[1];
        pc.correct_position(&mut ba, &mut bb, normal, pen, 0.016);
        // ba should move upward (positive Y), bb is static so no movement
        assert!(
            ba.position[1] >= orig_y,
            "body A should be corrected upward"
        );
    }

    #[test]
    fn position_correction_none_no_change() {
        let pc = PositionCorrection {
            method: PositionCorrectionMethod::None,
            ..Default::default()
        };
        let mut ba = SolverBody::dynamic([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        let mut bb = SolverBody::static_body([0.0, -1.0, 0.0]);
        let pos_before = ba.position;
        pc.correct_position(&mut ba, &mut bb, [0.0, 1.0, 0.0], 0.05, 0.016);
        assert_eq!(
            ba.position, pos_before,
            "None correction should not change position"
        );
    }

    // ---- BouncingModel ----

    #[test]
    fn bouncing_model_inelastic_below_min_speed() {
        let bm = BouncingModel::default();
        let e = bm.effective_e(0.1); // below min_speed = 0.2
        assert!(e.abs() < 1e-12, "e should be zero below min speed");
    }

    #[test]
    fn bouncing_model_newton_post_impact() {
        let bm = BouncingModel {
            model: RestitutionModel::Newton,
            coefficient: 0.5,
            min_speed: 0.0,
        };
        // v_n = -2.0 (approaching)
        let v_after = bm.post_impact_velocity(-2.0);
        assert!(
            (v_after - 1.0).abs() < 1e-10,
            "e=0.5, v_n=-2 → v_after=1.0, got {v_after}"
        );
    }

    #[test]
    fn bouncing_model_elastic() {
        let bm = BouncingModel {
            model: RestitutionModel::Energetic,
            coefficient: 1.0,
            min_speed: 0.0,
        };
        let v_after = bm.post_impact_velocity(-3.0);
        assert!(
            (v_after - 3.0).abs() < 1e-10,
            "perfectly elastic e=1 → v_after=3.0, got {v_after}"
        );
    }

    #[test]
    fn bouncing_model_impact_impulse() {
        let bm = BouncingModel {
            model: RestitutionModel::Newton,
            coefficient: 0.5,
            min_speed: 0.0,
        };
        // j = -(1+e)*v_n*eff_mass = -(1.5)*(-2)*0.5 = 1.5
        let j = bm.impact_impulse(-2.0, 0.5);
        assert!((j - 1.5).abs() < 1e-10, "impulse should be 1.5, got {j}");
    }

    // ---- SoftContact ----

    #[test]
    fn soft_contact_zero_penetration_zero_force() {
        let sc = SoftContact::default();
        let f = sc.normal_force(0.0, 0.0);
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn soft_contact_positive_penetration() {
        let sc = SoftContact {
            stiffness: 1000.0,
            damping: 0.0,
            tangential_stiffness: 500.0,
            max_force: 1e6,
        };
        let f = sc.normal_force(0.01, 0.0);
        assert!((f - 10.0).abs() < 1e-10, "1000*0.01=10, got {f}");
    }

    #[test]
    fn soft_contact_damping_reduces_force() {
        let sc = SoftContact {
            stiffness: 1000.0,
            damping: 100.0,
            tangential_stiffness: 0.0,
            max_force: 1e6,
        };
        // penetration=0.02, closing_vel=0.1 → f = 1000*0.02 - 100*0.1 = 20 - 10 = 10
        let f = sc.normal_force(0.02, 0.1);
        assert!(
            (f - 10.0).abs() < 1e-10,
            "force with damping should be 10, got {f}"
        );
    }

    #[test]
    fn soft_contact_hertz() {
        // Known: F = (4/3) * E* * R^0.5 * delta^1.5
        let e_star = 1.0;
        let r_star = 1.0;
        let delta = 1.0;
        let f = SoftContact::hertz_force(e_star, r_star, delta);
        assert!(
            (f - 4.0 / 3.0).abs() < 1e-12,
            "Hertz force should be 4/3 for unit inputs, got {f}"
        );
    }

    #[test]
    fn soft_contact_max_force_clamp() {
        let sc = SoftContact {
            stiffness: 1.0e10,
            damping: 0.0,
            tangential_stiffness: 0.0,
            max_force: 100.0,
        };
        let f = sc.normal_force(0.1, 0.0);
        assert!(
            (f - 100.0).abs() < 1e-10,
            "force should be clamped to max_force"
        );
    }

    // ---- RollingResistance ----

    #[test]
    fn rolling_resistance_torque_magnitude() {
        let rr = RollingResistance {
            mu_rolling: 0.01,
            mu_spinning: 0.005,
            contact_radius: 0.5,
        };
        // T = 0.01 * 100 * 0.5 = 0.5
        let t = rr.rolling_torque(100.0);
        assert!(
            (t - 0.5).abs() < 1e-12,
            "rolling torque should be 0.5, got {t}"
        );
    }

    #[test]
    fn rolling_resistance_torque_vector_opposing() {
        let rr = RollingResistance::default();
        let omega = [1.0, 0.0, 0.0];
        let tv = rr.torque_vector(omega, 100.0);
        // Torque should oppose omega (x-direction), so tv[0] < 0
        assert!(tv[0] < 0.0, "rolling torque should oppose angular velocity");
    }

    #[test]
    fn rolling_resistance_zero_omega_zero_torque() {
        let rr = RollingResistance::default();
        let tv = rr.torque_vector([0.0; 3], 100.0);
        let mag = (tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2]).sqrt();
        assert!(mag < 1e-12, "zero omega should produce zero torque");
    }

    // ---- ContactMaterial ----

    #[test]
    fn contact_material_register_and_get() {
        let mut cm = ContactMaterial::new();
        let mat = Material {
            name: "steel".into(),
            mu_static: 0.7,
            mu_dynamic: 0.5,
            ..Default::default()
        };
        let id = cm.register(mat);
        assert!(cm.get(id).is_some());
        assert_eq!(cm.get(id).unwrap().name, "steel");
    }

    #[test]
    fn contact_material_combined_geometric_mean() {
        let mut cm = ContactMaterial::new();
        let id_a = cm.register(Material {
            mu_static: 0.4,
            mu_dynamic: 0.3,
            ..Default::default()
        });
        let id_b = cm.register(Material {
            mu_static: 0.9,
            mu_dynamic: 0.6,
            ..Default::default()
        });
        let props = cm.combined(id_a, id_b);
        let expected_s = (0.4_f64 * 0.9).sqrt();
        let expected_d = (0.3_f64 * 0.6).sqrt();
        assert!((props.friction_static - expected_s).abs() < 1e-12);
        assert!((props.friction_dynamic - expected_d).abs() < 1e-12);
    }

    #[test]
    fn contact_material_restitution_minimum() {
        let mut cm = ContactMaterial::new();
        let id_a = cm.register(Material {
            restitution: 0.8,
            ..Default::default()
        });
        let id_b = cm.register(Material {
            restitution: 0.3,
            ..Default::default()
        });
        let props = cm.combined(id_a, id_b);
        assert!(
            (props.restitution - 0.3).abs() < 1e-12,
            "restitution should be min(0.8,0.3)=0.3"
        );
    }

    #[test]
    fn contact_material_override() {
        let mut cm = ContactMaterial::new();
        let id_a = cm.register(Material::default());
        let id_b = cm.register(Material::default());
        let override_props = PairedProperties {
            friction_static: 1.0,
            friction_dynamic: 0.9,
            restitution: 0.0,
            mu_rolling: 0.001,
        };
        cm.set_override(id_a, id_b, override_props.clone());
        let props = cm.combined(id_a, id_b);
        assert!(
            (props.friction_static - 1.0).abs() < 1e-12,
            "override should be used"
        );
    }

    #[test]
    fn contact_material_len_and_empty() {
        let mut cm = ContactMaterial::new();
        assert!(cm.is_empty());
        cm.register(Material::default());
        assert_eq!(cm.len(), 1);
        assert!(!cm.is_empty());
    }

    // ---- contact_tangent_frame ----

    #[test]
    fn tangent_frame_orthogonal() {
        let n = [0.0, 1.0, 0.0f64];
        let (t1, t2) = contact_tangent_frame(n);
        // t1 and t2 should be perpendicular to n and to each other
        let d_n_t1 = v3_dot(n, t1).abs();
        let d_n_t2 = v3_dot(n, t2).abs();
        let d_t1_t2 = v3_dot(t1, t2).abs();
        assert!(d_n_t1 < 1e-10, "t1 should be perpendicular to n");
        assert!(d_n_t2 < 1e-10, "t2 should be perpendicular to n");
        assert!(d_t1_t2 < 1e-10, "t1 should be perpendicular to t2");
    }

    #[test]
    fn tangent_frame_unit_length() {
        let n = [1.0f64, 0.0, 0.0];
        let (t1, t2) = contact_tangent_frame(n);
        let l1 = v3_norm(t1);
        let l2 = v3_norm(t2);
        assert!((l1 - 1.0).abs() < 1e-10, "t1 should be unit length");
        assert!((l2 - 1.0).abs() < 1e-10, "t2 should be unit length");
    }
}
