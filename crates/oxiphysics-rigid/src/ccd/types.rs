//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    ConvexSupport, add, dot, norm, normalize, scale, speculative_contact_candidate, sub,
};

/// A rigid body represented for CCD swept-volume tests.
#[derive(Debug, Clone, PartialEq)]
pub struct SweptBody {
    /// World-space position (metres).
    pub pos: [f64; 3],
    /// Linear velocity (m/s).
    pub vel: [f64; 3],
    /// Bounding sphere radius (metres).
    pub radius: f64,
    /// Inverse mass (1/kg); zero for static/infinite-mass bodies.
    pub inv_mass: f64,
}
impl SweptBody {
    /// Construct a new [`SweptBody`].
    pub fn new(pos: [f64; 3], vel: [f64; 3], radius: f64, inv_mass: f64) -> Self {
        Self {
            pos,
            vel,
            radius,
            inv_mass,
        }
    }
    /// Position advanced by `dt` seconds (linear only).
    pub fn position_at(&self, t: f64) -> [f64; 3] {
        add(self.pos, scale(self.vel, t))
    }
}
/// Swept collision detection utilities for sphere primitives and tunneling
/// analysis.
///
/// Provides TOI computation, substep count selection, and velocity-based
/// tunneling detection as a self-contained, allocation-free utility type.
pub struct CcdSweep;
impl CcdSweep {
    /// Compute the Time Of Impact (TOI) for two moving spheres.
    ///
    /// Solves the quadratic |r(t)|² = (r_a + r_b)² where
    /// r(t) = (pos_a - pos_b) + (vel_a - vel_b)·t.
    ///
    /// Returns the earliest positive root ∈ `[0, dt]`, or `None` if the
    /// spheres do not collide within the interval.
    ///
    /// - `pos_a`, `pos_b` : initial world positions (m)
    /// - `vel_a`, `vel_b` : linear velocities (m/s)
    /// - `radius_a`, `radius_b` : bounding sphere radii (m)
    /// - `dt` : time step to search within (s)
    pub fn compute_collision_time_sphere(
        pos_a: [f64; 3],
        vel_a: [f64; 3],
        radius_a: f64,
        pos_b: [f64; 3],
        vel_b: [f64; 3],
        radius_b: f64,
        dt: f64,
    ) -> Option<f64> {
        let r = sub(pos_a, pos_b);
        let v = sub(vel_a, vel_b);
        let sum_r = radius_a + radius_b;
        let aa = dot(v, v);
        let bb = 2.0 * dot(r, v);
        let cc = dot(r, r) - sum_r * sum_r;
        if cc <= 0.0 {
            return Some(0.0);
        }
        if aa < 1e-30 {
            return None;
        }
        let disc: f64 = bb * bb - 4.0 * aa * cc;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t0 = (-bb - sqrt_disc) / (2.0 * aa);
        let t1 = (-bb + sqrt_disc) / (2.0 * aa);
        [t0, t1].iter().find(|&&t| t >= 0.0 && t <= dt).copied()
    }
    /// Compute how many CCD substeps are needed so that no sphere travels
    /// more than a fraction `safety_fraction` of its radius per substep.
    ///
    /// `max_speed` is the maximum relative speed of any body pair (m/s),
    /// `min_radius` is the smallest bounding radius in the scene (m),
    /// `dt` is the full time step (s).
    ///
    /// The result is clamped to `[1, max_substeps]`.
    pub fn compute_sub_step_count(
        max_speed: f64,
        min_radius: f64,
        dt: f64,
        safety_fraction: f64,
        max_substeps: usize,
    ) -> usize {
        if max_speed < 1e-30 || min_radius < 1e-30 || dt < 1e-30 {
            return 1;
        }
        let safe_frac = safety_fraction.max(1e-6);
        let n_f = (max_speed * dt) / (safe_frac * min_radius);
        let n = n_f.ceil() as usize;
        n.clamp(1, max_substeps)
    }
    /// Detect whether any body in `bodies` is tunneling based on its velocity.
    ///
    /// A body is considered to be tunneling if its speed over `dt` exceeds
    /// `threshold` times its radius (i.e., it would travel more than
    /// `threshold` radii in one step).
    ///
    /// Returns a list of indices of bodies exhibiting potential tunneling.
    pub fn detect_tunneling(bodies: &[SweptBody], dt: f64, threshold: f64) -> Vec<usize> {
        let mut tunneling = Vec::new();
        for (i, body) in bodies.iter().enumerate() {
            if body.radius < 1e-30 {
                continue;
            }
            let speed = norm(body.vel);
            let travel = speed * dt;
            if travel > threshold * body.radius {
                tunneling.push(i);
            }
        }
        tunneling
    }
}
/// Swept body with angular velocity.
///
/// The bounding sphere grows conservatively by the tangential speed at the
/// outer radius: `|omega| * radius`.
#[derive(Debug, Clone)]
pub struct RotatingSweptBody {
    /// Linear position.
    pub pos: [f64; 3],
    /// Linear velocity.
    pub vel: [f64; 3],
    /// Angular velocity (axis-angle, rad/s).
    pub omega: [f64; 3],
    /// Bounding sphere radius.
    pub radius: f64,
    /// Inverse mass.
    pub inv_mass: f64,
}
impl RotatingSweptBody {
    /// Create a new rotating swept body.
    pub fn new(pos: [f64; 3], vel: [f64; 3], omega: [f64; 3], radius: f64, inv_mass: f64) -> Self {
        Self {
            pos,
            vel,
            omega,
            radius,
            inv_mass,
        }
    }
    /// Conservative bounding sphere radius accounting for rotation over `dt`.
    ///
    /// The surface tangential speed is `|omega| * radius`, so the effective
    /// radius grows by that amount per second.
    pub fn conservative_radius(&self, dt: f64) -> f64 {
        self.radius + norm(self.omega) * self.radius * dt
    }
    /// Linear position at time `t`.
    pub fn position_at(&self, t: f64) -> [f64; 3] {
        add(self.pos, scale(self.vel, t))
    }
}
/// Accumulated statistics for a CCD simulation pass.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct CcdStats {
    /// Number of body pairs evaluated.
    pub pairs_tested: usize,
    /// Number of TOI events found.
    pub events_found: usize,
    /// Number of impulses applied.
    pub impulses_applied: usize,
    /// Earliest TOI found (or f64::INFINITY if none).
    pub earliest_toi: f64,
    /// Total simulation time covered by CCD sub-steps.
    pub time_covered: f64,
}
impl CcdStats {
    /// Create fresh stats.
    pub fn new() -> Self {
        Self {
            earliest_toi: f64::INFINITY,
            ..Default::default()
        }
    }
    /// Record a new TOI event.
    pub fn record_event(&mut self, toi: f64) {
        self.events_found += 1;
        if toi < self.earliest_toi {
            self.earliest_toi = toi;
        }
    }
    /// Record a pair evaluation.
    pub fn record_pair_test(&mut self) {
        self.pairs_tested += 1;
    }
    /// Record an impulse application.
    pub fn record_impulse(&mut self) {
        self.impulses_applied += 1;
    }
    /// Record time advanced.
    pub fn record_time(&mut self, dt: f64) {
        self.time_covered += dt;
    }
    /// Reset stats to zero.
    pub fn reset(&mut self) {
        *self = Self::new();
    }
}
/// A capsule defined by two endpoint positions and a radius.
///
/// The capsule is the Minkowski sum of a line segment with a sphere of the
/// given radius.
#[derive(Debug, Clone, PartialEq)]
pub struct SweptCapsule {
    /// Start endpoint of the capsule axis in world space.
    pub p0: [f64; 3],
    /// End endpoint of the capsule axis in world space.
    pub p1: [f64; 3],
    /// Capsule radius.
    pub radius: f64,
    /// Linear velocity applied to both endpoints equally.
    pub vel: [f64; 3],
    /// Inverse mass.
    pub inv_mass: f64,
}
impl SweptCapsule {
    /// Construct a new swept capsule.
    pub fn new(p0: [f64; 3], p1: [f64; 3], radius: f64, vel: [f64; 3], inv_mass: f64) -> Self {
        Self {
            p0,
            p1,
            radius,
            vel,
            inv_mass,
        }
    }
    /// Centre position of the capsule (midpoint of the axis).
    pub fn center(&self) -> [f64; 3] {
        [
            (self.p0[0] + self.p1[0]) * 0.5,
            (self.p0[1] + self.p1[1]) * 0.5,
            (self.p0[2] + self.p1[2]) * 0.5,
        ]
    }
    /// Advance both endpoints by velocity × `t`.
    pub fn advance(&self, t: f64) -> SweptCapsule {
        SweptCapsule {
            p0: add(self.p0, scale(self.vel, t)),
            p1: add(self.p1, scale(self.vel, t)),
            radius: self.radius,
            vel: self.vel,
            inv_mass: self.inv_mass,
        }
    }
    /// Length of the capsule axis.
    pub fn length(&self) -> f64 {
        norm(sub(self.p1, self.p0))
    }
    /// Axis direction (unit vector from p0 to p1), or zero if degenerate.
    pub fn axis(&self) -> [f64; 3] {
        normalize(sub(self.p1, self.p0))
    }
    /// Bounding sphere that conservatively encloses the capsule.
    ///
    /// The sphere is centred at the capsule midpoint with radius
    /// `half_length + self.radius`.
    pub fn bounding_sphere_radius(&self) -> f64 {
        self.length() * 0.5 + self.radius
    }
}
/// A body slot in the CCD broadphase, combining swept geometry and filter.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CcdBodySlot {
    /// The swept body.
    pub body: SweptBody,
    /// ID assigned by the broadphase.
    pub id: usize,
    /// CCD filter.
    pub filter: CcdFilter,
    /// Whether CCD is enabled for this body.
    pub ccd_enabled: bool,
}
impl CcdBodySlot {
    /// Create a new CCD body slot.
    pub fn new(body: SweptBody, id: usize, filter: CcdFilter) -> Self {
        Self {
            body,
            id,
            filter,
            ccd_enabled: true,
        }
    }
    /// Disable CCD for this body (treated as discrete).
    pub fn disable_ccd(mut self) -> Self {
        self.ccd_enabled = false;
        self
    }
}
/// An AABB as a convex support.
#[derive(Debug, Clone)]
pub struct AabbSupport {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
/// Result of the GJK distance query.
#[derive(Debug, Clone)]
pub struct GjkResult {
    /// Minimum squared distance between the two shapes.
    pub sq_dist: f64,
    /// Closest point on shape A.
    pub closest_a: [f64; 3],
    /// Closest point on shape B.
    pub closest_b: [f64; 3],
    /// Whether the shapes overlap (penetrating).
    pub overlapping: bool,
}
/// A CCD broadphase that produces filtered candidate pairs.
#[allow(dead_code)]
pub struct CcdBroadphase {
    /// The time step.
    pub dt: f64,
}
impl CcdBroadphase {
    /// Create a broadphase with the given timestep.
    pub fn new(dt: f64) -> Self {
        Self { dt }
    }
    /// Collect candidate pairs from a list of body slots.
    ///
    /// Returns `(slot_a_id, slot_b_id)` pairs where:
    /// - Both slots have CCD enabled.
    /// - Their filters pass.
    /// - Their swept volumes overlap (speculative contact check).
    pub fn candidate_pairs(&self, slots: &[CcdBodySlot]) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..slots.len() {
            for j in (i + 1)..slots.len() {
                let sa = &slots[i];
                let sb = &slots[j];
                if !sa.ccd_enabled && !sb.ccd_enabled {
                    continue;
                }
                if !CcdFilter::should_check(&sa.filter, &sb.filter) {
                    continue;
                }
                if speculative_contact_candidate(&sa.body, &sb.body, self.dt) {
                    pairs.push((sa.id, sb.id));
                }
            }
        }
        pairs
    }
}
/// Broad + narrow CCD pipeline operating over a slice of [`SweptBody`]s.
pub struct CcdPipeline {
    /// Simulation sub-step size (seconds).
    pub dt: f64,
}
impl CcdPipeline {
    /// Create a new pipeline with the given time step.
    pub fn new(dt: f64) -> Self {
        Self { dt }
    }
    /// O(n²) all-pairs tunnelling detection.
    ///
    /// Returns a list of [`CcdEvent`]s sorted by ascending TOI.
    /// `ids` maps slice indices to user-defined body IDs stored in the event.
    pub fn detect_tunneling(&self, bodies: &[SweptBody], ids: &[usize]) -> Vec<CcdEvent> {
        let mut events = Vec::new();
        let n = bodies.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(toi) = TunnelTest::check_sphere_sphere(&bodies[i], &bodies[j], self.dt)
                {
                    let pa = bodies[i].position_at(toi);
                    let pb = bodies[j].position_at(toi);
                    let normal = normalize(sub(pa, pb));
                    let contact_point = [
                        0.5 * (pa[0] + pb[0]),
                        0.5 * (pa[1] + pb[1]),
                        0.5 * (pa[2] + pb[2]),
                    ];
                    events.push(CcdEvent {
                        toi,
                        body_a: ids[i],
                        body_b: ids[j],
                        contact_normal: normal,
                        contact_point,
                    });
                }
            }
        }
        events.sort_by(|x, y| {
            x.toi
                .partial_cmp(&y.toi)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        events
    }
    /// Return only the earliest [`CcdEvent`], if any.
    pub fn earliest_event(&self, bodies: &[SweptBody], ids: &[usize]) -> Option<CcdEvent> {
        self.detect_tunneling(bodies, ids).into_iter().next()
    }
}
impl CcdPipeline {
    /// Run a full CCD pass: detect, resolve, and return statistics.
    ///
    /// Advances all bodies to the earliest TOI, resolves the contact, then
    /// continues until `dt` is consumed or no more events occur.
    pub fn run_with_stats(
        &self,
        bodies: &mut Vec<SweptBody>,
        ids: &[usize],
        restitution: f64,
    ) -> CcdStats {
        let mut stats = CcdStats::new();
        let integrator = SubstepIntegrator::new(1);
        let mut remaining = self.dt;
        let mut local_bodies = bodies.clone();
        let local_ids = ids.to_vec();
        let n = local_ids.len();
        stats.pairs_tested += n * n.saturating_sub(1) / 2;
        let mut iter_cap = 64usize;
        while remaining > 1e-10 && iter_cap > 0 {
            iter_cap -= 1;
            let sub_pipeline = CcdPipeline::new(remaining);
            let events = sub_pipeline.detect_tunneling(&local_bodies, &local_ids);
            if events.is_empty() {
                integrator.step_to_toi(&mut local_bodies, remaining);
                stats.record_time(remaining);
                break;
            }
            let ev = events[0].clone();
            stats.record_event(ev.toi);
            let advance = ev.toi.max(1e-8);
            integrator.step_to_toi(&mut local_bodies, advance);
            integrator.resolve_event(&mut local_bodies, &ev, restitution);
            stats.record_impulse();
            stats.record_time(advance);
            remaining -= advance;
        }
        *bodies = local_bodies;
        let _ = local_ids;
        stats
    }
}
/// A sphere as a convex support.
#[derive(Debug, Clone)]
pub struct SphereSupport {
    /// Centre of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere.
    pub radius: f64,
}
/// Helper wrapper that translates a support function by an offset.
pub struct TranslatedSupport<'a, S: ConvexSupport> {
    pub(super) shape: &'a S,
    pub(super) offset: [f64; 3],
}
/// A detected continuous-collision event.
#[derive(Debug, Clone, PartialEq)]
pub struct CcdEvent {
    /// Time of impact (seconds from current time).
    pub toi: f64,
    /// Index into the body slice for the first body.
    pub body_a: usize,
    /// Index into the body slice for the second body.
    pub body_b: usize,
    /// Contact normal pointing from B toward A.
    pub contact_normal: [f64; 3],
    /// World-space contact point.
    pub contact_point: [f64; 3],
}
/// Filter flags for CCD broadphase pairs.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CcdFilter {
    /// Category bitmask of the body.
    pub category: u32,
    /// Collision mask (collides with bodies whose category ∩ mask ≠ 0).
    pub mask: u32,
}
impl CcdFilter {
    /// Create a new filter.
    pub fn new(category: u32, mask: u32) -> Self {
        Self { category, mask }
    }
    /// Test whether two filters would produce a CCD pair.
    pub fn should_check(a: &CcdFilter, b: &CcdFilter) -> bool {
        (a.category & b.mask) != 0 || (b.category & a.mask) != 0
    }
}
/// Result of a ray–sphere intersection test.
#[derive(Debug, Clone, PartialEq)]
pub struct RaySphereHit {
    /// Parameter `t` along the ray at which the first intersection occurs.
    pub t: f64,
    /// World-space hit point.
    pub point: [f64; 3],
    /// Outward surface normal at the hit point.
    pub normal: [f64; 3],
}
/// Coefficient-of-restitution models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RestitutionModel {
    /// Constant restitution coefficient (Newton's model).
    Constant(f64),
    /// Restitution varies linearly with impact speed:
    /// `e = e_max * (1 - v / v_max)`, clamped to `[0, e_max]`.
    SpeedDependent {
        /// Maximum restitution at zero impact speed.
        e_max: f64,
        /// Speed at which restitution falls to zero.
        v_max: f64,
    },
}
impl RestitutionModel {
    /// Evaluate the restitution coefficient for the given impact speed.
    pub fn evaluate(&self, impact_speed: f64) -> f64 {
        match *self {
            RestitutionModel::Constant(e) => e.clamp(0.0, 1.0),
            RestitutionModel::SpeedDependent { e_max, v_max } => {
                let e = e_max * (1.0 - impact_speed / v_max);
                e.clamp(0.0, e_max)
            }
        }
    }
}
/// Utility methods for detecting tunnelling (CCD narrow-phase).
pub struct TunnelTest;
impl TunnelTest {
    /// Sphere–sphere continuous collision test.
    ///
    /// Returns the time of impact `t ∈ [0, dt]` when the two spheres first
    /// touch, or `None` if they do not collide within `[0, dt]`.
    pub fn check_sphere_sphere(a: &SweptBody, b: &SweptBody, dt: f64) -> Option<f64> {
        let r = sub(a.pos, b.pos);
        let v = sub(a.vel, b.vel);
        let sum_r = a.radius + b.radius;
        let aa = dot(v, v);
        let bb = 2.0 * dot(r, v);
        let cc = dot(r, r) - sum_r * sum_r;
        if cc <= 0.0 {
            return Some(0.0);
        }
        if aa < 1e-30 {
            return None;
        }
        let disc: f64 = bb * bb - 4.0 * aa * cc;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t = (-bb - sqrt_disc) / (2.0 * aa);
        if t >= 0.0 && t <= dt { Some(t) } else { None }
    }
    /// Sphere vs. static infinite plane (`n·x = d`).
    ///
    /// Returns the TOI in `[0, dt]` when the sphere surface first touches the
    /// plane, or `None` if no contact occurs in that interval.
    pub fn check_sphere_plane(
        sphere: &SweptBody,
        plane_normal: [f64; 3],
        plane_d: f64,
        dt: f64,
    ) -> Option<f64> {
        let n = normalize(plane_normal);
        let dist = dot(n, sphere.pos) - plane_d - sphere.radius;
        if dist <= 0.0 {
            return Some(0.0);
        }
        let closing = -dot(n, sphere.vel);
        if closing <= 1e-12 {
            return None;
        }
        let t = dist / closing;
        if t <= dt { Some(t) } else { None }
    }
}
/// Substep integrator that advances bodies to a TOI and resolves collisions.
pub struct SubstepIntegrator {
    /// Number of substeps per full time-step.
    pub n_substeps: usize,
}
impl SubstepIntegrator {
    /// Create an integrator with the given number of substeps.
    pub fn new(n_substeps: usize) -> Self {
        Self { n_substeps }
    }
    /// Advance all bodies linearly to time `toi`.
    pub fn step_to_toi(&self, bodies: &mut Vec<SweptBody>, toi: f64) {
        for b in bodies.iter_mut() {
            b.pos = b.position_at(toi);
        }
    }
    /// Apply an impulsive response at the contact described by `event`.
    ///
    /// Uses a simple coefficient-of-restitution impulse along the contact
    /// normal.  Bodies with `inv_mass == 0` (static) are unaffected.
    pub fn resolve_event(&self, bodies: &mut Vec<SweptBody>, event: &CcdEvent, restitution: f64) {
        let idx_a = bodies.iter().position(|_| true).map(|_| event.body_a);
        let idx_b = Some(event.body_b);
        let (ia, ib) = match (idx_a, idx_b) {
            (Some(a), Some(b)) if a < bodies.len() && b < bodies.len() => (a, b),
            _ => return,
        };
        let n = event.contact_normal;
        let rel_vel = sub(bodies[ia].vel, bodies[ib].vel);
        let vn = dot(rel_vel, n);
        if vn >= 0.0 {
            return;
        }
        let inv_a = bodies[ia].inv_mass;
        let inv_b = bodies[ib].inv_mass;
        let denom = inv_a + inv_b;
        if denom < 1e-30 {
            return;
        }
        let j = -(1.0 + restitution) * vn / denom;
        let impulse = scale(n, j);
        bodies[ia].vel = add(bodies[ia].vel, scale(impulse, inv_a));
        bodies[ib].vel = sub(bodies[ib].vel, scale(impulse, inv_b));
    }
}
