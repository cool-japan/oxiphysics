//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
use crate::body::{BodyState, BodyType, RigidBody};
use crate::collider::Collider;
use crate::impulse::SubStepRestitutionCorrector;
use crate::pipeline::BodySnapshot;
use crate::sets::{ColliderSet, RigidBodySet};
use crate::solver::soft::SoftParams;
use oxiphysics_core::math::{Mat3, Vec3};
use oxiphysics_core::{BodyHandle, ColliderHandle};

/// A simple rigid body world that owns a flat list of bodies and steps them
/// forward with semi-implicit Euler integration under constant gravity.
pub struct RigidWorld {
    /// All bodies currently in the simulation.
    pub bodies: Vec<BodySnapshot>,
    /// Gravitational acceleration \[x, y, z\] in m/s².
    pub gravity: [f64; 3],
    /// Fixed time-step in seconds.
    pub dt: f64,
    /// Accumulated simulation time in seconds.
    pub time: f64,
    /// Monotonically-increasing ID counter.
    pub(super) next_id: u64,
}
impl RigidWorld {
    /// Creates a new `RigidWorld` with the given gravity and time-step.
    pub fn new(gravity: [f64; 3], dt: f64) -> Self {
        Self {
            bodies: Vec::new(),
            gravity,
            dt,
            time: 0.0,
            next_id: 0,
        }
    }
    /// Adds a body with the given position, velocity, and inverse mass.
    ///
    /// Returns the unique ID assigned to the new body.
    pub fn add_body(&mut self, pos: [f64; 3], vel: [f64; 3], inv_mass: f64) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        self.bodies.push(BodySnapshot {
            id,
            position: pos,
            velocity: vel,
            ang_velocity: [0.0, 0.0, 0.0],
            inv_mass,
            active: true,
        });
        id
    }
    /// Returns an immutable reference to the body with the given ID, if it exists.
    pub fn get_body(&self, id: u64) -> Option<&BodySnapshot> {
        self.bodies.iter().find(|b| b.id == id)
    }
    /// Returns a mutable reference to the body with the given ID, if it exists.
    pub fn get_body_mut(&mut self, id: u64) -> Option<&mut BodySnapshot> {
        self.bodies.iter_mut().find(|b| b.id == id)
    }
    /// Removes the body with the given ID.
    ///
    /// Returns `true` if a body was found and removed, `false` otherwise.
    pub fn remove_body(&mut self, id: u64) -> bool {
        if let Some(pos) = self.bodies.iter().position(|b| b.id == id) {
            self.bodies.swap_remove(pos);
            true
        } else {
            false
        }
    }
    /// Advances the simulation by one time-step.
    ///
    /// For every active dynamic body (inv_mass > 0) gravity is applied first,
    /// then positions are integrated with semi-implicit Euler.
    pub fn step(&mut self) {
        let g = self.gravity;
        let dt = self.dt;
        for body in self.bodies.iter_mut() {
            if !body.active {
                continue;
            }
            if body.inv_mass > 0.0 {
                body.velocity[0] += g[0] * dt;
                body.velocity[1] += g[1] * dt;
                body.velocity[2] += g[2] * dt;
            }
            body.position[0] += body.velocity[0] * dt;
            body.position[1] += body.velocity[1] * dt;
            body.position[2] += body.velocity[2] * dt;
        }
        self.time += dt;
    }
    /// Applies an impulse to the body with the given ID.
    ///
    /// The velocity change is `impulse * inv_mass`.
    pub fn apply_impulse_to_body(&mut self, id: u64, impulse: [f64; 3]) {
        if let Some(body) = self.get_body_mut(id) {
            let m = body.inv_mass;
            body.velocity[0] += impulse[0] * m;
            body.velocity[1] += impulse[1] * m;
            body.velocity[2] += impulse[2] * m;
        }
    }
    /// Returns the number of bodies currently in the world.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Computes the total kinetic energy: Σ 0.5 * m * |v|².
    ///
    /// Bodies with `inv_mass == 0` (static) contribute zero.
    pub fn kinetic_energy(&self) -> f64 {
        self.bodies
            .iter()
            .filter(|b| b.inv_mass > 0.0)
            .map(|b| {
                let v = b.velocity;
                let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                let mass = 1.0 / b.inv_mass;
                0.5 * mass * v_sq
            })
            .sum()
    }
}
impl RigidWorld {
    /// Export the world to a snapshot.
    pub fn to_snapshot(&self) -> RigidWorldSnapshot {
        RigidWorldSnapshot {
            bodies: self.bodies.clone(),
            gravity: self.gravity,
            dt: self.dt,
            time: self.time,
        }
    }
    /// Restore a world from a snapshot.
    pub fn from_snapshot(snap: RigidWorldSnapshot) -> Self {
        let next_id = snap.bodies.iter().map(|b| b.id + 1).max().unwrap_or(0);
        Self {
            bodies: snap.bodies,
            gravity: snap.gravity,
            dt: snap.dt,
            time: snap.time,
            next_id,
        }
    }
    /// Number of active (non-sleeping) bodies.
    pub fn active_body_count(&self) -> usize {
        self.bodies.iter().filter(|b| b.active).count()
    }
    /// Deactivate a body, preventing it from being stepped.
    pub fn deactivate_body(&mut self, id: u64) {
        if let Some(b) = self.bodies.iter_mut().find(|b| b.id == id) {
            b.active = false;
        }
    }
    /// Reactivate a body.
    pub fn activate_body(&mut self, id: u64) {
        if let Some(b) = self.bodies.iter_mut().find(|b| b.id == id) {
            b.active = true;
        }
    }
    /// Apply the same force (as an impulse this step) to all active dynamic bodies.
    pub fn apply_global_impulse(&mut self, impulse: [f64; 3]) {
        for body in self.bodies.iter_mut() {
            if body.active && body.inv_mass > 0.0 {
                body.velocity[0] += impulse[0] * body.inv_mass;
                body.velocity[1] += impulse[1] * body.inv_mass;
                body.velocity[2] += impulse[2] * body.inv_mass;
            }
        }
    }
    /// Returns the body with the highest kinetic energy, or `None` if empty.
    pub fn most_energetic_body(&self) -> Option<&BodySnapshot> {
        self.bodies
            .iter()
            .filter(|b| b.inv_mass > 0.0)
            .max_by(|a, b| {
                let ke = |body: &&BodySnapshot| {
                    let v = body.velocity;
                    v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
                };
                ke(a)
                    .partial_cmp(&ke(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }
    /// Step multiple times.
    pub fn step_n(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }
}
impl RigidWorld {
    /// Cast a ray through the world and collect **all** hit bodies, sorted by
    /// ascending `t`.
    ///
    /// Bodies are treated as spheres; the radius is `1.0 / inv_mass.cbrt() * 0.3`
    /// (the same heuristic used in the broadphase), clamped to at least 0.1 m.
    pub fn raycast_all(&self, ray: &Ray) -> Vec<RaycastResult> {
        let mut hits = Vec::new();
        for body in &self.bodies {
            if !body.active {
                continue;
            }
            let r = if body.inv_mass > 0.0 {
                (1.0 / body.inv_mass).cbrt() * 0.3
            } else {
                0.5_f64
            }
            .max(0.1);
            if let Some(t) = ray_sphere_intersect(ray, body.position, r) {
                let hit_point = ray.at(t);
                let nx = (hit_point[0] - body.position[0]) / r;
                let ny = (hit_point[1] - body.position[1]) / r;
                let nz = (hit_point[2] - body.position[2]) / r;
                hits.push(RaycastResult {
                    body_id: body.id,
                    hit_point,
                    normal: [nx, ny, nz],
                    t,
                });
            }
        }
        hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        hits
    }
    /// Cast a ray and return the **nearest** hit, or `None`.
    pub fn raycast_nearest(&self, ray: &Ray) -> Option<RaycastResult> {
        self.raycast_all(ray).into_iter().next()
    }
}
impl RigidWorld {
    /// Sweep a sphere of `query_radius` along `ray` and return all bodies
    /// whose bounding spheres it first touches, sorted by ascending `t`.
    ///
    /// Implemented by expanding each body's bounding sphere by `query_radius`
    /// and performing a standard ray–sphere test.
    pub fn sphere_cast(&self, ray: &Ray, query_radius: f64) -> Vec<SphereCastResult> {
        let mut results = Vec::new();
        for body in &self.bodies {
            if !body.active {
                continue;
            }
            let r = if body.inv_mass > 0.0 {
                (1.0 / body.inv_mass).cbrt() * 0.3
            } else {
                0.5_f64
            }
            .max(0.1);
            let expanded = r + query_radius;
            if let Some(t) = ray_sphere_intersect(ray, body.position, expanded) {
                let hit_centre = ray.at(t);
                results.push(SphereCastResult {
                    body_id: body.id,
                    hit_centre,
                    t,
                });
            }
        }
        results.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        results
    }
}
impl RigidWorld {
    /// Return the IDs of all active bodies whose centres lie inside `query`.
    pub fn aabb_overlap_query(&self, query: &AabbQuery) -> Vec<u64> {
        self.bodies
            .iter()
            .filter(|b| b.active && query.contains_point(b.position))
            .map(|b| b.id)
            .collect()
    }
    /// Return the IDs of all active bodies whose bounding spheres overlap `query`.
    pub fn aabb_bounding_sphere_query(&self, query: &AabbQuery) -> Vec<u64> {
        self.bodies
            .iter()
            .filter(|b| {
                if !b.active {
                    return false;
                }
                let r = if b.inv_mass > 0.0 {
                    (1.0 / b.inv_mass).cbrt() * 0.3
                } else {
                    0.5
                }
                .max(0.1);
                let sphere_aabb = AabbQuery::new(
                    [b.position[0] - r, b.position[1] - r, b.position[2] - r],
                    [b.position[0] + r, b.position[1] + r, b.position[2] + r],
                );
                query.overlaps(&sphere_aabb)
            })
            .map(|b| b.id)
            .collect()
    }
}
impl RigidWorld {
    /// Add multiple bodies from a slice of [`BodyConfig`]s.
    ///
    /// Returns the IDs assigned to the new bodies in the same order as the
    /// input configs.
    pub fn add_bodies_from_configs(&mut self, configs: &[BodyConfig]) -> Vec<u64> {
        configs
            .iter()
            .map(|cfg| {
                let id = self.next_id;
                self.next_id += 1;
                self.bodies.push(BodySnapshot {
                    id,
                    position: cfg.position,
                    velocity: cfg.velocity,
                    ang_velocity: [0.0; 3],
                    inv_mass: cfg.inv_mass,
                    active: cfg.active,
                });
                id
            })
            .collect()
    }
}
impl RigidWorld {
    /// Serialise the world state to a flat `Vec`f64`.
    ///
    /// Layout per body (9 values):
    /// `\[id_as_f64, px, py, pz, vx, vy, vz, inv_mass, active_as_f64\]`
    ///
    /// Followed by a 3-element gravity vector and 1-element time:
    /// `\[gx, gy, gz, time\]`
    ///
    /// Total: `9 * body_count + 4` values.
    pub fn to_state_vector(&self) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.bodies.len() * 9 + 4);
        for b in &self.bodies {
            v.push(b.id as f64);
            v.push(b.position[0]);
            v.push(b.position[1]);
            v.push(b.position[2]);
            v.push(b.velocity[0]);
            v.push(b.velocity[1]);
            v.push(b.velocity[2]);
            v.push(b.inv_mass);
            v.push(if b.active { 1.0 } else { 0.0 });
        }
        v.push(self.gravity[0]);
        v.push(self.gravity[1]);
        v.push(self.gravity[2]);
        v.push(self.time);
        v
    }
    /// Restore a world from a state vector produced by [`RigidWorld::to_state_vector`].
    ///
    /// # Panics
    /// Panics if the vector length is not consistent with `9 * N + 4` for some N.
    pub fn from_state_vector(v: &[f64], dt: f64) -> Self {
        assert!(v.len() >= 4, "state vector too short");
        let trailer_start = v.len() - 4;
        let body_count = trailer_start / 9;
        assert_eq!(body_count * 9 + 4, v.len(), "state vector length mismatch");
        let gravity = [v[trailer_start], v[trailer_start + 1], v[trailer_start + 2]];
        let time = v[trailer_start + 3];
        let mut bodies = Vec::with_capacity(body_count);
        let mut max_id = 0u64;
        for i in 0..body_count {
            let base = i * 9;
            let id = v[base] as u64;
            if id > max_id {
                max_id = id;
            }
            bodies.push(BodySnapshot {
                id,
                position: [v[base + 1], v[base + 2], v[base + 3]],
                velocity: [v[base + 4], v[base + 5], v[base + 6]],
                ang_velocity: [0.0; 3],
                inv_mass: v[base + 7],
                active: v[base + 8] > 0.5,
            });
        }
        Self {
            bodies,
            gravity,
            dt,
            time,
            next_id: max_id + 1,
        }
    }
}
impl RigidWorld {
    /// Apply all gravity regions to the world's bodies for one time-step.
    ///
    /// For each active dynamic body, the contributions from all regions whose
    /// extent includes the body are accumulated and applied to velocity.
    /// The world's own global gravity is *not* re-applied here.
    pub fn apply_gravity_regions(&mut self, regions: &[GravityRegion]) {
        let dt = self.dt;
        for body in self.bodies.iter_mut() {
            if !body.active || body.inv_mass == 0.0 {
                continue;
            }
            for region in regions {
                if let Some(g) = region.gravity_at(body.position) {
                    body.velocity[0] += g[0] * dt;
                    body.velocity[1] += g[1] * dt;
                    body.velocity[2] += g[2] * dt;
                }
            }
        }
    }
}
impl RigidWorld {
    /// Collect statistics about the current world state.
    pub fn stats(&self) -> RigidWorldStats {
        let total_bodies = self.bodies.len();
        let active_bodies = self
            .bodies
            .iter()
            .filter(|b| b.active && b.inv_mass > 0.0)
            .count();
        let inactive_bodies = self
            .bodies
            .iter()
            .filter(|b| !b.active && b.inv_mass > 0.0)
            .count();
        let static_bodies = self.bodies.iter().filter(|b| b.inv_mass == 0.0).count();
        RigidWorldStats {
            total_bodies,
            active_bodies,
            inactive_bodies,
            static_bodies,
            kinetic_energy: self.kinetic_energy(),
            sim_time: self.time,
        }
    }
}
/// Basic statistics for a `RigidWorld`.
#[derive(Debug, Clone)]
pub struct RigidWorldStats {
    /// Total number of bodies.
    pub total_bodies: usize,
    /// Number of active bodies.
    pub active_bodies: usize,
    /// Number of inactive (sleeping/deactivated) bodies.
    pub inactive_bodies: usize,
    /// Number of static bodies (inv_mass == 0).
    pub static_bodies: usize,
    /// Total kinetic energy.
    pub kinetic_energy: f64,
    /// Accumulated simulation time.
    pub sim_time: f64,
}
/// A serialisable snapshot of the rigid world state.
#[derive(Debug, Clone)]
pub struct RigidWorldSnapshot {
    /// Body snapshots.
    pub bodies: Vec<crate::pipeline::BodySnapshot>,
    /// Gravity.
    pub gravity: [f64; 3],
    /// Fixed time step.
    pub dt: f64,
    /// Accumulated simulation time.
    pub time: f64,
}
/// High-level configuration for a physics world.
#[derive(Debug, Clone)]
pub struct WorldConfig {
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
    /// Fixed time step in seconds.
    pub dt: f64,
    /// Linear sleep threshold (m/s).
    pub linear_sleep_threshold: f64,
    /// Angular sleep threshold (rad/s).
    pub angular_sleep_threshold: f64,
    /// Time before a body is put to sleep (seconds).
    pub time_before_sleep: f64,
    /// Solver configuration.
    pub solver: SolverConfig,
}
/// Result of a sphere-cast (moving sphere) query.
#[derive(Debug, Clone)]
pub struct SphereCastResult {
    /// Body ID of the first hit.
    pub body_id: u64,
    /// Centre of the query sphere at the moment of first contact.
    pub hit_centre: [f64; 3],
    /// Parametric distance `t` (in ray-space units).
    pub t: f64,
}
/// A detected contact between two bodies stored by their [`BodyHandle`]s.
#[derive(Debug, Clone)]
pub struct ContactPair {
    /// Handle of the first body.
    pub body_a: BodyHandle,
    /// Handle of the second body.
    pub body_b: BodyHandle,
    /// World-space contact normal (points from B toward A).
    pub normal: [f64; 3],
    /// Penetration depth in metres (positive = overlap).
    pub depth: f64,
    /// World-space contact point.
    pub contact_point: [f64; 3],
}
/// A contact cached at frame start for the small-steps solver. Stores the body
/// handles, the frame-start normal (from B toward A), and the two sphere radii
/// so penetration depth can be RE-PROJECTED from current body transforms each
/// substep without re-running broad/narrowphase.
#[derive(Debug, Clone)]
struct CachedContact {
    body_a: BodyHandle,
    body_b: BodyHandle,
    normal: [f64; 3],
    radius_a: f64,
    radius_b: f64,
}
/// Result of a raycast query against the world.
#[derive(Debug, Clone)]
pub struct RaycastResult {
    /// The body ID (for `RigidWorld`) or handle generation+index encoding.
    pub body_id: u64,
    /// The hit position in world space.
    pub hit_point: [f64; 3],
    /// The outward surface normal at the hit point.
    pub normal: [f64; 3],
    /// Parametric distance along the ray: `hit_point = origin + t * dir`.
    pub t: f64,
}
/// An axis-aligned bounding box query region.
#[derive(Debug, Clone, Copy)]
pub struct AabbQuery {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl AabbQuery {
    /// Create from two corners (order-independent).
    pub fn new(a: [f64; 3], b: [f64; 3]) -> Self {
        Self {
            min: [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])],
            max: [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])],
        }
    }
    /// Returns `true` if point `p` is inside (inclusive).
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
    /// Returns `true` if this AABB overlaps another AABB.
    pub fn overlaps(&self, other: &AabbQuery) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}
/// Configuration for the sequential impulse constraint solver.
#[derive(Debug, Clone)]
pub struct SolverConfig {
    /// Number of velocity solver iterations per step.
    pub velocity_iterations: usize,
    /// Number of position correction iterations per step.
    pub position_iterations: usize,
    /// Baumgarte stabilization factor (0..1).
    pub baumgarte_factor: f64,
    /// Minimum restitution threshold — below this value restitution is treated
    /// as zero to avoid jitter.
    pub restitution_threshold: f64,
    /// Number of velocity sub-steps per frame for the small-steps solver (default 1 = legacy single-step behaviour).
    pub substeps: usize,
    /// Coefficient of restitution used by the small-steps solver normal response.
    pub restitution: f64,
    /// Natural frequency (Hz) for soft contacts. `0.0` selects a rigid contact
    /// (the soft path then reduces to a Baumgarte-free accumulated-impulse
    /// projection). Only consulted when [`use_soft_contacts`](Self::use_soft_contacts)
    /// is `true`. Per Catto the effective frequency is internally capped to
    /// `0.25 / h` (a quarter of the sub-step rate) for stability.
    pub contact_hertz: f64,
    /// Damping ratio ζ for soft contacts. `1.0` = critically damped, `< 1` =
    /// under-damped (springy), `> 1` = over-damped. Only consulted when
    /// [`use_soft_contacts`](Self::use_soft_contacts) is `true`.
    pub contact_damping_ratio: f64,
    /// Enable soft (frequency / damping-ratio) contact response. When `false`
    /// (the default) the solver uses the existing rigid impulse path, which is
    /// byte-identical to prior behaviour.
    pub use_soft_contacts: bool,
}
/// A ray defined by an origin and a direction (need not be normalised; `t` is
/// in *ray-space* units i.e. multiples of the direction length).
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f64; 3],
    /// Ray direction (does not need to be unit length).
    pub direction: [f64; 3],
}
impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self { origin, direction }
    }
    /// Evaluate the ray at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}
/// Snapshot of physics-world runtime statistics.
#[derive(Debug, Clone)]
pub struct WorldStatistics {
    /// Total bodies (active + sleeping + static).
    pub total_bodies: usize,
    /// Number of active dynamic bodies.
    pub active_bodies: usize,
    /// Number of sleeping dynamic bodies.
    pub sleeping_bodies: usize,
    /// Number of static bodies.
    pub static_bodies: usize,
    /// Number of colliders.
    pub collider_count: usize,
    /// Number of contact pairs from the last step.
    pub contact_count: usize,
    /// Accumulated simulation time in seconds.
    pub sim_time: f64,
    /// Total kinetic energy of all dynamic bodies (J).
    pub kinetic_energy: f64,
}
/// Configuration for a single body added in batch.
#[derive(Debug, Clone)]
pub struct BodyConfig {
    /// Initial world-space position.
    pub position: [f64; 3],
    /// Initial linear velocity.
    pub velocity: [f64; 3],
    /// Inverse mass (0 = static).
    pub inv_mass: f64,
    /// Whether the body starts active.
    pub active: bool,
}
impl BodyConfig {
    /// Dynamic body at rest.
    pub fn dynamic_at(pos: [f64; 3]) -> Self {
        Self {
            position: pos,
            velocity: [0.0; 3],
            inv_mass: 1.0,
            active: true,
        }
    }
    /// Static body.
    pub fn static_at(pos: [f64; 3]) -> Self {
        Self {
            position: pos,
            velocity: [0.0; 3],
            inv_mass: 0.0,
            active: true,
        }
    }
}
/// Defines a localised gravity field that applies only to bodies inside a
/// spherical region.
#[derive(Debug, Clone)]
pub struct GravityRegion {
    /// Centre of the region in world space.
    pub centre: [f64; 3],
    /// Radius of influence (metres).
    pub radius: f64,
    /// Gravity acceleration applied inside the region.
    pub gravity: [f64; 3],
    /// Blend at the boundary: if `true`, gravity is linearly interpolated
    /// from zero at `radius` to full at `0` (falloff).
    pub falloff: bool,
}
impl GravityRegion {
    /// Create a new gravity region.
    pub fn new(centre: [f64; 3], radius: f64, gravity: [f64; 3]) -> Self {
        Self {
            centre,
            radius,
            gravity,
            falloff: false,
        }
    }
    /// Enable linear falloff.
    pub fn with_falloff(mut self) -> Self {
        self.falloff = true;
        self
    }
    /// Gravity acceleration for a body at `pos`.
    ///
    /// Returns `Some(g)` if `pos` is inside the region, `None` otherwise.
    pub fn gravity_at(&self, pos: [f64; 3]) -> Option<[f64; 3]> {
        let dx = pos[0] - self.centre[0];
        let dy = pos[1] - self.centre[1];
        let dz = pos[2] - self.centre[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist > self.radius {
            return None;
        }
        if self.falloff {
            let t = 1.0 - (dist / self.radius).min(1.0);
            Some([
                self.gravity[0] * t,
                self.gravity[1] * t,
                self.gravity[2] * t,
            ])
        } else {
            Some(self.gravity)
        }
    }
}
/// Full rigid body world: owns body storage, colliders, and runs the complete
/// simulation pipeline on each [`PhysicsWorld::step`] call.
///
/// Pipeline order per step:
/// 1. Gravity integration (velocity += g * dt)
/// 2. Force/torque accumulator integration
/// 3. Broadphase AABB overlap detection
/// 4. Narrowphase contact generation (sphere–sphere)
/// 5. Sequential-impulse constraint solver
/// 6. Position integration
/// 7. Sleep check / island deactivation
pub struct PhysicsWorld {
    /// Arena storage for all rigid bodies.
    pub bodies: RigidBodySet,
    /// Arena storage for all colliders.
    pub colliders: ColliderSet,
    /// Gravitational acceleration [x, y, z] in m/s².
    pub gravity: [f64; 3],
    /// Fixed simulation time-step in seconds.
    pub dt: f64,
    /// Constraint solver configuration.
    pub solver_config: SolverConfig,
    /// Accumulated simulation time in seconds.
    pub time: f64,
    /// Linear velocity threshold for sleep (m/s).
    pub linear_sleep_threshold: f64,
    /// Angular velocity threshold for sleep (rad/s).
    pub angular_sleep_threshold: f64,
    /// Time a body must be below thresholds before it sleeps (seconds).
    pub time_before_sleep: f64,
    /// Contact pairs from the last step.
    pub last_contacts: Vec<ContactPair>,
}
impl PhysicsWorld {
    /// Creates a new [`PhysicsWorld`] with the given gravity and fixed time-step.
    pub fn new(gravity: [f64; 3], dt: f64) -> Self {
        Self {
            bodies: RigidBodySet::new(),
            colliders: ColliderSet::new(),
            gravity,
            dt,
            solver_config: SolverConfig::default(),
            time: 0.0,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            last_contacts: Vec::new(),
        }
    }
    /// Creates a new [`PhysicsWorld`] with standard Earth gravity (0, -9.81, 0)
    /// and a 60 Hz time-step.
    pub fn default_earth() -> Self {
        Self::new([0.0, -9.81, 0.0], 1.0 / 60.0)
    }
    /// Sets the gravitational acceleration vector.
    pub fn set_gravity(&mut self, g: [f64; 3]) {
        self.gravity = g;
    }
    /// Inserts a [`RigidBody`] and returns a [`BodyHandle`] to it.
    pub fn add_rigid_body(&mut self, body: RigidBody) -> BodyHandle {
        self.bodies.insert(body)
    }
    /// Removes the body referenced by `handle`.
    ///
    /// Returns the removed body if the handle was valid.
    pub fn remove_rigid_body(&mut self, handle: BodyHandle) -> Option<RigidBody> {
        self.bodies.remove(handle)
    }
    /// Inserts a [`Collider`] with an optional parent body and returns a
    /// [`ColliderHandle`].
    pub fn add_collider(
        &mut self,
        mut collider: Collider,
        parent: Option<BodyHandle>,
    ) -> ColliderHandle {
        if let Some(h) = parent {
            collider = collider.with_body(h);
        }
        self.colliders.insert(collider)
    }
    /// Returns an immutable reference to the body referenced by `handle`.
    pub fn get_body(&self, handle: BodyHandle) -> Option<&RigidBody> {
        self.bodies.get(handle)
    }
    /// Returns a mutable reference to the body referenced by `handle`.
    pub fn get_body_mut(&mut self, handle: BodyHandle) -> Option<&mut RigidBody> {
        self.bodies.get_mut(handle)
    }
    /// Applies a persistent force (N) to the body's centre of mass.
    ///
    /// The force is added to the force accumulator and cleared after the next
    /// `step`.
    pub fn apply_force(&mut self, handle: BodyHandle, force: [f64; 3]) {
        if let Some(body) = self.bodies.get_mut(handle) {
            use oxiphysics_core::math::Vec3;
            body.apply_force(Vec3::new(force[0], force[1], force[2]));
        }
    }
    /// Applies an instantaneous linear impulse (N·s) to the body.
    pub fn apply_impulse(&mut self, handle: BodyHandle, impulse: [f64; 3]) {
        if let Some(body) = self.bodies.get_mut(handle) {
            body.apply_impulse(Vec3::new(impulse[0], impulse[1], impulse[2]));
        }
    }
    /// Total number of bodies currently in the world.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Number of bodies that are currently active (not sleeping).
    pub fn active_body_count(&self) -> usize {
        self.bodies
            .iter()
            .filter(|(_, b)| b.state == BodyState::Active)
            .count()
    }
    /// Number of bodies that are currently sleeping.
    pub fn sleeping_body_count(&self) -> usize {
        self.bodies
            .iter()
            .filter(|(_, b)| b.state == BodyState::Sleeping)
            .count()
    }
    /// Computes total kinetic energy (translational) of all dynamic bodies.
    pub fn kinetic_energy(&self) -> f64 {
        self.bodies
            .iter()
            .filter(|(_, b)| b.body_type == BodyType::Dynamic)
            .map(|(_, b)| {
                let v = &b.velocity;
                let v_sq = v.x * v.x + v.y * v.y + v.z * v.z;
                if b.inverse_mass > 0.0 {
                    let m = 1.0 / b.inverse_mass;
                    0.5 * m * v_sq
                } else {
                    0.0
                }
            })
            .sum()
    }
    /// Advances the simulation by one fixed time-step `dt`.
    ///
    /// Full pipeline:
    /// 1. Apply gravity + force accumulators → update velocities.
    /// 2. Broadphase: collect AABB overlap pairs.
    /// 3. Narrowphase: generate contacts for overlapping pairs.
    /// 4. Sequential-impulse solver: resolve velocity constraints.
    /// 5. Integrate positions.
    /// 6. Sleep check.
    pub fn step(&mut self) {
        let dt = self.dt;
        let gravity = self.gravity;
        let g_vec = Vec3::new(gravity[0], gravity[1], gravity[2]);
        for (_, body) in self.bodies.iter_mut() {
            body.integrate_forces(dt, &g_vec);
        }
        let candidate_pairs = self.broadphase_pairs();
        let contacts = self.narrowphase_contacts(&candidate_pairs);
        self.last_contacts = contacts.clone();
        let iters = self.solver_config.velocity_iterations;
        if self.solver_config.use_soft_contacts {
            // Soft (TGS-Soft / Catto) contact response. Build the soft
            // coefficients once for this step (`h == dt`, no sub-stepping in the
            // legacy `step` entry point), warm-start a per-contact accumulated
            // normal impulse across the velocity iterations, and route every
            // contact through the soft solver.
            let h = dt;
            let hz_eff = self.soft_contact_hertz(h);
            let zeta = self.solver_config.contact_damping_ratio;
            let soft = SoftParams::from_frequency(hz_eff, zeta, h);
            let mut lambda_acc = vec![0.0_f64; contacts.len()];
            for _ in 0..iters {
                for (i, cp) in contacts.iter().enumerate() {
                    self.solve_contact_velocity_soft(cp, &soft, &mut lambda_acc[i]);
                }
            }
        } else {
            for _ in 0..iters {
                for cp in &contacts {
                    self.solve_contact_velocity(cp);
                }
            }
        }
        for (_, body) in self.bodies.iter_mut() {
            body.integrate_velocity(dt);
        }
        let lin_thr = self.linear_sleep_threshold;
        let ang_thr = self.angular_sleep_threshold;
        let tbs = self.time_before_sleep;
        for (_, body) in self.bodies.iter_mut() {
            body.check_sleep(dt, lin_thr, ang_thr, tbs);
        }
        self.time += dt;
    }
    /// Returns candidate overlapping body pairs using AABB sphere-bound tests.
    ///
    /// Each body is represented by a sphere of radius `body.mass.cbrt() * 0.5`
    /// (heuristic) centred at `body.transform.position`.  Pairs whose bounding
    /// spheres overlap are returned.
    fn broadphase_pairs(&self) -> Vec<(BodyHandle, BodyHandle)> {
        let entries: Vec<(BodyHandle, [f64; 3], f64)> = self
            .bodies
            .iter()
            .map(|(h, b)| {
                let pos = [
                    b.transform.position.x,
                    b.transform.position.y,
                    b.transform.position.z,
                ];
                let r = if b.mass > 0.0 {
                    b.mass.cbrt() * 0.5
                } else {
                    0.5_f64
                };
                let r = r.max(0.1);
                (h, pos, r)
            })
            .collect();
        let mut pairs = Vec::new();
        let n = entries.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let (ha, pa, ra) = entries[i];
                let (hb, pb, rb) = entries[j];
                let dx = pa[0] - pb[0];
                let dy = pa[1] - pb[1];
                let dz = pa[2] - pb[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                let sum_r = ra + rb;
                if dist_sq <= sum_r * sum_r {
                    pairs.push((ha, hb));
                }
            }
        }
        pairs
    }
    /// Generates [`ContactPair`]s from the given broadphase candidate pairs.
    fn narrowphase_contacts(&self, pairs: &[(BodyHandle, BodyHandle)]) -> Vec<ContactPair> {
        let mut contacts = Vec::new();
        for &(ha, hb) in pairs {
            if let (Some(ba), Some(bb)) = (self.bodies.get(ha), self.bodies.get(hb)) {
                let pa = [
                    ba.transform.position.x,
                    ba.transform.position.y,
                    ba.transform.position.z,
                ];
                let pb = [
                    bb.transform.position.x,
                    bb.transform.position.y,
                    bb.transform.position.z,
                ];
                let ra = if ba.mass > 0.0 {
                    ba.mass.cbrt() * 0.5
                } else {
                    0.5_f64
                }
                .max(0.1);
                let rb = if bb.mass > 0.0 {
                    bb.mass.cbrt() * 0.5
                } else {
                    0.5_f64
                }
                .max(0.1);
                let dx = pa[0] - pb[0];
                let dy = pa[1] - pb[1];
                let dz = pa[2] - pb[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                let sum_r = ra + rb;
                if dist_sq < sum_r * sum_r {
                    let dist = dist_sq.sqrt();
                    let depth = sum_r - dist;
                    let (nx, ny, nz) = if dist > 1e-12 {
                        (dx / dist, dy / dist, dz / dist)
                    } else {
                        (0.0, 1.0, 0.0)
                    };
                    let contact_point = [pb[0] + nx * rb, pb[1] + ny * rb, pb[2] + nz * rb];
                    contacts.push(ContactPair {
                        body_a: ha,
                        body_b: hb,
                        normal: [nx, ny, nz],
                        depth,
                        contact_point,
                    });
                }
            }
        }
        contacts
    }
    /// Applies a single velocity-level impulse for one contact pair.
    fn solve_contact_velocity(&mut self, cp: &ContactPair) {
        let (va, inv_ma, ang_va) = match self.bodies.get(cp.body_a) {
            Some(b) if b.body_type == BodyType::Dynamic => {
                let v = [b.velocity.x, b.velocity.y, b.velocity.z];
                let av = [
                    b.angular_velocity.x,
                    b.angular_velocity.y,
                    b.angular_velocity.z,
                ];
                (v, b.inverse_mass, av)
            }
            _ => return,
        };
        let (vb, inv_mb, ang_vb) = match self.bodies.get(cp.body_b) {
            Some(b) => {
                let v = [b.velocity.x, b.velocity.y, b.velocity.z];
                let av = [
                    b.angular_velocity.x,
                    b.angular_velocity.y,
                    b.angular_velocity.z,
                ];
                let im = if b.body_type == BodyType::Dynamic {
                    b.inverse_mass
                } else {
                    0.0
                };
                (v, im, av)
            }
            _ => return,
        };
        let _ = ang_va;
        let _ = ang_vb;
        let n = cp.normal;
        let rel_vn = (va[0] - vb[0]) * n[0] + (va[1] - vb[1]) * n[1] + (va[2] - vb[2]) * n[2];
        if rel_vn <= 0.0 {
            return;
        }
        let denom = inv_ma + inv_mb;
        if denom < 1e-30 {
            return;
        }
        let restitution = 0.3_f64;
        let j = -(1.0 + restitution) * rel_vn / denom;
        if let Some(ba) = self.bodies.get_mut(cp.body_a)
            && ba.body_type == BodyType::Dynamic
        {
            ba.velocity.x += j * n[0] * inv_ma;
            ba.velocity.y += j * n[1] * inv_ma;
            ba.velocity.z += j * n[2] * inv_ma;
        }
        if let Some(bb) = self.bodies.get_mut(cp.body_b)
            && bb.body_type == BodyType::Dynamic
        {
            bb.velocity.x -= j * n[0] * inv_mb;
            bb.velocity.y -= j * n[1] * inv_mb;
            bb.velocity.z -= j * n[2] * inv_mb;
        }
    }
    /// Effective soft-contact frequency for sub-step size `h`, capped to Catto's
    /// stability limit of a quarter of the sub-step rate (`0.25 / h`). A
    /// `contact_hertz` of `0.0` is preserved (selects rigid soft params).
    fn soft_contact_hertz(&self, h: f64) -> f64 {
        let cap = if h > 1e-9 { 0.25 / h } else { f64::INFINITY };
        self.solver_config.contact_hertz.min(cap)
    }
    /// Core soft normal-contact solve shared by the per-step ([`step`](Self::step))
    /// and small-steps ([`step_small_steps`](Self::step_small_steps)) paths.
    ///
    /// Builds the full contact Jacobian `J = [n, r_a×n, -n, -(r_b×n)]` (linear +
    /// angular), computes the effective mass `J·M⁻¹·Jᵀ`, and performs one
    /// accumulated-impulse TGS-Soft iteration (Catto / Box2D v3):
    ///
    /// ```text
    /// jv        = J · [v_a; ω_a; v_b; ω_b]              (separation speed, >0 apart)
    /// target    = max(bias_rate · depth, restitution)  (never their sum)
    /// numerator = mass_scale · (jv − target) + impulse_scale · λ
    /// Δλ_raw    = −numerator / eff_mass
    /// λ         = max(λ + Δλ_raw, 0)                    (unilateral)
    /// ```
    ///
    /// The signed gap fed into the soft bias is `c = −depth` (negative when
    /// penetrating), so `bias_rate · depth` is the desired separation speed; the
    /// Baumgarte position push is folded into `target` with `max` (not a sum) per
    /// the P2 restitution lesson, so penetration recovery never inflates a
    /// bounce. The accumulated impulse `λ` (`lambda_acc`) is warm-started by the
    /// caller across iterations / sub-steps and clamped non-negative
    /// (unilateral). With [`SoftParams::rigid`] and `restitution == 0` this
    /// reduces to a plain non-penetration velocity projection.
    ///
    /// A body that is not [`BodyType::Dynamic`] contributes zero inverse mass /
    /// inertia and is never written, so the solve is symmetric in `body_a` /
    /// `body_b` (either may be static). Restitution and its threshold are passed
    /// explicitly so the small-steps path can supply a per-sub-step value.
    fn apply_soft_normal_constraint(
        &mut self,
        cp: &ContactPair,
        soft: &SoftParams,
        restitution: f64,
        restitution_threshold: f64,
        lambda_acc: &mut f64,
    ) {
        let (com_a, va, wa, inv_ma, inv_i_a, a_dynamic) = match self.bodies.get(cp.body_a) {
            Some(b) => {
                let dyn_a = b.body_type == BodyType::Dynamic;
                (
                    b.transform.position,
                    b.velocity,
                    b.angular_velocity,
                    if dyn_a { b.inverse_mass } else { 0.0 },
                    if dyn_a {
                        b.world_inverse_inertia
                    } else {
                        Mat3::zeros()
                    },
                    dyn_a,
                )
            }
            None => return,
        };
        let (com_b, vb, wb, inv_mb, inv_i_b, b_dynamic) = match self.bodies.get(cp.body_b) {
            Some(b) => {
                let dyn_b = b.body_type == BodyType::Dynamic;
                (
                    b.transform.position,
                    b.velocity,
                    b.angular_velocity,
                    if dyn_b { b.inverse_mass } else { 0.0 },
                    if dyn_b {
                        b.world_inverse_inertia
                    } else {
                        Mat3::zeros()
                    },
                    dyn_b,
                )
            }
            None => return,
        };
        if !a_dynamic && !b_dynamic {
            return;
        }
        let n = Vec3::new(cp.normal[0], cp.normal[1], cp.normal[2]);
        let point = Vec3::new(
            cp.contact_point[0],
            cp.contact_point[1],
            cp.contact_point[2],
        );
        let r_a = point - com_a;
        let r_b = point - com_b;
        let ra_x_n = r_a.cross(&n);
        let rb_x_n = r_b.cross(&n);
        // Effective mass J·M⁻¹·Jᵀ. The angular terms vanish for sphere contacts
        // (contact point on the line of centres ⇒ r×n = 0) but are retained for
        // correctness with off-centre contacts.
        let ang_a = ra_x_n.dot(&(inv_i_a * ra_x_n));
        let ang_b = rb_x_n.dot(&(inv_i_b * rb_x_n));
        let eff_mass = inv_ma + inv_mb + ang_a + ang_b;
        if eff_mass < 1e-30 {
            return;
        }
        // Relative normal velocity at the contact point (separating > 0,
        // approaching < 0); normal points from B toward A.
        let jv = n.dot(&(va - vb)) + ra_x_n.dot(&wa) - rb_x_n.dot(&wb);
        let depth = cp.depth.max(0.0);
        let soft_bias = soft.bias_rate * depth;
        let restitution_term = if jv < -restitution_threshold {
            -restitution * jv
        } else {
            0.0
        };
        let target = soft_bias.max(restitution_term);
        let numerator = soft.mass_scale * (jv - target) + soft.impulse_scale * *lambda_acc;
        let impulse_raw = -numerator / eff_mass;
        let lambda_new = (*lambda_acc + impulse_raw).max(0.0);
        let d_lambda = lambda_new - *lambda_acc;
        *lambda_acc = lambda_new;
        if a_dynamic && let Some(ba) = self.bodies.get_mut(cp.body_a) {
            ba.velocity += n * (d_lambda * inv_ma);
            ba.angular_velocity += inv_i_a * (ra_x_n * d_lambda);
        }
        if b_dynamic && let Some(bb) = self.bodies.get_mut(cp.body_b) {
            bb.velocity -= n * (d_lambda * inv_mb);
            bb.angular_velocity -= inv_i_b * (rb_x_n * d_lambda);
        }
    }
    /// Soft (TGS-Soft) velocity solve for one contact pair, used by the per-step
    /// [`step`](Self::step) entry point. Restitution and its threshold come from
    /// [`SolverConfig`] (the threshold is clamped to `0.5` to match the
    /// small-steps rule); `lambda_acc` is the caller-owned accumulated normal
    /// impulse for this contact, warm-started across velocity iterations.
    fn solve_contact_velocity_soft(
        &mut self,
        cp: &ContactPair,
        soft: &SoftParams,
        lambda_acc: &mut f64,
    ) {
        let restitution = self.solver_config.restitution;
        let restitution_threshold = self.solver_config.restitution_threshold.min(0.5);
        self.apply_soft_normal_constraint(cp, soft, restitution, restitution_threshold, lambda_acc);
    }
    /// Soft (TGS-Soft) velocity solve for one cached small-steps contact.
    ///
    /// Re-projects penetration depth from the CURRENT body transforms (the
    /// normal is held at its cached frame-start value, exactly like
    /// [`solve_cached_contact_velocity`](Self::solve_cached_contact_velocity));
    /// separated contacts (negative re-projected depth) are skipped. The
    /// sphere-sphere contact point lies on the line of centres, so the angular
    /// Jacobian terms vanish and the solve reduces to the linear soft normal
    /// response. `e_sub` is the per-sub-step restitution; `lambda_acc` is
    /// warm-started across the frame.
    fn solve_cached_contact_velocity_soft(
        &mut self,
        c: &CachedContact,
        soft: &SoftParams,
        e_sub: f64,
        lambda_acc: &mut f64,
    ) {
        let pa = match self.bodies.get(c.body_a) {
            Some(b) => b.transform.position,
            None => return,
        };
        let pb = match self.bodies.get(c.body_b) {
            Some(b) => b.transform.position,
            None => return,
        };
        let delta = pa - pb;
        let dist = delta.norm();
        let depth = (c.radius_a + c.radius_b) - dist;
        if depth < 0.0 {
            return;
        }
        // Sphere-sphere contact point: on A's surface toward B. The normal
        // points from B toward A, so step back along -n by radius_a from A.
        let contact_point = [
            pa.x - c.normal[0] * c.radius_a,
            pa.y - c.normal[1] * c.radius_a,
            pa.z - c.normal[2] * c.radius_a,
        ];
        let cp = ContactPair {
            body_a: c.body_a,
            body_b: c.body_b,
            normal: c.normal,
            depth,
            contact_point,
        };
        let restitution_threshold = self.solver_config.restitution_threshold.min(0.5);
        self.apply_soft_normal_constraint(&cp, soft, e_sub, restitution_threshold, lambda_acc);
    }
}
impl PhysicsWorld {
    /// Apply a [`WorldConfig`] to this world.
    pub fn apply_config(&mut self, cfg: &WorldConfig) {
        self.gravity = cfg.gravity;
        self.dt = cfg.dt;
        self.linear_sleep_threshold = cfg.linear_sleep_threshold;
        self.angular_sleep_threshold = cfg.angular_sleep_threshold;
        self.time_before_sleep = cfg.time_before_sleep;
        self.solver_config = cfg.solver.clone();
    }
    /// Snapshot the current configuration.
    pub fn get_config(&self) -> WorldConfig {
        WorldConfig {
            gravity: self.gravity,
            dt: self.dt,
            linear_sleep_threshold: self.linear_sleep_threshold,
            angular_sleep_threshold: self.angular_sleep_threshold,
            time_before_sleep: self.time_before_sleep,
            solver: self.solver_config.clone(),
        }
    }
}
impl PhysicsWorld {
    /// Collect world statistics.
    pub fn statistics(&self) -> WorldStatistics {
        use crate::body::BodyState;
        use crate::body::BodyType;
        let total_bodies = self.bodies.len();
        let active_bodies = self
            .bodies
            .iter()
            .filter(|(_, b)| b.state == BodyState::Active && b.body_type == BodyType::Dynamic)
            .count();
        let sleeping_bodies = self
            .bodies
            .iter()
            .filter(|(_, b)| b.state == BodyState::Sleeping)
            .count();
        let static_bodies = self
            .bodies
            .iter()
            .filter(|(_, b)| b.body_type == BodyType::Static)
            .count();
        WorldStatistics {
            total_bodies,
            active_bodies,
            sleeping_bodies,
            static_bodies,
            collider_count: self.colliders.len(),
            contact_count: self.last_contacts.len(),
            sim_time: self.time,
            kinetic_energy: self.kinetic_energy(),
        }
    }
}
impl PhysicsWorld {
    /// Return handles of all bodies whose centre is within `radius` of `point`.
    pub fn bodies_within_radius(
        &self,
        point: [f64; 3],
        radius: f64,
    ) -> Vec<oxiphysics_core::BodyHandle> {
        self.bodies
            .iter()
            .filter_map(|(h, b)| {
                let p = [
                    b.transform.position.x,
                    b.transform.position.y,
                    b.transform.position.z,
                ];
                let dx = p[0] - point[0];
                let dy = p[1] - point[1];
                let dz = p[2] - point[2];
                if dx * dx + dy * dy + dz * dz <= radius * radius {
                    Some(h)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Return handles of all bodies above the given Y-coordinate.
    pub fn bodies_above_y(&self, y: f64) -> Vec<oxiphysics_core::BodyHandle> {
        self.bodies
            .iter()
            .filter_map(|(h, b)| {
                if b.transform.position.y > y {
                    Some(h)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Return handles of all bodies with kinetic energy exceeding `threshold`.
    pub fn bodies_exceeding_ke(&self, threshold: f64) -> Vec<oxiphysics_core::BodyHandle> {
        self.bodies
            .iter()
            .filter_map(|(h, b)| {
                if b.body_type != BodyType::Dynamic || b.inverse_mass <= 0.0 {
                    return None;
                }
                let v = &b.velocity;
                let v_sq = v.x * v.x + v.y * v.y + v.z * v.z;
                let ke = 0.5 / b.inverse_mass * v_sq;
                if ke > threshold { Some(h) } else { None }
            })
            .collect()
    }
}
impl PhysicsWorld {
    /// Apply the same force vector to a list of bodies.
    pub fn apply_force_batch(&mut self, handles: &[oxiphysics_core::BodyHandle], force: [f64; 3]) {
        let f = Vec3::new(force[0], force[1], force[2]);
        for &h in handles {
            if let Some(b) = self.bodies.get_mut(h) {
                b.apply_force(f);
            }
        }
    }
    /// Apply an impulse to a list of bodies.
    pub fn apply_impulse_batch(
        &mut self,
        handles: &[oxiphysics_core::BodyHandle],
        impulse: [f64; 3],
    ) {
        let imp = Vec3::new(impulse[0], impulse[1], impulse[2]);
        for &h in handles {
            if let Some(b) = self.bodies.get_mut(h) {
                b.apply_impulse(imp);
            }
        }
    }
    /// Wake up all sleeping bodies (set state back to Active).
    pub fn wake_all(&mut self) {
        for (_, b) in self.bodies.iter_mut() {
            if b.state == BodyState::Sleeping {
                b.state = BodyState::Active;
                b.sleep_timer = 0.0;
            }
        }
    }
    /// Set gravity scale on all bodies matching `body_type`.
    pub fn set_gravity_scale_all(&mut self, body_type: crate::body::BodyType, scale: f64) {
        for (_, b) in self.bodies.iter_mut() {
            if b.body_type == body_type {
                b.gravity_scale = scale;
            }
        }
    }
}
impl PhysicsWorld {
    /// Advance the simulation by `total_dt` seconds using fixed sub-steps of
    /// size `self.dt`.
    ///
    /// The number of sub-steps is `ceil(total_dt / dt)`.  The last sub-step
    /// may use a smaller `dt` if `total_dt` is not an exact multiple.
    ///
    /// This is useful for decoupling the render frame rate from the physics
    /// integration frequency.
    pub fn step_substep(&mut self, total_dt: f64) {
        if total_dt <= 0.0 || self.dt <= 0.0 {
            return;
        }
        let mut remaining = total_dt;
        while remaining > 1e-15 {
            let sub_dt = remaining.min(self.dt);
            let saved_dt = self.dt;
            self.dt = sub_dt;
            self.step();
            self.dt = saved_dt;
            remaining -= sub_dt;
        }
    }
    /// Compute the energy conservation error relative to a reference energy.
    ///
    /// Returns `|E_current - E_reference| / (|E_reference| + 1.0)` so that
    /// the result is always finite even when `E_reference` is zero.
    ///
    /// A small value (≪ 1) indicates good energy conservation; a large value
    /// suggests numerical drift or constraint errors.
    pub fn compute_energy_error(&self, reference_energy: f64) -> f64 {
        let ke = self.kinetic_energy();
        let diff = (ke - reference_energy).abs();
        diff / (reference_energy.abs() + 1.0)
    }
    /// Apply a uniform wind force to all active dynamic bodies.
    ///
    /// The wind force on each body is modelled as aerodynamic drag:
    ///
    /// ```text
    /// F_wind = 0.5 * rho_air * Cd * A * |v_rel|² * v_rel_hat
    /// ```
    ///
    /// where `v_rel = wind_velocity - body_velocity`.
    ///
    /// For simplicity the product `Cd * A` (drag coefficient × reference area)
    /// is assumed to be `mass^(2/3) * 0.1` (a heuristic for a unit-density
    /// sphere).  This gives a physically plausible force without requiring per-
    /// body aerodynamic properties.
    ///
    /// * `wind_velocity` — world-space wind velocity vector (m/s).
    /// * `air_density`   — air density (kg/m³), e.g. 1.2 at sea level.
    pub fn apply_wind_force(&mut self, wind_velocity: [f64; 3], air_density: f64) {
        let rho = air_density.max(0.0);
        for (_, body) in self.bodies.iter_mut() {
            if body.inverse_mass <= 0.0 || body.body_type != BodyType::Dynamic {
                continue;
            }
            let m = if body.inverse_mass > 1e-30 {
                1.0 / body.inverse_mass
            } else {
                0.0
            };
            if m < 1e-30 {
                continue;
            }
            let vx = wind_velocity[0] - body.velocity.x;
            let vy = wind_velocity[1] - body.velocity.y;
            let vz = wind_velocity[2] - body.velocity.z;
            let v_sq = vx * vx + vy * vy + vz * vz;
            if v_sq < 1e-30 {
                continue;
            }
            let v_mag = v_sq.sqrt();
            let cd_a = m.powf(2.0 / 3.0) * 0.1;
            let force_mag = 0.5 * rho * cd_a * v_sq;
            let f = Vec3::new(
                force_mag * vx / v_mag,
                force_mag * vy / v_mag,
                force_mag * vz / v_mag,
            );
            body.apply_force(f);
        }
    }
}
impl PhysicsWorld {
    /// Advance the simulation by one frame using the PhysX-5 / Macklin et al.
    /// 2019 "Small Steps" sub-stepping scheme.
    ///
    /// # The small-steps contract
    ///
    /// Unlike [`step`](Self::step) (which runs the whole pipeline once per
    /// frame) and unlike the naive [`step_substep`](Self::step_substep) (which
    /// re-runs the *entire* broad/narrowphase + solve pipeline for every
    /// sub-`dt`), this method performs collision **detection exactly ONCE per
    /// frame on the frame-start state**. The detected contacts are cached
    /// together with their frame-start normal (which points from B toward A)
    /// and the two sphere radii. The frame is then advanced over `N =
    /// self.solver_config.substeps` sub-steps of size `h = dt / N`, where
    /// `dt = self.dt`.
    ///
    /// Each sub-step:
    /// 1. integrates gravity into velocity over `h` (`integrate_forces`);
    /// 2. runs **one** velocity solve over the cached contacts. The penetration
    ///    depth is **RE-PROJECTED from the CURRENT body transforms along the
    ///    cached normal** — broad/narrowphase is **NOT** re-run, and the normal
    ///    is held fixed at its frame-start value (this is the defining
    ///    "re-project, do not re-detect" contract of small steps). Any contact
    ///    whose re-projected depth has become negative (the bodies have
    ///    separated) is **SKIPPED** for that sub-step;
    /// 3. integrates position over `h` (`integrate_velocity`).
    ///
    /// After the sub-step loop a short final **relax pass** (two iterations
    /// with Baumgarte bias and restitution disabled) removes residual approach
    /// velocity for stack stability. Simulation time is advanced by the full
    /// `dt` once at the end (NOT per sub-step).
    ///
    /// This scheme supersedes [`step_substep`](Self::step_substep) on the hot
    /// path because detecting once and re-projecting is dramatically cheaper
    /// and far more stable for stacks than re-running the full pipeline per
    /// sub-`dt`. The legacy [`step`](Self::step) and
    /// [`step_substep`](Self::step_substep) entry points are left unchanged;
    /// `self.dt` controls the frame step and `self.solver_config.substeps`
    /// controls the sub-step count (default `1`, i.e. a single small step).
    ///
    /// # Restitution threshold clamping
    ///
    /// The [`SubStepRestitutionCorrector`] is constructed with its velocity
    /// threshold clamped to `restitution_threshold.min(0.5)`. The default
    /// `restitution_threshold` is `1.0` m/s, which would erroneously suppress a
    /// genuine `0.5` m/s bounce; clamping to `0.5` keeps resting-jitter
    /// suppression while still allowing real bounces through. (The per-sub-step
    /// restitution `e_sub` does not currently consume the velocity threshold —
    /// the corrector is constructed exactly as specified for documentation and
    /// forward-compatibility, and `e_sub` is its
    /// [`per_substep_restitution`](SubStepRestitutionCorrector::per_substep_restitution)
    /// output.)
    pub fn step_small_steps(&mut self) {
        let dt = self.dt;
        let substeps = self.solver_config.substeps.max(1);
        let h = dt / substeps as f64;
        let gravity = self.gravity;
        let g_vec = Vec3::new(gravity[0], gravity[1], gravity[2]);

        // Detect ONCE on frame-start state.
        let pairs = self.broadphase_pairs();
        let contacts = self.narrowphase_contacts(&pairs);
        self.last_contacts = contacts.clone();

        let mut cached: Vec<CachedContact> = Vec::with_capacity(contacts.len());
        for cp in &contacts {
            let (ra, rb) = {
                let ba = match self.bodies.get(cp.body_a) {
                    Some(b) => b,
                    None => continue,
                };
                let bb = match self.bodies.get(cp.body_b) {
                    Some(b) => b,
                    None => continue,
                };
                let ra = if ba.mass > 0.0 {
                    ba.mass.cbrt() * 0.5
                } else {
                    0.5_f64
                };
                let ra = ra.max(0.1);
                let rb = if bb.mass > 0.0 {
                    bb.mass.cbrt() * 0.5
                } else {
                    0.5_f64
                };
                let rb = rb.max(0.1);
                (ra, rb)
            };
            cached.push(CachedContact {
                body_a: cp.body_a,
                body_b: cp.body_b,
                normal: cp.normal,
                radius_a: ra,
                radius_b: rb,
            });
        }

        let corrector = SubStepRestitutionCorrector::new(
            self.solver_config.restitution,
            self.solver_config.restitution_threshold.min(0.5),
            substeps,
        );
        let e_sub = corrector.per_substep_restitution();

        // Move (not clone) into the vec iterated by the solve loops so we never
        // hold a borrow of `self.bodies` across the position/velocity writes.
        let cached_clone = cached;

        if self.solver_config.use_soft_contacts {
            // Soft (TGS-Soft / Catto) small-steps path. The soft coefficients
            // are built once for the sub-step size `h = dt / substeps`; a
            // per-contact accumulated normal impulse is warm-started across the
            // whole frame (sub-steps + relax pass). Restitution stays
            // per-sub-step (`e_sub`) so the total restitution is unchanged.
            let hz_eff = self.soft_contact_hertz(h);
            let zeta = self.solver_config.contact_damping_ratio;
            let soft = SoftParams::from_frequency(hz_eff, zeta, h);
            let mut lambda_acc = vec![0.0_f64; cached_clone.len()];
            for _ in 0..substeps {
                for (_, body) in self.bodies.iter_mut() {
                    body.integrate_forces(h, &g_vec);
                }
                for (i, c) in cached_clone.iter().enumerate() {
                    self.solve_cached_contact_velocity_soft(c, &soft, e_sub, &mut lambda_acc[i]);
                }
                for (_, body) in self.bodies.iter_mut() {
                    body.integrate_velocity(h);
                }
            }
            // Final relax pass: pure non-penetration projection (rigid soft
            // params disable both bias and restitution), warm-started.
            let relax = SoftParams::rigid();
            for _ in 0..2 {
                for (i, c) in cached_clone.iter().enumerate() {
                    self.solve_cached_contact_velocity_soft(c, &relax, 0.0, &mut lambda_acc[i]);
                }
            }
        } else {
            for _ in 0..substeps {
                for (_, body) in self.bodies.iter_mut() {
                    body.integrate_forces(h, &g_vec);
                }
                for c in &cached_clone {
                    self.solve_cached_contact_velocity(c, e_sub, h, true);
                }
                for (_, body) in self.bodies.iter_mut() {
                    body.integrate_velocity(h);
                }
            }

            // Final relax pass: settle residual approach velocity for stack
            // stability (Baumgarte bias and restitution disabled).
            for _ in 0..2 {
                for c in &cached_clone {
                    self.solve_cached_contact_velocity(c, 0.0, h, false);
                }
            }
        }

        self.time += dt;
    }
    /// One velocity-level solve for a single cached small-steps contact.
    ///
    /// The penetration depth is **re-projected from the CURRENT body
    /// transforms** (`depth = (radius_a + radius_b) - |pa - pb|`) while the
    /// contact normal is held **fixed** at its cached frame-start value; if the
    /// re-projected depth is negative the bodies have separated and the contact
    /// is skipped.
    ///
    /// Sign convention: the normal points from B toward A. A resting /
    /// penetrating contact has `rel_vn` slightly negative and (in bias mode)
    /// `bias > 0`, so the solved impulse `j > 0` pushes A along `+n` and B along
    /// `-n` (apart). A fast impact has `rel_vn` strongly negative and
    /// `restitution_term = -e_sub * rel_vn > 0`, producing a large `j` and hence
    /// a bounce. `j` is clamped `>= 0` so bodies are never pulled together. The
    /// solve target is `max(bias, restitution_term)` — the *larger* of the
    /// Baumgarte position-correction bias and the restitution bias, never their
    /// sum — so the penetration push cannot inflate the bounce. Energy stays
    /// bounded because `j` is bounded by `|rel_vn|` plus the larger of
    /// `e_sub * |rel_vn|` (with `e_sub < 1`) and the position-correction bias.
    ///
    /// When `with_bias` is `false` (the relax/settle pass) both the Baumgarte
    /// position-correction bias and the restitution term are disabled, leaving a
    /// pure non-penetration velocity projection.
    fn solve_cached_contact_velocity(
        &mut self,
        c: &CachedContact,
        e_sub: f64,
        h: f64,
        with_bias: bool,
    ) {
        let pa = match self.bodies.get(c.body_a) {
            Some(b) => [
                b.transform.position.x,
                b.transform.position.y,
                b.transform.position.z,
            ],
            None => return,
        };
        let pb = match self.bodies.get(c.body_b) {
            Some(b) => [
                b.transform.position.x,
                b.transform.position.y,
                b.transform.position.z,
            ],
            None => return,
        };
        let dx = pa[0] - pb[0];
        let dy = pa[1] - pb[1];
        let dz = pa[2] - pb[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let depth = (c.radius_a + c.radius_b) - dist;
        if depth < 0.0 {
            return;
        }
        let (va, inv_ma) = match self.bodies.get(c.body_a) {
            Some(b) => {
                let im = if b.body_type == BodyType::Dynamic {
                    b.inverse_mass
                } else {
                    0.0
                };
                ([b.velocity.x, b.velocity.y, b.velocity.z], im)
            }
            None => return,
        };
        let (vb, inv_mb) = match self.bodies.get(c.body_b) {
            Some(b) => {
                let im = if b.body_type == BodyType::Dynamic {
                    b.inverse_mass
                } else {
                    0.0
                };
                ([b.velocity.x, b.velocity.y, b.velocity.z], im)
            }
            None => return,
        };
        let n = c.normal;
        let rel_vn = (va[0] - vb[0]) * n[0] + (va[1] - vb[1]) * n[1] + (va[2] - vb[2]) * n[2];
        let denom = inv_ma + inv_mb;
        if denom < 1e-30 {
            return;
        }
        let bias = if with_bias {
            let slop = 0.005;
            let beta = self.solver_config.baumgarte_factor;
            (beta / h) * (depth - slop).max(0.0)
        } else {
            0.0
        };
        let restitution_term = if rel_vn < 0.0 { -e_sub * rel_vn } else { 0.0 };
        // Post-solve target separation speed: the LARGER of the Baumgarte
        // position-correction bias and the restitution bias, NOT their sum.
        // Summing double-counts (the penetration push would inflate the bounce
        // and inject energy); taking the maximum is the standard Catto/Box2D
        // rule and keeps a real bounce while bounding energy.
        let target = bias.max(restitution_term);
        let mut j = (target - rel_vn) / denom;
        if j < 0.0 {
            j = 0.0;
        }
        if let Some(ba) = self.bodies.get_mut(c.body_a)
            && ba.body_type == BodyType::Dynamic
        {
            ba.velocity.x += j * n[0] * inv_ma;
            ba.velocity.y += j * n[1] * inv_ma;
            ba.velocity.z += j * n[2] * inv_ma;
        }
        if let Some(bb) = self.bodies.get_mut(c.body_b)
            && bb.body_type == BodyType::Dynamic
        {
            bb.velocity.x -= j * n[0] * inv_mb;
            bb.velocity.y -= j * n[1] * inv_mb;
            bb.velocity.z -= j * n[2] * inv_mb;
        }
    }
}
