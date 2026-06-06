// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::math::Real;
use crate::math::{Quat, Vec3};
use crate::types::{Aabb, BodyHandle, PhysicsConfig, TimeStep, Transform};
use std::collections::VecDeque;

/// A force field applied to bodies inside a given AABB.
#[derive(Debug, Clone)]
pub struct ForceField {
    /// The kind of force applied.
    pub kind: ForceFieldKind,
    /// Region of influence. `None` = global.
    pub region: Option<crate::types::Aabb>,
    /// Whether this force field is active.
    pub active: bool,
}
impl ForceField {
    /// Create a constant (wind-like) force field.
    pub fn constant(force: [Real; 3]) -> Self {
        Self {
            kind: ForceFieldKind::Constant { force },
            region: None,
            active: true,
        }
    }
    /// Create a point-gravity force field.
    pub fn point_gravity(position: [Real; 3], gm: Real, epsilon: Real) -> Self {
        Self {
            kind: ForceFieldKind::PointGravity {
                position,
                gm,
                epsilon,
            },
            region: None,
            active: true,
        }
    }
    /// Create a quadratic drag field.
    pub fn drag(cd: Real) -> Self {
        Self {
            kind: ForceFieldKind::QuadraticDrag { cd },
            region: None,
            active: true,
        }
    }
    /// Compute the force on a body with given position, velocity, and mass.
    pub fn force_on(&self, pos: [Real; 3], vel: [Real; 3], mass: Real) -> [Real; 3] {
        if !self.active {
            return [0.0; 3];
        }
        if let Some(region) = &self.region {
            let p = crate::math::Vec3::new(pos[0], pos[1], pos[2]);
            if !region.contains_point(&p) {
                return [0.0; 3];
            }
        }
        match &self.kind {
            ForceFieldKind::Constant { force } => *force,
            ForceFieldKind::PointGravity {
                position,
                gm,
                epsilon,
            } => {
                let dx = position[0] - pos[0];
                let dy = position[1] - pos[1];
                let dz = position[2] - pos[2];
                let r2 = dx * dx + dy * dy + dz * dz + epsilon * epsilon;
                let r = r2.sqrt();
                let f = mass * gm / (r2 * r);
                [f * dx, f * dy, f * dz]
            }
            ForceFieldKind::QuadraticDrag { cd } => {
                let v2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
                let v = v2.sqrt();
                if v < 1e-12 {
                    return [0.0; 3];
                }
                let f = -cd * v;
                [f * vel[0], f * vel[1], f * vel[2]]
            }
            ForceFieldKind::Vortex {
                center,
                axis,
                omega,
                decay,
            } => {
                let rx = pos[0] - center[0];
                let ry = pos[1] - center[1];
                let rz = pos[2] - center[2];
                let rr = rx * axis[0] + ry * axis[1] + rz * axis[2];
                let rx_rad = rx - rr * axis[0];
                let ry_rad = ry - rr * axis[1];
                let rz_rad = rz - rr * axis[2];
                let rad = (rx_rad * rx_rad + ry_rad * ry_rad + rz_rad * rz_rad).sqrt();
                if rad < 1e-12 {
                    return [0.0; 3];
                }
                let scale = omega * (-rad / decay).exp();
                let tx = axis[1] * rz_rad - axis[2] * ry_rad;
                let ty = axis[2] * rx_rad - axis[0] * rz_rad;
                let tz = axis[0] * ry_rad - axis[1] * rx_rad;
                let t_mag = (tx * tx + ty * ty + tz * tz).sqrt().max(1e-20);
                [
                    mass * scale * tx / t_mag,
                    mass * scale * ty / t_mag,
                    mass * scale * tz / t_mag,
                ]
            }
        }
    }
}
/// All simulation state for a single rigid body.
#[derive(Debug, Clone)]
pub struct Body {
    /// Current world-space transform.
    pub transform: Transform,
    /// Previous-frame transform (for CCD and interpolation).
    pub prev_transform: Transform,
    /// Current velocities.
    pub velocity: BodyVelocity,
    /// Accumulated force for this step (Newtons).
    pub force: [Real; 3],
    /// Accumulated torque for this step (N·m).
    pub torque: [Real; 3],
    /// Body mass (kg). 0 means static.
    pub mass: Real,
    /// Inverse mass. 0 means static/kinematic.
    pub inv_mass: Real,
    /// Isotropic moment of inertia scalar.
    pub inertia: Real,
    /// Linear damping coefficient.
    pub linear_damping: Real,
    /// Angular damping coefficient.
    pub angular_damping: Real,
    /// Restitution coefficient (bounciness).
    pub restitution: Real,
    /// Simulation mode.
    pub mode: BodyMode,
    /// Whether the body is currently active (participating in simulation).
    pub active: bool,
    /// How many consecutive steps the body has been below sleep thresholds.
    pub sleep_timer: Real,
    /// Generation counter (for handle validation).
    pub generation: u32,
}
impl Body {
    /// Create a dynamic body at the given position with given mass.
    pub fn new_dynamic(position: [Real; 3], mass: Real) -> Self {
        use crate::math::{Quat, Vec3};
        let pos = Vec3::new(position[0], position[1], position[2]);
        Body {
            transform: Transform::new(pos, Quat::identity()),
            prev_transform: Transform::new(pos, Quat::identity()),
            mass,
            inv_mass: if mass > 0.0 { 1.0 / mass } else { 0.0 },
            inertia: mass * 0.1,
            ..Default::default()
        }
    }
    /// Create a static body at the given position.
    pub fn new_static(position: [Real; 3]) -> Self {
        let pos = Vec3::new(position[0], position[1], position[2]);
        Body {
            transform: Transform::new(pos, Quat::identity()),
            prev_transform: Transform::new(pos, Quat::identity()),
            mass: 0.0,
            inv_mass: 0.0,
            inertia: 0.0,
            mode: BodyMode::Static,
            ..Default::default()
        }
    }
    /// Apply a force to the body (accumulates until cleared).
    pub fn apply_force(&mut self, force: [Real; 3]) {
        if self.mode == BodyMode::Dynamic {
            self.force[0] += force[0];
            self.force[1] += force[1];
            self.force[2] += force[2];
        }
    }
    /// Apply a torque to the body.
    pub fn apply_torque(&mut self, torque: [Real; 3]) {
        if self.mode == BodyMode::Dynamic {
            self.torque[0] += torque[0];
            self.torque[1] += torque[1];
            self.torque[2] += torque[2];
        }
    }
    /// Apply an impulse at a world-space offset (cross-product contributes to angular).
    pub fn apply_impulse_at(&mut self, impulse: [Real; 3], offset: [Real; 3]) {
        if self.mode != BodyMode::Dynamic {
            return;
        }
        self.velocity.apply_linear_impulse(impulse, self.inv_mass);
        let tx = offset[1] * impulse[2] - offset[2] * impulse[1];
        let ty = offset[2] * impulse[0] - offset[0] * impulse[2];
        let tz = offset[0] * impulse[1] - offset[1] * impulse[0];
        let inv_i = if self.inertia > 0.0 {
            1.0 / self.inertia
        } else {
            0.0
        };
        self.velocity.apply_angular_impulse([tx, ty, tz], inv_i);
    }
    /// Clear accumulated forces and torques.
    pub fn clear_forces(&mut self) {
        self.force = [0.0; 3];
        self.torque = [0.0; 3];
    }
    /// Integrate velocities from forces over dt (semi-implicit Euler).
    pub fn integrate_velocities(&mut self, gravity: [Real; 3], dt: Real) {
        if self.mode != BodyMode::Dynamic {
            return;
        }
        let inv_m = self.inv_mass;
        for (k, vl) in self.velocity.linear.iter_mut().enumerate() {
            *vl += (self.force[k] * inv_m + gravity[k]) * dt;
        }
        let inv_i = if self.inertia > 0.0 {
            1.0 / self.inertia
        } else {
            0.0
        };
        for k in 0..3 {
            self.velocity.angular[k] += self.torque[k] * inv_i * dt;
        }
        self.velocity
            .damp(self.linear_damping, self.angular_damping, dt);
    }
    /// Integrate positions from velocities over dt.
    pub fn integrate_positions(&mut self, dt: Real) {
        if self.mode == BodyMode::Static {
            return;
        }
        self.prev_transform = self.transform.clone();
        self.transform.position.x += self.velocity.linear[0] * dt;
        self.transform.position.y += self.velocity.linear[1] * dt;
        self.transform.position.z += self.velocity.linear[2] * dt;
        let wx = self.velocity.angular[0] * dt * 0.5;
        let wy = self.velocity.angular[1] * dt * 0.5;
        let wz = self.velocity.angular[2] * dt * 0.5;
        let q = self.transform.rotation;
        use nalgebra::Quaternion;
        let dq = Quaternion::new(
            -q.i * wx - q.j * wy - q.k * wz,
            q.w * wx + q.j * wz - q.k * wy,
            q.w * wy + q.k * wx - q.i * wz,
            q.w * wz + q.i * wy - q.j * wx,
        );
        let new_q = Quaternion::new(q.w + dq.w, q.i + dq.i, q.j + dq.j, q.k + dq.k);
        let norm =
            (new_q.w * new_q.w + new_q.i * new_q.i + new_q.j * new_q.j + new_q.k * new_q.k).sqrt();
        if norm > 1e-10 {
            use crate::math::Quat;
            self.transform.rotation = Quat::from_quaternion(Quaternion::new(
                new_q.w / norm,
                new_q.i / norm,
                new_q.j / norm,
                new_q.k / norm,
            ));
        }
    }
    /// Return an interpolated transform at sub-step fraction α ∈ \[0,1\].
    pub fn interpolated_transform(&self, alpha: Real) -> Transform {
        self.prev_transform.lerp(&self.transform, alpha)
    }
}
/// Six planes defining a view frustum (near, far, left, right, top, bottom).
///
/// Each plane is stored as `[nx, ny, nz, d]` where `dot(n, p) + d >= 0` means
/// the point `p` is inside that half-space.
#[derive(Debug, Clone)]
pub struct Frustum {
    /// Six frustum planes: `[nx, ny, nz, d]` each.
    pub planes: [[Real; 4]; 6],
}
impl Frustum {
    /// Create a frustum from six planes.
    pub fn new(planes: [[Real; 4]; 6]) -> Self {
        Self { planes }
    }
    /// Test whether the given AABB (min/max corners) is at least partially
    /// inside the frustum (conservative inclusion test).
    pub fn intersects_aabb(&self, min: [Real; 3], max: [Real; 3]) -> bool {
        for plane in &self.planes {
            let [nx, ny, nz, d] = *plane;
            let px = if nx >= 0.0 { max[0] } else { min[0] };
            let py = if ny >= 0.0 { max[1] } else { min[1] };
            let pz = if nz >= 0.0 { max[2] } else { min[2] };
            if nx * px + ny * py + nz * pz + d < 0.0 {
                return false;
            }
        }
        true
    }
}
/// Type of constraint between two bodies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintKind {
    /// Fixes the relative position of two anchor points.
    BallSocket,
    /// Fixes relative position and one rotational axis.
    Hinge,
    /// Maintains a minimum distance between bodies.
    DistanceMin,
    /// Maintains a maximum distance between bodies.
    DistanceMax,
    /// Maintains an exact distance between bodies.
    DistanceFixed,
    /// Contact constraint (collision response).
    Contact,
}
/// A constraint between two bodies (or one body and the world).
#[derive(Debug, Clone)]
pub struct Constraint {
    /// First body handle. Required.
    pub body_a: BodyHandle,
    /// Second body handle. `None` means the world.
    pub body_b: Option<BodyHandle>,
    /// Constraint type.
    pub kind: ConstraintKind,
    /// Anchor point in body A's local space.
    pub anchor_a: [Real; 3],
    /// Anchor point in body B's local space (or world space if `body_b` is None).
    pub anchor_b: [Real; 3],
    /// Target value (e.g., distance for distance constraints).
    pub target: Real,
    /// Compliance (softness) parameter; 0 = rigid.
    pub compliance: Real,
    /// Whether this constraint is currently active.
    pub active: bool,
}
impl Constraint {
    /// Create a ball-socket constraint between two bodies.
    pub fn ball_socket(
        body_a: BodyHandle,
        body_b: BodyHandle,
        anchor_a: [Real; 3],
        anchor_b: [Real; 3],
    ) -> Self {
        Self {
            body_a,
            body_b: Some(body_b),
            kind: ConstraintKind::BallSocket,
            anchor_a,
            anchor_b,
            target: 0.0,
            compliance: 0.0,
            active: true,
        }
    }
    /// Create a fixed-distance constraint between two bodies.
    pub fn distance(
        body_a: BodyHandle,
        body_b: BodyHandle,
        anchor_a: [Real; 3],
        anchor_b: [Real; 3],
        distance: Real,
    ) -> Self {
        Self {
            body_a,
            body_b: Some(body_b),
            kind: ConstraintKind::DistanceFixed,
            anchor_a,
            anchor_b,
            target: distance,
            compliance: 0.0,
            active: true,
        }
    }
}
/// Simulation mode for a body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyMode {
    /// Full rigid-body dynamics.
    Dynamic,
    /// Infinite-mass body that is moved by user code; affects dynamic bodies.
    Kinematic,
    /// Immovable body; participates in collision but never moves.
    Static,
    /// Body that is asleep and not integrated until woken.
    Sleeping,
}
/// Represents a single contact between two bodies.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B (`None` = static world).
    pub body_b: Option<BodyHandle>,
    /// Contact point in world space.
    pub point: [Real; 3],
    /// Contact normal pointing from B toward A.
    pub normal: [Real; 3],
    /// Penetration depth.
    pub depth: Real,
    /// Combined restitution coefficient.
    pub restitution: Real,
    /// Combined friction coefficient.
    pub friction: Real,
}
impl Contact {
    /// Create a new contact.
    pub fn new(
        body_a: BodyHandle,
        body_b: Option<BodyHandle>,
        point: [Real; 3],
        normal: [Real; 3],
        depth: Real,
        restitution: Real,
        friction: Real,
    ) -> Self {
        Self {
            body_a,
            body_b,
            point,
            normal,
            depth,
            restitution,
            friction,
        }
    }
}
/// Snapshot of the complete state of a single body.
#[derive(Debug, Clone)]
pub struct BodySnapshot {
    /// Body handle.
    pub handle: BodyHandle,
    /// Transform at snapshot time.
    pub transform: crate::types::Transform,
    /// Velocity at snapshot time.
    pub velocity: BodyVelocity,
    /// Accumulated forces.
    pub force: [Real; 3],
    /// Accumulated torques.
    pub torque: [Real; 3],
    /// Sleep timer.
    pub sleep_timer: Real,
    /// Body mode.
    pub mode: BodyMode,
}
impl BodySnapshot {
    /// Take a snapshot of a body.
    pub fn from_body(handle: BodyHandle, body: &Body) -> Self {
        Self {
            handle,
            transform: body.transform.clone(),
            velocity: body.velocity.clone(),
            force: body.force,
            torque: body.torque,
            sleep_timer: body.sleep_timer,
            mode: body.mode,
        }
    }
    /// Restore a snapshot into the given body.
    pub fn restore(&self, body: &mut Body) {
        body.transform = self.transform.clone();
        body.velocity = self.velocity.clone();
        body.force = self.force;
        body.torque = self.torque;
        body.sleep_timer = self.sleep_timer;
        body.mode = self.mode;
    }
}
/// A spatial hash grid for broad-phase collision detection.
///
/// The grid divides world space into cells of size `cell_size`.  Each body
/// occupies one or more cells; candidate pairs are extracted by looking up
/// all bodies in the same cell.
#[derive(Debug, Default)]
pub struct SpatialHashGrid {
    /// Cell size in world units.
    pub cell_size: Real,
    /// Map from cell (ix, iy, iz) to list of body slot indices.
    pub(super) cells: std::collections::HashMap<(i32, i32, i32), Vec<u64>>,
}
impl SpatialHashGrid {
    /// Create a new grid with the given cell size.
    pub fn new(cell_size: Real) -> Self {
        Self {
            cell_size,
            cells: std::collections::HashMap::new(),
        }
    }
    /// Quantise a world coordinate to a cell index.
    fn cell_coord(&self, x: Real) -> i32 {
        (x / self.cell_size).floor() as i32
    }
    /// Insert body `id` into all cells covered by its AABB.
    pub fn insert(&mut self, id: u64, min: [Real; 3], max: [Real; 3]) {
        let ix0 = self.cell_coord(min[0]);
        let iy0 = self.cell_coord(min[1]);
        let iz0 = self.cell_coord(min[2]);
        let ix1 = self.cell_coord(max[0]);
        let iy1 = self.cell_coord(max[1]);
        let iz1 = self.cell_coord(max[2]);
        for ix in ix0..=ix1 {
            for iy in iy0..=iy1 {
                for iz in iz0..=iz1 {
                    self.cells.entry((ix, iy, iz)).or_default().push(id);
                }
            }
        }
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.cells.clear();
    }
    /// Return candidate overlapping pairs (may contain duplicates).
    pub fn candidate_pairs(&self) -> Vec<(u64, u64)> {
        let mut pairs: std::collections::HashSet<(u64, u64)> = std::collections::HashSet::new();
        for bodies in self.cells.values() {
            for i in 0..bodies.len() {
                for j in i + 1..bodies.len() {
                    let a = bodies[i].min(bodies[j]);
                    let b = bodies[i].max(bodies[j]);
                    pairs.insert((a, b));
                }
            }
        }
        pairs.into_iter().collect()
    }
}
/// A human-readable snapshot of world state as a vector of lines.
#[derive(Debug, Clone)]
pub struct WorldSerialSnapshot {
    /// Lines of the snapshot.
    pub lines: Vec<String>,
}
/// A node in a static binary BVH built over triangle soup or body AABBs.
#[derive(Debug, Clone)]
pub enum BvhNode {
    /// Leaf node holding a list of body/primitive indices.
    Leaf {
        /// Bounding box of this leaf.
        aabb: crate::types::Aabb,
        /// Primitive indices contained here.
        indices: Vec<usize>,
    },
    /// Internal node with two children.
    Internal {
        /// Bounding box encompassing both children.
        aabb: crate::types::Aabb,
        /// Left child.
        left: Box<BvhNode>,
        /// Right child.
        right: Box<BvhNode>,
    },
}
impl BvhNode {
    /// Borrow the AABB of this node.
    pub fn aabb(&self) -> &crate::types::Aabb {
        match self {
            BvhNode::Leaf { aabb, .. } => aabb,
            BvhNode::Internal { aabb, .. } => aabb,
        }
    }
    /// Build a BVH over a list of (aabb, index) pairs using the SAH heuristic.
    /// Falls back to median split when SAH gives no benefit.
    pub fn build(items: &mut [(crate::types::Aabb, usize)]) -> Option<Self> {
        if items.is_empty() {
            return None;
        }
        let mut bound = items[0].0.clone();
        for (a, _) in items.iter().skip(1) {
            bound = bound.merge(a);
        }
        if items.len() <= 4 {
            return Some(BvhNode::Leaf {
                aabb: bound,
                indices: items.iter().map(|(_, i)| *i).collect(),
            });
        }
        let axis = bound.longest_axis();
        items.sort_by(|(a, _), (b, _)| {
            let ca = match axis {
                0 => (a.min.x + a.max.x) * 0.5,
                1 => (a.min.y + a.max.y) * 0.5,
                _ => (a.min.z + a.max.z) * 0.5,
            };
            let cb = match axis {
                0 => (b.min.x + b.max.x) * 0.5,
                1 => (b.min.y + b.max.y) * 0.5,
                _ => (b.min.z + b.max.z) * 0.5,
            };
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = items.len() / 2;
        let left = BvhNode::build(&mut items[..mid])?;
        let right = BvhNode::build(&mut items[mid..])?;
        Some(BvhNode::Internal {
            aabb: bound,
            left: Box::new(left),
            right: Box::new(right),
        })
    }
    /// Query all leaf indices whose AABB overlaps `query`.
    pub fn query_aabb(&self, query: &crate::types::Aabb, out: &mut Vec<usize>) {
        if !self.aabb().intersects(query) {
            return;
        }
        match self {
            BvhNode::Leaf { indices, .. } => {
                out.extend_from_slice(indices);
            }
            BvhNode::Internal { left, right, .. } => {
                left.query_aabb(query, out);
                right.query_aabb(query, out);
            }
        }
    }
    /// Query all leaf indices whose AABB is hit by `ray`.
    pub fn query_ray(&self, ray: &crate::types::Ray, out: &mut Vec<usize>) {
        if self.aabb().ray_intersect(ray).is_none() {
            return;
        }
        match self {
            BvhNode::Leaf { indices, .. } => {
                out.extend_from_slice(indices);
            }
            BvhNode::Internal { left, right, .. } => {
                left.query_ray(ray, out);
                right.query_ray(ray, out);
            }
        }
    }
}
/// Events emitted during simulation.
#[derive(Debug, Clone)]
pub enum PhysicsEvent {
    /// A collision began between two bodies.
    CollisionStarted {
        /// Handle of the first colliding body.
        body_a: BodyHandle,
        /// Handle of the second colliding body.
        body_b: BodyHandle,
        /// Contact normal (from A to B).
        normal: [Real; 3],
        /// Penetration depth.
        depth: Real,
    },
    /// An ongoing collision ended.
    CollisionEnded {
        /// Handle of the first body.
        body_a: BodyHandle,
        /// Handle of the second body.
        body_b: BodyHandle,
    },
    /// A body's speed fell below the sleep threshold.
    BodySlept {
        /// Handle of the sleeping body.
        body: BodyHandle,
    },
    /// A previously sleeping body was woken.
    BodyWoke {
        /// Handle of the woken body.
        body: BodyHandle,
    },
    /// A constraint was violated beyond its compliance range.
    ConstraintViolated {
        /// Index of the violated constraint.
        constraint_index: usize,
        /// Magnitude of the violation.
        violation: Real,
    },
}
/// Compact statistics snapshot for the world.
#[derive(Debug, Clone)]
pub struct WorldStats {
    /// Number of live bodies.
    pub body_count: usize,
    /// Number of active constraints.
    pub constraint_count: usize,
    /// Number of islands.
    pub island_count: usize,
    /// Approximate memory used by body arena (in bytes).
    pub body_arena_bytes: usize,
    /// Total simulation time elapsed.
    pub simulation_time: Real,
}
/// Linear/angular velocity state for a rigid body.
#[derive(Debug, Clone, Default)]
pub struct BodyVelocity {
    /// Linear velocity (m/s).
    pub linear: [Real; 3],
    /// Angular velocity (rad/s).
    pub angular: [Real; 3],
}
impl BodyVelocity {
    /// Create a new velocity state.
    pub fn new(linear: [Real; 3], angular: [Real; 3]) -> Self {
        Self { linear, angular }
    }
    /// Return the squared magnitude of the linear speed.
    pub fn linear_speed_sq(&self) -> Real {
        self.linear[0] * self.linear[0]
            + self.linear[1] * self.linear[1]
            + self.linear[2] * self.linear[2]
    }
    /// Return the squared magnitude of the angular speed.
    pub fn angular_speed_sq(&self) -> Real {
        self.angular[0] * self.angular[0]
            + self.angular[1] * self.angular[1]
            + self.angular[2] * self.angular[2]
    }
    /// Apply a linear impulse (change in linear velocity) given inverse mass.
    pub fn apply_linear_impulse(&mut self, impulse: [Real; 3], inv_mass: Real) {
        self.linear[0] += impulse[0] * inv_mass;
        self.linear[1] += impulse[1] * inv_mass;
        self.linear[2] += impulse[2] * inv_mass;
    }
    /// Apply an angular impulse (change in angular velocity) given inverse inertia scalar.
    pub fn apply_angular_impulse(&mut self, torque_impulse: [Real; 3], inv_inertia: Real) {
        self.angular[0] += torque_impulse[0] * inv_inertia;
        self.angular[1] += torque_impulse[1] * inv_inertia;
        self.angular[2] += torque_impulse[2] * inv_inertia;
    }
    /// Damp velocities by the given linear/angular damping factors.
    pub fn damp(&mut self, linear_damping: Real, angular_damping: Real, dt: Real) {
        let lf = (1.0 - linear_damping * dt).max(0.0);
        let af = (1.0 - angular_damping * dt).max(0.0);
        for c in &mut self.linear {
            *c *= lf;
        }
        for c in &mut self.angular {
            *c *= af;
        }
    }
    /// Zero all velocities (put body at rest).
    pub fn zero(&mut self) {
        self.linear = [0.0; 3];
        self.angular = [0.0; 3];
    }
}
/// An island is a group of dynamically connected bodies that can be solved
/// independently from the rest of the world.
#[derive(Debug, Clone, Default)]
pub struct Island {
    /// Indices into the body storage that belong to this island.
    pub body_indices: Vec<usize>,
    /// Whether all bodies in this island are sleeping.
    pub all_sleeping: bool,
}
impl Island {
    /// Create an empty island.
    pub fn new() -> Self {
        Self::default()
    }
    /// Return the number of bodies in this island.
    pub fn size(&self) -> usize {
        self.body_indices.len()
    }
}
/// A generational arena storing bodies.
#[derive(Debug, Default)]
pub struct BodyArena {
    pub(super) slots: Vec<BodySlot>,
    pub(super) free_head: Option<usize>,
    pub(super) count: usize,
}
impl BodyArena {
    /// Create an empty arena.
    pub fn new() -> Self {
        Self::default()
    }
    /// Insert a body and return its handle.
    pub fn insert(&mut self, body: Body) -> BodyHandle {
        if let Some(idx) = self.free_head {
            let (prev_gen, next) = match &self.slots[idx] {
                BodySlot::Free {
                    generation,
                    next_free,
                } => (*generation, *next_free),
                _ => panic!("free_head pointed to occupied slot"),
            };
            let new_gen = prev_gen.wrapping_add(1);
            self.slots[idx] = BodySlot::Occupied {
                body: Box::new(body),
                generation: new_gen,
            };
            self.free_head = next;
            self.count += 1;
            BodyHandle::new(idx as u32, new_gen)
        } else {
            let idx = self.slots.len();
            let generation = 0u32;
            self.slots.push(BodySlot::Occupied {
                body: Box::new(body),
                generation,
            });
            self.count += 1;
            BodyHandle::new(idx as u32, generation)
        }
    }
    /// Remove a body and free its slot.
    pub fn remove(&mut self, handle: BodyHandle) -> Option<Body> {
        let idx = handle.index as usize;
        if idx >= self.slots.len() {
            return None;
        }
        let slot = &self.slots[idx];
        match slot {
            BodySlot::Occupied { generation, .. } if *generation == handle.generation => {}
            _ => return None,
        }
        let old = std::mem::replace(
            &mut self.slots[idx],
            BodySlot::Free {
                next_free: self.free_head,
                generation: handle.generation,
            },
        );
        self.free_head = Some(idx);
        self.count -= 1;
        match old {
            BodySlot::Occupied { body, .. } => Some(*body),
            _ => unreachable!(),
        }
    }
    /// Get an immutable reference to a body by handle.
    pub fn get(&self, handle: BodyHandle) -> Option<&Body> {
        let idx = handle.index as usize;
        match self.slots.get(idx)? {
            BodySlot::Occupied { body, generation } if *generation == handle.generation => {
                Some(body.as_ref())
            }
            _ => None,
        }
    }
    /// Get a mutable reference to a body by handle.
    pub fn get_mut(&mut self, handle: BodyHandle) -> Option<&mut Body> {
        let idx = handle.index as usize;
        match self.slots.get_mut(idx)? {
            BodySlot::Occupied { body, generation } if *generation == handle.generation => {
                Some(body.as_mut())
            }
            _ => None,
        }
    }
    /// Return the number of live bodies.
    pub fn len(&self) -> usize {
        self.count
    }
    /// Whether the arena is empty.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
    /// Iterate over all live body indices.
    pub fn live_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.slots.iter().enumerate().filter_map(|(i, s)| match s {
            BodySlot::Occupied { .. } => Some(i),
            _ => None,
        })
    }
    /// Iterate over all live bodies (index, handle, body).
    pub fn iter(&self) -> impl Iterator<Item = (usize, BodyHandle, &Body)> + '_ {
        self.slots.iter().enumerate().filter_map(|(i, s)| match s {
            BodySlot::Occupied { body, generation } => {
                Some((i, BodyHandle::new(i as u32, *generation), body.as_ref()))
            }
            _ => None,
        })
    }
}
/// Top-level physics world that owns all simulation state.
#[derive(Debug)]
pub struct PhysicsWorld {
    /// Simulation time accumulated so far.
    pub time: Real,
    /// Physics configuration.
    pub config: PhysicsConfig,
    /// Fixed time step used for stepping.
    pub fixed_dt: Real,
    /// Accumulated time not yet simulated.
    pub accumulator: Real,
    /// All registered bodies.
    pub bodies: BodyArena,
    /// All constraints.
    pub constraints: Vec<Constraint>,
    /// Event queue (pending events from last step).
    pub events: VecDeque<PhysicsEvent>,
    /// Islands computed at last step.
    pub islands: Vec<Island>,
    /// Maximum events held in the queue before oldest are dropped.
    pub max_events: usize,
}
impl PhysicsWorld {
    /// Create a new empty physics world with default config.
    pub fn new() -> Self {
        Self::default()
    }
    /// Create a new physics world with the given configuration.
    pub fn with_config(config: PhysicsConfig) -> Self {
        Self {
            config,
            ..Self::default()
        }
    }
    /// Add a body to the world and return its handle.
    pub fn add_body(&mut self, body: Body) -> BodyHandle {
        self.bodies.insert(body)
    }
    /// Remove a body from the world. Returns the body if the handle was valid.
    pub fn remove_body(&mut self, handle: BodyHandle) -> Option<Body> {
        self.constraints
            .retain(|c| c.body_a != handle && c.body_b != Some(handle));
        self.bodies.remove(handle)
    }
    /// Get an immutable reference to a body.
    pub fn body(&self, handle: BodyHandle) -> Option<&Body> {
        self.bodies.get(handle)
    }
    /// Get a mutable reference to a body.
    pub fn body_mut(&mut self, handle: BodyHandle) -> Option<&mut Body> {
        self.bodies.get_mut(handle)
    }
    /// Return the number of bodies currently in the world.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Add a constraint between bodies and return its index.
    pub fn add_constraint(&mut self, constraint: Constraint) -> usize {
        let idx = self.constraints.len();
        self.constraints.push(constraint);
        idx
    }
    /// Remove a constraint by index. Swaps with the last to avoid shifting.
    pub fn remove_constraint(&mut self, index: usize) {
        if index < self.constraints.len() {
            self.constraints.swap_remove(index);
        }
    }
    /// Push an event onto the queue, dropping the oldest if over capacity.
    pub fn push_event(&mut self, event: PhysicsEvent) {
        if self.events.len() >= self.max_events {
            self.events.pop_front();
        }
        self.events.push_back(event);
    }
    /// Drain all pending events from the queue.
    pub fn drain_events(&mut self) -> Vec<PhysicsEvent> {
        self.events.drain(..).collect()
    }
    /// Set the world gravity vector.
    pub fn set_gravity(&mut self, gx: Real, gy: Real, gz: Real) {
        use crate::math::Vec3;
        self.config.gravity = Vec3::new(gx, gy, gz);
    }
    /// Step the world forward by a fixed `TimeStep`.
    pub fn step(&mut self, time_step: &TimeStep) {
        self.step_dt(time_step.dt);
    }
    /// Step the world forward by the given delta time, using fixed substeps.
    /// Returns the number of substeps taken.
    pub fn step_with_accumulator(&mut self, dt: Real) -> u32 {
        self.accumulator += dt;
        let mut steps = 0u32;
        while self.accumulator >= self.fixed_dt {
            self.accumulator -= self.fixed_dt;
            self.step_dt(self.fixed_dt);
            steps += 1;
        }
        steps
    }
    /// Perform one full simulation substep of duration `dt`.
    fn step_dt(&mut self, dt: Real) {
        let gravity = [
            self.config.gravity.x,
            self.config.gravity.y,
            self.config.gravity.z,
        ];
        let lin_sleep = self.config.linear_sleep_threshold;
        let ang_sleep = self.config.angular_sleep_threshold;
        let t_sleep = self.config.time_before_sleep;
        let indices: Vec<usize> = self.bodies.live_indices().collect();
        for &idx in &indices {
            let slot = &mut self.bodies.slots[idx];
            if let BodySlot::Occupied { body, generation } = slot {
                let handle = BodyHandle::new(idx as u32, *generation);
                body.integrate_velocities(gravity, dt);
                body.integrate_positions(dt);
                body.clear_forces();
                if body.mode == BodyMode::Dynamic {
                    let fast = body.velocity.linear_speed_sq() > lin_sleep * lin_sleep
                        || body.velocity.angular_speed_sq() > ang_sleep * ang_sleep;
                    if fast {
                        if body.mode == BodyMode::Sleeping {
                            body.mode = BodyMode::Dynamic;
                            self.events
                                .push_back(PhysicsEvent::BodyWoke { body: handle });
                        }
                        body.sleep_timer = 0.0;
                    } else {
                        body.sleep_timer += dt;
                        if body.sleep_timer >= t_sleep {
                            body.mode = BodyMode::Sleeping;
                            body.velocity.zero();
                            self.events
                                .push_back(PhysicsEvent::BodySlept { body: handle });
                        }
                    }
                }
            }
        }
        self.time += dt;
    }
    /// Recompute simulation islands using union-find on constraint graph.
    /// Islands group bodies that interact so they can be solved independently.
    pub fn compute_islands(&mut self) {
        let n = self.bodies.slots.len();
        if n == 0 {
            self.islands.clear();
            return;
        }
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }
        fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
        for c in &self.constraints {
            if !c.active {
                continue;
            }
            let a_idx = c.body_a.index as usize;
            if let Some(b_handle) = c.body_b {
                let b_idx = b_handle.index as usize;
                if a_idx < n && b_idx < n {
                    union(&mut parent, a_idx, b_idx);
                }
            }
        }
        let mut island_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for idx in self.bodies.live_indices() {
            let root = find(&mut parent, idx);
            island_map.entry(root).or_default().push(idx);
        }
        self.islands = island_map
            .into_values()
            .map(|body_indices| {
                let all_sleeping = body_indices.iter().all(|&i| {
                    matches!(
                        & self.bodies.slots[i], BodySlot::Occupied { body, .. } if
                        body.mode == BodyMode::Sleeping || body.mode ==
                        BodyMode::Static
                    )
                });
                Island {
                    body_indices,
                    all_sleeping,
                }
            })
            .collect();
    }
    /// Return handles of all dynamic bodies whose axis-aligned bounding box
    /// (approximated as a sphere of radius `query_radius` around `center`)
    /// overlaps the given AABB.
    pub fn bodies_in_aabb(&self, query: &Aabb) -> Vec<BodyHandle> {
        let mut result = Vec::new();
        for (idx, handle, body) in self.bodies.iter() {
            let _ = idx;
            let p = &body.transform.position;
            if p.x >= query.min.x
                && p.x <= query.max.x
                && p.y >= query.min.y
                && p.y <= query.max.y
                && p.z >= query.min.z
                && p.z <= query.max.z
            {
                result.push(handle);
            }
        }
        result
    }
    /// Return handles of all bodies within `radius` of `center`.
    pub fn bodies_in_radius(&self, center: [Real; 3], radius: Real) -> Vec<BodyHandle> {
        let r2 = radius * radius;
        let mut result = Vec::new();
        for (_, handle, body) in self.bodies.iter() {
            let p = &body.transform.position;
            let dx = p.x - center[0];
            let dy = p.y - center[1];
            let dz = p.z - center[2];
            if dx * dx + dy * dy + dz * dz <= r2 {
                result.push(handle);
            }
        }
        result
    }
    /// Total kinetic energy (translational + rotational) in the world.
    pub fn total_kinetic_energy(&self) -> Real {
        let mut ke = 0.0;
        for (_, _, body) in self.bodies.iter() {
            if body.mass > 0.0 {
                let v2 = body.velocity.linear_speed_sq();
                let w2 = body.velocity.angular_speed_sq();
                ke += 0.5 * body.mass * v2 + 0.5 * body.inertia * w2;
            }
        }
        ke
    }
    /// Total linear momentum of the world.
    pub fn total_momentum(&self) -> [Real; 3] {
        let mut p = [0.0; 3];
        for (_, _, body) in self.bodies.iter() {
            for (k, pk) in p.iter_mut().enumerate() {
                *pk += body.mass * body.velocity.linear[k];
            }
        }
        p
    }
}
impl PhysicsWorld {
    /// Take a complete snapshot of the world state.
    pub fn snapshot(&self) -> WorldSnapshot {
        let bodies = self
            .bodies
            .iter()
            .map(|(_, handle, body)| BodySnapshot::from_body(handle, body))
            .collect();
        WorldSnapshot {
            time: self.time,
            accumulator: self.accumulator,
            bodies,
        }
    }
    /// Restore the world to a previously taken snapshot.
    pub fn restore_snapshot(&mut self, snap: &WorldSnapshot) {
        self.time = snap.time;
        self.accumulator = snap.accumulator;
        for bs in &snap.bodies {
            if let Some(b) = self.bodies.get_mut(bs.handle) {
                bs.restore(b);
            }
        }
    }
    /// Apply all force fields to all dynamic bodies.
    pub fn apply_force_fields(&mut self, fields: &[ForceField]) {
        let live: Vec<(usize, BodyHandle)> = self.bodies.iter().map(|(i, h, _)| (i, h)).collect();
        for (idx, _handle) in live {
            let (pos, vel, mass, is_dynamic) = {
                let slot = &self.bodies.slots[idx];
                match slot {
                    BodySlot::Occupied { body, .. } => {
                        let p = body.transform.position;
                        let pos = [p.x, p.y, p.z];
                        let vel = body.velocity.linear;
                        (pos, vel, body.mass, body.mode == BodyMode::Dynamic)
                    }
                    _ => continue,
                }
            };
            if !is_dynamic {
                continue;
            }
            let mut total_force = [0.0_f64; 3];
            for field in fields {
                let f = field.force_on(pos, vel, mass);
                total_force[0] += f[0];
                total_force[1] += f[1];
                total_force[2] += f[2];
            }
            if let BodySlot::Occupied { body, .. } = &mut self.bodies.slots[idx] {
                body.force[0] += total_force[0];
                body.force[1] += total_force[1];
                body.force[2] += total_force[2];
            }
        }
    }
    /// Resolve a list of contacts (sequential impulse, one pass).
    pub fn resolve_contacts(&mut self, contacts: &[Contact]) {
        for c in contacts {
            resolve_contact(&mut self.bodies, c);
        }
    }
    /// Compute total angular momentum of the world.
    pub fn total_angular_momentum(&self) -> [Real; 3] {
        let mut l = [0.0_f64; 3];
        for (_, _, body) in self.bodies.iter() {
            if body.mass <= 0.0 {
                continue;
            }
            let p = &body.transform.position;
            let v = &body.velocity.linear;
            l[0] += body.mass * (p.y * v[2] - p.z * v[1]);
            l[1] += body.mass * (p.z * v[0] - p.x * v[2]);
            l[2] += body.mass * (p.x * v[1] - p.y * v[0]);
            l[0] += body.inertia * body.velocity.angular[0];
            l[1] += body.inertia * body.velocity.angular[1];
            l[2] += body.inertia * body.velocity.angular[2];
        }
        l
    }
    /// Apply a damping "blast" to all bodies within `radius` of `center`,
    /// scaling their velocities by `factor`.
    pub fn damp_region(&mut self, center: [Real; 3], radius: Real, factor: Real) {
        let r2 = radius * radius;
        for slot in &mut self.bodies.slots {
            if let BodySlot::Occupied { body, .. } = slot {
                let p = &body.transform.position;
                let dx = p.x - center[0];
                let dy = p.y - center[1];
                let dz = p.z - center[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    for k in 0..3 {
                        body.velocity.linear[k] *= factor;
                        body.velocity.angular[k] *= factor;
                    }
                }
            }
        }
    }
    /// Return the body with the highest kinetic energy, or `None` if empty.
    pub fn most_energetic_body(&self) -> Option<BodyHandle> {
        let mut best_ke = -1.0_f64;
        let mut best = None;
        for (_, handle, body) in self.bodies.iter() {
            let ke = 0.5 * body.mass * body.velocity.linear_speed_sq()
                + 0.5 * body.inertia * body.velocity.angular_speed_sq();
            if ke > best_ke {
                best_ke = ke;
                best = Some(handle);
            }
        }
        best
    }
    /// Apply a position correction to a body (used after PBD/XPBD solve).
    pub fn apply_position_correction(&mut self, handle: BodyHandle, delta: [Real; 3]) {
        if let Some(b) = self.bodies.get_mut(handle)
            && b.mode == BodyMode::Dynamic
        {
            b.transform.position.x += delta[0];
            b.transform.position.y += delta[1];
            b.transform.position.z += delta[2];
        }
    }
    /// Set the velocity of a body by handle.
    pub fn set_velocity(&mut self, handle: BodyHandle, linear: [Real; 3], angular: [Real; 3]) {
        if let Some(b) = self.bodies.get_mut(handle) {
            b.velocity.linear = linear;
            b.velocity.angular = angular;
        }
    }
    /// Teleport a body to a new position without affecting velocity.
    pub fn set_position(&mut self, handle: BodyHandle, position: [Real; 3]) {
        if let Some(b) = self.bodies.get_mut(handle) {
            b.transform.position = crate::math::Vec3::new(position[0], position[1], position[2]);
        }
    }
    /// Wake all sleeping bodies in the world.
    pub fn wake_all(&mut self) {
        for slot in &mut self.bodies.slots {
            if let BodySlot::Occupied { body, .. } = slot
                && body.mode == BodyMode::Sleeping
            {
                body.mode = BodyMode::Dynamic;
                body.sleep_timer = 0.0;
            }
        }
    }
    /// Return a summary of the world state for diagnostics.
    pub fn diagnostics(&self) -> WorldDiagnostics {
        let mut n_dynamic = 0usize;
        let mut n_static = 0usize;
        let mut n_sleeping = 0usize;
        let mut n_kinematic = 0usize;
        for (_, _, b) in self.bodies.iter() {
            match b.mode {
                BodyMode::Dynamic => n_dynamic += 1,
                BodyMode::Static => n_static += 1,
                BodyMode::Sleeping => n_sleeping += 1,
                BodyMode::Kinematic => n_kinematic += 1,
            }
        }
        WorldDiagnostics {
            total_bodies: self.bodies.len(),
            n_dynamic,
            n_static,
            n_sleeping,
            n_kinematic,
            n_constraints: self.constraints.len(),
            n_islands: self.islands.len(),
            simulation_time: self.time,
            total_kinetic_energy: self.total_kinetic_energy(),
        }
    }
}
impl PhysicsWorld {
    /// Return all body handles whose position is within `radius` of `center`.
    pub fn query_sphere(&self, center: [Real; 3], radius: Real) -> Vec<BodyHandle> {
        let r2 = radius * radius;
        self.bodies
            .iter()
            .filter_map(|(_, handle, body)| {
                let p = body.transform.position;
                let dx = p.x - center[0];
                let dy = p.y - center[1];
                let dz = p.z - center[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    Some(handle)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Return all body handles whose AABB is at least partially inside `frustum`.
    ///
    /// Bodies without a stored `aabb` use a small unit-box centred at their position.
    pub fn query_frustum(&self, frustum: &Frustum) -> Vec<BodyHandle> {
        self.bodies
            .iter()
            .filter_map(|(_, handle, body)| {
                let p = body.transform.position;
                let half = 0.5_f64;
                let min = [p.x - half, p.y - half, p.z - half];
                let max = [p.x + half, p.y + half, p.z + half];
                if frustum.intersects_aabb(min, max) {
                    Some(handle)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Compute a compact statistics snapshot of the current world state.
    pub fn compute_stats(&self) -> WorldStats {
        WorldStats {
            body_count: self.bodies.len(),
            constraint_count: self.constraints.len(),
            island_count: self.islands.len(),
            body_arena_bytes: self.bodies.slots.len() * 256,
            simulation_time: self.time,
        }
    }
    /// Serialise the world state to a list of human-readable text lines.
    ///
    /// Useful for debugging, logging, and lightweight persistence.
    pub fn serialize_snapshot(&self) -> WorldSerialSnapshot {
        let mut lines = Vec::new();
        lines.push(format!("time={:.6}", self.time));
        lines.push(format!("body_count={}", self.bodies.len()));
        lines.push(format!("constraint_count={}", self.constraints.len()));
        lines.push(format!("island_count={}", self.islands.len()));
        for (_, handle, body) in self.bodies.iter() {
            let p = body.transform.position;
            let v = &body.velocity;
            lines
                .push(
                    format!(
                        "body idx={} gen={} pos=[{:.4},{:.4},{:.4}] vel=[{:.4},{:.4},{:.4}] mass={:.4} mode={:?}",
                        handle.index, handle.generation, p.x, p.y, p.z, v.linear[0], v
                        .linear[1], v.linear[2], body.mass, body.mode,
                    ),
                );
        }
        WorldSerialSnapshot { lines }
    }
}
/// Type of continuous force field applied to bodies.
#[derive(Debug, Clone)]
pub enum ForceFieldKind {
    /// Constant directional force (e.g., wind).
    Constant {
        /// Force vector applied to all matching bodies.
        force: [Real; 3],
    },
    /// Point attractor / repulsor.
    PointGravity {
        /// Position of the source.
        position: [Real; 3],
        /// Gravitational constant * source mass (G*M). Negative = repulsion.
        gm: Real,
        /// Softening parameter to avoid singularity at r=0.
        epsilon: Real,
    },
    /// Drag force proportional to velocity squared.
    QuadraticDrag {
        /// Drag coefficient.
        cd: Real,
    },
    /// Vortex field (rotation about an axis).
    Vortex {
        /// Center of the vortex.
        center: [Real; 3],
        /// Axis of rotation (unit vector).
        axis: [Real; 3],
        /// Angular velocity of the vortex.
        omega: Real,
        /// Radial falloff (force = omega*r * exp(-r/decay)).
        decay: Real,
    },
}
/// A world-level snapshot for full rollback.
#[derive(Debug, Clone)]
pub struct WorldSnapshot {
    /// Time at snapshot.
    pub time: Real,
    /// Accumulator at snapshot.
    pub accumulator: Real,
    /// Per-body snapshots.
    pub bodies: Vec<BodySnapshot>,
}
/// Summary diagnostic information about the physics world.
#[derive(Debug, Clone)]
pub struct WorldDiagnostics {
    /// Total number of bodies.
    pub total_bodies: usize,
    /// Number of dynamic bodies.
    pub n_dynamic: usize,
    /// Number of static bodies.
    pub n_static: usize,
    /// Number of sleeping bodies.
    pub n_sleeping: usize,
    /// Number of kinematic bodies.
    pub n_kinematic: usize,
    /// Number of active constraints.
    pub n_constraints: usize,
    /// Number of simulation islands.
    pub n_islands: usize,
    /// Total elapsed simulation time.
    pub simulation_time: Real,
    /// Current total kinetic energy.
    pub total_kinetic_energy: Real,
}
/// Slot in the generational body arena.
///
/// The `Occupied` variant boxes the body to avoid the `large_enum_variant`
/// penalty: `Body` is ~200 bytes while `Free` is only ~16 bytes.
#[derive(Debug, Clone)]
pub(super) enum BodySlot {
    /// A live body occupying this slot.
    Occupied {
        /// The boxed body data.
        body: Box<Body>,
        /// Generation counter for handle validation.
        generation: u32,
    },
    /// A free slot in the generational arena.
    Free {
        /// Next free slot index, or `None` if this is the last.
        next_free: Option<usize>,
        /// Generation counter carried across free/occupy cycles.
        generation: u32,
    },
}
