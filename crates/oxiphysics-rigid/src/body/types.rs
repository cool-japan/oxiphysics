//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::{Quat, Real, Vec3};

use super::rigidbody_type::RigidBody;

/// A fixed-size ring buffer holding body states for sub-step interpolation.
#[derive(Debug, Clone)]
pub struct StateInterpolator {
    /// Stored (time, position, rotation) samples.
    pub(super) samples: Vec<(Real, Vec3, Quat)>,
    /// Maximum number of samples to keep.
    pub(super) max_samples: usize,
}
impl StateInterpolator {
    /// Create a new interpolator with the given sample capacity.
    #[allow(dead_code)]
    pub fn new(max_samples: usize) -> Self {
        Self {
            samples: Vec::new(),
            max_samples: max_samples.max(2),
        }
    }
    /// Record a body state sample.
    #[allow(dead_code)]
    pub fn push(&mut self, time: Real, position: Vec3, rotation: Quat) {
        if self.samples.len() >= self.max_samples {
            self.samples.remove(0);
        }
        self.samples.push((time, position, rotation));
    }
    /// Push from a `RigidBody` snapshot.
    #[allow(dead_code)]
    pub fn push_from_body(&mut self, time: Real, body: &RigidBody) {
        self.push(time, body.transform.position, body.transform.rotation);
    }
    /// Linearly interpolate position at an arbitrary time `t`.
    /// Returns `None` if fewer than 2 samples are available or `t` is out of range.
    #[allow(dead_code)]
    pub fn interpolate_position(&self, t: Real) -> Option<Vec3> {
        if self.samples.len() < 2 {
            return None;
        }
        for i in 0..self.samples.len() - 1 {
            let (t0, p0, _) = self.samples[i];
            let (t1, p1, _) = self.samples[i + 1];
            if t >= t0 && t <= t1 {
                let dt = t1 - t0;
                if dt < 1e-30 {
                    return Some(p0);
                }
                let alpha = (t - t0) / dt;
                return Some(p0 + (p1 - p0) * alpha);
            }
        }
        None
    }
    /// Spherical linear interpolation (slerp) of rotation at time `t`.
    #[allow(dead_code)]
    pub fn interpolate_rotation(&self, t: Real) -> Option<Quat> {
        if self.samples.len() < 2 {
            return None;
        }
        for i in 0..self.samples.len() - 1 {
            let (t0, _, r0) = self.samples[i];
            let (t1, _, r1) = self.samples[i + 1];
            if t >= t0 && t <= t1 {
                let dt = t1 - t0;
                if dt < 1e-30 {
                    return Some(r0);
                }
                let alpha = ((t - t0) / dt) as f32;
                return Some(r0.slerp(&r1, alpha as Real));
            }
        }
        None
    }
    /// Number of samples stored.
    #[allow(dead_code)]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }
    /// Time range of stored samples `(t_min, t_max)`.
    #[allow(dead_code)]
    pub fn time_range(&self) -> Option<(Real, Real)> {
        if self.samples.is_empty() {
            return None;
        }
        Some((
            self.samples
                .first()
                .expect("collection should not be empty")
                .0,
            self.samples
                .last()
                .expect("collection should not be empty")
                .0,
        ))
    }
    /// Clear all samples.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.samples.clear();
    }
}
/// Tracks the position of the center of mass over time for analysis.
#[derive(Debug, Clone, Default)]
pub struct CenterOfMassTracker {
    /// History of center-of-mass positions (time, position).
    pub history: Vec<(Real, Vec3)>,
    /// Maximum history length before oldest entries are dropped.
    pub max_history: usize,
}
impl CenterOfMassTracker {
    /// Create a tracker with a specified history length.
    #[allow(dead_code)]
    pub fn new(max_history: usize) -> Self {
        Self {
            history: Vec::new(),
            max_history,
        }
    }
    /// Record the current center of mass.
    #[allow(dead_code)]
    pub fn record(&mut self, time: Real, position: Vec3) {
        if self.max_history > 0 && self.history.len() >= self.max_history {
            self.history.remove(0);
        }
        self.history.push((time, position));
    }
    /// Most recent recorded position.
    #[allow(dead_code)]
    pub fn latest(&self) -> Option<Vec3> {
        self.history.last().map(|(_, p)| *p)
    }
    /// Total distance traveled (sum of step-to-step displacements).
    #[allow(dead_code)]
    pub fn total_path_length(&self) -> Real {
        if self.history.len() < 2 {
            return 0.0;
        }
        self.history
            .windows(2)
            .map(|w| (w[1].1 - w[0].1).norm())
            .sum()
    }
    /// Average speed over the recorded history (m/s).
    #[allow(dead_code)]
    pub fn average_speed(&self) -> Real {
        if self.history.len() < 2 {
            return 0.0;
        }
        let first_t = self.history.first().map(|(t, _)| *t).unwrap_or(0.0);
        let last_t = self.history.last().map(|(t, _)| *t).unwrap_or(0.0);
        let dt = last_t - first_t;
        if dt < 1e-30 {
            return 0.0;
        }
        self.total_path_length() / dt
    }
    /// Number of recorded samples.
    #[allow(dead_code)]
    pub fn sample_count(&self) -> usize {
        self.history.len()
    }
    /// Clear all history.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.history.clear();
    }
}
/// Rolling without slipping constraint for a rigid body on a surface.
///
/// Enforces the kinematic constraint: v_contact = 0, i.e.,
///   v_cm + omega × r = 0  (at the contact point)
/// Applied as a velocity-level correction each frame.
#[derive(Debug, Clone, Copy)]
pub struct RollingConstraint {
    /// Contact point radius (distance from CoM to contact).
    pub radius: Real,
    /// Surface normal (unit vector, pointing away from surface).
    pub surface_normal: [Real; 3],
    /// Rolling friction coefficient (dimensionless, for energy dissipation).
    pub rolling_friction: Real,
}
impl RollingConstraint {
    /// Create a rolling constraint for a sphere of given radius.
    pub fn sphere(radius: Real) -> Self {
        Self {
            radius,
            surface_normal: [0.0, 1.0, 0.0],
            rolling_friction: 0.01,
        }
    }
    /// Slip velocity at the contact point: v_slip = v_cm + omega × (-r * n).
    pub fn slip_velocity(&self, body: &RigidBody) -> [Real; 3] {
        let n = self.surface_normal;
        let r_vec = [
            -self.radius * n[0],
            -self.radius * n[1],
            -self.radius * n[2],
        ];
        let [ox, oy, oz] = [
            body.angular_velocity.x,
            body.angular_velocity.y,
            body.angular_velocity.z,
        ];
        let [rx, ry, rz] = r_vec;
        let cross = [oy * rz - oz * ry, oz * rx - ox * rz, ox * ry - oy * rx];
        [
            body.velocity.x + cross[0],
            body.velocity.y + cross[1],
            body.velocity.z + cross[2],
        ]
    }
    /// Apply rolling constraint (zero slip) to body velocity.
    ///
    /// Projects out the tangential velocity component and adjusts omega.
    pub fn enforce(&self, body: &mut RigidBody) {
        if !body.is_dynamic() {
            return;
        }
        let slip = self.slip_velocity(body);
        let n = self.surface_normal;
        let slip_n = slip[0] * n[0] + slip[1] * n[1] + slip[2] * n[2];
        let slip_tan = [
            slip[0] - slip_n * n[0],
            slip[1] - slip_n * n[1],
            slip[2] - slip_n * n[2],
        ];
        let slip_speed =
            (slip_tan[0] * slip_tan[0] + slip_tan[1] * slip_tan[1] + slip_tan[2] * slip_tan[2])
                .sqrt();
        if slip_speed < 1e-12 {
            return;
        }
        let alpha = 2.0 / 3.0;
        body.velocity.x -= alpha * slip_tan[0];
        body.velocity.y -= alpha * slip_tan[1];
        body.velocity.z -= alpha * slip_tan[2];
        let inv_r = if self.radius > 1e-12 {
            1.0 / self.radius
        } else {
            0.0
        };
        body.angular_velocity.x += alpha * slip_tan[2] * inv_r;
        body.angular_velocity.y -= alpha * slip_tan[0] * inv_r;
        body.angular_velocity.z -= alpha * slip_tan[1] * inv_r;
    }
    /// Rolling friction torque magnitude: tau = mu_r * m * g * r.
    pub fn rolling_friction_torque(&self, body_mass: Real, gravity: Real) -> Real {
        self.rolling_friction * body_mass * gravity * self.radius
    }
}
/// Serializable snapshot of a rigid body's state.
#[derive(Debug, Clone)]
pub struct BodySnapshot {
    /// Position \[x, y, z\].
    pub position: [Real; 3],
    /// Rotation quaternion \[x, y, z, w\].
    pub rotation: [Real; 4],
    /// Linear velocity \[vx, vy, vz\].
    pub velocity: [Real; 3],
    /// Angular velocity \[wx, wy, wz\].
    pub angular_velocity: [Real; 3],
    /// Mass.
    pub mass: Real,
    /// Body type as string tag.
    pub body_type: &'static str,
    /// Activity state as string tag.
    pub state: &'static str,
    /// Linear damping.
    pub linear_damping: Real,
    /// Angular damping.
    pub angular_damping: Real,
    /// Gravity scale.
    pub gravity_scale: Real,
}
/// Compound rigid body assembled from multiple primitives.
///
/// The compound body computes its centre of mass and inertia tensor using
/// the parallel axis theorem.
#[derive(Debug, Clone)]
pub struct CompoundRigidBody {
    /// Component primitives.
    pub primitives: Vec<CompoundPrimitive>,
    /// Total mass.
    pub total_mass: Real,
    /// Centre of mass in local body space.
    pub local_com: [Real; 3],
    /// Combined inertia tensor diagonal (principal axes) at the CoM.
    pub inertia_diag: [Real; 3],
}
impl CompoundRigidBody {
    /// Build a compound body from the given primitives.
    pub fn from_primitives(primitives: Vec<CompoundPrimitive>) -> Self {
        let total_mass: Real = primitives.iter().map(|p| p.mass).sum();
        let com = if total_mass > 1e-30 {
            let cx: Real = primitives
                .iter()
                .map(|p| p.mass * p.local_position[0])
                .sum::<Real>()
                / total_mass;
            let cy: Real = primitives
                .iter()
                .map(|p| p.mass * p.local_position[1])
                .sum::<Real>()
                / total_mass;
            let cz: Real = primitives
                .iter()
                .map(|p| p.mass * p.local_position[2])
                .sum::<Real>()
                / total_mass;
            [cx, cy, cz]
        } else {
            [0.0; 3]
        };
        let mut ixx = 0.0;
        let mut iyy = 0.0;
        let mut izz = 0.0;
        for prim in &primitives {
            let dx = prim.local_position[0] - com[0];
            let dy = prim.local_position[1] - com[1];
            let dz = prim.local_position[2] - com[2];
            let m = prim.mass;
            ixx += prim.local_inertia_diag[0] + m * (dy * dy + dz * dz);
            iyy += prim.local_inertia_diag[1] + m * (dx * dx + dz * dz);
            izz += prim.local_inertia_diag[2] + m * (dx * dx + dy * dy);
        }
        Self {
            primitives,
            total_mass,
            local_com: com,
            inertia_diag: [ixx, iyy, izz],
        }
    }
    /// Build a RigidBody from this compound body.
    pub fn to_rigid_body(&self) -> RigidBody {
        use oxiphysics_core::math::Mat3;
        let mut body = RigidBody::new(self.total_mass);
        let [ixx, iyy, izz] = self.inertia_diag;
        let inertia = Mat3::from_diagonal(&oxiphysics_core::math::Vec3::new(ixx, iyy, izz));
        body.local_inertia = inertia;
        body.update_world_inertia();
        body
    }
    /// Number of primitives.
    pub fn primitive_count(&self) -> usize {
        self.primitives.len()
    }
    /// True if the compound body is valid (mass > 0).
    pub fn is_valid(&self) -> bool {
        self.total_mass > 0.0 && self.inertia_diag.iter().all(|&i| i >= 0.0)
    }
}
/// Kinematic target for scripted body motion.
#[derive(Debug, Clone)]
pub struct KinematicTarget {
    /// Target position to move towards.
    pub target_position: Vec3,
    /// Target rotation to rotate towards.
    pub target_rotation: Quat,
}
/// A named persistent force that is re-applied every integration step.
///
/// Use this for e.g. wind forces, thrusters, or spring forces that should
/// persist across frames without needing to re-apply them manually each step.
#[derive(Debug, Clone)]
pub struct PersistentForce {
    /// Human-readable label.
    pub label: String,
    /// Force vector in world space.
    pub force: Vec3,
    /// Optional world-space application point (None = center of mass).
    pub point: Option<Vec3>,
    /// Whether this force is currently enabled.
    pub enabled: bool,
}
impl PersistentForce {
    /// Create a new enabled persistent force at the center of mass.
    #[allow(dead_code)]
    pub fn new(label: &str, force: Vec3) -> Self {
        Self {
            label: label.to_string(),
            force,
            point: None,
            enabled: true,
        }
    }
    /// Create a persistent force applied at a world-space point.
    #[allow(dead_code)]
    pub fn at_point(label: &str, force: Vec3, point: Vec3) -> Self {
        Self {
            label: label.to_string(),
            force,
            point: Some(point),
            enabled: true,
        }
    }
    /// Disable this force.
    #[allow(dead_code)]
    pub fn disable(&mut self) {
        self.enabled = false;
    }
    /// Enable this force.
    #[allow(dead_code)]
    pub fn enable(&mut self) {
        self.enabled = true;
    }
}
/// A single contact event recorded in a body's history.
#[derive(Debug, Clone)]
pub struct ContactHistoryEntry {
    /// Simulation time of the contact.
    pub time: Real,
    /// Handle of the other body involved.
    pub other_handle: u32,
    /// Contact point in world space.
    pub point: Vec3,
    /// Contact normal (pointing from other to this body).
    pub normal: Vec3,
    /// Impulse magnitude applied.
    pub impulse: Real,
}
/// A constraint record attached to a body.
#[derive(Debug, Clone)]
pub struct BodyConstraint {
    /// Unique constraint identifier.
    pub id: u32,
    /// Constraint type.
    pub constraint_type: ConstraintType,
    /// Human-readable description.
    pub description: String,
    /// Whether the constraint is currently active.
    pub active: bool,
}
impl BodyConstraint {
    /// Create a new active constraint.
    #[allow(dead_code)]
    pub fn new(id: u32, constraint_type: ConstraintType, description: &str) -> Self {
        Self {
            id,
            constraint_type,
            description: description.to_string(),
            active: true,
        }
    }
    /// Deactivate this constraint.
    #[allow(dead_code)]
    pub fn deactivate(&mut self) {
        self.active = false;
    }
    /// Activate this constraint.
    #[allow(dead_code)]
    pub fn activate(&mut self) {
        self.active = true;
    }
}
/// List of constraints attached to a body.
#[derive(Debug, Clone, Default)]
pub struct ConstraintList {
    pub(super) constraints: Vec<BodyConstraint>,
    pub(super) next_id: u32,
}
impl ConstraintList {
    /// Create a new empty constraint list.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a constraint, returning its ID.
    #[allow(dead_code)]
    pub fn add(&mut self, constraint_type: ConstraintType, description: &str) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        self.constraints
            .push(BodyConstraint::new(id, constraint_type, description));
        id
    }
    /// Remove constraint by ID.
    #[allow(dead_code)]
    pub fn remove(&mut self, id: u32) -> bool {
        let before = self.constraints.len();
        self.constraints.retain(|c| c.id != id);
        self.constraints.len() < before
    }
    /// Get a constraint by ID.
    #[allow(dead_code)]
    pub fn get(&self, id: u32) -> Option<&BodyConstraint> {
        self.constraints.iter().find(|c| c.id == id)
    }
    /// Get all active constraints.
    #[allow(dead_code)]
    pub fn active_constraints(&self) -> Vec<&BodyConstraint> {
        self.constraints.iter().filter(|c| c.active).collect()
    }
    /// Number of constraints.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.constraints.len()
    }
    /// Whether there are no constraints.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.constraints.is_empty()
    }
    /// Apply active constraints to a body (modifies velocity/angular_velocity).
    #[allow(dead_code)]
    pub fn apply_to(&self, body: &mut RigidBody) {
        for c in self.constraints.iter().filter(|c| c.active) {
            match c.constraint_type {
                ConstraintType::Fixed => {
                    body.velocity = Vec3::zeros();
                    body.angular_velocity = Vec3::zeros();
                }
                ConstraintType::PositionLock => {
                    body.velocity = Vec3::zeros();
                }
                ConstraintType::RotationLock => {
                    body.angular_velocity = Vec3::zeros();
                }
                ConstraintType::AxisLock(axis) => {
                    let ax = axis.min(2) as usize;
                    body.velocity[ax] = 0.0;
                }
                ConstraintType::AngularAxisLock(axis) => {
                    let ax = axis.min(2) as usize;
                    body.angular_velocity[ax] = 0.0;
                }
                ConstraintType::Custom => {}
            }
        }
    }
}
/// Archimedes buoyancy parameters.
#[derive(Debug, Clone, Copy)]
pub struct BuoyancyParams {
    /// Fluid density (kg/m³), e.g. 1000 for water.
    pub fluid_density: Real,
    /// Gravitational acceleration magnitude (m/s²).
    pub gravity: Real,
    /// Submerged volume fraction in \[0, 1\] (estimated or exact).
    pub submerged_fraction: Real,
}
impl BuoyancyParams {
    /// Standard water buoyancy.
    pub fn water() -> Self {
        Self {
            fluid_density: 1000.0,
            gravity: 9.81,
            submerged_fraction: 1.0,
        }
    }
    /// Compute buoyant force (upward, +Y) given object volume.
    pub fn buoyant_force(&self, volume: Real) -> [Real; 3] {
        let f = self.fluid_density * self.gravity * volume * self.submerged_fraction;
        [0.0, f, 0.0]
    }
    /// Apply buoyancy to a rigid body given its volume.
    pub fn apply_to_body(&self, body: &mut RigidBody, volume: Real) {
        use oxiphysics_core::math::Vec3;
        let [_, fy, _] = self.buoyant_force(volume);
        body.apply_force(Vec3::new(0.0, fy, 0.0));
    }
    /// Net vertical force (buoyancy − weight).
    pub fn net_vertical_force(&self, body_mass: Real, volume: Real) -> Real {
        let fb = self.fluid_density * self.gravity * volume * self.submerged_fraction;
        let fg = body_mass * self.gravity;
        fb - fg
    }
}
/// Aerodynamic drag and wind force model.
#[derive(Debug, Clone, Copy)]
pub struct AerodynamicParams {
    /// Air density (kg/m³), standard is 1.225 at sea level.
    pub air_density: Real,
    /// Drag coefficient (dimensionless).
    pub drag_coefficient: Real,
    /// Reference cross-sectional area (m²).
    pub reference_area: Real,
    /// Ambient wind velocity in world space.
    pub wind_velocity: [Real; 3],
}
impl AerodynamicParams {
    /// Standard air parameters with no wind.
    pub fn standard() -> Self {
        Self {
            air_density: 1.225,
            drag_coefficient: 0.47,
            reference_area: 1.0,
            wind_velocity: [0.0; 3],
        }
    }
    /// Compute aerodynamic drag force on a body with the given velocity.
    /// F_drag = -0.5 * rho * Cd * A * |v_rel|^2 * v_rel_hat
    pub fn drag_force(&self, body_velocity: [Real; 3]) -> [Real; 3] {
        let vrel = [
            body_velocity[0] - self.wind_velocity[0],
            body_velocity[1] - self.wind_velocity[1],
            body_velocity[2] - self.wind_velocity[2],
        ];
        let speed_sq = vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2];
        let speed = speed_sq.sqrt();
        if speed < 1e-12 {
            return [0.0; 3];
        }
        let coeff = -0.5 * self.air_density * self.drag_coefficient * self.reference_area * speed;
        [coeff * vrel[0], coeff * vrel[1], coeff * vrel[2]]
    }
    /// Apply drag to a rigid body.
    pub fn apply_drag(&self, body: &mut RigidBody) {
        let v = [body.velocity.x, body.velocity.y, body.velocity.z];
        let [fx, fy, fz] = self.drag_force(v);
        body.apply_force(Vec3::new(fx, fy, fz));
    }
    /// Wind dynamic pressure: q = 0.5 * rho * v_wind^2.
    pub fn dynamic_pressure(&self) -> Real {
        let speed_sq = self.wind_velocity.iter().map(|&v| v * v).sum::<Real>();
        0.5 * self.air_density * speed_sq
    }
}
/// Force accumulation manager for a rigid body.
///
/// Tracks named persistent forces and one-shot impulses applied in a
/// single step, with per-step history for debugging.
#[derive(Debug, Clone, Default)]
pub struct ForceAccumulator {
    /// Named persistent forces re-applied every step.
    pub persistent: Vec<PersistentForce>,
    /// Accumulated one-shot force for the current step.
    pub current_force: Vec3,
    /// Accumulated one-shot torque for the current step.
    pub current_torque: Vec3,
    /// Peak force magnitude observed across all steps.
    pub peak_force_magnitude: Real,
    /// Total number of force application calls since creation.
    pub total_force_applications: u64,
}
impl ForceAccumulator {
    /// Create a new empty accumulator.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a named persistent force.
    #[allow(dead_code)]
    pub fn add_persistent(&mut self, force: PersistentForce) {
        self.persistent.push(force);
    }
    /// Remove a persistent force by label.
    #[allow(dead_code)]
    pub fn remove_persistent(&mut self, label: &str) {
        self.persistent.retain(|f| f.label != label);
    }
    /// Apply a one-shot force for the current step.
    #[allow(dead_code)]
    pub fn apply_force(&mut self, force: Vec3) {
        self.current_force += force;
        self.total_force_applications += 1;
        let mag = force.norm();
        if mag > self.peak_force_magnitude {
            self.peak_force_magnitude = mag;
        }
    }
    /// Apply a one-shot torque for the current step.
    #[allow(dead_code)]
    pub fn apply_torque(&mut self, torque: Vec3) {
        self.current_torque += torque;
    }
    /// Collect all active forces (persistent + current) into the body.
    /// Clears one-shot force/torque after accumulation.
    #[allow(dead_code)]
    pub fn flush(&mut self, body: &mut RigidBody) {
        for pf in &self.persistent {
            if !pf.enabled {
                continue;
            }
            if let Some(pt) = pf.point {
                body.apply_force_at_point(pf.force, pt);
            } else {
                body.apply_force(pf.force);
            }
        }
        body.apply_force(self.current_force);
        body.apply_torque(self.current_torque);
        self.current_force = Vec3::zeros();
        self.current_torque = Vec3::zeros();
    }
    /// Number of persistent forces.
    #[allow(dead_code)]
    pub fn persistent_count(&self) -> usize {
        self.persistent.len()
    }
    /// Number of enabled persistent forces.
    #[allow(dead_code)]
    pub fn enabled_count(&self) -> usize {
        self.persistent.iter().filter(|f| f.enabled).count()
    }
}
/// Rolling history of recent contacts for a body.
#[derive(Debug, Clone)]
pub struct ContactHistory {
    pub(super) entries: Vec<ContactHistoryEntry>,
    pub(super) max_entries: usize,
}
impl ContactHistory {
    /// Create a new contact history with the given capacity.
    #[allow(dead_code)]
    pub fn new(max_entries: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_entries: max_entries.max(1),
        }
    }
    /// Record a new contact event.
    #[allow(dead_code)]
    pub fn record(
        &mut self,
        time: Real,
        other_handle: u32,
        point: Vec3,
        normal: Vec3,
        impulse: Real,
    ) {
        if self.entries.len() >= self.max_entries {
            self.entries.remove(0);
        }
        self.entries.push(ContactHistoryEntry {
            time,
            other_handle,
            point,
            normal,
            impulse,
        });
    }
    /// Most recent contact, if any.
    #[allow(dead_code)]
    pub fn latest(&self) -> Option<&ContactHistoryEntry> {
        self.entries.last()
    }
    /// All contacts with a given body handle.
    #[allow(dead_code)]
    pub fn contacts_with(&self, handle: u32) -> Vec<&ContactHistoryEntry> {
        self.entries
            .iter()
            .filter(|e| e.other_handle == handle)
            .collect()
    }
    /// Peak impulse magnitude over all recorded contacts.
    #[allow(dead_code)]
    pub fn peak_impulse(&self) -> Real {
        self.entries.iter().map(|e| e.impulse).fold(0.0, f64::max)
    }
    /// Total number of recorded contacts.
    #[allow(dead_code)]
    pub fn count(&self) -> usize {
        self.entries.len()
    }
    /// Whether the history is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear all history.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Average impulse across all recorded contacts.
    #[allow(dead_code)]
    pub fn average_impulse(&self) -> Real {
        if self.entries.is_empty() {
            return 0.0;
        }
        let sum: Real = self.entries.iter().map(|e| e.impulse).sum();
        sum / self.entries.len() as Real
    }
}
/// Magnus effect force for spinning bodies (balls, bullets, etc.).
///
/// F_magnus = rho * V_ball * (omega × v)
/// where V_ball is the volume of the ball.
#[derive(Debug, Clone, Copy)]
pub struct MagnusEffect {
    /// Fluid density (kg/m³).
    pub fluid_density: Real,
    /// Volume of the spinning body (m³).
    pub volume: Real,
    /// Magnus coefficient (dimensionless scaling factor).
    pub coefficient: Real,
}
impl MagnusEffect {
    /// Create with default coefficient 1.0.
    pub fn new(fluid_density: Real, volume: Real) -> Self {
        Self {
            fluid_density,
            volume,
            coefficient: 1.0,
        }
    }
    /// Sphere magnus effect.
    pub fn sphere(radius: Real, fluid_density: Real) -> Self {
        let vol = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        Self::new(fluid_density, vol)
    }
    /// Compute Magnus force: F = coefficient * rho * V * (omega × v).
    pub fn force(&self, angular_velocity: [Real; 3], linear_velocity: [Real; 3]) -> [Real; 3] {
        let [ox, oy, oz] = angular_velocity;
        let [vx, vy, vz] = linear_velocity;
        let cx = oy * vz - oz * vy;
        let cy = oz * vx - ox * vz;
        let cz = ox * vy - oy * vx;
        let scale = self.coefficient * self.fluid_density * self.volume;
        [scale * cx, scale * cy, scale * cz]
    }
    /// Apply Magnus force to a rigid body.
    pub fn apply_to_body(&self, body: &mut RigidBody) {
        let omega = [
            body.angular_velocity.x,
            body.angular_velocity.y,
            body.angular_velocity.z,
        ];
        let vel = [body.velocity.x, body.velocity.y, body.velocity.z];
        let [fx, fy, fz] = self.force(omega, vel);
        body.apply_force(Vec3::new(fx, fy, fz));
    }
}
/// Type of constraint applied to a body.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstraintType {
    /// Lock all degrees of freedom (position + rotation).
    Fixed,
    /// Lock position, allow rotation.
    PositionLock,
    /// Lock rotation, allow translation.
    RotationLock,
    /// Lock a single axis of translation.
    AxisLock(u8),
    /// Lock a single axis of rotation.
    AngularAxisLock(u8),
    /// Custom user-defined constraint label.
    Custom,
}
/// Activity state of a rigid body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyState {
    /// Body is actively simulated.
    Active,
    /// Body is sleeping (not simulated until disturbed).
    Sleeping,
}
/// A primitive shape with local offset and density for compound body building.
#[derive(Debug, Clone)]
pub struct CompoundPrimitive {
    /// Local position offset.
    pub local_position: [Real; 3],
    /// Local rotation (quaternion \[x, y, z, w\]).
    pub local_rotation: [Real; 4],
    /// Mass of this primitive.
    pub mass: Real,
    /// Local inertia diagonal (principal axes assumed aligned).
    pub local_inertia_diag: [Real; 3],
    /// Human-readable tag.
    pub tag: &'static str,
}
impl CompoundPrimitive {
    /// Box primitive.
    pub fn new_box(local_position: [Real; 3], half_extents: [Real; 3], density: Real) -> Self {
        let vol = 8.0 * half_extents[0] * half_extents[1] * half_extents[2];
        let mass = density * vol;
        let [hx, hy, hz] = half_extents;
        let ixx = mass / 3.0 * (hy * hy + hz * hz);
        let iyy = mass / 3.0 * (hx * hx + hz * hz);
        let izz = mass / 3.0 * (hx * hx + hy * hy);
        Self {
            local_position,
            local_rotation: [0.0, 0.0, 0.0, 1.0],
            mass,
            local_inertia_diag: [ixx, iyy, izz],
            tag: "box",
        }
    }
    /// Sphere primitive.
    pub fn new_sphere(local_position: [Real; 3], radius: Real, density: Real) -> Self {
        let mass = density * (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        let i = 0.4 * mass * radius * radius;
        Self {
            local_position,
            local_rotation: [0.0, 0.0, 0.0, 1.0],
            mass,
            local_inertia_diag: [i, i, i],
            tag: "sphere",
        }
    }
    /// Cylinder primitive (axis along Y).
    pub fn new_cylinder(
        local_position: [Real; 3],
        radius: Real,
        height: Real,
        density: Real,
    ) -> Self {
        let mass = density * std::f64::consts::PI * radius * radius * height;
        let ixx = mass * (3.0 * radius * radius + height * height) / 12.0;
        let iyy = 0.5 * mass * radius * radius;
        Self {
            local_position,
            local_rotation: [0.0, 0.0, 0.0, 1.0],
            mass,
            local_inertia_diag: [ixx, iyy, ixx],
            tag: "cylinder",
        }
    }
}
/// Type of rigid body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyType {
    /// Fully simulated body affected by forces and collisions.
    Dynamic,
    /// Body moved only by user code, not affected by forces.
    Kinematic,
    /// Immovable body (infinite mass).
    Static,
}
/// A rigid body with variable mass, modelling propellant consumption.
///
/// Implements the Meshchersky equation:
///   m(t) * dv/dt = F_ext + v_rel * dm/dt
/// where v_rel is the exhaust velocity relative to the body.
#[derive(Debug, Clone)]
pub struct VariableMassBody {
    /// Underlying rigid body (mass is updated each step).
    pub body: RigidBody,
    /// Initial propellant mass.
    pub propellant_mass: Real,
    /// Propellant already consumed.
    pub consumed_mass: Real,
    /// Specific impulse (Isp in seconds).
    pub isp: Real,
    /// Current thrust level \[0, 1\].
    pub throttle: Real,
    /// Maximum thrust force (N).
    pub max_thrust: Real,
    /// Thrust direction in body local frame.
    pub thrust_direction: [Real; 3],
}
impl VariableMassBody {
    /// Create a new variable-mass body.
    pub fn new(
        dry_mass: Real,
        propellant_mass: Real,
        isp: Real,
        max_thrust: Real,
        thrust_direction: [Real; 3],
    ) -> Self {
        let total_mass = dry_mass + propellant_mass;
        Self {
            body: RigidBody::new(total_mass),
            propellant_mass,
            consumed_mass: 0.0,
            isp,
            throttle: 0.0,
            max_thrust,
            thrust_direction,
        }
    }
    /// Remaining propellant mass.
    pub fn remaining_propellant(&self) -> Real {
        (self.propellant_mass - self.consumed_mass).max(0.0)
    }
    /// Current total mass.
    pub fn current_mass(&self) -> Real {
        self.body.mass
    }
    /// Mass flow rate (kg/s) at current throttle.
    pub fn mass_flow_rate(&self) -> Real {
        let g0 = 9.80665;
        if self.isp > 0.0 {
            self.max_thrust * self.throttle / (self.isp * g0)
        } else {
            0.0
        }
    }
    /// Integrate one step: consume propellant, apply thrust, update mass.
    pub fn step(&mut self, dt: Real, gravity: &[Real; 3]) {
        if self.remaining_propellant() > 0.0 && self.throttle > 0.0 {
            let mdot = self.mass_flow_rate();
            let dm = mdot * dt;
            let actual_dm = dm.min(self.remaining_propellant());
            self.consumed_mass += actual_dm;
            let new_mass = (self.body.mass - actual_dm).max(1e-6);
            self.body.mass = new_mass;
            self.body.inverse_mass = 1.0 / new_mass;
            let td = self.thrust_direction;
            let thrust_mag = self.max_thrust * self.throttle;
            let thrust = Vec3::new(td[0] * thrust_mag, td[1] * thrust_mag, td[2] * thrust_mag);
            self.body.apply_force(thrust);
        }
        let grav = Vec3::new(gravity[0], gravity[1], gravity[2]);
        self.body.integrate_forces(dt, &grav);
        self.body.integrate_velocity(dt);
    }
    /// Tsiolkovsky delta-v: Δv = Isp * g0 * ln(m0 / mf).
    pub fn tsiolkovsky_delta_v(&self, dry_mass: Real) -> Real {
        let g0 = 9.80665;
        let m0 = dry_mass + self.propellant_mass;
        let mf = dry_mass;
        if mf > 0.0 && m0 > mf {
            self.isp * g0 * (m0 / mf).ln()
        } else {
            0.0
        }
    }
}
