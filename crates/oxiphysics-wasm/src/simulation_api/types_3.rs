//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use super::types::{
    WasmColliderShape, WasmContactEvent, WasmJointHandle, WasmRaycastResult, WasmRigidBody,
};
use super::types_4::{
    WasmBodyState, WasmBodyType, WasmCollider, WasmColliderHandle, WasmJoint, WasmOverlapResult,
    WasmRigidBodyHandle, WasmSimulationConfig,
};
use crate::wasm_helpers::to_js_value;

/// A complete WASM physics simulation world.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct WasmWorld {
    /// Simulation configuration.
    #[wasm_bindgen(skip)]
    pub config: WasmSimulationConfig,
    /// Rigid bodies stored by handle.
    pub(super) bodies: HashMap<u32, WasmRigidBody>,
    /// Colliders stored by handle.
    pub(super) colliders: HashMap<u32, WasmCollider>,
    /// Joints stored by handle.
    joints: HashMap<u32, WasmJoint>,
    /// Next body handle ID.
    next_body_id: u32,
    /// Next collider handle ID.
    next_collider_id: u32,
    /// Next joint handle ID.
    next_joint_id: u32,
    /// Pending physics events.
    events: Vec<WasmPhysicsEvent>,
    /// Contact events from last step.
    contact_events: Vec<WasmContactEvent>,
    /// Simulation time (s).
    pub time: f64,
    /// Step count.
    pub step_count: u64,
}
impl WasmWorld {
    /// Create a new world with the given configuration.
    pub fn new(config: WasmSimulationConfig) -> Self {
        Self {
            config,
            bodies: HashMap::new(),
            colliders: HashMap::new(),
            joints: HashMap::new(),
            next_body_id: 0,
            next_collider_id: 0,
            next_joint_id: 0,
            events: Vec::new(),
            contact_events: Vec::new(),
            time: 0.0,
            step_count: 0,
        }
    }
    /// Create a world with Earth gravity.
    pub fn earth() -> Self {
        Self::new(WasmSimulationConfig::earth())
    }
    /// Add a dynamic rigid body. Returns its handle.
    pub fn add_dynamic_body(
        &mut self,
        mass: f64,
        px: f64,
        py: f64,
        pz: f64,
    ) -> WasmRigidBodyHandle {
        let handle = WasmRigidBodyHandle(self.next_body_id);
        self.next_body_id += 1;
        let body = WasmRigidBody::new_dynamic(handle, [px, py, pz], mass);
        self.bodies.insert(handle.0, body);
        handle
    }
    /// Add a static rigid body. Returns its handle.
    pub fn add_static_body(&mut self, px: f64, py: f64, pz: f64) -> WasmRigidBodyHandle {
        let handle = WasmRigidBodyHandle(self.next_body_id);
        self.next_body_id += 1;
        let body = WasmRigidBody::new_static(handle, [px, py, pz]);
        self.bodies.insert(handle.0, body);
        handle
    }
    /// Remove a rigid body and its attached colliders.
    pub fn remove_body(&mut self, handle: WasmRigidBodyHandle) -> bool {
        let removed = self.bodies.remove(&handle.0).is_some();
        if removed {
            self.colliders.retain(|_, c| c.body != Some(handle));
        }
        removed
    }
    /// Get the position of a body.
    pub fn get_position(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 3]> {
        self.bodies.get(&handle.0).map(|b| b.position)
    }
    /// Set the position of a body.
    pub fn set_position(&mut self, handle: WasmRigidBodyHandle, pos: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.position = pos;
            true
        } else {
            false
        }
    }
    /// Get the velocity of a body.
    pub fn get_velocity(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 3]> {
        self.bodies.get(&handle.0).map(|b| b.linear_vel)
    }
    /// Set the velocity of a body.
    pub fn set_velocity(&mut self, handle: WasmRigidBodyHandle, vel: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.linear_vel = vel;
            true
        } else {
            false
        }
    }
    /// Apply a force to a body at its center of mass.
    pub fn apply_force(&mut self, handle: WasmRigidBodyHandle, force: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.force_accum[0] += force[0];
            body.force_accum[1] += force[1];
            body.force_accum[2] += force[2];
            true
        } else {
            false
        }
    }
    /// Apply an impulse to a body at its center of mass.
    pub fn apply_impulse(&mut self, handle: WasmRigidBodyHandle, impulse: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            if body.inv_mass > 1e-15 {
                body.linear_vel[0] += impulse[0] * body.inv_mass;
                body.linear_vel[1] += impulse[1] * body.inv_mass;
                body.linear_vel[2] += impulse[2] * body.inv_mass;
            }
            true
        } else {
            false
        }
    }
    /// Apply torque to a body.
    pub fn apply_torque(&mut self, handle: WasmRigidBodyHandle, torque: [f64; 3]) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.torque_accum[0] += torque[0];
            body.torque_accum[1] += torque[1];
            body.torque_accum[2] += torque[2];
            true
        } else {
            false
        }
    }
    /// Get the rotation quaternion of a body (x,y,z,w).
    pub fn get_rotation(&self, handle: WasmRigidBodyHandle) -> Option<[f64; 4]> {
        self.bodies.get(&handle.0).map(|b| b.rotation)
    }
    /// Set the gravity scale of a body.
    pub fn set_gravity_scale(&mut self, handle: WasmRigidBodyHandle, scale: f64) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            body.gravity_scale = scale;
            true
        } else {
            false
        }
    }
    /// Wake a sleeping body.
    pub fn wake_body(&mut self, handle: WasmRigidBodyHandle) -> bool {
        if let Some(body) = self.bodies.get_mut(&handle.0) {
            if body.sleeping {
                body.sleeping = false;
                body.sleep_time = 0.0;
                self.events
                    .push(WasmPhysicsEvent::BodyWoke { body: handle });
            }
            true
        } else {
            false
        }
    }
    /// Get body state.
    pub fn get_body_state(&self, handle: WasmRigidBodyHandle) -> Option<WasmBodyState> {
        self.bodies.get(&handle.0).map(|b| b.state())
    }
    /// Number of rigid bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Attach a sphere collider to a body.
    pub fn add_sphere_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        radius: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(handle, WasmColliderShape::Sphere { radius });
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }
    /// Attach a box collider to a body.
    pub fn add_box_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        hx: f64,
        hy: f64,
        hz: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(
            handle,
            WasmColliderShape::Box {
                half_extents: [hx, hy, hz],
            },
        );
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }
    /// Attach a plane collider to a body.
    pub fn add_plane_collider(
        &mut self,
        body: WasmRigidBodyHandle,
        nx: f64,
        ny: f64,
        nz: f64,
        d: f64,
    ) -> WasmColliderHandle {
        let handle = WasmColliderHandle(self.next_collider_id);
        self.next_collider_id += 1;
        let mut coll = WasmCollider::new(
            handle,
            WasmColliderShape::Plane {
                normal: [nx, ny, nz],
                offset: d,
            },
        );
        coll.body = Some(body);
        self.colliders.insert(handle.0, coll);
        handle
    }
    /// Set collider friction.
    pub fn set_friction(&mut self, handle: WasmColliderHandle, friction: f64) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.friction = friction.clamp(0.0, f64::INFINITY);
            true
        } else {
            false
        }
    }
    /// Set collider restitution (bounciness).
    pub fn set_restitution(&mut self, handle: WasmColliderHandle, restitution: f64) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.restitution = restitution.clamp(0.0, 1.0);
            true
        } else {
            false
        }
    }
    /// Make a collider a sensor (no collision response).
    pub fn set_sensor(&mut self, handle: WasmColliderHandle, is_sensor: bool) -> bool {
        if let Some(coll) = self.colliders.get_mut(&handle.0) {
            coll.is_sensor = is_sensor;
            true
        } else {
            false
        }
    }
    /// Number of colliders.
    pub fn collider_count(&self) -> usize {
        self.colliders.len()
    }
    /// Add a ball joint between two bodies.
    pub fn add_ball_joint(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> WasmJointHandle {
        let handle = WasmJointHandle(self.next_joint_id);
        self.next_joint_id += 1;
        let joint = WasmJoint::ball(handle, body_a, body_b, anchor_a, anchor_b);
        self.joints.insert(handle.0, joint);
        handle
    }
    /// Add a distance constraint.
    pub fn add_distance_joint(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        target_dist: f64,
    ) -> WasmJointHandle {
        let handle = WasmJointHandle(self.next_joint_id);
        self.next_joint_id += 1;
        let joint = WasmJoint::distance(handle, body_a, body_b, target_dist);
        self.joints.insert(handle.0, joint);
        handle
    }
    /// Remove a joint.
    pub fn remove_joint(&mut self, handle: WasmJointHandle) -> bool {
        self.joints.remove(&handle.0).is_some()
    }
    /// Number of joints.
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }
    /// Advance the simulation by the configured dt.
    pub fn step(&mut self) {
        let dt = self.config.dt;
        let gravity = self.config.gravity;
        self.events.clear();
        self.contact_events.clear();
        let handles: Vec<u32> = self.bodies.keys().copied().collect();
        for id in handles {
            if let Some(body) = self.bodies.get_mut(&id) {
                body.integrate(dt, gravity);
            }
        }
        self.detect_sphere_plane_collisions(dt);
        self.update_sleep(dt);
        self.time += dt;
        self.step_count += 1;
    }
    /// Advance the simulation by a custom delta time.
    pub fn step_with_dt(&mut self, dt: f64) {
        let old_dt = self.config.dt;
        self.config.dt = dt;
        self.step();
        self.config.dt = old_dt;
    }
    fn detect_sphere_plane_collisions(&mut self, dt: f64) {
        let sphere_data: Vec<(WasmColliderHandle, Option<WasmRigidBodyHandle>, f64)> = self
            .colliders
            .values()
            .filter_map(|c| {
                if let WasmColliderShape::Sphere { radius } = c.shape {
                    Some((c.handle, c.body, radius))
                } else {
                    None
                }
            })
            .collect();
        let plane_data: Vec<(WasmColliderHandle, [f64; 3], f64)> = self
            .colliders
            .values()
            .filter_map(|c| {
                if let WasmColliderShape::Plane { normal, offset } = c.shape {
                    Some((c.handle, normal, offset))
                } else {
                    None
                }
            })
            .collect();
        for (sc_handle, body_opt, radius) in &sphere_data {
            let Some(body_handle) = body_opt else {
                continue;
            };
            let body_pos = match self.bodies.get(&body_handle.0) {
                Some(b) => b.position,
                None => continue,
            };
            for (pc_handle, plane_normal, plane_offset) in &plane_data {
                let dist = body_pos[0] * plane_normal[0]
                    + body_pos[1] * plane_normal[1]
                    + body_pos[2] * plane_normal[2]
                    - plane_offset;
                let penetration = radius - dist;
                if penetration > 0.0 {
                    if let Some(body) = self.bodies.get_mut(&body_handle.0) {
                        body.position[0] += penetration * plane_normal[0];
                        body.position[1] += penetration * plane_normal[1];
                        body.position[2] += penetration * plane_normal[2];
                        let restitution = self
                            .colliders
                            .get(&sc_handle.0)
                            .map(|c| c.restitution)
                            .unwrap_or(0.3);
                        let vn = body.linear_vel[0] * plane_normal[0]
                            + body.linear_vel[1] * plane_normal[1]
                            + body.linear_vel[2] * plane_normal[2];
                        if vn < 0.0 {
                            let impulse = -(1.0 + restitution) * vn;
                            body.linear_vel[0] += impulse * plane_normal[0];
                            body.linear_vel[1] += impulse * plane_normal[1];
                            body.linear_vel[2] += impulse * plane_normal[2];
                        }
                        let friction = self
                            .colliders
                            .get(&sc_handle.0)
                            .map(|c| c.friction)
                            .unwrap_or(0.5);
                        body.linear_vel[0] *= (1.0 - friction * dt).max(0.0);
                        body.linear_vel[2] *= (1.0 - friction * dt).max(0.0);
                    }
                    self.contact_events.push(WasmContactEvent {
                        collider_a: *sc_handle,
                        collider_b: *pc_handle,
                        normal: *plane_normal,
                        contact_point: [
                            body_pos[0] - radius * plane_normal[0],
                            body_pos[1] - radius * plane_normal[1],
                            body_pos[2] - radius * plane_normal[2],
                        ],
                        depth: penetration,
                        impulse: penetration,
                        friction_impulse: [0.0; 3],
                    });
                }
            }
        }
    }
    fn update_sleep(&mut self, dt: f64) {
        if !self.config.allow_sleeping {
            return;
        }
        let lin_thresh = self.config.sleep_linear_threshold;
        let ang_thresh = self.config.sleep_angular_threshold;
        let sleep_time = self.config.sleep_time_threshold;
        let handles: Vec<u32> = self.bodies.keys().copied().collect();
        for id in handles {
            if let Some(body) = self.bodies.get_mut(&id) {
                if body.body_type != WasmBodyType::Dynamic {
                    continue;
                }
                let lin_speed2 = body.linear_vel[0] * body.linear_vel[0]
                    + body.linear_vel[1] * body.linear_vel[1]
                    + body.linear_vel[2] * body.linear_vel[2];
                let ang_speed2 = body.angular_vel[0] * body.angular_vel[0]
                    + body.angular_vel[1] * body.angular_vel[1]
                    + body.angular_vel[2] * body.angular_vel[2];
                if lin_speed2 < lin_thresh * lin_thresh && ang_speed2 < ang_thresh * ang_thresh {
                    body.sleep_time += dt;
                    if body.sleep_time >= sleep_time && !body.sleeping {
                        body.sleeping = true;
                        let h = WasmRigidBodyHandle(id);
                        self.events.push(WasmPhysicsEvent::BodySlept { body: h });
                    }
                } else {
                    body.sleep_time = 0.0;
                    if body.sleeping {
                        body.sleeping = false;
                        let h = WasmRigidBodyHandle(id);
                        self.events.push(WasmPhysicsEvent::BodyWoke { body: h });
                    }
                }
            }
        }
    }
    /// Cast a ray and return the first hit.
    pub fn raycast(
        &self,
        origin: [f64; 3],
        direction: [f64; 3],
        max_dist: f64,
    ) -> WasmRaycastResult {
        let dir_len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        if dir_len < 1e-15 {
            return WasmRaycastResult::miss();
        }
        let dir = [
            direction[0] / dir_len,
            direction[1] / dir_len,
            direction[2] / dir_len,
        ];
        let mut best_toi = max_dist;
        let mut best_result = WasmRaycastResult::miss();
        for coll in self.colliders.values() {
            let body_pos = coll
                .body
                .and_then(|h| self.bodies.get(&h.0))
                .map(|b| b.position)
                .unwrap_or([0.0; 3]);
            let hit = match &coll.shape {
                WasmColliderShape::Sphere { radius } => {
                    let oc = [
                        origin[0] - body_pos[0],
                        origin[1] - body_pos[1],
                        origin[2] - body_pos[2],
                    ];
                    let b = oc[0] * dir[0] + oc[1] * dir[1] + oc[2] * dir[2];
                    let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
                    let discriminant = b * b - c;
                    if discriminant < 0.0 {
                        None
                    } else {
                        let t = -b - discriminant.sqrt();
                        if t > 1e-4 && t < best_toi {
                            let hit_pos = [
                                origin[0] + t * dir[0],
                                origin[1] + t * dir[1],
                                origin[2] + t * dir[2],
                            ];
                            let normal = [
                                (hit_pos[0] - body_pos[0]) / radius,
                                (hit_pos[1] - body_pos[1]) / radius,
                                (hit_pos[2] - body_pos[2]) / radius,
                            ];
                            Some((t, hit_pos, normal))
                        } else {
                            None
                        }
                    }
                }
                WasmColliderShape::Plane { normal, offset } => {
                    let denom = dir[0] * normal[0] + dir[1] * normal[1] + dir[2] * normal[2];
                    if denom.abs() < 1e-10 {
                        None
                    } else {
                        let t = (offset
                            - (origin[0] * normal[0]
                                + origin[1] * normal[1]
                                + origin[2] * normal[2]))
                            / denom;
                        if t > 1e-4 && t < best_toi {
                            let hit_pos = [
                                origin[0] + t * dir[0],
                                origin[1] + t * dir[1],
                                origin[2] + t * dir[2],
                            ];
                            Some((t, hit_pos, *normal))
                        } else {
                            None
                        }
                    }
                }
                _ => None,
            };
            if let Some((t, pos, norm)) = hit {
                best_toi = t;
                best_result = WasmRaycastResult {
                    hit: true,
                    collider: Some(coll.handle),
                    body: coll.body,
                    point: pos,
                    normal: norm,
                    toi: t,
                    feature_id: 0,
                };
            }
        }
        best_result
    }
    /// Test all colliders against an AABB and return overlapping ones.
    pub fn aabb_overlap(&self, min: [f64; 3], max: [f64; 3]) -> WasmOverlapResult {
        let mut result = WasmOverlapResult::empty();
        for coll in self.colliders.values() {
            let body_pos = coll
                .body
                .and_then(|h| self.bodies.get(&h.0))
                .map(|b| b.position)
                .unwrap_or([0.0; 3]);
            let r = coll.shape.bounding_radius();
            let aabb_min = [body_pos[0] - r, body_pos[1] - r, body_pos[2] - r];
            let aabb_max = [body_pos[0] + r, body_pos[1] + r, body_pos[2] + r];
            if aabb_min[0] <= max[0]
                && aabb_max[0] >= min[0]
                && aabb_min[1] <= max[1]
                && aabb_max[1] >= min[1]
                && aabb_min[2] <= max[2]
                && aabb_max[2] >= min[2]
            {
                result.colliders.push(coll.handle);
                result.bodies.push(coll.body);
            }
        }
        result
    }
    /// Drain pending physics events.
    pub fn drain_events(&mut self) -> Vec<WasmPhysicsEvent> {
        std::mem::take(&mut self.events)
    }
    /// Drain pending contact events.
    pub fn drain_contact_events(&mut self) -> Vec<WasmContactEvent> {
        std::mem::take(&mut self.contact_events)
    }
    /// Remove all bodies, colliders, and joints.
    pub fn reset(&mut self) {
        self.bodies.clear();
        self.colliders.clear();
        self.joints.clear();
        self.events.clear();
        self.contact_events.clear();
        self.next_body_id = 0;
        self.next_collider_id = 0;
        self.next_joint_id = 0;
        self.time = 0.0;
        self.step_count = 0;
    }
}
#[wasm_bindgen]
impl WasmWorld {
    /// Construct a new world with default Earth-gravity configuration.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> Self {
        Self::new(WasmSimulationConfig::default())
    }
    /// Build a world from an explicit configuration.
    #[wasm_bindgen(js_name = "withConfig")]
    pub fn with_config_js(config: WasmSimulationConfig) -> Self {
        Self::new(config)
    }
    /// Build a world preconfigured for Earth-standard gravity.
    #[wasm_bindgen(js_name = "earth")]
    pub fn earth_js() -> Self {
        Self::earth()
    }
    /// Get a clone of the simulation configuration.
    #[wasm_bindgen(js_name = "getConfig")]
    pub fn get_config_js(&self) -> WasmSimulationConfig {
        self.config.clone()
    }
    /// Replace the simulation configuration.
    #[wasm_bindgen(js_name = "setConfig")]
    pub fn set_config_js(&mut self, config: WasmSimulationConfig) {
        self.config = config;
    }
    /// Add a dynamic rigid body. Returns its handle.
    #[wasm_bindgen(js_name = "addDynamicBody")]
    pub fn add_dynamic_body_js(
        &mut self,
        mass: f64,
        px: f64,
        py: f64,
        pz: f64,
    ) -> WasmRigidBodyHandle {
        self.add_dynamic_body(mass, px, py, pz)
    }
    /// Add a static rigid body. Returns its handle.
    #[wasm_bindgen(js_name = "addStaticBody")]
    pub fn add_static_body_js(&mut self, px: f64, py: f64, pz: f64) -> WasmRigidBodyHandle {
        self.add_static_body(px, py, pz)
    }
    /// Remove a rigid body and any attached colliders.
    #[wasm_bindgen(js_name = "removeBody")]
    pub fn remove_body_js(&mut self, handle: WasmRigidBodyHandle) -> bool {
        self.remove_body(handle)
    }
    /// Returns `true` if the given body exists.
    #[wasm_bindgen(js_name = "hasBody")]
    pub fn has_body_js(&self, handle: WasmRigidBodyHandle) -> bool {
        self.bodies.contains_key(&handle.0)
    }
    /// Get position as `[x, y, z]`. Returns an empty vector if the body is missing.
    #[wasm_bindgen(js_name = "getPosition")]
    pub fn get_position_js(&self, handle: WasmRigidBodyHandle) -> Vec<f64> {
        self.get_position(handle)
            .map(|p| p.to_vec())
            .unwrap_or_default()
    }
    /// Set body position. Returns `true` on success.
    #[wasm_bindgen(js_name = "setPosition")]
    pub fn set_position_js(&mut self, handle: WasmRigidBodyHandle, x: f64, y: f64, z: f64) -> bool {
        self.set_position(handle, [x, y, z])
    }
    /// Get linear velocity as `[x, y, z]`. Returns an empty vector when missing.
    #[wasm_bindgen(js_name = "getVelocity")]
    pub fn get_velocity_js(&self, handle: WasmRigidBodyHandle) -> Vec<f64> {
        self.get_velocity(handle)
            .map(|v| v.to_vec())
            .unwrap_or_default()
    }
    /// Set linear velocity. Returns `true` on success.
    #[wasm_bindgen(js_name = "setVelocity")]
    pub fn set_velocity_js(&mut self, handle: WasmRigidBodyHandle, x: f64, y: f64, z: f64) -> bool {
        self.set_velocity(handle, [x, y, z])
    }
    /// Apply a force at the body's centre of mass.
    #[wasm_bindgen(js_name = "applyForce")]
    pub fn apply_force_js(
        &mut self,
        handle: WasmRigidBodyHandle,
        fx: f64,
        fy: f64,
        fz: f64,
    ) -> bool {
        self.apply_force(handle, [fx, fy, fz])
    }
    /// Apply an impulse at the body's centre of mass.
    #[wasm_bindgen(js_name = "applyImpulse")]
    pub fn apply_impulse_js(
        &mut self,
        handle: WasmRigidBodyHandle,
        ix: f64,
        iy: f64,
        iz: f64,
    ) -> bool {
        self.apply_impulse(handle, [ix, iy, iz])
    }
    /// Apply a torque to the body.
    #[wasm_bindgen(js_name = "applyTorque")]
    pub fn apply_torque_js(
        &mut self,
        handle: WasmRigidBodyHandle,
        tx: f64,
        ty: f64,
        tz: f64,
    ) -> bool {
        self.apply_torque(handle, [tx, ty, tz])
    }
    /// Get rotation quaternion `[x, y, z, w]`. Empty when missing.
    #[wasm_bindgen(js_name = "getRotation")]
    pub fn get_rotation_js(&self, handle: WasmRigidBodyHandle) -> Vec<f64> {
        self.get_rotation(handle)
            .map(|q| q.to_vec())
            .unwrap_or_default()
    }
    /// Set per-body gravity scale. Returns `true` on success.
    #[wasm_bindgen(js_name = "setGravityScale")]
    pub fn set_gravity_scale_js(&mut self, handle: WasmRigidBodyHandle, scale: f64) -> bool {
        self.set_gravity_scale(handle, scale)
    }
    /// Wake a sleeping body.
    #[wasm_bindgen(js_name = "wakeBody")]
    pub fn wake_body_js(&mut self, handle: WasmRigidBodyHandle) -> bool {
        self.wake_body(handle)
    }
    /// Get a body state as a JS object via `serde-wasm-bindgen`.
    ///
    /// Returns `null` (a `JsValue::NULL`) when the handle is unknown.
    #[wasm_bindgen(js_name = "getBodyState")]
    pub fn get_body_state_js(&self, handle: WasmRigidBodyHandle) -> Result<JsValue, JsValue> {
        match self.get_body_state(handle) {
            Some(state) => to_js_value(&state),
            None => Ok(JsValue::NULL),
        }
    }
    /// Number of rigid bodies.
    #[wasm_bindgen(js_name = "bodyCount")]
    pub fn body_count_js(&self) -> usize {
        self.body_count()
    }
    /// Attach a sphere collider to a body.
    #[wasm_bindgen(js_name = "addSphereCollider")]
    pub fn add_sphere_collider_js(
        &mut self,
        body: WasmRigidBodyHandle,
        radius: f64,
    ) -> WasmColliderHandle {
        self.add_sphere_collider(body, radius)
    }
    /// Attach a box collider with the given half-extents.
    #[wasm_bindgen(js_name = "addBoxCollider")]
    pub fn add_box_collider_js(
        &mut self,
        body: WasmRigidBodyHandle,
        hx: f64,
        hy: f64,
        hz: f64,
    ) -> WasmColliderHandle {
        self.add_box_collider(body, hx, hy, hz)
    }
    /// Attach a plane collider (normal + offset).
    #[wasm_bindgen(js_name = "addPlaneCollider")]
    pub fn add_plane_collider_js(
        &mut self,
        body: WasmRigidBodyHandle,
        nx: f64,
        ny: f64,
        nz: f64,
        d: f64,
    ) -> WasmColliderHandle {
        self.add_plane_collider(body, nx, ny, nz, d)
    }
    /// Set the friction coefficient on a collider.
    #[wasm_bindgen(js_name = "setFriction")]
    pub fn set_friction_js(&mut self, handle: WasmColliderHandle, friction: f64) -> bool {
        self.set_friction(handle, friction)
    }
    /// Set the restitution (bounciness) coefficient on a collider.
    #[wasm_bindgen(js_name = "setRestitution")]
    pub fn set_restitution_js(&mut self, handle: WasmColliderHandle, restitution: f64) -> bool {
        self.set_restitution(handle, restitution)
    }
    /// Mark or unmark a collider as a sensor (no collision response).
    #[wasm_bindgen(js_name = "setSensor")]
    pub fn set_sensor_js(&mut self, handle: WasmColliderHandle, is_sensor: bool) -> bool {
        self.set_sensor(handle, is_sensor)
    }
    /// Number of colliders in the world.
    #[wasm_bindgen(js_name = "colliderCount")]
    pub fn collider_count_js(&self) -> usize {
        self.collider_count()
    }
    /// Add a ball joint between two bodies.
    #[wasm_bindgen(js_name = "addBallJoint")]
    pub fn add_ball_joint_js(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        ax: f64,
        ay: f64,
        az: f64,
        bx: f64,
        by: f64,
        bz: f64,
    ) -> WasmJointHandle {
        self.add_ball_joint(body_a, body_b, [ax, ay, az], [bx, by, bz])
    }
    /// Add a distance constraint between two bodies.
    #[wasm_bindgen(js_name = "addDistanceJoint")]
    pub fn add_distance_joint_js(
        &mut self,
        body_a: WasmRigidBodyHandle,
        body_b: WasmRigidBodyHandle,
        target_dist: f64,
    ) -> WasmJointHandle {
        self.add_distance_joint(body_a, body_b, target_dist)
    }
    /// Remove a joint.
    #[wasm_bindgen(js_name = "removeJoint")]
    pub fn remove_joint_js(&mut self, handle: WasmJointHandle) -> bool {
        self.remove_joint(handle)
    }
    /// Number of joints.
    #[wasm_bindgen(js_name = "jointCount")]
    pub fn joint_count_js(&self) -> usize {
        self.joint_count()
    }
    /// Advance the simulation by the configured `dt`.
    #[wasm_bindgen(js_name = "step")]
    pub fn step_js(&mut self) {
        self.step();
    }
    /// Advance the simulation by a custom delta time.
    #[wasm_bindgen(js_name = "stepWithDt")]
    pub fn step_with_dt_js(&mut self, dt: f64) {
        self.step_with_dt(dt);
    }
    /// Cast a ray and return the first hit as a JS object.
    #[wasm_bindgen(js_name = "raycast")]
    pub fn raycast_js(
        &self,
        ox: f64,
        oy: f64,
        oz: f64,
        dx: f64,
        dy: f64,
        dz: f64,
        max_dist: f64,
    ) -> Result<JsValue, JsValue> {
        let hit = self.raycast([ox, oy, oz], [dx, dy, dz], max_dist);
        to_js_value(&hit)
    }
    /// Test all colliders against an AABB and return overlapping ones.
    #[wasm_bindgen(js_name = "aabbOverlap")]
    pub fn aabb_overlap_js(
        &self,
        min_x: f64,
        min_y: f64,
        min_z: f64,
        max_x: f64,
        max_y: f64,
        max_z: f64,
    ) -> Result<JsValue, JsValue> {
        let result = self.aabb_overlap([min_x, min_y, min_z], [max_x, max_y, max_z]);
        to_js_value(&result)
    }
    /// Drain pending physics events as a JS array.
    #[wasm_bindgen(js_name = "drainEvents")]
    pub fn drain_events_js(&mut self) -> Result<JsValue, JsValue> {
        let events = self.drain_events();
        to_js_value(&events)
    }
    /// Drain pending contact events as a JS array.
    #[wasm_bindgen(js_name = "drainContactEvents")]
    pub fn drain_contact_events_js(&mut self) -> Result<JsValue, JsValue> {
        let events = self.drain_contact_events();
        to_js_value(&events)
    }
    /// Remove all bodies, colliders, joints, and reset time.
    #[wasm_bindgen(js_name = "reset")]
    pub fn reset_js(&mut self) {
        self.reset();
    }
    /// Current simulated time in seconds.
    #[wasm_bindgen(js_name = "currentTime")]
    pub fn current_time_js(&self) -> f64 {
        self.time
    }
    /// Number of completed simulation steps.
    #[wasm_bindgen(js_name = "currentStepCount")]
    pub fn current_step_count_js(&self) -> u64 {
        self.step_count
    }
}
/// A physics event type.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum WasmPhysicsEvent {
    /// Two colliders started touching.
    CollisionStarted {
        collider_a: WasmColliderHandle,
        collider_b: WasmColliderHandle,
    },
    /// Two colliders stopped touching.
    CollisionEnded {
        collider_a: WasmColliderHandle,
        collider_b: WasmColliderHandle,
    },
    /// A sensor collider was entered.
    SensorEntered {
        sensor: WasmColliderHandle,
        other: WasmColliderHandle,
    },
    /// A sensor collider was exited.
    SensorExited {
        sensor: WasmColliderHandle,
        other: WasmColliderHandle,
    },
    /// A body went to sleep.
    BodySlept { body: WasmRigidBodyHandle },
    /// A body woke up.
    BodyWoke { body: WasmRigidBodyHandle },
}
/// Serialized state of a single body for JSON output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializedBody {
    /// Handle ID.
    pub id: u32,
    /// Body type as string.
    pub body_type: String,
    /// Position.
    pub position: [f64; 3],
    /// Rotation (x, y, z, w).
    pub rotation: [f64; 4],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular velocity.
    pub angular_vel: [f64; 3],
    /// Mass.
    pub mass: f64,
    /// Sleeping state.
    pub sleeping: bool,
    /// User data.
    pub user_data: u64,
}
/// Serialized state of a single collider.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializedCollider {
    /// Handle ID.
    pub id: u32,
    /// Associated body ID (if any).
    pub body_id: Option<u32>,
    /// Shape description.
    pub shape_type: String,
    /// Shape parameters (radius, half-extents, etc.).
    pub shape_params: Vec<f64>,
    /// Friction.
    pub friction: f64,
    /// Restitution.
    pub restitution: f64,
    /// Is sensor.
    pub is_sensor: bool,
}
