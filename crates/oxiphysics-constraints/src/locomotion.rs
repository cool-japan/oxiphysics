// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Locomotion controllers and character physics.
//!
//! Provides controllers for walking, jumping, climbing, swimming, procedural IK,
//! balance, vaulting, crowd simulation, and steering behaviours.

use std::f64::consts::{PI, TAU};

// ---------------------------------------------------------------------------
// Vector helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length of a 3-vector.
#[inline]
pub fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalize a 3-vector; returns zero vector if degenerate.
#[inline]
pub fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-12 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// 2-vector dot product (xy plane).
#[inline]
pub fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// 2-vector length.
#[inline]
pub fn len2(v: [f64; 2]) -> f64 {
    dot2(v, v).sqrt()
}

/// 2-vector normalization.
#[inline]
pub fn norm2(v: [f64; 2]) -> [f64; 2] {
    let l = len2(v);
    if l < 1e-12 {
        [0.0; 2]
    } else {
        [v[0] / l, v[1] / l]
    }
}

// ---------------------------------------------------------------------------
// CharacterController
// ---------------------------------------------------------------------------

/// Capsule-based character controller.
///
/// Stores the body dimensions, movement limits, and current state for a
/// physics-driven character.
pub struct CharacterController {
    /// Capsule radius (metres).
    pub radius: f64,
    /// Capsule total height (metres).
    pub height: f64,
    /// Maximum step height the character can auto-step over (metres).
    pub step_height: f64,
    /// Maximum slope angle in radians the character can walk on.
    pub slope_limit: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: f64,
    /// Current position (world-space).
    pub position: [f64; 3],
    /// Current velocity.
    pub velocity: [f64; 3],
    /// True when the character is on the ground.
    pub grounded: bool,
    /// Normal of the ground surface beneath the character.
    pub ground_normal: [f64; 3],
}

impl CharacterController {
    /// Create a new character controller with default parameters.
    pub fn new(radius: f64, height: f64) -> Self {
        Self {
            radius,
            height,
            step_height: 0.3,
            slope_limit: 45.0_f64.to_radians(),
            gravity: -9.81,
            position: [0.0; 3],
            velocity: [0.0; 3],
            grounded: false,
            ground_normal: [0.0, 1.0, 0.0],
        }
    }

    /// Perform a ground check given a height-field query result.
    ///
    /// `ground_y` is the y-coordinate of the ground directly below the
    /// controller.  Returns `true` when contact is detected.
    pub fn check_ground(&mut self, ground_y: f64) -> bool {
        let foot_y = self.position[1] - self.height * 0.5;
        let penetration = ground_y - foot_y;
        if penetration >= -0.01 {
            self.grounded = true;
            if penetration > 0.0 {
                self.position[1] += penetration;
                if self.velocity[1] < 0.0 {
                    self.velocity[1] = 0.0;
                }
            }
        } else {
            self.grounded = false;
        }
        self.grounded
    }

    /// Integrate position by `dt` seconds applying gravity when airborne.
    pub fn integrate(&mut self, dt: f64) {
        if !self.grounded {
            self.velocity[1] += self.gravity * dt;
        }
        self.position = add3(self.position, scale3(self.velocity, dt));
    }

    /// Return the slope angle (radians) relative to the ground normal.
    pub fn slope_angle(&self) -> f64 {
        let up = [0.0, 1.0, 0.0_f64];
        dot3(self.ground_normal, up).clamp(-1.0, 1.0).acos()
    }

    /// Return whether the current ground slope exceeds the slope limit.
    pub fn is_too_steep(&self) -> bool {
        self.slope_angle() > self.slope_limit
    }
}

// ---------------------------------------------------------------------------
// WalkController
// ---------------------------------------------------------------------------

/// Walking locomotion controller.
///
/// Models target velocity, turning, acceleration, and deceleration on slopes.
pub struct WalkController {
    /// Target walk speed (m/s).
    pub target_speed: f64,
    /// Maximum angular turn rate (rad/s).
    pub turn_rate: f64,
    /// Acceleration (m/s²).
    pub acceleration: f64,
    /// Deceleration when stopping (m/s²).
    pub deceleration: f64,
    /// Slope deceleration factor (dimensionless, applied per radian of slope).
    pub slope_decel_factor: f64,
    /// Current facing direction angle (radians, XZ plane).
    pub facing_angle: f64,
    /// Current horizontal speed (m/s).
    pub current_speed: f64,
}

impl WalkController {
    /// Create a walk controller with default parameters.
    pub fn new(target_speed: f64) -> Self {
        Self {
            target_speed,
            turn_rate: PI,
            acceleration: 10.0,
            deceleration: 20.0,
            slope_decel_factor: 0.5,
            facing_angle: 0.0,
            current_speed: 0.0,
        }
    }

    /// Update the controller given `input_dir` (unit XZ vector) and `slope_angle` (radians).
    ///
    /// Returns the new horizontal velocity vector.
    pub fn update(&mut self, input_dir: [f64; 2], slope_angle: f64, dt: f64) -> [f64; 2] {
        let input_len = len2(input_dir);
        let desired_speed = self.target_speed * input_len.min(1.0);

        // Adjust desired speed based on slope
        let effective_desired =
            desired_speed * (1.0 - slope_decel(slope_angle, self.slope_decel_factor));

        // Accelerate / decelerate
        if self.current_speed < effective_desired {
            self.current_speed =
                (self.current_speed + self.acceleration * dt).min(effective_desired);
        } else {
            self.current_speed =
                (self.current_speed - self.deceleration * dt).max(effective_desired);
        }

        // Turn towards input direction
        if input_len > 0.01 {
            let target_angle = input_dir[1].atan2(input_dir[0]);
            let delta = angle_diff(target_angle, self.facing_angle);
            let max_turn = self.turn_rate * dt;
            self.facing_angle += delta.clamp(-max_turn, max_turn);
        }

        let (s, c) = self.facing_angle.sin_cos();
        [c * self.current_speed, s * self.current_speed]
    }
}

/// Compute slope deceleration fraction (0..1).
fn slope_decel(slope_angle: f64, factor: f64) -> f64 {
    (slope_angle * factor).clamp(0.0, 1.0)
}

/// Compute signed shortest angular difference from `from` to `to`.
fn angle_diff(to: f64, from: f64) -> f64 {
    let mut d = to - from;
    while d > PI {
        d -= TAU;
    }
    while d < -PI {
        d += TAU;
    }
    d
}

// ---------------------------------------------------------------------------
// JumpController
// ---------------------------------------------------------------------------

/// Jump controller with coyote time and jump buffering.
pub struct JumpController {
    /// Vertical impulse applied at jump (m/s).
    pub jump_impulse: f64,
    /// Maximum jump height (metres) — used for validation only.
    pub max_height: f64,
    /// Coyote time window (seconds) after leaving a ledge during which a jump is still allowed.
    pub coyote_time: f64,
    /// Jump-buffer window (seconds): how early before landing the player can press jump.
    pub jump_buffer: f64,
    /// Remaining coyote time (seconds).
    pub coyote_timer: f64,
    /// Remaining jump-buffer time (seconds).
    pub buffer_timer: f64,
    /// True when currently airborne from a jump.
    pub is_jumping: bool,
}

impl JumpController {
    /// Create a jump controller.
    pub fn new(jump_impulse: f64) -> Self {
        let max_height = jump_impulse * jump_impulse / (2.0 * 9.81);
        Self {
            jump_impulse,
            max_height,
            coyote_time: 0.1,
            jump_buffer: 0.1,
            coyote_timer: 0.0,
            buffer_timer: 0.0,
            is_jumping: false,
        }
    }

    /// Tick timers.  `grounded` is the current grounded state.
    pub fn tick(&mut self, grounded: bool, dt: f64) {
        if grounded {
            self.coyote_timer = self.coyote_time;
            self.is_jumping = false;
        } else {
            self.coyote_timer = (self.coyote_timer - dt).max(0.0);
        }
        self.buffer_timer = (self.buffer_timer - dt).max(0.0);
    }

    /// Register a jump button press (fills the buffer).
    pub fn press_jump(&mut self) {
        self.buffer_timer = self.jump_buffer;
    }

    /// Attempt to execute a jump.  Returns the vertical impulse to apply, or 0 if not possible.
    pub fn try_jump(&mut self) -> f64 {
        if self.buffer_timer > 0.0 && (self.coyote_timer > 0.0 || !self.is_jumping) {
            self.buffer_timer = 0.0;
            self.coyote_timer = 0.0;
            self.is_jumping = true;
            self.jump_impulse
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// ClimbController
// ---------------------------------------------------------------------------

/// Ladder / wall climbing controller.
pub struct ClimbController {
    /// Speed along the climbing surface (m/s).
    pub climb_speed: f64,
    /// Snap distance for attaching to a ladder/wall (metres).
    pub snap_distance: f64,
    /// Whether the character is currently climbing.
    pub is_climbing: bool,
    /// Surface normal of the climbing surface (points away from surface).
    pub surface_normal: [f64; 3],
    /// Position snapped to the surface.
    pub snapped_position: [f64; 3],
    /// Ceiling height at the current climb position.
    pub ceiling_height: f64,
}

impl ClimbController {
    /// Create a climb controller.
    pub fn new(climb_speed: f64) -> Self {
        Self {
            climb_speed,
            snap_distance: 0.5,
            is_climbing: false,
            surface_normal: [0.0, 0.0, 1.0],
            snapped_position: [0.0; 3],
            ceiling_height: f64::MAX,
        }
    }

    /// Attempt to attach to a surface at `surface_point` with `normal`.
    ///
    /// Returns `true` if snap was successful.
    pub fn try_attach(
        &mut self,
        character_pos: [f64; 3],
        surface_point: [f64; 3],
        normal: [f64; 3],
    ) -> bool {
        let dist = len3(sub3(character_pos, surface_point));
        if dist <= self.snap_distance {
            self.is_climbing = true;
            self.surface_normal = normal;
            // Snap position to surface
            let offset = scale3(normal, self.snap_distance * 0.5);
            self.snapped_position = add3(surface_point, offset);
            true
        } else {
            false
        }
    }

    /// Detach from the climbing surface.
    pub fn detach(&mut self) {
        self.is_climbing = false;
    }

    /// Compute the climb velocity given input (vertical + lateral).
    ///
    /// Returns velocity delta to apply this frame.
    pub fn climb_velocity(&self, vertical_input: f64, lateral_input: f64) -> [f64; 3] {
        if !self.is_climbing {
            return [0.0; 3];
        }
        let up = [0.0, 1.0, 0.0_f64];
        let right = cross3(self.surface_normal, up);
        let up_scaled = scale3(up, vertical_input * self.climb_speed);
        let right_scaled = scale3(right, lateral_input * self.climb_speed);
        add3(up_scaled, right_scaled)
    }

    /// Check if the character has reached the ceiling (cannot climb further up).
    pub fn ceiling_blocked(&self, current_y: f64) -> bool {
        current_y >= self.ceiling_height - 0.1
    }
}

// ---------------------------------------------------------------------------
// SwimController
// ---------------------------------------------------------------------------

/// Underwater swimming controller.
pub struct SwimController {
    /// Buoyancy force (N/kg — upward acceleration equivalent).
    pub buoyancy: f64,
    /// Drag coefficient in water.
    pub drag: f64,
    /// Maximum swim speed (m/s).
    pub swim_speed: f64,
    /// Water surface Y coordinate.
    pub water_surface_y: f64,
    /// Whether the character is currently submerged.
    pub is_submerged: bool,
    /// Current velocity while swimming.
    pub velocity: [f64; 3],
}

impl SwimController {
    /// Create a swim controller.
    pub fn new(swim_speed: f64, water_surface_y: f64) -> Self {
        Self {
            buoyancy: 4.0,
            drag: 2.0,
            swim_speed,
            water_surface_y,
            is_submerged: false,
            velocity: [0.0; 3],
        }
    }

    /// Update submersion state given character Y position.
    pub fn update_submersion(&mut self, character_y: f64) {
        self.is_submerged = character_y < self.water_surface_y;
    }

    /// Apply buoyancy and drag; returns net acceleration vector.
    pub fn apply_forces(&self, gravity: f64) -> [f64; 3] {
        if !self.is_submerged {
            return [0.0, gravity, 0.0];
        }
        let net_vertical = self.buoyancy + gravity; // gravity is negative
        let drag_force = scale3(self.velocity, -self.drag);
        add3([0.0, net_vertical, 0.0], drag_force)
    }

    /// Integrate swimming given direction input and dt.
    pub fn integrate(&mut self, swim_dir: [f64; 3], dt: f64, gravity: f64) {
        if self.is_submerged {
            let desired_vel = scale3(norm3(swim_dir), self.swim_speed);
            // Blend towards desired
            for (vi, dv) in self.velocity.iter_mut().zip(desired_vel.iter()) {
                *vi += (dv - *vi) * (1.0 - (-self.drag * dt).exp());
            }
        } else {
            self.velocity[1] += gravity * dt;
        }
    }

    /// Compute surface transition: returns extra upward impulse when near surface.
    pub fn surface_transition_impulse(&self, character_y: f64) -> f64 {
        let depth = self.water_surface_y - character_y;
        if depth < 0.5 && depth > 0.0 {
            (0.5 - depth) * 2.0 * self.buoyancy
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// ProceduralFootIk
// ---------------------------------------------------------------------------

/// Two-bone IK for procedural foot placement on terrain.
pub struct ProceduralFootIk {
    /// Upper leg length (metres).
    pub upper_length: f64,
    /// Lower leg length (metres).
    pub lower_length: f64,
    /// Current hip position.
    pub hip: [f64; 3],
    /// Current foot target position.
    pub foot_target: [f64; 3],
    /// Actual foot position (IK solution).
    pub foot_pos: [f64; 3],
    /// Knee pole vector (world space).
    pub pole: [f64; 3],
    /// Step timing interval (seconds).
    pub step_interval: f64,
    /// Time since last step.
    pub step_timer: f64,
}

impl ProceduralFootIk {
    /// Create a foot IK solver.
    pub fn new(upper_length: f64, lower_length: f64) -> Self {
        Self {
            upper_length,
            lower_length,
            hip: [0.0; 3],
            foot_target: [0.0; 3],
            foot_pos: [0.0; 3],
            pole: [0.0, 0.0, 1.0],
            step_interval: 0.4,
            step_timer: 0.0,
        }
    }

    /// Solve the two-bone IK chain given `hip` and `target`.
    ///
    /// Returns `(upper_joint, lower_joint)` world positions, or `None` if unreachable.
    pub fn solve(&self, hip: [f64; 3], target: [f64; 3]) -> Option<([f64; 3], [f64; 3])> {
        let to_target = sub3(target, hip);
        let dist = len3(to_target);
        let max_reach = self.upper_length + self.lower_length;

        if dist > max_reach {
            // Over-extension: stretch toward target
            let dir = norm3(to_target);
            let knee = add3(hip, scale3(dir, self.upper_length));
            return Some((knee, target));
        }
        if dist < (self.upper_length - self.lower_length).abs() {
            return None; // too close — singularity
        }

        // Law of cosines for knee angle
        let cos_upper = (dist * dist + self.upper_length * self.upper_length
            - self.lower_length * self.lower_length)
            / (2.0 * dist * self.upper_length);
        let upper_angle = cos_upper.clamp(-1.0, 1.0).acos();

        // Forward direction
        let forward = norm3(to_target);
        // Bend axis using pole vector
        let pole_offset = sub3(self.pole, hip);
        let bend_axis = norm3(cross3(pole_offset, forward));
        let side = cross3(forward, bend_axis);

        let knee = add3(
            hip,
            add3(
                scale3(forward, self.upper_length * upper_angle.cos()),
                scale3(side, self.upper_length * upper_angle.sin()),
            ),
        );
        Some((knee, target))
    }

    /// Update step timer; when elapsed, snap foot target to terrain height.
    pub fn tick_step(&mut self, dt: f64, terrain_y: f64) {
        self.step_timer += dt;
        if self.step_timer >= self.step_interval {
            self.step_timer = 0.0;
            self.foot_target[1] = terrain_y;
        }
    }
}

// ---------------------------------------------------------------------------
// BalanceController
// ---------------------------------------------------------------------------

/// Inverted-pendulum balance controller with ZMP.
pub struct BalanceController {
    /// Height of the centre of mass above the ground (metres).
    pub com_height: f64,
    /// Proportional gain for balance recovery.
    pub kp: f64,
    /// Derivative gain for balance recovery.
    pub kd: f64,
    /// Current centre-of-mass position.
    pub com_pos: [f64; 3],
    /// Current centre-of-mass velocity.
    pub com_vel: [f64; 3],
    /// Support polygon vertices (XZ plane), convex hull.
    pub support_polygon: Vec<[f64; 2]>,
    /// Tipping threshold: ZMP distance outside polygon (metres).
    pub tipping_threshold: f64,
}

impl BalanceController {
    /// Create a balance controller.
    pub fn new(com_height: f64) -> Self {
        Self {
            com_height,
            kp: 40.0,
            kd: 8.0,
            com_pos: [0.0; 3],
            com_vel: [0.0; 3],
            support_polygon: vec![[-0.1, -0.1], [0.1, -0.1], [0.1, 0.1], [-0.1, 0.1]],
            tipping_threshold: 0.05,
        }
    }

    /// Compute the Zero-Moment Point (ZMP) in the XZ plane.
    pub fn zmp(&self, gravity: f64) -> [f64; 2] {
        // ZMP = CoM_xz - h/g * CoM_accel_xz  (linearised pendulum)
        // Without external accel, ZMP ≈ CoM projection
        let g = gravity.abs();
        [
            self.com_pos[0] - (self.com_height / g) * self.com_vel[0],
            self.com_pos[2] - (self.com_height / g) * self.com_vel[2],
        ]
    }

    /// Returns `true` when the ZMP is outside the support polygon (tipping).
    pub fn is_tipping(&self, gravity: f64) -> bool {
        let zmp = self.zmp(gravity);
        !point_in_polygon_xz(zmp, &self.support_polygon)
    }

    /// Compute a corrective ankle torque to recover balance (PD control).
    pub fn corrective_torque(&self, gravity: f64) -> [f64; 3] {
        let zmp = self.zmp(gravity);
        let err_x = self.com_pos[0] - zmp[0];
        let err_z = self.com_pos[2] - zmp[1];
        [
            -(self.kp * err_x + self.kd * self.com_vel[0]),
            0.0,
            -(self.kp * err_z + self.kd * self.com_vel[2]),
        ]
    }
}

/// Test whether `point` (XZ) is inside a convex polygon.
fn point_in_polygon_xz(point: [f64; 2], polygon: &[[f64; 2]]) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    let mut inside = true;
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let edge = [b[0] - a[0], b[1] - a[1]];
        let to_point = [point[0] - a[0], point[1] - a[1]];
        let cross = edge[0] * to_point[1] - edge[1] * to_point[0];
        if cross < 0.0 {
            inside = false;
            break;
        }
    }
    inside
}

// ---------------------------------------------------------------------------
// VaultController
// ---------------------------------------------------------------------------

/// Vault controller: detects a low obstacle and computes the vault arc.
pub struct VaultController {
    /// Maximum obstacle height to vault (metres).
    pub max_vault_height: f64,
    /// Maximum obstacle depth (metres) the character can vault.
    pub max_vault_depth: f64,
    /// Duration of the vault animation (seconds).
    pub vault_duration: f64,
    /// Whether a vault is in progress.
    pub vaulting: bool,
    /// Progress \[0..1\] of the current vault.
    pub vault_progress: f64,
    /// Start position of the vault.
    pub vault_start: [f64; 3],
    /// End position of the vault.
    pub vault_end: [f64; 3],
    /// Arc apex height (metres above start/end midpoint).
    pub arc_height: f64,
}

impl VaultController {
    /// Create a vault controller.
    pub fn new() -> Self {
        Self {
            max_vault_height: 1.2,
            max_vault_depth: 0.8,
            vault_duration: 0.5,
            vaulting: false,
            vault_progress: 0.0,
            vault_start: [0.0; 3],
            vault_end: [0.0; 3],
            arc_height: 0.3,
        }
    }

    /// Detect whether an obstacle at `obstacle_top_y` and depth `depth` can be vaulted.
    pub fn can_vault(&self, character_y: f64, obstacle_top_y: f64, depth: f64) -> bool {
        let height = obstacle_top_y - character_y;
        height >= 0.0 && height <= self.max_vault_height && depth <= self.max_vault_depth
    }

    /// Begin a vault from `start` to `end`.
    pub fn begin_vault(&mut self, start: [f64; 3], end: [f64; 3]) {
        self.vaulting = true;
        self.vault_progress = 0.0;
        self.vault_start = start;
        self.vault_end = end;
    }

    /// Advance the vault and return the interpolated position.
    pub fn tick(&mut self, dt: f64) -> [f64; 3] {
        if !self.vaulting {
            return self.vault_start;
        }
        self.vault_progress = (self.vault_progress + dt / self.vault_duration).min(1.0);
        if self.vault_progress >= 1.0 {
            self.vaulting = false;
        }
        vault_arc_position(
            self.vault_start,
            self.vault_end,
            self.arc_height,
            self.vault_progress,
        )
    }
}

impl Default for VaultController {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute position on a parabolic vault arc at normalised time `t` ∈ \[0, 1\].
pub fn vault_arc_position(start: [f64; 3], end: [f64; 3], arc_height: f64, t: f64) -> [f64; 3] {
    let base = add3(scale3(start, 1.0 - t), scale3(end, t));
    let arc_y = 4.0 * arc_height * t * (1.0 - t);
    [base[0], base[1] + arc_y, base[2]]
}

// ---------------------------------------------------------------------------
// CrowdAgent (ORCA)
// ---------------------------------------------------------------------------

/// A single agent in an ORCA-based crowd simulation.
pub struct CrowdAgent {
    /// Agent radius (metres).
    pub radius: f64,
    /// Maximum speed (m/s).
    pub max_speed: f64,
    /// Preferred (goal) velocity.
    pub preferred_velocity: [f64; 2],
    /// Current velocity.
    pub velocity: [f64; 2],
    /// Position.
    pub position: [f64; 2],
    /// Look-ahead time for ORCA (seconds).
    pub time_horizon: f64,
}

impl CrowdAgent {
    /// Create a crowd agent.
    pub fn new(radius: f64, max_speed: f64) -> Self {
        Self {
            radius,
            max_speed,
            preferred_velocity: [0.0; 2],
            velocity: [0.0; 2],
            position: [0.0; 2],
            time_horizon: 2.0,
        }
    }

    /// Compute new velocity respecting ORCA half-planes from `neighbors`.
    pub fn orca_step(&mut self, neighbors: &[CrowdAgent], dt: f64) {
        let half_planes: Vec<OrcaHalfPlane> = neighbors
            .iter()
            .filter_map(|other| {
                velocity_obstacle(
                    self.position,
                    self.velocity,
                    self.radius,
                    other.position,
                    other.velocity,
                    other.radius,
                    self.time_horizon,
                )
            })
            .collect();

        // Project preferred velocity onto feasible region (simplified LP)
        let mut new_vel = self.preferred_velocity;
        for hp in &half_planes {
            let pen = dot2(new_vel, hp.normal) - hp.offset;
            if pen < 0.0 {
                new_vel = add2(new_vel, scale2(hp.normal, -pen));
            }
        }

        // Clamp to max speed
        let spd = len2(new_vel);
        if spd > self.max_speed {
            new_vel = scale2(norm2(new_vel), self.max_speed);
        }
        self.velocity = new_vel;
        self.position = add2(self.position, scale2(self.velocity, dt));
    }
}

/// An ORCA half-plane constraint: `dot(v, normal) >= offset`.
pub struct OrcaHalfPlane {
    /// Outward normal of the half-plane.
    pub normal: [f64; 2],
    /// Offset along the normal.
    pub offset: f64,
}

/// Compute the velocity obstacle and return an ORCA half-plane for `agent` w.r.t. `other`.
pub fn velocity_obstacle(
    pos_a: [f64; 2],
    vel_a: [f64; 2],
    rad_a: f64,
    pos_b: [f64; 2],
    vel_b: [f64; 2],
    rad_b: f64,
    time_horizon: f64,
) -> Option<OrcaHalfPlane> {
    let rel_pos = sub2(pos_b, pos_a);
    let rel_vel = sub2(vel_a, vel_b);
    let combined_radius = rad_a + rad_b;
    let dist_sq = dot2(rel_pos, rel_pos);
    let dist = dist_sq.sqrt();

    if dist < 1e-6 {
        return None;
    }

    // Truncated velocity obstacle
    let inv_tau = 1.0 / time_horizon;
    let w = sub2(rel_vel, scale2(rel_pos, inv_tau));
    let w_len = len2(w);

    if w_len < 1e-9 {
        return None;
    }

    let dot_w_rel = dot2(w, rel_pos);
    let leg = (dist_sq - combined_radius * combined_radius * inv_tau * inv_tau).sqrt();
    let _ = leg; // used in full ORCA but omitted in simplified version

    let normal = if dot_w_rel < 0.0 {
        norm2(scale2(w, -1.0))
    } else {
        let perp = [-rel_pos[1] / dist, rel_pos[0] / dist];
        norm2(perp)
    };

    let offset = dot2(normal, rel_vel) - combined_radius * inv_tau;
    Some(OrcaHalfPlane { normal, offset })
}

/// Compute an ORCA half-plane from a velocity obstacle cone.
pub fn orca_half_plane(
    rel_pos: [f64; 2],
    rel_vel: [f64; 2],
    combined_radius: f64,
    time_horizon: f64,
) -> OrcaHalfPlane {
    let inv_tau = 1.0 / time_horizon;
    let w = sub2(rel_vel, scale2(rel_pos, inv_tau));
    let w_len = len2(w);
    let normal = if w_len > 1e-9 { norm2(w) } else { [1.0, 0.0] };
    let offset = dot2(normal, rel_vel) - combined_radius * inv_tau;
    OrcaHalfPlane { normal, offset }
}

/// Check capsule–ground contact: returns penetration depth (positive = penetrating).
pub fn capsule_ground_contact(capsule_base_y: f64, ground_y: f64) -> f64 {
    ground_y - capsule_base_y
}

/// Detect whether the character is stepping up (positive) or down (negative) at a ledge.
pub fn step_detection(current_y: f64, ahead_y: f64, step_height: f64) -> Option<f64> {
    let delta = ahead_y - current_y;
    if delta.abs() <= step_height {
        Some(delta)
    } else {
        None
    }
}

// 2D vector helpers (not exported separately)
fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale2(v: [f64; 2], s: f64) -> [f64; 2] {
    [v[0] * s, v[1] * s]
}

// ---------------------------------------------------------------------------
// SteeringBehavior
// ---------------------------------------------------------------------------

/// Classic steering behaviours for autonomous agents.
pub struct SteeringBehavior {
    /// Maximum force (N equivalent).
    pub max_force: f64,
    /// Maximum speed (m/s).
    pub max_speed: f64,
    /// Arrival slowing radius (metres).
    pub arrival_radius: f64,
    /// Wander circle radius.
    pub wander_radius: f64,
    /// Wander circle distance ahead.
    pub wander_distance: f64,
    /// Current wander angle.
    pub wander_angle: f64,
    /// Wander angle jitter per frame (radians).
    pub wander_jitter: f64,
    /// Separation neighbourhood radius (metres).
    pub separation_radius: f64,
    /// Cohesion neighbourhood radius (metres).
    pub cohesion_radius: f64,
    /// Alignment neighbourhood radius (metres).
    pub alignment_radius: f64,
}

impl SteeringBehavior {
    /// Create with defaults.
    pub fn new(max_force: f64, max_speed: f64) -> Self {
        Self {
            max_force,
            max_speed,
            arrival_radius: 2.0,
            wander_radius: 1.0,
            wander_distance: 2.0,
            wander_angle: 0.0,
            wander_jitter: 0.3,
            separation_radius: 1.5,
            cohesion_radius: 5.0,
            alignment_radius: 5.0,
        }
    }

    /// Seek: steer towards `target`.
    pub fn seek(&self, pos: [f64; 2], vel: [f64; 2], target: [f64; 2]) -> [f64; 2] {
        let desired = scale2(norm2(sub2(target, pos)), self.max_speed);
        truncate2(sub2(desired, vel), self.max_force)
    }

    /// Flee: steer away from `threat`.
    pub fn flee(&self, pos: [f64; 2], vel: [f64; 2], threat: [f64; 2]) -> [f64; 2] {
        let desired = scale2(norm2(sub2(pos, threat)), self.max_speed);
        truncate2(sub2(desired, vel), self.max_force)
    }

    /// Arrive: seek with slow-down near target.
    pub fn arrive(&self, pos: [f64; 2], vel: [f64; 2], target: [f64; 2]) -> [f64; 2] {
        let to_target = sub2(target, pos);
        let dist = len2(to_target);
        if dist < 0.01 {
            return [0.0; 2];
        }
        let speed = if dist < self.arrival_radius {
            self.max_speed * dist / self.arrival_radius
        } else {
            self.max_speed
        };
        let desired = scale2(norm2(to_target), speed);
        truncate2(sub2(desired, vel), self.max_force)
    }

    /// Pursue: predict `quarry_pos` based on `quarry_vel` and seek.
    pub fn pursue(
        &self,
        pos: [f64; 2],
        vel: [f64; 2],
        quarry_pos: [f64; 2],
        quarry_vel: [f64; 2],
    ) -> [f64; 2] {
        let to_quarry = sub2(quarry_pos, pos);
        let dist = len2(to_quarry);
        let t = dist / (self.max_speed + 0.001);
        let predicted = add2(quarry_pos, scale2(quarry_vel, t));
        self.seek(pos, vel, predicted)
    }

    /// Evade: predict `threat` and flee.
    pub fn evade(
        &self,
        pos: [f64; 2],
        vel: [f64; 2],
        threat_pos: [f64; 2],
        threat_vel: [f64; 2],
    ) -> [f64; 2] {
        let to_threat = sub2(threat_pos, pos);
        let dist = len2(to_threat);
        let t = dist / (self.max_speed + 0.001);
        let predicted = add2(threat_pos, scale2(threat_vel, t));
        self.flee(pos, vel, predicted)
    }

    /// Wander: add small random jitter to a circle ahead of the agent.
    pub fn wander(&mut self, pos: [f64; 2], vel: [f64; 2], jitter: f64) -> [f64; 2] {
        self.wander_angle += jitter * self.wander_jitter;
        let heading = if len2(vel) > 0.01 {
            norm2(vel)
        } else {
            [1.0, 0.0]
        };
        let circle_centre = add2(pos, scale2(heading, self.wander_distance));
        let wander_offset = [
            self.wander_radius * self.wander_angle.cos(),
            self.wander_radius * self.wander_angle.sin(),
        ];
        let target = add2(circle_centre, wander_offset);
        self.seek(pos, vel, target)
    }

    /// Separation: steer away from nearby neighbours.
    pub fn separation(&self, pos: [f64; 2], vel: [f64; 2], neighbors: &[[f64; 2]]) -> [f64; 2] {
        let mut steering = [0.0_f64; 2];
        let mut count = 0;
        for &n in neighbors {
            let d = len2(sub2(pos, n));
            if d > 0.001 && d < self.separation_radius {
                let away = scale2(norm2(sub2(pos, n)), 1.0 / d);
                steering = add2(steering, away);
                count += 1;
            }
        }
        if count > 0 {
            steering = scale2(steering, 1.0 / count as f64);
            steering = sub2(scale2(norm2(steering), self.max_speed), vel);
        }
        truncate2(steering, self.max_force)
    }

    /// Cohesion: steer towards average position of neighbours.
    pub fn cohesion(&self, pos: [f64; 2], vel: [f64; 2], neighbors: &[[f64; 2]]) -> [f64; 2] {
        let mut centre = [0.0_f64; 2];
        let mut count = 0;
        for &n in neighbors {
            if len2(sub2(pos, n)) < self.cohesion_radius {
                centre = add2(centre, n);
                count += 1;
            }
        }
        if count == 0 {
            return [0.0; 2];
        }
        let target = scale2(centre, 1.0 / count as f64);
        self.arrive(pos, vel, target)
    }

    /// Alignment: steer to match average velocity of neighbours.
    pub fn alignment(&self, vel: [f64; 2], neighbor_vels: &[[f64; 2]]) -> [f64; 2] {
        let mut avg = [0.0_f64; 2];
        let n = neighbor_vels.len();
        if n == 0 {
            return [0.0; 2];
        }
        for &v in neighbor_vels {
            avg = add2(avg, v);
        }
        avg = scale2(avg, 1.0 / n as f64);
        let desired = scale2(norm2(avg), self.max_speed);
        truncate2(sub2(desired, vel), self.max_force)
    }
}

/// Truncate a 2D vector to maximum magnitude.
fn truncate2(v: [f64; 2], max_len: f64) -> [f64; 2] {
    let l = len2(v);
    if l > max_len {
        scale2(v, max_len / l)
    } else {
        v
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. CharacterController: initial state
    #[test]
    fn test_char_controller_initial() {
        let c = CharacterController::new(0.3, 1.8);
        assert!(!c.grounded);
        assert_eq!(c.position, [0.0; 3]);
    }

    // 2. CharacterController: ground check snaps position
    #[test]
    fn test_char_controller_ground_check() {
        let mut c = CharacterController::new(0.3, 1.8);
        c.position = [0.0, 0.0, 0.0];
        let grounded = c.check_ground(0.9); // foot at -0.9, ground at 0.9 → pen = 1.8
        assert!(grounded);
        assert!(c.position[1] > 0.0);
    }

    // 3. CharacterController: above ground → not grounded
    #[test]
    fn test_char_controller_airborne() {
        let mut c = CharacterController::new(0.3, 1.8);
        c.position = [0.0, 5.0, 0.0];
        let grounded = c.check_ground(0.0);
        assert!(!grounded);
    }

    // 4. CharacterController: gravity integration
    #[test]
    fn test_char_controller_gravity() {
        let mut c = CharacterController::new(0.3, 1.8);
        c.position = [0.0, 10.0, 0.0];
        c.grounded = false;
        c.integrate(0.1);
        assert!(c.velocity[1] < 0.0);
    }

    // 5. CharacterController: slope angle calculation
    #[test]
    fn test_char_controller_slope_angle() {
        let mut c = CharacterController::new(0.3, 1.8);
        c.ground_normal = [0.0, 1.0, 0.0];
        assert!((c.slope_angle()).abs() < 1e-9);
    }

    // 6. WalkController: accelerates towards target speed
    #[test]
    fn test_walk_controller_accelerate() {
        let mut wc = WalkController::new(5.0);
        wc.update([1.0, 0.0], 0.0, 0.1);
        assert!(wc.current_speed > 0.0);
    }

    // 7. WalkController: no input → decelerates
    #[test]
    fn test_walk_controller_decelerate() {
        let mut wc = WalkController::new(5.0);
        wc.current_speed = 3.0;
        wc.update([0.0, 0.0], 0.0, 0.1);
        assert!(wc.current_speed < 3.0);
    }

    // 8. WalkController: slope reduces speed
    #[test]
    fn test_walk_controller_slope() {
        let mut wc = WalkController::new(5.0);
        wc.current_speed = 0.0;
        let v1 = wc.update([1.0, 0.0], 0.0, 1.0);
        let mut wc2 = WalkController::new(5.0);
        wc2.current_speed = 0.0;
        let v2 = wc2.update([1.0, 0.0], 0.5, 1.0);
        assert!(len2(v1) >= len2(v2));
    }

    // 9. JumpController: pressing jump fills buffer
    #[test]
    fn test_jump_buffer_fill() {
        let mut jc = JumpController::new(6.0);
        jc.press_jump();
        assert!(jc.buffer_timer > 0.0);
    }

    // 10. JumpController: jump executes with coyote time
    #[test]
    fn test_jump_coyote() {
        let mut jc = JumpController::new(6.0);
        jc.coyote_timer = 0.05;
        jc.press_jump();
        let impulse = jc.try_jump();
        assert!((impulse - 6.0).abs() < 1e-9);
    }

    // 11. JumpController: no jump without coyote and buffer
    #[test]
    fn test_jump_no_coyote() {
        let mut jc = JumpController::new(6.0);
        let impulse = jc.try_jump();
        assert_eq!(impulse, 0.0);
    }

    // 12. ClimbController: attach within snap distance
    #[test]
    fn test_climb_attach() {
        let mut cc = ClimbController::new(2.0);
        let attached = cc.try_attach([0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [0.0, 0.0, -1.0]);
        assert!(attached);
        assert!(cc.is_climbing);
    }

    // 13. ClimbController: detach
    #[test]
    fn test_climb_detach() {
        let mut cc = ClimbController::new(2.0);
        cc.is_climbing = true;
        cc.detach();
        assert!(!cc.is_climbing);
    }

    // 14. ClimbController: velocity zero when not climbing
    #[test]
    fn test_climb_velocity_not_climbing() {
        let cc = ClimbController::new(2.0);
        let v = cc.climb_velocity(1.0, 0.0);
        assert_eq!(v, [0.0; 3]);
    }

    // 15. SwimController: submerged when below surface
    #[test]
    fn test_swim_submerged() {
        let mut sc = SwimController::new(3.0, 0.0);
        sc.update_submersion(-1.0);
        assert!(sc.is_submerged);
    }

    // 16. SwimController: not submerged above surface
    #[test]
    fn test_swim_not_submerged() {
        let mut sc = SwimController::new(3.0, 0.0);
        sc.update_submersion(1.0);
        assert!(!sc.is_submerged);
    }

    // 17. SwimController: buoyancy partially cancels gravity
    #[test]
    fn test_swim_buoyancy() {
        let mut sc = SwimController::new(3.0, 0.0);
        sc.is_submerged = true;
        let f = sc.apply_forces(-9.81);
        // net vertical should be smaller magnitude than raw gravity
        assert!(f[1].abs() < 9.81);
    }

    // 18. ProceduralFootIk: IK solution within reach
    #[test]
    fn test_foot_ik_solution() {
        let ik = ProceduralFootIk::new(0.5, 0.5);
        let result = ik.solve([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(result.is_some());
    }

    // 19. ProceduralFootIk: over-extension returns stretched result
    #[test]
    fn test_foot_ik_overreach() {
        let ik = ProceduralFootIk::new(0.5, 0.5);
        let result = ik.solve([0.0, 0.0, 0.0], [0.0, 5.0, 0.0]);
        assert!(result.is_some());
    }

    // 20. BalanceController: ZMP at rest is CoM projection
    #[test]
    fn test_zmp_at_rest() {
        let mut bc = BalanceController::new(1.0);
        bc.com_pos = [0.1, 1.0, 0.0];
        bc.com_vel = [0.0; 3];
        let zmp = bc.zmp(9.81);
        assert!((zmp[0] - 0.1).abs() < 1e-6);
    }

    // 21. VaultController: can_vault check
    #[test]
    fn test_vault_can_vault() {
        let vc = VaultController::new();
        assert!(vc.can_vault(0.0, 1.0, 0.5));
        assert!(!vc.can_vault(0.0, 2.0, 0.5)); // too high
    }

    // 22. VaultController: arc position at t=0.5 has apex
    #[test]
    fn test_vault_arc_apex() {
        let mid = vault_arc_position([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.5);
        assert!((mid[1] - 1.0).abs() < 1e-9); // 4*1*0.5*0.5 = 1.0
    }

    // 23. CrowdAgent: initial velocity is zero
    #[test]
    fn test_crowd_agent_initial() {
        let a = CrowdAgent::new(0.3, 1.5);
        assert_eq!(a.velocity, [0.0; 2]);
    }

    // 24. SteeringBehavior: seek moves towards target
    #[test]
    fn test_steering_seek() {
        let sb = SteeringBehavior::new(10.0, 5.0);
        let f = sb.seek([0.0, 0.0], [0.0, 0.0], [3.0, 0.0]);
        assert!(f[0] > 0.0);
    }

    // 25. SteeringBehavior: flee moves away
    #[test]
    fn test_steering_flee() {
        let sb = SteeringBehavior::new(10.0, 5.0);
        let f = sb.flee([0.0, 0.0], [0.0, 0.0], [3.0, 0.0]);
        assert!(f[0] < 0.0);
    }

    // 26. SteeringBehavior: arrive returns zero at target
    #[test]
    fn test_steering_arrive_at_target() {
        let sb = SteeringBehavior::new(10.0, 5.0);
        let f = sb.arrive([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]);
        assert_eq!(f, [0.0; 2]);
    }

    // 27. velocity_obstacle: returns half-plane
    #[test]
    fn test_velocity_obstacle() {
        let hp = velocity_obstacle(
            [0.0, 0.0],
            [1.0, 0.0],
            0.3,
            [3.0, 0.0],
            [0.0, 0.0],
            0.3,
            2.0,
        );
        assert!(hp.is_some());
    }

    // 28. capsule_ground_contact: positive when penetrating
    #[test]
    fn test_capsule_ground_contact_positive() {
        let pen = capsule_ground_contact(0.5, 1.0);
        assert!(pen > 0.0);
    }

    // 29. step_detection: step within limit
    #[test]
    fn test_step_detection_within() {
        let s = step_detection(0.0, 0.2, 0.3);
        assert!(s.is_some());
        assert!((s.unwrap() - 0.2).abs() < 1e-9);
    }

    // 30. step_detection: step exceeds limit
    #[test]
    fn test_step_detection_exceed() {
        let s = step_detection(0.0, 1.0, 0.3);
        assert!(s.is_none());
    }

    // 31. angle_diff wrapping
    #[test]
    fn test_angle_diff_wrap() {
        let d = angle_diff(0.1, 6.0);
        assert!(d.abs() < PI + 0.001);
    }

    // 32. point_in_polygon_xz: inside unit square
    #[test]
    fn test_point_in_polygon_inside() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(point_in_polygon_xz([0.5, 0.5], &poly));
    }

    // 33. orca_half_plane: normal is unit
    #[test]
    fn test_orca_half_plane_unit_normal() {
        let hp = orca_half_plane([1.0, 0.0], [0.5, 0.0], 0.5, 2.0);
        let l = len2(hp.normal);
        assert!((l - 1.0).abs() < 1e-9);
    }
}
