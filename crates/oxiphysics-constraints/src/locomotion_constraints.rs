// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Locomotion constraints for bipedal and legged robotics.
//!
//! This module provides constraint-based formulations for bipedal walking,
//! balance control, and gait planning. Key components include:
//!
//! - **Inverted pendulum model** (linear / 3D): simplified CoM dynamics for
//!   walking analysis and control.
//! - **Zero Moment Point (ZMP)**: trajectory planning and stability checking.
//! - **Footstep planning**: discrete step placement with timing.
//! - **Swing leg trajectory**: Bezier-based foot trajectories.
//! - **Contact schedule**: time-indexed contact state for each foot.
//! - **Gait phase detection**: stance, swing, double-support identification.
//! - **Dynamic balance**: centroidal momentum and angular momentum regulation.
//! - **Capture point / DCM**: instantaneous capture point and divergent
//!   component of motion for push-recovery.
//!
//! All vectors use `[f64; 3]` arrays (no nalgebra dependency).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Add two 3-vectors.
#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors (a - b).
#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn v3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Length of a 3-vector.
#[inline]
fn v3_len(v: [f64; 3]) -> f64 {
    v3_dot(v, v).sqrt()
}

/// Normalize a 3-vector; returns zero if degenerate.
#[inline]
fn v3_norm(v: [f64; 3]) -> [f64; 3] {
    let l = v3_len(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        v3_scale(v, 1.0 / l)
    }
}

/// Linear interpolation between two 3-vectors.
#[inline]
fn v3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Length of a 2D vector (x, y).
#[inline]
fn v2_len(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Inverted Pendulum Model
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for a Linear Inverted Pendulum Model (LIPM).
///
/// The LIPM approximates bipedal walking by constraining the centre-of-mass
/// (CoM) to move at a constant height above the ground while the support leg
/// acts as an inverted pendulum.
#[derive(Debug, Clone, Copy)]
pub struct LipmParams {
    /// CoM height above ground \[m\].
    pub com_height: f64,
    /// Gravitational acceleration \[m/s^2\].
    pub gravity: f64,
    /// Total body mass \[kg\].
    pub mass: f64,
}

impl LipmParams {
    /// Create new LIPM parameters.
    pub fn new(com_height: f64, gravity: f64, mass: f64) -> Self {
        Self {
            com_height,
            gravity,
            mass,
        }
    }

    /// Natural frequency of the LIPM: omega = sqrt(g / z_c).
    pub fn omega(&self) -> f64 {
        if self.com_height.abs() < 1e-15 {
            return 0.0;
        }
        (self.gravity / self.com_height).sqrt()
    }

    /// Time constant: 1 / omega.
    pub fn time_constant(&self) -> f64 {
        let w = self.omega();
        if w.abs() < 1e-15 { 0.0 } else { 1.0 / w }
    }
}

/// State of a 2D LIPM in the sagittal or lateral plane.
///
/// Horizontal position and velocity relative to the support foot.
#[derive(Debug, Clone, Copy)]
pub struct LipmState2d {
    /// Horizontal position \[m\].
    pub x: f64,
    /// Horizontal velocity \[m/s\].
    pub xdot: f64,
}

impl LipmState2d {
    /// Create a new 2D LIPM state.
    pub fn new(x: f64, xdot: f64) -> Self {
        Self { x, xdot }
    }

    /// Integrate the LIPM forward by `dt` seconds using the closed-form
    /// hyperbolic solution.
    ///
    /// x(t) = x0 * cosh(wt) + (xdot0 / w) * sinh(wt)
    /// xdot(t) = x0 * w * sinh(wt) + xdot0 * cosh(wt)
    pub fn step(&self, params: &LipmParams, dt: f64) -> Self {
        let w = params.omega();
        if w.abs() < 1e-15 {
            return Self {
                x: self.x + self.xdot * dt,
                xdot: self.xdot,
            };
        }
        let wt = w * dt;
        let ch = wt.cosh();
        let sh = wt.sinh();
        Self {
            x: self.x * ch + (self.xdot / w) * sh,
            xdot: self.x * w * sh + self.xdot * ch,
        }
    }

    /// Compute the orbital energy of the LIPM.
    ///
    /// E = 0.5 * (xdot^2 - omega^2 * x^2)
    pub fn orbital_energy(&self, params: &LipmParams) -> f64 {
        let w = params.omega();
        0.5 * (self.xdot * self.xdot - w * w * self.x * self.x)
    }
}

/// 3D LIPM state (position + velocity in horizontal plane).
#[derive(Debug, Clone, Copy)]
pub struct LipmState3d {
    /// CoM position (x, y, z_c).
    pub pos: [f64; 3],
    /// CoM velocity.
    pub vel: [f64; 3],
}

impl LipmState3d {
    /// Create a new 3D LIPM state.
    pub fn new(pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self { pos, vel }
    }

    /// Step forward by `dt` using independent x/y LIPM integration.
    pub fn step(&self, params: &LipmParams, dt: f64) -> Self {
        let sx = LipmState2d::new(self.pos[0], self.vel[0]).step(params, dt);
        let sy = LipmState2d::new(self.pos[1], self.vel[1]).step(params, dt);
        Self {
            pos: [sx.x, sy.x, params.com_height],
            vel: [sx.xdot, sy.xdot, 0.0],
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Zero Moment Point (ZMP)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Zero Moment Point from CoM position and acceleration.
///
/// ZMP_x = x_com - (z_com / g) * x_ddot
/// ZMP_y = y_com - (z_com / g) * y_ddot
///
/// Returns the ZMP position as `[x, y, 0]`.
pub fn compute_zmp(com_pos: [f64; 3], com_accel: [f64; 3], gravity: f64) -> [f64; 3] {
    if gravity.abs() < 1e-15 {
        return [com_pos[0], com_pos[1], 0.0];
    }
    let ratio = com_pos[2] / gravity;
    [
        com_pos[0] - ratio * com_accel[0],
        com_pos[1] - ratio * com_accel[1],
        0.0,
    ]
}

/// Check if a ZMP lies within a rectangular support polygon.
///
/// The support polygon is defined by its centre `[cx, cy]` and half-widths
/// `[hx, hy]`.
pub fn zmp_in_support_rect(zmp: [f64; 3], centre: [f64; 2], half_widths: [f64; 2]) -> bool {
    (zmp[0] - centre[0]).abs() <= half_widths[0] && (zmp[1] - centre[1]).abs() <= half_widths[1]
}

/// Plan a ZMP trajectory as a piecewise-linear path between waypoints.
///
/// Each waypoint is `(position_xy, duration)`. Returns the ZMP position at
/// `time` by linear interpolation.
pub fn zmp_trajectory_at(waypoints: &[([f64; 2], f64)], time: f64) -> [f64; 2] {
    if waypoints.is_empty() {
        return [0.0, 0.0];
    }
    if waypoints.len() == 1 || time <= 0.0 {
        return waypoints[0].0;
    }
    let mut t_acc = 0.0;
    for i in 0..waypoints.len() - 1 {
        let dur = waypoints[i].1;
        if time <= t_acc + dur {
            let frac = if dur.abs() < 1e-15 {
                1.0
            } else {
                (time - t_acc) / dur
            };
            let a = waypoints[i].0;
            let b = waypoints[i + 1].0;
            return [a[0] + (b[0] - a[0]) * frac, a[1] + (b[1] - a[1]) * frac];
        }
        t_acc += dur;
    }
    waypoints.last().expect("collection should not be empty").0
}

/// ZMP stability margin: minimum distance from ZMP to the edges of a
/// rectangular support polygon.
pub fn zmp_stability_margin(zmp: [f64; 3], centre: [f64; 2], half_widths: [f64; 2]) -> f64 {
    let dx = half_widths[0] - (zmp[0] - centre[0]).abs();
    let dy = half_widths[1] - (zmp[1] - centre[1]).abs();
    dx.min(dy)
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Footstep Planning
// ─────────────────────────────────────────────────────────────────────────────

/// Which foot is in contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Foot {
    /// Left foot.
    Left,
    /// Right foot.
    Right,
}

impl Foot {
    /// Return the opposite foot.
    pub fn opposite(self) -> Self {
        match self {
            Foot::Left => Foot::Right,
            Foot::Right => Foot::Left,
        }
    }
}

/// A planned footstep.
#[derive(Debug, Clone, Copy)]
pub struct Footstep {
    /// Which foot.
    pub foot: Foot,
    /// Position of the foot centre on the ground plane \[x, y, 0\].
    pub position: [f64; 3],
    /// Heading angle of the foot \[rad\].
    pub heading: f64,
    /// Timing: start of this step \[s\].
    pub t_start: f64,
    /// Duration of the swing phase \[s\].
    pub swing_duration: f64,
    /// Duration of the double-support phase after landing \[s\].
    pub double_support_duration: f64,
}

impl Footstep {
    /// Create a new footstep.
    pub fn new(
        foot: Foot,
        position: [f64; 3],
        heading: f64,
        t_start: f64,
        swing_duration: f64,
        double_support_duration: f64,
    ) -> Self {
        Self {
            foot,
            position,
            heading,
            t_start,
            swing_duration,
            double_support_duration,
        }
    }

    /// Total duration (swing + double-support).
    pub fn total_duration(&self) -> f64 {
        self.swing_duration + self.double_support_duration
    }

    /// End time of the swing phase.
    pub fn swing_end_time(&self) -> f64 {
        self.t_start + self.swing_duration
    }

    /// End time of this entire footstep phase.
    pub fn end_time(&self) -> f64 {
        self.t_start + self.total_duration()
    }
}

/// Generate a straight-line footstep plan.
///
/// `n_steps` alternating footsteps along the forward direction, starting
/// with `start_foot`. Step length and lateral offset are configurable.
pub fn plan_straight_footsteps(
    start_pos: [f64; 3],
    step_length: f64,
    lateral_offset: f64,
    swing_dur: f64,
    ds_dur: f64,
    n_steps: usize,
    start_foot: Foot,
) -> Vec<Footstep> {
    let mut steps = Vec::with_capacity(n_steps);
    let mut current_foot = start_foot;
    let mut t = 0.0;
    for i in 0..n_steps {
        let x = start_pos[0] + step_length * (i as f64 + 1.0);
        let y_sign = match current_foot {
            Foot::Left => 1.0,
            Foot::Right => -1.0,
        };
        let y = start_pos[1] + lateral_offset * y_sign;
        steps.push(Footstep::new(
            current_foot,
            [x, y, 0.0],
            0.0,
            t,
            swing_dur,
            ds_dur,
        ));
        t += swing_dur + ds_dur;
        current_foot = current_foot.opposite();
    }
    steps
}

/// Compute the centre of the support polygon for a double-support phase
/// given left and right foot positions.
pub fn double_support_centre(left_pos: [f64; 3], right_pos: [f64; 3]) -> [f64; 3] {
    [
        0.5 * (left_pos[0] + right_pos[0]),
        0.5 * (left_pos[1] + right_pos[1]),
        0.0,
    ]
}

/// Compute the maximum step length based on kinematic reach.
///
/// `leg_length` is the total leg length and `max_reach_fraction` is the
/// fraction of leg length available as step length (typically 0.3 - 0.6).
pub fn max_step_length(leg_length: f64, max_reach_fraction: f64) -> f64 {
    leg_length * max_reach_fraction.clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Swing Leg Trajectory
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a swing foot trajectory using a cubic Bezier curve in 3D.
///
/// The trajectory goes from `start` to `end` with an apex height of
/// `step_height` above the midpoint. Returns the position at parameter `t`
/// in `[0, 1]`.
pub fn swing_bezier_trajectory(
    start: [f64; 3],
    end: [f64; 3],
    step_height: f64,
    t: f64,
) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    // Control points: P0=start, P1=start+lift, P2=end+lift, P3=end
    let mid_z = start[2].max(end[2]) + step_height;
    let p1 = [start[0], start[1], mid_z];
    let p2 = [end[0], end[1], mid_z];

    let mt = 1.0 - t;
    let mt2 = mt * mt;
    let mt3 = mt2 * mt;
    let t2 = t * t;
    let t3 = t2 * t;

    [
        mt3 * start[0] + 3.0 * mt2 * t * p1[0] + 3.0 * mt * t2 * p2[0] + t3 * end[0],
        mt3 * start[1] + 3.0 * mt2 * t * p1[1] + 3.0 * mt * t2 * p2[1] + t3 * end[1],
        mt3 * start[2] + 3.0 * mt2 * t * p1[2] + 3.0 * mt * t2 * p2[2] + t3 * end[2],
    ]
}

/// Compute the velocity of the swing bezier trajectory at parameter `t`.
///
/// Derivative of the cubic Bezier: 3*(1-t)^2*(P1-P0) + 6*(1-t)*t*(P2-P1)
/// + 3*t^2*(P3-P2).
pub fn swing_bezier_velocity(start: [f64; 3], end: [f64; 3], step_height: f64, t: f64) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let mid_z = start[2].max(end[2]) + step_height;
    let p1 = [start[0], start[1], mid_z];
    let p2 = [end[0], end[1], mid_z];

    let mt = 1.0 - t;
    let d01 = v3_sub(p1, start);
    let d12 = v3_sub(p2, p1);
    let d23 = v3_sub(end, p2);

    v3_add(
        v3_add(v3_scale(d01, 3.0 * mt * mt), v3_scale(d12, 6.0 * mt * t)),
        v3_scale(d23, 3.0 * t * t),
    )
}

/// Generate a parabolic swing trajectory (simpler alternative to Bezier).
///
/// The foot moves linearly in x-y and follows a parabola in z.
pub fn swing_parabolic_trajectory(
    start: [f64; 3],
    end: [f64; 3],
    step_height: f64,
    t: f64,
) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let xy = v3_lerp(start, end, t);
    // Parabola: z = z_interp + 4*h*t*(1-t)
    let z_base = start[2] + (end[2] - start[2]) * t;
    let z = z_base + 4.0 * step_height * t * (1.0 - t);
    [xy[0], xy[1], z]
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Contact Schedule
// ─────────────────────────────────────────────────────────────────────────────

/// Contact state of a single foot at a given instant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactState {
    /// Foot is in contact with the ground.
    InContact,
    /// Foot is in the air (swing phase).
    Swing,
}

/// Entry in a contact schedule.
#[derive(Debug, Clone, Copy)]
pub struct ContactEvent {
    /// The foot this event applies to.
    pub foot: Foot,
    /// Contact state.
    pub state: ContactState,
    /// Start time of this event.
    pub t_start: f64,
    /// Duration of this event.
    pub duration: f64,
}

impl ContactEvent {
    /// Create a new contact event.
    pub fn new(foot: Foot, state: ContactState, t_start: f64, duration: f64) -> Self {
        Self {
            foot,
            state,
            t_start,
            duration,
        }
    }

    /// End time of this event.
    pub fn t_end(&self) -> f64 {
        self.t_start + self.duration
    }
}

/// A contact schedule is a time-ordered list of contact events.
#[derive(Debug, Clone)]
pub struct ContactSchedule {
    /// Ordered list of events.
    pub events: Vec<ContactEvent>,
}

impl ContactSchedule {
    /// Create an empty schedule.
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }

    /// Add an event.
    pub fn push(&mut self, event: ContactEvent) {
        self.events.push(event);
    }

    /// Query the contact state of a foot at a given time.
    pub fn state_at(&self, foot: Foot, time: f64) -> Option<ContactState> {
        for ev in &self.events {
            if ev.foot == foot && time >= ev.t_start && time < ev.t_end() {
                return Some(ev.state);
            }
        }
        None
    }

    /// Build a contact schedule from a footstep plan.
    ///
    /// Each footstep produces a swing event for the stepping foot, followed
    /// by an in-contact event. The stance foot is in contact throughout.
    pub fn from_footsteps(steps: &[Footstep]) -> Self {
        let mut sched = Self::new();
        for step in steps {
            sched.push(ContactEvent::new(
                step.foot,
                ContactState::Swing,
                step.t_start,
                step.swing_duration,
            ));
            sched.push(ContactEvent::new(
                step.foot,
                ContactState::InContact,
                step.swing_end_time(),
                step.double_support_duration,
            ));
        }
        sched
    }

    /// Total duration covered by the schedule.
    pub fn total_duration(&self) -> f64 {
        self.events
            .iter()
            .map(|e| e.t_end())
            .fold(0.0_f64, f64::max)
    }
}

impl Default for ContactSchedule {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Gait Phase Detection
// ─────────────────────────────────────────────────────────────────────────────

/// High-level gait phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaitPhase {
    /// Left single-support (right foot swinging).
    LeftStance,
    /// Right single-support (left foot swinging).
    RightStance,
    /// Both feet on the ground.
    DoubleSupport,
    /// Neither foot on the ground (flight, e.g. running).
    Flight,
}

/// Detect the gait phase from left and right foot contact states.
pub fn detect_gait_phase(left: ContactState, right: ContactState) -> GaitPhase {
    match (left, right) {
        (ContactState::InContact, ContactState::InContact) => GaitPhase::DoubleSupport,
        (ContactState::InContact, ContactState::Swing) => GaitPhase::LeftStance,
        (ContactState::Swing, ContactState::InContact) => GaitPhase::RightStance,
        (ContactState::Swing, ContactState::Swing) => GaitPhase::Flight,
    }
}

/// Compute the gait phase fraction within a periodic gait cycle.
///
/// `time` is the current simulation time and `cycle_duration` is the full
/// gait cycle period. Returns a value in `[0, 1)`.
pub fn gait_phase_fraction(time: f64, cycle_duration: f64) -> f64 {
    if cycle_duration.abs() < 1e-15 {
        return 0.0;
    }
    let phase = (time % cycle_duration) / cycle_duration;
    phase.clamp(0.0, 1.0 - 1e-15)
}

/// Determine left/right swing timing from a duty factor and phase offset.
///
/// Returns `(left_is_swing, right_is_swing)` at the given phase fraction.
/// `duty_factor` is the fraction of the cycle that each foot spends on the
/// ground (0.5 = walking, < 0.5 = running).
pub fn bilateral_swing_timing(
    phase: f64,
    duty_factor: f64,
    right_phase_offset: f64,
) -> (bool, bool) {
    let swing_frac = 1.0 - duty_factor;
    let left_swing = phase >= duty_factor && phase < 1.0;
    let right_phase = (phase + right_phase_offset) % 1.0;
    let right_swing = right_phase >= duty_factor && right_phase < duty_factor + swing_frac;
    (left_swing, right_swing)
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Dynamic Balance
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the centre of mass of a set of point masses.
///
/// Each entry is `(position, mass)`.
pub fn compute_com(bodies: &[([f64; 3], f64)]) -> [f64; 3] {
    let mut total_mass = 0.0;
    let mut weighted = [0.0; 3];
    for &(pos, m) in bodies {
        total_mass += m;
        weighted[0] += pos[0] * m;
        weighted[1] += pos[1] * m;
        weighted[2] += pos[2] * m;
    }
    if total_mass.abs() < 1e-15 {
        return [0.0; 3];
    }
    v3_scale(weighted, 1.0 / total_mass)
}

/// Compute the total linear momentum.
///
/// Each entry is `(velocity, mass)`.
pub fn compute_linear_momentum(bodies: &[([f64; 3], f64)]) -> [f64; 3] {
    let mut p = [0.0; 3];
    for &(vel, m) in bodies {
        p[0] += vel[0] * m;
        p[1] += vel[1] * m;
        p[2] += vel[2] * m;
    }
    p
}

/// Compute the centroidal angular momentum about a reference point.
///
/// Each body contributes `(r_i - ref) x (m_i * v_i)`.
pub fn compute_centroidal_angular_momentum(
    bodies: &[([f64; 3], [f64; 3], f64)], // (position, velocity, mass)
    reference: [f64; 3],
) -> [f64; 3] {
    let mut l = [0.0; 3];
    for &(pos, vel, m) in bodies {
        let r = v3_sub(pos, reference);
        let mv = v3_scale(vel, m);
        let c = v3_cross(r, mv);
        l[0] += c[0];
        l[1] += c[1];
        l[2] += c[2];
    }
    l
}

/// Rate of change of centroidal angular momentum (net torque about CoM).
///
/// `forces` is a list of `(application_point, force_vector)`.
pub fn centroidal_angular_momentum_rate(
    forces: &[([f64; 3], [f64; 3])],
    com: [f64; 3],
) -> [f64; 3] {
    let mut tau = [0.0; 3];
    for &(point, force) in forces {
        let r = v3_sub(point, com);
        let c = v3_cross(r, force);
        tau[0] += c[0];
        tau[1] += c[1];
        tau[2] += c[2];
    }
    tau
}

/// Check dynamic balance: the net moment about the CoM projected onto
/// the ground must be zero for static equilibrium.
///
/// Returns the magnitude of the horizontal torque components.
pub fn dynamic_balance_error(forces: &[([f64; 3], [f64; 3])], com: [f64; 3]) -> f64 {
    let tau = centroidal_angular_momentum_rate(forces, com);
    // Only horizontal components affect tipping
    v2_len(tau[0], tau[1])
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. Centroidal Momentum Matrix (simplified)
// ─────────────────────────────────────────────────────────────────────────────

/// Centroidal momentum vector: `[linear_x, linear_y, linear_z, angular_x, angular_y, angular_z]`.
pub type CentroidalMomentum = [f64; 6];

/// Compute the full 6D centroidal momentum.
pub fn compute_centroidal_momentum(
    bodies: &[([f64; 3], [f64; 3], f64)], // (pos, vel, mass)
    com: [f64; 3],
) -> CentroidalMomentum {
    let mut h = [0.0; 6];
    for &(pos, vel, m) in bodies {
        // Linear momentum
        h[0] += m * vel[0];
        h[1] += m * vel[1];
        h[2] += m * vel[2];
        // Angular momentum about CoM
        let r = v3_sub(pos, com);
        let mv = v3_scale(vel, m);
        let c = v3_cross(r, mv);
        h[3] += c[0];
        h[4] += c[1];
        h[5] += c[2];
    }
    h
}

/// Desired centroidal momentum rate: compute the required net wrench
/// to achieve a target momentum rate.
///
/// Returns the error between desired and current momentum rates.
pub fn centroidal_momentum_error(
    current: CentroidalMomentum,
    desired: CentroidalMomentum,
) -> CentroidalMomentum {
    let mut err = [0.0; 6];
    for i in 0..6 {
        err[i] = desired[i] - current[i];
    }
    err
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. Capture Point
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Instantaneous Capture Point (ICP) / extrapolated CoM.
///
/// ICP = com_xy + (1/omega) * com_vel_xy
///
/// The capture point is where the robot must place its foot to come to a
/// stop on one foot.
pub fn capture_point(com_pos: [f64; 3], com_vel: [f64; 3], omega: f64) -> [f64; 2] {
    if omega.abs() < 1e-15 {
        return [com_pos[0], com_pos[1]];
    }
    let inv_w = 1.0 / omega;
    [
        com_pos[0] + inv_w * com_vel[0],
        com_pos[1] + inv_w * com_vel[1],
    ]
}

/// Compute the capture point from LIPM parameters.
pub fn capture_point_lipm(com_pos: [f64; 3], com_vel: [f64; 3], params: &LipmParams) -> [f64; 2] {
    capture_point(com_pos, com_vel, params.omega())
}

/// Check if the capture point is within a rectangular support region.
pub fn capture_point_feasible(
    cp: [f64; 2],
    support_centre: [f64; 2],
    half_widths: [f64; 2],
) -> bool {
    (cp[0] - support_centre[0]).abs() <= half_widths[0]
        && (cp[1] - support_centre[1]).abs() <= half_widths[1]
}

/// Distance from the capture point to the nearest edge of the support
/// polygon (rectangular). Positive means inside.
pub fn capture_point_margin(cp: [f64; 2], support_centre: [f64; 2], half_widths: [f64; 2]) -> f64 {
    let dx = half_widths[0] - (cp[0] - support_centre[0]).abs();
    let dy = half_widths[1] - (cp[1] - support_centre[1]).abs();
    dx.min(dy)
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. Divergent Component of Motion (DCM)
// ─────────────────────────────────────────────────────────────────────────────

/// The Divergent Component of Motion in the horizontal plane.
///
/// DCM = com_xy + (1/omega) * com_vel_xy
///
/// (Same as the capture point formula but used in DCM-based controllers.)
#[derive(Debug, Clone, Copy)]
pub struct DcmState {
    /// DCM position in the horizontal plane \[x, y\].
    pub xi: [f64; 2],
}

impl DcmState {
    /// Compute the DCM from CoM state.
    pub fn from_com(com_pos: [f64; 3], com_vel: [f64; 3], omega: f64) -> Self {
        let cp = capture_point(com_pos, com_vel, omega);
        Self { xi: cp }
    }

    /// Propagate DCM forward in time given a constant virtual repellent
    /// point (VRP / eCMP).
    ///
    /// xi(t) = vrp + exp(omega * t) * (xi_0 - vrp)
    pub fn propagate(&self, vrp: [f64; 2], omega: f64, dt: f64) -> Self {
        let exp_wt = (omega * dt).exp();
        Self {
            xi: [
                vrp[0] + exp_wt * (self.xi[0] - vrp[0]),
                vrp[1] + exp_wt * (self.xi[1] - vrp[1]),
            ],
        }
    }

    /// Compute the required VRP to steer the DCM towards a reference in
    /// one step.
    ///
    /// vrp = xi - (1/omega) * (xi_dot_des)
    ///
    /// Using simple proportional control:
    /// vrp = xi_ref + (1 / (omega * dt_remaining)) * (xi_ref - xi)
    /// (simplified: first-order DCM tracking.)
    pub fn required_vrp(&self, xi_ref: [f64; 2], _omega: f64, _kp: f64) -> [f64; 2] {
        // DCM tracking: vrp = xi + (1/omega) * omega * (xi - xi_ref)
        // = xi + (xi - xi_ref) = 2*xi - xi_ref (proportional deadbeat)
        // More stable: vrp = xi_ref - (1/omega) * kp * (xi_ref - xi)
        // For deadbeat: kp = omega -> vrp = xi
        [
            self.xi[0] + (self.xi[0] - xi_ref[0]),
            self.xi[1] + (self.xi[1] - xi_ref[1]),
        ]
    }
}

/// Plan a DCM trajectory for a sequence of footstep positions.
///
/// Uses backward recursion: the final DCM target equals the last footstep
/// position, and each preceding DCM waypoint is computed from the next via
/// the exponential decay formula.
///
/// Returns DCM waypoints (one per footstep) in forward order.
pub fn plan_dcm_trajectory(
    footstep_positions: &[[f64; 2]],
    step_durations: &[f64],
    omega: f64,
) -> Vec<[f64; 2]> {
    let n = footstep_positions.len();
    if n == 0 {
        return Vec::new();
    }
    let mut dcm = vec![[0.0; 2]; n];
    // Terminal DCM = last footstep
    dcm[n - 1] = footstep_positions[n - 1];
    // Backward recursion
    for i in (0..n - 1).rev() {
        let dt = if i < step_durations.len() {
            step_durations[i]
        } else {
            0.5
        };
        let exp_neg_wt = (-omega * dt).exp();
        // xi_i = foot_i + exp(-w*T) * (xi_{i+1} - foot_i)
        dcm[i] = [
            footstep_positions[i][0] + exp_neg_wt * (dcm[i + 1][0] - footstep_positions[i][0]),
            footstep_positions[i][1] + exp_neg_wt * (dcm[i + 1][1] - footstep_positions[i][1]),
        ];
    }
    dcm
}

/// DCM error: Euclidean distance between current and reference DCM.
pub fn dcm_error(current: [f64; 2], reference: [f64; 2]) -> f64 {
    v2_len(current[0] - reference[0], current[1] - reference[1])
}

// ─────────────────────────────────────────────────────────────────────────────
// 11. Push Recovery
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the required stepping location for push recovery.
///
/// If the capture point exits the support polygon, the robot must step to
/// a new location near the capture point. This function offsets the capture
/// point by a margin towards the robot's centre.
pub fn push_recovery_step_target(capture_pt: [f64; 2], com_xy: [f64; 2], margin: f64) -> [f64; 2] {
    let dx = capture_pt[0] - com_xy[0];
    let dy = capture_pt[1] - com_xy[1];
    let dist = v2_len(dx, dy);
    if dist < 1e-15 {
        return capture_pt;
    }
    // Step slightly beyond the capture point in the push direction
    let scale = (dist + margin) / dist;
    [com_xy[0] + dx * scale, com_xy[1] + dy * scale]
}

/// Evaluate whether a push is recoverable without stepping.
///
/// Returns `true` if the capture point falls within the support polygon.
pub fn is_push_recoverable_in_place(
    com_pos: [f64; 3],
    com_vel: [f64; 3],
    omega: f64,
    support_centre: [f64; 2],
    half_widths: [f64; 2],
) -> bool {
    let cp = capture_point(com_pos, com_vel, omega);
    capture_point_feasible(cp, support_centre, half_widths)
}

// ─────────────────────────────────────────────────────────────────────────────
// 12. Gait Pattern Generator
// ─────────────────────────────────────────────────────────────────────────────

/// Gait pattern type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaitPattern {
    /// Walking gait (duty factor > 0.5).
    Walk,
    /// Trotting gait (duty factor ~ 0.5, diagonal pairs).
    Trot,
    /// Running gait (duty factor < 0.5, flight phase).
    Run,
}

/// Parameters for a periodic bipedal gait.
#[derive(Debug, Clone, Copy)]
pub struct BipedalGaitParams {
    /// Full gait cycle duration \[s\].
    pub cycle_duration: f64,
    /// Fraction of the cycle each foot is on the ground.
    pub duty_factor: f64,
    /// Phase offset of the right foot (0.5 for alternating walk).
    pub right_phase_offset: f64,
}

impl BipedalGaitParams {
    /// Create standard walking parameters.
    pub fn walking(cycle_duration: f64) -> Self {
        Self {
            cycle_duration,
            duty_factor: 0.6,
            right_phase_offset: 0.5,
        }
    }

    /// Create running parameters.
    pub fn running(cycle_duration: f64) -> Self {
        Self {
            cycle_duration,
            duty_factor: 0.35,
            right_phase_offset: 0.5,
        }
    }

    /// Detect the gait pattern from the duty factor.
    pub fn pattern(&self) -> GaitPattern {
        if self.duty_factor > 0.5 {
            GaitPattern::Walk
        } else if self.duty_factor > 0.35 {
            GaitPattern::Trot
        } else {
            GaitPattern::Run
        }
    }

    /// Stride frequency in Hz.
    pub fn stride_frequency(&self) -> f64 {
        if self.cycle_duration.abs() < 1e-15 {
            0.0
        } else {
            1.0 / self.cycle_duration
        }
    }

    /// Query contact state at a given time.
    pub fn contact_at(&self, time: f64) -> (ContactState, ContactState) {
        let phase = gait_phase_fraction(time, self.cycle_duration);
        let (l_sw, r_sw) = bilateral_swing_timing(phase, self.duty_factor, self.right_phase_offset);
        let left = if l_sw {
            ContactState::Swing
        } else {
            ContactState::InContact
        };
        let right = if r_sw {
            ContactState::Swing
        } else {
            ContactState::InContact
        };
        (left, right)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 13. CoP / Wrench Distribution
// ─────────────────────────────────────────────────────────────────────────────

/// Centre of Pressure computation from a set of vertical contact forces.
///
/// Each force is `(position_xy, normal_force)`. Returns the CoP in the
/// horizontal plane.
pub fn compute_cop(forces: &[([f64; 2], f64)]) -> [f64; 2] {
    let mut total_fz = 0.0;
    let mut wx = 0.0;
    let mut wy = 0.0;
    for &(pos, fz) in forces {
        total_fz += fz;
        wx += pos[0] * fz;
        wy += pos[1] * fz;
    }
    if total_fz.abs() < 1e-15 {
        return [0.0, 0.0];
    }
    [wx / total_fz, wy / total_fz]
}

/// Distribute a desired net vertical force and CoP between two feet.
///
/// Returns `(left_fz, right_fz)` such that the net force and moment about
/// the y-axis are satisfied.
pub fn distribute_force_two_feet(
    desired_fz: f64,
    desired_cop_x: f64,
    left_x: f64,
    right_x: f64,
) -> (f64, f64) {
    let dx = right_x - left_x;
    if dx.abs() < 1e-15 {
        return (0.5 * desired_fz, 0.5 * desired_fz);
    }
    let right_fz = desired_fz * (desired_cop_x - left_x) / dx;
    let left_fz = desired_fz - right_fz;
    (left_fz.max(0.0), right_fz.max(0.0))
}

// ─────────────────────────────────────────────────────────────────────────────
// 14. Stepping Controller (simple)
// ─────────────────────────────────────────────────────────────────────────────

/// A simple stepping controller that decides when and where to step.
#[derive(Debug, Clone)]
pub struct SteppingController {
    /// Maximum allowed capture-point distance before a step is triggered.
    pub cp_threshold: f64,
    /// Nominal step duration \[s\].
    pub step_duration: f64,
    /// Foot lateral offset \[m\].
    pub lateral_offset: f64,
    /// Current stepping foot.
    pub next_foot: Foot,
    /// Time of last step.
    pub last_step_time: f64,
}

impl SteppingController {
    /// Create a new stepping controller.
    pub fn new(cp_threshold: f64, step_duration: f64, lateral_offset: f64) -> Self {
        Self {
            cp_threshold,
            step_duration,
            lateral_offset,
            next_foot: Foot::Right,
            last_step_time: 0.0,
        }
    }

    /// Decide whether a step is needed and return the target if so.
    pub fn update(
        &mut self,
        time: f64,
        com_pos: [f64; 3],
        com_vel: [f64; 3],
        omega: f64,
        support_centre: [f64; 2],
        half_widths: [f64; 2],
    ) -> Option<Footstep> {
        let cp = capture_point(com_pos, com_vel, omega);
        let margin = capture_point_margin(cp, support_centre, half_widths);

        if margin < self.cp_threshold && time - self.last_step_time > self.step_duration * 0.5 {
            let y_offset = match self.next_foot {
                Foot::Left => self.lateral_offset,
                Foot::Right => -self.lateral_offset,
            };
            let target = Footstep::new(
                self.next_foot,
                [cp[0], cp[1] + y_offset, 0.0],
                0.0,
                time,
                self.step_duration * 0.6,
                self.step_duration * 0.4,
            );
            self.last_step_time = time;
            self.next_foot = self.next_foot.opposite();
            Some(target)
        } else {
            None
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 15. Walking Speed Estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate comfortable walking speed from leg length using the Froude
/// number relation: v = sqrt(Fr * g * L).
///
/// For comfortable walking, Fr ~ 0.25.
pub fn comfortable_walking_speed(leg_length: f64, gravity: f64, froude: f64) -> f64 {
    (froude * gravity * leg_length).sqrt()
}

/// Estimate step length from walking speed and cadence.
///
/// step_length = speed / cadence  (cadence in steps/s).
pub fn step_length_from_cadence(speed: f64, cadence: f64) -> f64 {
    if cadence.abs() < 1e-15 {
        return 0.0;
    }
    speed / cadence
}

/// Estimate the cost of transport (dimensionless) for walking.
///
/// CoT = P / (m * g * v), where P is metabolic power.
pub fn cost_of_transport(power: f64, mass: f64, gravity: f64, speed: f64) -> f64 {
    let denom = mass * gravity * speed;
    if denom.abs() < 1e-15 {
        return 0.0;
    }
    power / denom
}

// ─────────────────────────────────────────────────────────────────────────────
// 16. Angular Momentum Regulation
// ─────────────────────────────────────────────────────────────────────────────

/// Desired angular momentum rate for balance: a PD controller on angular
/// momentum about the CoM.
///
/// tau_des = -kp * L - kd * L_dot
///
/// Where `L` is the current angular momentum and `L_dot` its rate.
pub fn angular_momentum_pd(
    ang_mom: [f64; 3],
    ang_mom_rate: [f64; 3],
    kp: f64,
    kd: f64,
) -> [f64; 3] {
    [
        -kp * ang_mom[0] - kd * ang_mom_rate[0],
        -kp * ang_mom[1] - kd * ang_mom_rate[1],
        -kp * ang_mom[2] - kd * ang_mom_rate[2],
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// 17. Terrain-Aware Stepping
// ─────────────────────────────────────────────────────────────────────────────

/// Adjust a footstep position for a sloped terrain.
///
/// Given a terrain gradient `[dz/dx, dz/dy]` and a nominal footstep
/// position, project the foot onto the terrain surface.
pub fn terrain_adjusted_footstep(
    nominal: [f64; 3],
    terrain_height_at_xy: f64,
    terrain_normal: [f64; 3],
) -> [f64; 3] {
    let _n = v3_norm(terrain_normal);
    [nominal[0], nominal[1], terrain_height_at_xy]
}

/// Compute the maximum safe step inclination angle from a terrain normal.
///
/// Returns the angle in radians between the terrain normal and vertical.
pub fn terrain_inclination(terrain_normal: [f64; 3]) -> f64 {
    let n = v3_norm(terrain_normal);
    let cos_angle = n[2].abs(); // dot with [0,0,1]
    cos_angle.clamp(-1.0, 1.0).acos()
}

/// Check if terrain is too steep to walk on.
pub fn terrain_walkable(terrain_normal: [f64; 3], max_slope_rad: f64) -> bool {
    terrain_inclination(terrain_normal) <= max_slope_rad
}

// ─────────────────────────────────────────────────────────────────────────────
// 18. Foot Force Distribution
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the ground reaction force for static standing.
///
/// For a biped standing with both feet on the ground, the total GRF must
/// equal `m*g` vertically and the ZMP must be at the CoM projection.
pub fn standing_grf(mass: f64, gravity: f64) -> f64 {
    mass * gravity
}

/// Friction cone check: ensure the tangential force does not exceed
/// mu * normal_force.
pub fn friction_cone_valid(tangential_force: [f64; 2], normal_force: f64, mu: f64) -> bool {
    let ft = v2_len(tangential_force[0], tangential_force[1]);
    ft <= mu * normal_force + 1e-10
}

// ─────────────────────────────────────────────────────────────────────────────
// 19. LIPM Preview Control
// ─────────────────────────────────────────────────────────────────────────────

/// Simple LIPM preview controller: compute CoM jerk to track a ZMP reference.
///
/// Uses the preview control relation:
/// u = -Kx * state + sum_i(Kp_i * zmp_ref_i)
///
/// This simplified version uses a proportional-derivative approach:
/// x_ddot = omega^2 * (x - zmp_ref) + kd * xdot
pub fn lipm_preview_accel(x: f64, xdot: f64, zmp_ref: f64, omega: f64, kd: f64) -> f64 {
    omega * omega * (x - zmp_ref) - kd * xdot
}

/// Compute the required ZMP to achieve a desired CoM acceleration.
///
/// zmp = x - x_ddot / omega^2
pub fn required_zmp_for_accel(x: f64, x_ddot: f64, omega: f64) -> f64 {
    if omega.abs() < 1e-15 {
        return x;
    }
    x - x_ddot / (omega * omega)
}

// ─────────────────────────────────────────────────────────────────────────────
// 20. Miscellaneous Utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the inverted pendulum energy at a given state.
///
/// E = 0.5 * m * v^2 - m * g * h * cos(theta)
///
/// For small angles: E ~ 0.5 * m * (v^2 - g * h * theta^2)
pub fn pendulum_energy(mass: f64, velocity: f64, gravity: f64, height: f64, theta: f64) -> f64 {
    0.5 * mass * velocity * velocity - mass * gravity * height * theta.cos()
}

/// Natural period of a simple pendulum: T = 2*pi*sqrt(L/g).
pub fn pendulum_period(length: f64, gravity: f64) -> f64 {
    if gravity.abs() < 1e-15 || length < 0.0 {
        return 0.0;
    }
    2.0 * PI * (length / gravity).sqrt()
}

/// Cadence (steps per minute) from cycle duration.
pub fn cadence_from_cycle(cycle_duration: f64) -> f64 {
    if cycle_duration.abs() < 1e-15 {
        return 0.0;
    }
    // Two steps per cycle for bipedal walking
    120.0 / cycle_duration
}

/// Stride length from step length (one full gait cycle = 2 steps).
pub fn stride_length(step_length: f64) -> f64 {
    2.0 * step_length
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // ── LIPM tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_lipm_omega() {
        let p = LipmParams::new(1.0, 9.81, 70.0);
        let w = p.omega();
        assert!((w - 9.81_f64.sqrt()).abs() < 1e-6, "omega={w}");
    }

    #[test]
    fn test_lipm_step_stationary() {
        let p = LipmParams::new(0.8, 9.81, 70.0);
        let s = LipmState2d::new(0.0, 0.0);
        let s2 = s.step(&p, 0.01);
        assert!(s2.x.abs() < EPS);
        assert!(s2.xdot.abs() < EPS);
    }

    #[test]
    fn test_lipm_step_forward() {
        let p = LipmParams::new(0.8, 9.81, 70.0);
        let s = LipmState2d::new(0.1, 0.5);
        let s2 = s.step(&p, 0.01);
        // Position should increase with positive velocity
        assert!(s2.x > s.x, "x should increase");
    }

    #[test]
    fn test_lipm_orbital_energy() {
        let p = LipmParams::new(1.0, 9.81, 70.0);
        let s = LipmState2d::new(0.0, 1.0);
        let e = s.orbital_energy(&p);
        // E = 0.5 * (1.0 - 0.0) = 0.5
        assert!((e - 0.5).abs() < EPS);
    }

    #[test]
    fn test_lipm_3d_step() {
        let p = LipmParams::new(0.8, 9.81, 70.0);
        let s = LipmState3d::new([0.1, 0.0, 0.8], [0.5, 0.0, 0.0]);
        let s2 = s.step(&p, 0.01);
        assert!(s2.pos[0] > s.pos[0]);
        assert!((s2.pos[2] - 0.8).abs() < EPS);
    }

    // ── ZMP tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_zmp_stationary() {
        let zmp = compute_zmp([0.0, 0.0, 0.8], [0.0, 0.0, 0.0], 9.81);
        assert!(zmp[0].abs() < EPS);
        assert!(zmp[1].abs() < EPS);
    }

    #[test]
    fn test_zmp_with_accel() {
        let zmp = compute_zmp([1.0, 0.0, 1.0], [2.0, 0.0, 0.0], 9.81);
        // zmp_x = 1.0 - (1.0/9.81)*2.0 ~ 0.796
        let expected = 1.0 - (1.0 / 9.81) * 2.0;
        assert!((zmp[0] - expected).abs() < 1e-6, "zmp_x={}", zmp[0]);
    }

    #[test]
    fn test_zmp_in_support() {
        let zmp = [0.05, 0.02, 0.0];
        assert!(zmp_in_support_rect(zmp, [0.0, 0.0], [0.1, 0.05]));
        assert!(!zmp_in_support_rect(zmp, [0.0, 0.0], [0.01, 0.01]));
    }

    #[test]
    fn test_zmp_trajectory() {
        let wps = vec![([0.0, 0.0], 1.0), ([1.0, 0.0], 1.0), ([2.0, 0.0], 1.0)];
        let p = zmp_trajectory_at(&wps, 0.5);
        assert!((p[0] - 0.5).abs() < 1e-6);
        let p2 = zmp_trajectory_at(&wps, 1.5);
        assert!((p2[0] - 1.5).abs() < 1e-6);
    }

    #[test]
    fn test_zmp_stability_margin() {
        let margin = zmp_stability_margin([0.0, 0.0, 0.0], [0.0, 0.0], [0.1, 0.05]);
        assert!((margin - 0.05).abs() < EPS);
    }

    // ── Footstep tests ──────────────────────────────────────────────────────

    #[test]
    fn test_footstep_planning() {
        let steps = plan_straight_footsteps([0.0, 0.0, 0.0], 0.3, 0.1, 0.4, 0.1, 4, Foot::Left);
        assert_eq!(steps.len(), 4);
        assert_eq!(steps[0].foot, Foot::Left);
        assert_eq!(steps[1].foot, Foot::Right);
        assert!((steps[0].position[0] - 0.3).abs() < EPS);
        assert!((steps[0].position[1] - 0.1).abs() < EPS);
        assert!((steps[1].position[1] + 0.1).abs() < EPS);
    }

    #[test]
    fn test_footstep_timing() {
        let fs = Footstep::new(Foot::Left, [0.0; 3], 0.0, 1.0, 0.4, 0.1);
        assert!((fs.total_duration() - 0.5).abs() < EPS);
        assert!((fs.swing_end_time() - 1.4).abs() < EPS);
        assert!((fs.end_time() - 1.5).abs() < EPS);
    }

    #[test]
    fn test_double_support_centre() {
        let c = double_support_centre([0.0, 0.1, 0.0], [0.0, -0.1, 0.0]);
        assert!(c[0].abs() < EPS);
        assert!(c[1].abs() < EPS);
    }

    // ── Swing trajectory tests ──────────────────────────────────────────────

    #[test]
    fn test_swing_bezier_endpoints() {
        let start = [0.0, 0.0, 0.0];
        let end = [0.3, 0.0, 0.0];
        let p0 = swing_bezier_trajectory(start, end, 0.05, 0.0);
        let p1 = swing_bezier_trajectory(start, end, 0.05, 1.0);
        for i in 0..3 {
            assert!((p0[i] - start[i]).abs() < 1e-6, "p0[{i}]={}", p0[i]);
            assert!((p1[i] - end[i]).abs() < 1e-6, "p1[{i}]={}", p1[i]);
        }
    }

    #[test]
    fn test_swing_bezier_apex() {
        let p = swing_bezier_trajectory([0.0; 3], [0.3, 0.0, 0.0], 0.05, 0.5);
        assert!(p[2] > 0.01, "mid-swing z={} should be elevated", p[2]);
    }

    #[test]
    fn test_swing_parabolic() {
        let p = swing_parabolic_trajectory([0.0; 3], [0.3, 0.0, 0.0], 0.05, 0.5);
        assert!((p[2] - 0.05).abs() < 1e-6, "peak z={}", p[2]);
        let p_end = swing_parabolic_trajectory([0.0; 3], [0.3, 0.0, 0.0], 0.05, 1.0);
        assert!(p_end[2].abs() < 1e-6);
    }

    // ── Contact schedule tests ──────────────────────────────────────────────

    #[test]
    fn test_contact_schedule_from_footsteps() {
        let steps = plan_straight_footsteps([0.0; 3], 0.3, 0.1, 0.4, 0.1, 2, Foot::Left);
        let sched = ContactSchedule::from_footsteps(&steps);
        // 2 steps -> 4 events (swing + contact for each)
        assert_eq!(sched.events.len(), 4);
        // At t=0.2, left foot should be swinging
        assert_eq!(sched.state_at(Foot::Left, 0.2), Some(ContactState::Swing));
    }

    #[test]
    fn test_contact_schedule_total_duration() {
        let mut sched = ContactSchedule::new();
        sched.push(ContactEvent::new(Foot::Left, ContactState::Swing, 0.0, 0.4));
        sched.push(ContactEvent::new(
            Foot::Left,
            ContactState::InContact,
            0.4,
            0.1,
        ));
        assert!((sched.total_duration() - 0.5).abs() < EPS);
    }

    // ── Gait phase tests ────────────────────────────────────────────────────

    #[test]
    fn test_gait_phase_detection() {
        assert_eq!(
            detect_gait_phase(ContactState::InContact, ContactState::InContact),
            GaitPhase::DoubleSupport
        );
        assert_eq!(
            detect_gait_phase(ContactState::InContact, ContactState::Swing),
            GaitPhase::LeftStance
        );
        assert_eq!(
            detect_gait_phase(ContactState::Swing, ContactState::Swing),
            GaitPhase::Flight
        );
    }

    #[test]
    fn test_gait_phase_fraction() {
        let f = gait_phase_fraction(0.5, 1.0);
        assert!((f - 0.5).abs() < EPS);
        let f2 = gait_phase_fraction(1.5, 1.0);
        assert!((f2 - 0.5).abs() < EPS);
    }

    // ── Dynamic balance tests ───────────────────────────────────────────────

    #[test]
    fn test_compute_com() {
        let bodies = vec![([0.0, 0.0, 1.0], 10.0), ([2.0, 0.0, 1.0], 10.0)];
        let com = compute_com(&bodies);
        assert!((com[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_linear_momentum() {
        let bodies = vec![([1.0, 0.0, 0.0], 5.0), ([-1.0, 0.0, 0.0], 5.0)];
        let p = compute_linear_momentum(&bodies);
        assert!(p[0].abs() < EPS);
    }

    #[test]
    fn test_centroidal_angular_momentum() {
        // Two masses moving in opposite directions equidistant from origin
        let bodies = vec![
            ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0),
            ([-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 1.0),
        ];
        let l = compute_centroidal_angular_momentum(&bodies, [0.0; 3]);
        // Both contribute positive angular momentum about z
        assert!(l[2] > 0.0, "Lz={}", l[2]);
    }

    // ── Capture point tests ─────────────────────────────────────────────────

    #[test]
    fn test_capture_point_stationary() {
        let cp = capture_point([0.0, 0.0, 0.8], [0.0; 3], 3.0);
        assert!(cp[0].abs() < EPS);
        assert!(cp[1].abs() < EPS);
    }

    #[test]
    fn test_capture_point_moving() {
        let omega = 3.0;
        let cp = capture_point([0.0, 0.0, 0.8], [1.0, 0.0, 0.0], omega);
        assert!((cp[0] - 1.0 / omega).abs() < EPS);
    }

    #[test]
    fn test_capture_point_feasibility() {
        assert!(capture_point_feasible(
            [0.05, 0.02],
            [0.0, 0.0],
            [0.1, 0.05]
        ));
        assert!(!capture_point_feasible([0.2, 0.0], [0.0, 0.0], [0.1, 0.05]));
    }

    // ── DCM tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_dcm_from_com() {
        let dcm = DcmState::from_com([0.0, 0.0, 0.8], [1.0, 0.0, 0.0], 3.0);
        assert!((dcm.xi[0] - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_dcm_propagate() {
        let dcm = DcmState { xi: [0.1, 0.0] };
        let vrp = [0.0, 0.0];
        let d2 = dcm.propagate(vrp, 3.0, 0.01);
        // Should diverge slightly from VRP
        assert!(d2.xi[0] > dcm.xi[0]);
    }

    #[test]
    fn test_dcm_plan_trajectory() {
        let feet = vec![[0.0, 0.0], [0.3, 0.1], [0.6, -0.1]];
        let durs = vec![0.5, 0.5, 0.5];
        let dcm_traj = plan_dcm_trajectory(&feet, &durs, 3.0);
        assert_eq!(dcm_traj.len(), 3);
        // Last DCM = last footstep
        assert!((dcm_traj[2][0] - 0.6).abs() < EPS);
    }

    #[test]
    fn test_dcm_error() {
        let e = dcm_error([0.1, 0.0], [0.0, 0.0]);
        assert!((e - 0.1).abs() < EPS);
    }

    // ── Gait pattern tests ──────────────────────────────────────────────────

    #[test]
    fn test_walking_gait_params() {
        let gp = BipedalGaitParams::walking(1.0);
        assert_eq!(gp.pattern(), GaitPattern::Walk);
        assert!((gp.stride_frequency() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_running_gait_params() {
        let gp = BipedalGaitParams::running(0.8);
        assert_eq!(gp.pattern(), GaitPattern::Run);
    }

    // ── CoP / force distribution tests ──────────────────────────────────────

    #[test]
    fn test_compute_cop() {
        let forces = vec![([0.0, 0.0], 50.0), ([1.0, 0.0], 50.0)];
        let cop = compute_cop(&forces);
        assert!((cop[0] - 0.5).abs() < EPS);
    }

    #[test]
    fn test_distribute_force() {
        let (lf, rf) = distribute_force_two_feet(100.0, 0.0, -0.1, 0.1);
        assert!((lf - 50.0).abs() < EPS);
        assert!((rf - 50.0).abs() < EPS);
    }

    // ── Misc tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_pendulum_period() {
        let t = pendulum_period(1.0, 9.81);
        let expected = 2.0 * PI * (1.0 / 9.81_f64).sqrt();
        assert!((t - expected).abs() < 1e-6);
    }

    #[test]
    fn test_comfortable_walking_speed() {
        let v = comfortable_walking_speed(0.9, 9.81, 0.25);
        // v ~ sqrt(0.25 * 9.81 * 0.9) ~ 1.485
        assert!(v > 1.0 && v < 2.0, "v={v}");
    }

    #[test]
    fn test_friction_cone() {
        assert!(friction_cone_valid([1.0, 0.0], 10.0, 0.5));
        assert!(!friction_cone_valid([6.0, 0.0], 10.0, 0.5));
    }

    #[test]
    fn test_terrain_inclination_flat() {
        let inc = terrain_inclination([0.0, 0.0, 1.0]);
        assert!(inc.abs() < 1e-6);
    }

    #[test]
    fn test_terrain_walkable() {
        assert!(terrain_walkable([0.0, 0.0, 1.0], 0.5));
        // 45-degree slope normal
        let n = [0.0, 0.707, 0.707];
        assert!(!terrain_walkable(n, 0.3));
    }

    #[test]
    fn test_cadence_from_cycle() {
        let c = cadence_from_cycle(1.0);
        assert!((c - 120.0).abs() < EPS);
    }

    #[test]
    fn test_max_step_length() {
        let ml = max_step_length(0.9, 0.5);
        assert!((ml - 0.45).abs() < EPS);
    }

    #[test]
    fn test_push_recovery_target() {
        let target = push_recovery_step_target([0.3, 0.0], [0.0, 0.0], 0.05);
        // Should be beyond the capture point
        assert!(target[0] > 0.3);
    }

    #[test]
    fn test_lipm_preview_accel() {
        let a = lipm_preview_accel(0.1, 0.0, 0.0, 3.0, 0.0);
        // omega^2 * (0.1 - 0.0) = 9 * 0.1 = 0.9
        assert!((a - 0.9).abs() < EPS);
    }

    #[test]
    fn test_required_zmp() {
        let zmp = required_zmp_for_accel(1.0, 0.0, 3.0);
        assert!((zmp - 1.0).abs() < EPS);
    }

    #[test]
    fn test_centroidal_momentum_error() {
        let current = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let desired = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let err = centroidal_momentum_error(current, desired);
        assert!((err[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_standing_grf() {
        let grf = standing_grf(70.0, 9.81);
        assert!((grf - 686.7).abs() < 0.1);
    }
}
