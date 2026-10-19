// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conservative advancement for bodies with bounded linear and angular motion.
//!
//! Conservative advancement (Mirtich 1996) finds the time of impact of two
//! moving bodies without solving their motion exactly: if the bodies are a
//! distance `gap` apart and no point of either body can approach a point of the
//! other faster than `closing_speed`, then they cannot touch within the next
//! `gap / closing_speed` seconds, so time can safely be advanced by that
//! amount. Repeating until the gap is below a tolerance converges on the first
//! contact from below and therefore never tunnels through it.

use super::{add, len, scale, sub};

/// A rigid body state for conservative advancement.
#[derive(Debug, Clone)]
pub struct RotatingBodyState {
    /// Center of mass position at the start of the advanced interval.
    pub position: [f64; 3],
    /// Linear velocity.
    pub linear_vel: [f64; 3],
    /// Angular speed (magnitude of angular velocity vector), radians/s.
    pub angular_speed: f64,
    /// Bounding sphere radius (enclosing the entire body).
    pub bound_radius: f64,
}

impl RotatingBodyState {
    /// Create a new rotating body state.
    pub fn new(
        position: [f64; 3],
        linear_vel: [f64; 3],
        angular_speed: f64,
        bound_radius: f64,
    ) -> Self {
        Self {
            position,
            linear_vel,
            angular_speed,
            bound_radius,
        }
    }

    /// Upper bound on the speed of any point on this body.
    ///
    /// A point at distance `r <= bound_radius` from the centre of mass moves
    /// with the centre of mass and turns around it, so its speed is at most
    /// `|v_linear| + |angular_speed| * bound_radius`.
    pub fn max_point_speed(&self) -> f64 {
        len(self.linear_vel) + self.angular_speed.abs() * self.bound_radius
    }

    /// Upper bound on how far any point on this body can move in time `dt`.
    ///
    /// This is a distance: `(|v_linear| + |angular_speed| * bound_radius) * dt`
    /// (see [`RotatingBodyState::max_point_speed`] for the rate).
    pub fn motion_bound(&self, dt: f64) -> f64 {
        self.max_point_speed() * dt
    }

    /// Advance position by `dt`.
    pub fn advance(&self, dt: f64) -> [f64; 3] {
        add(self.position, scale(self.linear_vel, dt))
    }
}

/// Upper bound on the rate at which the distance between any point of `a` and
/// any point of `b` can shrink.
///
/// The relative velocity of two points is the difference of the centres' linear
/// velocities plus the rotational surface velocities `|angular_speed| * r` of
/// both points, so its norm is at most the returned value.
fn closing_speed_bound(a: &RotatingBodyState, b: &RotatingBodyState) -> f64 {
    len(sub(a.linear_vel, b.linear_vel))
        + a.angular_speed.abs() * a.bound_radius
        + b.angular_speed.abs() * b.bound_radius
}

/// Distance between the bounding spheres of the two bodies `dt` after the start
/// of the interval (negative while the spheres overlap).
fn bounding_sphere_gap(a: &RotatingBodyState, b: &RotatingBodyState, dt: f64) -> f64 {
    len(sub(b.advance(dt), a.advance(dt))) - (a.bound_radius + b.bound_radius)
}

/// Conservative advancement loop shared by the public entry points.
///
/// `initial_gap` is the separation of the bodies at `t_start`; `gap_at(t)` is
/// their separation at the absolute time `t`. Both must not exceed the true
/// distance of the shapes, otherwise the advancement can overshoot the contact.
fn advance_to_contact<F>(
    closing_speed: f64,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
    max_iter: usize,
    initial_gap: f64,
    gap_at: F,
) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let tolerance = tolerance.max(0.0);
    let mut t = t_start;
    let mut gap = initial_gap;
    let mut steps = 0usize;
    loop {
        if gap.is_nan() {
            return None;
        }
        if gap <= tolerance {
            return Some(t);
        }
        let remaining = t_end - t;
        if remaining <= 0.0 || closing_speed <= 0.0 {
            // The interval is used up, or the bodies cannot approach each other.
            return None;
        }
        if steps >= max_iter {
            // Contact before `t` is excluded; whether the bodies touch later is
            // undecided, so report the time reached as a lower bound.
            return Some(t);
        }
        // Rate bound: no point pair closes `gap` faster than `closing_speed`
        // allows, so nothing can touch before `gap / closing_speed` has passed.
        let advance = (gap / closing_speed).min(remaining);
        t = if advance >= remaining {
            t_end
        } else {
            t + advance
        };
        gap = gap_at(t);
        steps += 1;
    }
}

/// Conservative advancement time of impact for two rotating bodies, using
/// their bounding spheres.
///
/// The positions of `a` and `b` are those at `t_start`; both bodies move with
/// their constant linear velocities (per unit of time) and turn with their
/// angular speeds, which only enter through the bound on how fast their
/// surfaces can approach each other. Time is advanced from `t_start` towards
/// `t_end` in steps of `gap / closing_speed`, where `gap` is the current
/// distance between the bounding spheres and `closing_speed` bounds how fast
/// any point of one body can approach any point of the other (the relative
/// linear speed plus the rotational surface speeds `|angular_speed| *
/// bound_radius` of both bodies). No point pair can close the gap within one
/// step, so the advancement never steps past a contact.
///
/// `separation` is the caller's distance between the actual shapes at
/// `t_start` (for example from a GJK query), negative when they penetrate. The
/// first step uses the smaller of it and the bounding-sphere gap.
///
/// Returns the absolute time `t` in `[t_start, t_end]` (not a fraction of the
/// interval) at which the gap between the bounding spheres has shrunk to at
/// most `tolerance`; the spheres do not touch before it. `Some(t_start)` is
/// returned when the bodies already overlap or `separation` is negative, and
/// `None` when the bodies cannot come within `tolerance` of each other during
/// the interval. If `max_iter` steps are not enough to converge (grazing
/// approaches converge slowly) the time reached so far is returned: contact
/// before it is excluded, so the result is a lower bound of the time of impact.
///
/// Bounding spheres overestimate the shapes, so the returned time is early for
/// shapes that do not fill their sphere; use
/// [`conservative_advancement_rotating_with_distance`] to advance on the true
/// separation of the shapes.
pub fn conservative_advancement_rotating(
    a: &RotatingBodyState,
    b: &RotatingBodyState,
    separation: f64,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
    max_iter: usize,
) -> Option<f64> {
    let total_dt = t_end - t_start;
    if total_dt.is_nan() || total_dt <= 0.0 {
        return None;
    }
    let sphere_gap = bounding_sphere_gap(a, b, 0.0);
    if sphere_gap <= 0.0 || separation < 0.0 {
        return Some(t_start);
    }
    advance_to_contact(
        closing_speed_bound(a, b),
        t_start,
        t_end,
        tolerance,
        max_iter,
        sphere_gap.min(separation),
        |t| bounding_sphere_gap(a, b, t - t_start),
    )
}

/// Conservative advancement time of impact for two rotating bodies, advancing
/// on the true separation of their shapes.
///
/// Like [`conservative_advancement_rotating`], but the gap is taken from
/// `shape_distance(t)` instead of the bounding spheres. `shape_distance` returns
/// the separation of the two actual shapes at the absolute time `t` in
/// `[t_start, t_end]` (the poses at `t` follow from the bodies' constant
/// velocities; the callback typically runs a GJK distance query on them) and is
/// negative or zero when they touch or penetrate. It must not overestimate the
/// distance, or the advancement can step past the contact; an exact distance
/// gives the tightest time of impact.
///
/// Returns the absolute time at which the separation has shrunk to at most
/// `tolerance` (`Some(t_start)` if it already is), or `None` when the shapes
/// cannot come that close during the interval. If `max_iter` steps do not
/// suffice, the time reached so far is returned as a lower bound of the time of
/// impact.
pub fn conservative_advancement_rotating_with_distance<F>(
    a: &RotatingBodyState,
    b: &RotatingBodyState,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
    max_iter: usize,
    shape_distance: F,
) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let total_dt = t_end - t_start;
    if total_dt.is_nan() || total_dt <= 0.0 {
        return None;
    }
    advance_to_contact(
        closing_speed_bound(a, b),
        t_start,
        t_end,
        tolerance,
        max_iter,
        shape_distance(t_start),
        shape_distance,
    )
}
