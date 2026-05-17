// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Motion interpolation — smooth damp, exponential decay, spring followers,
//! and scalar / vector utility functions.
//!
//! ## Functions
//!
//! - `smooth_damp` / `smooth_damp3` — Unity-style critically-damped
//!   smoothing with velocity output and optional `max_speed` clamp.
//! - `exp_decay` / `exp_decay3` — Simple exponential decay (no velocity
//!   tracking; cheaper than smooth damp).
//! - `lerp`, `lerp3`, `remap`, `smoothstep`, `smootherstep` —
//!   scalar and vector utilities.
//!
//! ## Types
//!
//! - `SpringFollower` — critically-damped spring that tracks a 1D target.
//! - `SpringFollower3` — same for a `[f64; 3]` target.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::interpolator::{smooth_damp, SpringFollower};
//!
//! let mut vel = 0.0_f64;
//! let out = smooth_damp(0.0, 10.0, &mut vel, 0.2, f64::INFINITY, 0.016);
//! assert!(out > 0.0 && out < 10.0);
//!
//! let mut spring = SpringFollower::new(100.0);
//! let pos = spring.update(5.0, 0.016);
//! assert!(pos > 0.0);
//! ```

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Scalar utilities
// ---------------------------------------------------------------------------

/// Linear interpolation: `a + t*(b−a)`.
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Linear interpolation for `[f64; 3]` vectors.
#[inline]
pub fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        lerp(a[0], b[0], t),
        lerp(a[1], b[1], t),
        lerp(a[2], b[2], t),
    ]
}

/// Remap `value` from range `[in_min, in_max]` to `[out_min, out_max]`.
///
/// Does **not** clamp the input, so values outside `[in_min, in_max]`
/// extrapolate linearly.
#[inline]
pub fn remap(value: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    let t = if in_max == in_min {
        0.0
    } else {
        (value - in_min) / (in_max - in_min)
    };
    lerp(out_min, out_max, t)
}

/// Cubic Hermite smoothstep: `3t² − 2t³`, clamped input.
///
/// Returns 0 for `x <= edge0`, 1 for `x >= edge1`.
#[inline]
pub fn smoothstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Quintic smootherstep: `6t⁵ − 15t⁴ + 10t³` (zero first and second
/// derivatives at 0 and 1).
#[inline]
pub fn smootherstep(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

// ---------------------------------------------------------------------------
// Exponential decay
// ---------------------------------------------------------------------------

/// Exponential decay toward `target`: `current + (target − current)(1 − e^(−rate·dt))`.
///
/// No velocity tracking — cheap and appropriate when velocity continuity is
/// not required.
#[inline]
pub fn exp_decay(current: f64, target: f64, rate: f64, dt: f64) -> f64 {
    let t = (-rate * dt).exp();
    target + (current - target) * t
}

/// Exponential decay for `[f64; 3]`.
#[inline]
pub fn exp_decay3(current: [f64; 3], target: [f64; 3], rate: f64, dt: f64) -> [f64; 3] {
    [
        exp_decay(current[0], target[0], rate, dt),
        exp_decay(current[1], target[1], rate, dt),
        exp_decay(current[2], target[2], rate, dt),
    ]
}

// ---------------------------------------------------------------------------
// Smooth damp (Unity-style critically-damped follower)
// ---------------------------------------------------------------------------

/// Unity-style smooth damp — critically-damped spring with velocity continuity.
///
/// Moves `current` toward `target` over approximately `smooth_time` seconds.
/// The caller must maintain `velocity` across frames.
///
/// `max_speed` clamps the maximum displacement per call.  Pass
/// `f64::INFINITY` for no limit.
///
/// Returns the new `current` value and writes the updated velocity to
/// `velocity`.
///
/// # Implementation
///
/// Uses a Padé approximant of the exponential `e^(−ω·dt)` so that it
/// degrades gracefully with large `dt` rather than exploding.
pub fn smooth_damp(
    current: f64,
    target: f64,
    velocity: &mut f64,
    smooth_time: f64,
    max_speed: f64,
    dt: f64,
) -> f64 {
    let smooth_time = smooth_time.max(1e-6);
    let omega = 2.0 / smooth_time;

    // Padé approximant: e^(-omega*dt) ≈ (1 - k*omega*dt) for stability
    let x = omega * dt;
    // exp approximation: 1 / (1 + x + 0.5x² + x³/6 + x⁴/24)
    let exp_approx = 1.0 / (1.0 + x + 0.48 * x * x + 0.235 * x * x * x);

    let original_to = target;
    let max_change = max_speed * smooth_time;
    let delta = (current - target).clamp(-max_change, max_change);
    let target_adj = current - delta;

    let temp = (*velocity + omega * delta) * dt;
    *velocity = (*velocity - omega * temp) * exp_approx;
    let output = target_adj + (delta + temp) * exp_approx;

    // Prevent overshoot
    if (original_to - current > 0.0) == (output > original_to) {
        *velocity = (output - original_to) / dt;
        return original_to;
    }

    output
}

/// Smooth damp for `[f64; 3]` vectors.
///
/// `velocity` is a mutable `[f64; 3]` that must be maintained by the caller.
pub fn smooth_damp3(
    current: [f64; 3],
    target: [f64; 3],
    velocity: &mut [f64; 3],
    smooth_time: f64,
    max_speed: f64,
    dt: f64,
) -> [f64; 3] {
    [
        smooth_damp(
            current[0],
            target[0],
            &mut velocity[0],
            smooth_time,
            max_speed,
            dt,
        ),
        smooth_damp(
            current[1],
            target[1],
            &mut velocity[1],
            smooth_time,
            max_speed,
            dt,
        ),
        smooth_damp(
            current[2],
            target[2],
            &mut velocity[2],
            smooth_time,
            max_speed,
            dt,
        ),
    ]
}

// ---------------------------------------------------------------------------
// SpringFollower
// ---------------------------------------------------------------------------

/// Critically-damped spring follower for a 1D scalar.
///
/// Maintains velocity state internally.  Call [`update`][Self::update] each
/// frame with the current target and `dt`.
///
/// `damping` is set automatically to `2 * sqrt(stiffness)` (critical damping)
/// so the system converges without oscillation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpringFollower {
    /// Spring stiffness coefficient (units: 1/s²).
    pub stiffness: f64,
    /// Damping coefficient (set to `2 * sqrt(stiffness)` by [`new`][Self::new]).
    pub damping: f64,
    /// Current simulated value.
    pub value: f64,
    velocity: f64,
}

impl SpringFollower {
    /// Create a critically-damped spring starting at `value = 0`.
    ///
    /// `stiffness` must be positive; larger values make the spring stiffer
    /// and faster to respond.
    pub fn new(stiffness: f64) -> Self {
        let stiffness = stiffness.max(0.0);
        Self {
            stiffness,
            damping: 2.0 * stiffness.sqrt(),
            value: 0.0,
            velocity: 0.0,
        }
    }

    /// Create with explicit stiffness and damping.
    pub fn with_damping(stiffness: f64, damping: f64) -> Self {
        Self {
            stiffness: stiffness.max(0.0),
            damping: damping.max(0.0),
            value: 0.0,
            velocity: 0.0,
        }
    }

    /// Set the current value without affecting velocity (teleport).
    pub fn set_value(&mut self, value: f64) {
        self.value = value;
    }

    /// Reset both value and velocity.
    pub fn reset(&mut self, value: f64) {
        self.value = value;
        self.velocity = 0.0;
    }

    /// Advance the spring toward `target` by `dt` seconds.
    ///
    /// Returns the new value.
    pub fn update(&mut self, target: f64, dt: f64) -> f64 {
        let error = self.value - target;
        let accel = -self.stiffness * error - self.damping * self.velocity;
        self.velocity += accel * dt;
        self.value += self.velocity * dt;
        self.value
    }
}

// ---------------------------------------------------------------------------
// SpringFollower3
// ---------------------------------------------------------------------------

/// Critically-damped spring follower for a `[f64; 3]` vector.
///
/// Each axis is independent.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpringFollower3 {
    /// Spring stiffness coefficient (units: 1/s²).
    pub stiffness: f64,
    /// Damping coefficient (`2 * sqrt(stiffness)` by default).
    pub damping: f64,
    /// Current simulated position.
    pub value: [f64; 3],
    velocity: [f64; 3],
}

impl SpringFollower3 {
    /// Create a critically-damped spring starting at the origin.
    pub fn new(stiffness: f64) -> Self {
        let stiffness = stiffness.max(0.0);
        Self {
            stiffness,
            damping: 2.0 * stiffness.sqrt(),
            value: [0.0; 3],
            velocity: [0.0; 3],
        }
    }

    /// Create with explicit stiffness and damping.
    pub fn with_damping(stiffness: f64, damping: f64) -> Self {
        Self {
            stiffness: stiffness.max(0.0),
            damping: damping.max(0.0),
            value: [0.0; 3],
            velocity: [0.0; 3],
        }
    }

    /// Set the current value without affecting velocity.
    pub fn set_value(&mut self, value: [f64; 3]) {
        self.value = value;
    }

    /// Reset value and velocity.
    pub fn reset(&mut self, value: [f64; 3]) {
        self.value = value;
        self.velocity = [0.0; 3];
    }

    /// Advance the spring toward `target` by `dt` seconds.
    ///
    /// Returns the new value.
    pub fn update(&mut self, target: [f64; 3], dt: f64) -> [f64; 3] {
        for ((val, vel), tgt) in self
            .value
            .iter_mut()
            .zip(self.velocity.iter_mut())
            .zip(target.iter())
        {
            let error = *val - tgt;
            let accel = -self.stiffness * error - self.damping * *vel;
            *vel += accel * dt;
            *val += *vel * dt;
        }
        self.value
    }
}
