// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Kinematic capsule character controller with sweep-and-slide, step-up, and slope handling.
//!
//! This module provides a self-contained character controller that uses a
//! pluggable sweep callback so callers can hook it into any collision system.
//! No sub-crate imports are required.
//!
//! ## Types
//!
//! - `CharacterShape` — capsule geometry (radius + half-height).
//! - `CharacterConfig` — tuning knobs (slope, step offset, skin width, up axis, iterations).
//! - `SweepHit` — result returned by the user-supplied sweep callback.
//! - `CharacterMove` — output of `CharacterController::move_and_slide`.
//! - `CharacterController` — the main controller struct.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::character::{CharacterController, CharacterConfig, CharacterShape, SweepHit};
//!
//! let shape = CharacterShape { radius: 0.4, half_height: 0.9 };
//! let config = CharacterConfig::default();
//!
//! let mut controller = CharacterController::new([0.0, 1.0, 0.0], shape, config);
//!
//! // Move down by 1 unit; caller-provided sweep returns no hit.
//! let result = controller.move_and_slide([0.0, -1.0, 0.0], |_origin, _dir| None);
//! let _ = result.translation;
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Inline math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Project vector `v` onto the plane defined by `normal` (must be unit length).
#[inline]
fn project_onto_plane(v: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let d = dot(v, normal);
    sub(v, scale(normal, d))
}

// ---------------------------------------------------------------------------
// Hit classification
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HitKind {
    Ground,
    Wall,
}

/// Classify a surface hit based on the angle between its normal and the up axis.
///
/// `dot(normal, up) > cos(max_slope_deg)` means the surface is shallow enough
/// to be treated as ground.
fn classify_hit(normal: [f64; 3], up: [f64; 3], max_slope_deg: f64) -> HitKind {
    let cos_threshold = max_slope_deg.to_radians().cos();
    if dot(normal, up) > cos_threshold {
        HitKind::Ground
    } else {
        HitKind::Wall
    }
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Capsule geometry used by the character controller.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CharacterShape {
    /// Capsule radius (metres).
    pub radius: f64,
    /// Half-height of the cylindrical portion (metres).
    pub half_height: f64,
}

/// Configuration knobs for the character controller.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CharacterConfig {
    /// Maximum walkable slope (degrees). Default: 45.0.
    pub max_slope_deg: f64,
    /// Maximum height the controller can step up over (metres). Default: 0.3.
    pub step_offset: f64,
    /// Collision skin distance (metres). Default: 1e-3.
    pub skin_width: f64,
    /// World-space up direction (should be unit length). Default: \[0, 1, 0\].
    pub up_axis: [f64; 3],
    /// Maximum number of slide iterations per call. Default: 4.
    pub max_iterations: usize,
}

impl Default for CharacterConfig {
    fn default() -> Self {
        Self {
            max_slope_deg: 45.0,
            step_offset: 0.3,
            skin_width: 1e-3,
            up_axis: [0.0, 1.0, 0.0],
            max_iterations: 4,
        }
    }
}

/// A surface hit returned by the caller-supplied sweep callback.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepHit {
    /// Time of impact in \[0, 1\] relative to the full sweep direction vector.
    pub toi: f64,
    /// Surface normal at the hit point, pointing away from the surface.
    pub normal: [f64; 3],
}

/// Output of [`CharacterController::move_and_slide`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CharacterMove {
    /// Final resolved displacement (may differ from the requested delta).
    pub translation: [f64; 3],
    /// Whether the controller is standing on a ground surface.
    pub is_grounded: bool,
    /// Normal of the ground surface if grounded.
    pub ground_normal: Option<[f64; 3]>,
    /// Number of sweep iterations performed.
    pub iteration_count: usize,
}

/// Kinematic capsule character controller.
///
/// Maintains position, velocity, and grounding state.  The actual collision
/// queries are delegated to a caller-supplied sweep closure so this module
/// remains independent of any particular collision backend.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CharacterController {
    /// World-space position of the capsule origin.
    pub position: [f64; 3],
    /// Caller-maintained world-space velocity (m/s).
    pub velocity: [f64; 3],
    /// Whether the controller is currently on a ground surface.
    pub is_grounded: bool,
    /// Normal of the ground surface (when grounded).
    pub ground_normal: Option<[f64; 3]>,
    /// Capsule shape.
    pub shape: CharacterShape,
    /// Configuration.
    pub config: CharacterConfig,
}

impl CharacterController {
    /// Create a new character controller at `position`.
    pub fn new(position: [f64; 3], shape: CharacterShape, config: CharacterConfig) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            is_grounded: false,
            ground_normal: None,
            shape,
            config,
        }
    }

    // -----------------------------------------------------------------------
    // Public motion API
    // -----------------------------------------------------------------------

    /// Move the character by `desired_delta` with collision response.
    ///
    /// `sweep(capsule_origin, direction_scaled_by_length) -> Option<SweepHit>`
    ///
    /// The `toi` in the returned hit is in \[0, 1\] relative to the full length
    /// of the direction vector.
    pub fn move_and_slide<F>(&mut self, desired_delta: [f64; 3], sweep: F) -> CharacterMove
    where
        F: Fn([f64; 3], [f64; 3]) -> Option<SweepHit>,
    {
        let initial_pos = self.position;
        // Reset grounding; re-detect from this move.
        self.is_grounded = false;
        self.ground_normal = None;

        let mut remaining = desired_delta;
        let mut iteration_count = 0_usize;

        for _ in 0..self.config.max_iterations {
            let mag = len(remaining);
            if mag < self.config.skin_width {
                break;
            }

            iteration_count += 1;

            match sweep(self.position, remaining) {
                None => {
                    // No obstacle — advance fully and we're done.
                    self.position = add(self.position, remaining);
                    break;
                }
                Some(hit) => {
                    let t = hit.toi.clamp(0.0, 1.0);
                    let normal = hit.normal;

                    // Advance to just before the surface.
                    let skin_ratio = (self.config.skin_width / mag).min(t);
                    let effective_t = (t - skin_ratio).max(0.0);
                    self.position = add(self.position, scale(remaining, effective_t));

                    let kind = classify_hit(normal, self.config.up_axis, self.config.max_slope_deg);
                    match kind {
                        HitKind::Ground => {
                            self.is_grounded = true;
                            self.ground_normal = Some(normal);
                            // Project remaining onto ground plane.
                            remaining = project_onto_plane(remaining, normal);
                        }
                        HitKind::Wall => {
                            // Try to step up over the obstacle.
                            let stepped =
                                try_step_up(self.position, remaining, &self.config, &sweep);
                            if let Some(new_pos) = stepped {
                                self.position = new_pos;
                                // After a successful step-up the remaining
                                // motion is consumed; break out of the loop.
                                break;
                            }
                            // Fall back: project remaining onto wall plane.
                            remaining = project_onto_plane(remaining, normal);
                        }
                    }
                }
            }
        }

        // Ground snap: try to stick to ground over small steps/slopes.
        let snap_result = ground_snap(self.position, &self.config, &sweep);
        if let Some((snapped_pos, snap_normal)) = snap_result {
            self.position = snapped_pos;
            self.is_grounded = true;
            self.ground_normal = Some(snap_normal);
        }

        let translation = sub(self.position, initial_pos);

        CharacterMove {
            translation,
            is_grounded: self.is_grounded,
            ground_normal: self.ground_normal,
            iteration_count,
        }
    }

    /// Set the velocity component along `up_axis` to `speed` and clear grounding.
    ///
    /// Use this to initiate a jump.  The horizontal velocity components are
    /// left unchanged.
    pub fn jump(&mut self, speed: f64) {
        let up = self.config.up_axis;
        let v_up = dot(self.velocity, up);
        // Remove current up component and replace with `speed`.
        self.velocity = add(self.velocity, scale(up, speed - v_up));
        self.is_grounded = false;
        self.ground_normal = None;
    }

    /// Apply gravity acceleration for one time step.
    ///
    /// Adds `gravity_accel * dt` downward (i.e., along `-up_axis`) to velocity.
    pub fn apply_gravity(&mut self, gravity_accel: f64, dt: f64) {
        let down = scale(self.config.up_axis, -1.0);
        self.velocity = add(self.velocity, scale(down, gravity_accel * dt));
    }
}

// ---------------------------------------------------------------------------
// Step-up helper
// ---------------------------------------------------------------------------

/// Attempt to step over an obstacle by sweeping upward, then forward, then
/// snapping back down.
///
/// Returns the new position if the step-up succeeds, or `None` if blocked.
fn try_step_up<F>(
    pos: [f64; 3],
    remaining: [f64; 3],
    config: &CharacterConfig,
    sweep: &F,
) -> Option<[f64; 3]>
where
    F: Fn([f64; 3], [f64; 3]) -> Option<SweepHit>,
{
    let up = config.up_axis;
    let step = config.step_offset;
    let skin = config.skin_width;

    // 1. Sweep upward by step_offset.
    let up_delta = scale(up, step);
    let up_pos = match sweep(pos, up_delta) {
        None => add(pos, up_delta),
        Some(hit) => {
            // Partially cleared the step height.
            let mag = len(up_delta);
            let skin_ratio = (skin / mag).min(hit.toi);
            let t = (hit.toi - skin_ratio).max(0.0);
            if t < 1e-6 {
                return None; // Fully blocked upward.
            }
            add(pos, scale(up_delta, t))
        }
    };

    // 2. Sweep forward (horizontal) from elevated position.
    let fwd_pos = match sweep(up_pos, remaining) {
        None => add(up_pos, remaining),
        Some(_) => return None, // Still blocked after lifting; give up.
    };

    // 3. Snap back down.
    let down_delta = scale(up, -(step + skin));
    match sweep(fwd_pos, down_delta) {
        None => {
            // No floor found; reject step-up (floating).
            None
        }
        Some(hit) => {
            let mag = len(down_delta);
            let skin_ratio = (skin / mag).min(hit.toi);
            let t = (hit.toi - skin_ratio).max(0.0);
            let snapped = add(fwd_pos, scale(down_delta, t));
            // Make sure we didn't land on too-steep a wall.
            if classify_hit(hit.normal, up, config.max_slope_deg) == HitKind::Ground {
                Some(snapped)
            } else {
                None
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Ground-snap helper
// ---------------------------------------------------------------------------

/// After the main slide loop, sweep downward by `step_offset` to stay
/// attached to ground over small slopes and ledge edges.
///
/// Returns `(new_position, ground_normal)` if a suitable ground surface is
/// found within range, or `None` otherwise.
fn ground_snap<F>(
    pos: [f64; 3],
    config: &CharacterConfig,
    sweep: &F,
) -> Option<([f64; 3], [f64; 3])>
where
    F: Fn([f64; 3], [f64; 3]) -> Option<SweepHit>,
{
    let up = config.up_axis;
    let down = scale(up, -1.0);
    let down_delta = scale(down, config.step_offset);
    let mag = len(down_delta);

    let hit = sweep(pos, down_delta)?;
    // Only snap if the surface qualifies as ground.
    if classify_hit(hit.normal, up, config.max_slope_deg) != HitKind::Ground {
        return None;
    }
    let skin_ratio = (config.skin_width / mag).min(hit.toi);
    let t = (hit.toi - skin_ratio).max(0.0);
    let snapped_pos = add(pos, scale(down_delta, t));
    Some((snapped_pos, hit.normal))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    // Helper: build a default controller at the origin.
    fn default_controller() -> CharacterController {
        CharacterController::new(
            [0.0, 0.0, 0.0],
            CharacterShape {
                radius: 0.4,
                half_height: 0.9,
            },
            CharacterConfig::default(),
        )
    }

    // ---------------------------------------------------------------------------
    // 1. Flat floor slide
    // ---------------------------------------------------------------------------
    /// Sweep hits at toi=0.5 with normal=[0,1,0].
    /// Desired delta has both a forward (x) and downward (y) component.
    /// Character should move forward and be detected as grounded.
    #[test]
    fn test_flat_floor_slide() {
        let mut ctrl = default_controller();

        // Desired: move right 1 and down 1.
        let desired = [1.0_f64, -1.0, 0.0];

        // Sweep: first call hits at toi=0.5 (floor), subsequent calls return None.
        let call_count = Cell::new(0_usize);
        let result = ctrl.move_and_slide(desired, |_origin, _dir| {
            let n = call_count.get();
            call_count.set(n + 1);
            if n == 0 {
                // First call: floor hit at toi=0.5
                Some(SweepHit {
                    toi: 0.5,
                    normal: [0.0, 1.0, 0.0],
                })
            } else {
                None
            }
        });

        assert!(result.is_grounded, "Should be grounded after floor hit");
        assert!(result.ground_normal.is_some());

        // After the floor hit (toi=0.5), remaining is projected onto the floor
        // plane (y=0 normal), so only x-motion continues.
        // The character should have moved roughly half the x-distance on the
        // first sub-step, then the remainder in the second.
        // Total x translation should be close to 1.0, y should be close to -0.5.
        let tx = result.translation[0];
        let ty = result.translation[1];
        assert!(tx > 0.4, "Should have moved forward (tx={tx})");
        assert!(ty < 0.0, "Should have moved down to the floor (ty={ty})");
        // After floor hit: advance = 0.5 * |[1,-1,0]| scaled fraction then slide.
        // ty should be around -0.5 (stopped at skin before floor).
        assert!(
            ty > -1.1,
            "Should not have sunk through the floor (ty={ty})"
        );
    }

    // ---------------------------------------------------------------------------
    // 2. Wall stop
    // ---------------------------------------------------------------------------
    /// Sweep hits wall at toi=0.1, normal=[1,0,0].
    /// After projection onto the wall plane, the x-component of remaining
    /// becomes zero (wall faces along x), so horizontal motion is killed.
    #[test]
    fn test_wall_stop() {
        let mut ctrl = default_controller();
        let desired = [1.0_f64, 0.0, 0.0];

        // Wall hit at toi=0.1; subsequent snap returns None.
        let call_count = Cell::new(0_usize);
        let result = ctrl.move_and_slide(desired, |_origin, _dir| {
            let n = call_count.get();
            call_count.set(n + 1);
            if n == 0 {
                Some(SweepHit {
                    toi: 0.1,
                    normal: [1.0, 0.0, 0.0],
                })
            } else {
                None
            }
        });

        // Not grounded (wall hit is not ground).
        assert!(!result.is_grounded, "Wall hit should not be grounded");

        // After a wall hit at toi=0.1 the remaining is projected onto the wall plane.
        // Wall normal [1,0,0] → projection removes x component of remaining.
        // Remaining after projection = [0,0,0], so total x translation is small.
        let tx = result.translation[0];
        assert!(
            tx.abs() < 0.2,
            "Wall hit: x translation should be small (tx={tx})"
        );
    }

    // ---------------------------------------------------------------------------
    // 3. Slope climb (40° slope — under 45° threshold → treated as ground)
    // ---------------------------------------------------------------------------
    #[test]
    fn test_slope_climb() {
        let mut ctrl = default_controller();
        // Desired: move forward and slightly down.
        let desired = [1.0_f64, -0.2, 0.0];

        // 40° slope: normal is [−sin40°, cos40°, 0]
        let angle = 40.0_f64.to_radians();
        let nx = -angle.sin();
        let ny = angle.cos();
        let slope_normal = [nx, ny, 0.0_f64];

        // dot(slope_normal, [0,1,0]) = cos(40°) ≈ 0.766 > cos(45°) ≈ 0.707 → ground
        let call_count = Cell::new(0_usize);
        let result = ctrl.move_and_slide(desired, |_origin, _dir| {
            let n = call_count.get();
            call_count.set(n + 1);
            if n == 0 {
                Some(SweepHit {
                    toi: 0.5,
                    normal: slope_normal,
                })
            } else {
                None
            }
        });

        assert!(result.is_grounded, "40° slope should be treated as ground");
        assert!(
            result.ground_normal.is_some(),
            "ground_normal should be set for 40° slope"
        );
    }

    // ---------------------------------------------------------------------------
    // 4. Slope reject (60° — over 45° threshold → treated as wall)
    // ---------------------------------------------------------------------------
    #[test]
    fn test_slope_reject() {
        let mut ctrl = default_controller();
        let desired = [1.0_f64, 0.0, 0.0];

        // 60° slope: dot with [0,1,0] = cos(60°) = 0.5 < cos(45°) ≈ 0.707 → wall
        let angle = 60.0_f64.to_radians();
        let nx = -angle.sin();
        let ny = angle.cos();
        let slope_normal = [nx, ny, 0.0_f64];

        let call_count = Cell::new(0_usize);
        let result = ctrl.move_and_slide(desired, |_origin, _dir| {
            let n = call_count.get();
            call_count.set(n + 1);
            if n == 0 {
                Some(SweepHit {
                    toi: 0.5,
                    normal: slope_normal,
                })
            } else {
                None
            }
        });

        assert!(
            !result.is_grounded,
            "60° slope should be treated as wall, not ground"
        );
        // x translation is significantly reduced after wall projection.
        let tx = result.translation[0];
        assert!(
            tx < 0.8,
            "Steep slope: x translation should be reduced (tx={tx})"
        );
    }

    // ---------------------------------------------------------------------------
    // 5. Step-up
    // ---------------------------------------------------------------------------
    /// On an initial wall hit, the controller tries to step up.
    /// Upward sweep clears; forward sweep clears; down-snap finds ground.
    /// Controller should advance to the stepped-up position.
    #[test]
    fn test_step_up() {
        let mut ctrl = default_controller();
        // Position slightly below a small step.
        ctrl.position = [0.0, 0.0, 0.0];
        let desired = [0.5_f64, 0.0, 0.0];

        // Call sequence:
        //   0: main sweep → wall hit (step ahead)
        //   1: step-up upward sweep → None (clears)
        //   2: step-up forward sweep → None (clears)
        //   3: step-up down-snap   → ground hit
        //   4: ground snap (post-loop) → None (already snapped in step-up)
        let call_count = Cell::new(0_usize);
        let result = ctrl.move_and_slide(desired, |_origin, dir| {
            let n = call_count.get();
            call_count.set(n + 1);
            match n {
                0 => {
                    // Wall directly ahead at toi=0.1
                    Some(SweepHit {
                        toi: 0.1,
                        normal: [1.0, 0.0, 0.0],
                    })
                }
                1 => None, // upward: clears
                2 => None, // forward (after lift): clears
                3 => {
                    // down-snap: ground below
                    // dir is downward; toi=0.5 → lands at mid-point
                    let _ = dir;
                    Some(SweepHit {
                        toi: 0.5,
                        normal: [0.0, 1.0, 0.0],
                    })
                }
                _ => None,
            }
        });

        // After step-up the character should have moved forward.
        let tx = result.translation[0];
        assert!(
            tx > 0.1,
            "Step-up: controller should have advanced forward (tx={tx})"
        );
    }

    // ---------------------------------------------------------------------------
    // 6. Jump
    // ---------------------------------------------------------------------------
    #[test]
    fn test_jump() {
        let mut ctrl = default_controller();
        ctrl.is_grounded = true;
        ctrl.velocity = [2.0, 0.0, 0.0]; // horizontal velocity before jump

        ctrl.jump(5.0);

        // Up component (y) should be 5.0.
        assert!(
            (ctrl.velocity[1] - 5.0).abs() < 1e-10,
            "After jump, velocity.y should be 5.0, got {}",
            ctrl.velocity[1]
        );
        // Horizontal component preserved.
        assert!(
            (ctrl.velocity[0] - 2.0).abs() < 1e-10,
            "Jump should preserve horizontal velocity, got {}",
            ctrl.velocity[0]
        );
        assert!(!ctrl.is_grounded, "Should not be grounded after jump");
        assert!(ctrl.ground_normal.is_none());
    }

    // ---------------------------------------------------------------------------
    // 7. Serde round-trip
    // ---------------------------------------------------------------------------
    #[test]
    fn test_serde_round_trip() {
        let original = CharacterController {
            position: [1.0, 2.5, -3.0],
            velocity: [0.5, 0.0, 1.2],
            is_grounded: true,
            ground_normal: Some([0.0, 1.0, 0.0]),
            shape: CharacterShape {
                radius: 0.4,
                half_height: 0.9,
            },
            config: CharacterConfig {
                max_slope_deg: 45.0,
                step_offset: 0.3,
                skin_width: 1e-3,
                up_axis: [0.0, 1.0, 0.0],
                max_iterations: 4,
            },
        };

        let json = serde_json::to_string(&original).expect("serialize should succeed");
        let restored: CharacterController =
            serde_json::from_str(&json).expect("deserialize should succeed");

        assert!((restored.position[0] - original.position[0]).abs() < 1e-12);
        assert!((restored.position[1] - original.position[1]).abs() < 1e-12);
        assert!((restored.position[2] - original.position[2]).abs() < 1e-12);
        assert!((restored.velocity[0] - original.velocity[0]).abs() < 1e-12);
        assert_eq!(restored.is_grounded, original.is_grounded);
        assert!(restored.ground_normal.is_some());
        assert!((restored.shape.radius - original.shape.radius).abs() < 1e-12);
        assert!((restored.config.max_slope_deg - original.config.max_slope_deg).abs() < 1e-12);
        assert_eq!(
            restored.config.max_iterations,
            original.config.max_iterations
        );
    }

    // ---------------------------------------------------------------------------
    // 8. No-hit path
    // ---------------------------------------------------------------------------
    #[test]
    fn test_no_hit_path() {
        let mut ctrl = default_controller();
        let desired = [3.0_f64, 0.0, 0.0];

        let result = ctrl.move_and_slide(desired, |_origin, _dir| None);

        // With no hits, character moves the full requested delta.
        assert!(
            (result.translation[0] - 3.0).abs() < 1e-10,
            "No-hit: should move full delta x (tx={})",
            result.translation[0]
        );
        assert!(
            result.translation[1].abs() < 1e-10,
            "No-hit: y translation should be zero"
        );
        assert!(
            result.translation[2].abs() < 1e-10,
            "No-hit: z translation should be zero"
        );
        // Not grounded (no sweep hit anything).
        // (ground snap also returns None when sweep returns None)
        assert!(!result.is_grounded, "No-hit: should not be grounded");
    }

    // ---------------------------------------------------------------------------
    // Additional: apply_gravity
    // ---------------------------------------------------------------------------
    #[test]
    fn test_apply_gravity() {
        let mut ctrl = default_controller();
        ctrl.velocity = [0.0, 0.0, 0.0];
        ctrl.apply_gravity(9.81, 1.0);
        // Down = -up = [0, -1, 0]; gravity_accel * dt = 9.81
        assert!(
            (ctrl.velocity[1] - (-9.81)).abs() < 1e-10,
            "apply_gravity: velocity.y should be -9.81, got {}",
            ctrl.velocity[1]
        );
    }

    // ---------------------------------------------------------------------------
    // Additional: CharacterConfig default values
    // ---------------------------------------------------------------------------
    #[test]
    fn test_config_defaults() {
        let cfg = CharacterConfig::default();
        assert!((cfg.max_slope_deg - 45.0).abs() < 1e-10);
        assert!((cfg.step_offset - 0.3).abs() < 1e-10);
        assert!((cfg.skin_width - 1e-3).abs() < 1e-15);
        assert_eq!(cfg.up_axis, [0.0, 1.0, 0.0]);
        assert_eq!(cfg.max_iterations, 4);
    }
}
