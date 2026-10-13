// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Continuous collision detection (CCD) by time-of-impact sub-stepping.
//!
//! The contacts of a step are detected and solved before anything moves, with
//! or without CCD (see the pipeline documentation), so resting contacts
//! elsewhere in the scene are handled in every frame whether or not a fast body
//! produces a CCD event, and the bodies are moved with their *solved*
//! velocities. A CCD step then does not integrate the world in one go. It
//! repeatedly
//!
//! 1. predicts where every collider will be at the end of the *remaining*
//!    interval (position **and** orientation),
//! 2. looks for the earliest time of impact (TOI) among the fast pairs,
//! 3. advances the whole world to that time and resolves the impact with a
//!    normal impulse,
//!
//! until no event is left or the sub-step budget of the pipeline
//! ([`PhysicsPipeline::max_ccd_substeps`]) is used up. The caller integrates
//! what remains of the step.
//!
//! Only pairs that are still approaching at their time of impact count as
//! events, so a pair that was just resolved does not trigger again, while two
//! walls a bullet ping-pongs between can each be hit several times per step.
//! A pair is "fast" if its relative displacement over the *whole* step reaches
//! [`CCD_SPEED_THRESHOLD`], however short the rest of the step has become: a
//! bullet that has just bounced must still be caught by thin geometry in what
//! is left of the step.
//!
//! Bodies that follow a scripted kinematic target are treated as stationary
//! while sub-stepping and reach their target in the final integration of the
//! step.

use std::cmp::Ordering;
use std::collections::HashSet;
use std::time::{Duration, Instant};

use oxiphysics_collision::BroadPhase;
use oxiphysics_collision::ccd::{ToiResult, time_of_impact};
use oxiphysics_core::math::{Quat, Real, Vec3};
use oxiphysics_core::{Aabb, Transform};
use oxiphysics_rigid::{BodyState, BodyType, Collider, RigidBody, RigidBodySet};

use super::PhysicsPipeline;
use super::collision::{ColliderEntry, pair_needs_resolution, shares_body, world_aabb};
use super::solver::{Impact, resolve_impact};

/// Minimum relative displacement (m per step) for CCD to be attempted on a pair.
const CCD_SPEED_THRESHOLD: Real = 0.5;

/// Remaining step length (s) below which no further CCD sub-step is attempted.
const MIN_CCD_INTERVAL: Real = 1e-9;

/// Approach speed (m/s) a pair must exceed to count as closing. It keeps the
/// round-off of a contact the solver has just brought to rest (a fast body
/// sliding along a floor) from producing events.
const MIN_APPROACH_SPEED: Real = 1e-6;

/// Default number of time-of-impact events resolved per step.
pub(super) const DEFAULT_MAX_CCD_SUBSTEPS: usize = 8;

/// A time-of-impact event between two colliders.
struct CcdEvent {
    /// Time of impact within the interval and contact geometry.
    toi: ToiResult,
    /// Index of collider A in the entry list.
    a: usize,
    /// Index of collider B in the entry list.
    b: usize,
}

/// Result of one swept broad phase and time-of-impact search.
struct CcdDetection {
    events: Vec<CcdEvent>,
    broadphase_time: Duration,
    toi_time: Duration,
}

/// Whether the body is driven by a scripted kinematic target.
pub(super) fn follows_target(body: &RigidBody) -> bool {
    body.body_type == BodyType::Kinematic && body.kinematic_target.is_some()
}

/// Linear and angular velocity CCD predicts for a body; bodies that do not move
/// by their velocity (static, sleeping, scripted) are predicted to stand still.
fn swept_motion(body: Option<&RigidBody>) -> (Vec3, Vec3) {
    match body {
        Some(b)
            if b.body_type != BodyType::Static
                && b.state != BodyState::Sleeping
                && !follows_target(b) =>
        {
            (b.velocity, b.angular_velocity)
        }
        _ => (Vec3::zeros(), Vec3::zeros()),
    }
}

/// Speed (m/s) at which A approaches B along `normal` (pointing from B toward
/// A; negative while approaching) at `point`, `elapsed` seconds from now.
fn approach_speed(
    body_a: Option<&RigidBody>,
    body_b: Option<&RigidBody>,
    point: &Vec3,
    normal: &Vec3,
    elapsed: Real,
) -> Real {
    let point_velocity = |body: Option<&RigidBody>| -> Vec3 {
        let (velocity, angular_velocity) = swept_motion(body);
        match body {
            Some(b) => {
                let arm = point - (b.transform.position + velocity * elapsed);
                velocity + angular_velocity.cross(&arm)
            }
            None => velocity,
        }
    };
    (point_velocity(body_a) - point_velocity(body_b)).dot(normal)
}

/// Transform of `body` after `dt` of constant linear and angular velocity.
///
/// The orientation is advanced exactly like
/// [`RigidBody::integrate_velocity`] does it (rotation by `|w| * dt` about the
/// axis of the angular velocity), so the returned pose is the pose the
/// integrator will actually reach after `dt`.
pub(super) fn predicted_transform(body: &RigidBody, dt: Real) -> Transform {
    if body.body_type == BodyType::Static
        || body.state == BodyState::Sleeping
        || follows_target(body)
    {
        return body.transform.clone();
    }
    let position = body.transform.position + body.velocity * dt;
    let rotation = if body.angular_velocity.norm() > 1e-10 {
        Quat::from_scaled_axis(body.angular_velocity * dt) * body.transform.rotation
    } else {
        body.transform.rotation
    };
    Transform::new(position, rotation)
}

/// World-space bound of a collider for the swept broad phase: the union of its
/// (rotation-aware) boxes at `start` and at `end`.
///
/// This bounds a translation exactly, and a rotation as long as the shape is
/// not wider in between than at either end, which holds for the few degrees per
/// step of a normal simulation. A shape that turns by a large part of a quarter
/// turn within one step can reach slightly beyond the bound.
pub(super) fn swept_aabb(collider: &Collider, start: &Transform, end: &Transform) -> Aabb {
    world_aabb(collider, start).merge(&world_aabb(collider, end))
}

/// Integrate every body that moves by its velocity over `dt`.
fn advance_bodies(bodies: &mut RigidBodySet, dt: Real) {
    for (_, body) in bodies.iter_mut() {
        if follows_target(body) {
            continue;
        }
        body.integrate_velocity(dt);
    }
}

impl PhysicsPipeline {
    /// Swept broad phase plus time-of-impact search over the next `dt`.
    ///
    /// `dt` is the part of the step that is still to be integrated and
    /// `step_dt` the length of the whole step. A pair is only worth a
    /// time-of-impact search if it moves fast over the *whole* step
    /// ([`CCD_SPEED_THRESHOLD`]); judging it by what is left of the step would
    /// let a bullet that has just bounced tunnel through thin geometry in the
    /// remainder, because the remainder is too short to count as "fast".
    fn detect_ccd_events(
        &self,
        dt: Real,
        step_dt: Real,
        bodies: &RigidBodySet,
        entries: &[ColliderEntry<'_>],
        skipped: &HashSet<(usize, usize)>,
    ) -> CcdDetection {
        let t_broadphase = Instant::now();

        // World transform of every collider at the start and at the end of the
        // interval.
        let motion: Vec<(Transform, Transform)> = entries
            .iter()
            .map(|entry| {
                let start = entry.world_transform(bodies);
                let end = match entry.body(bodies) {
                    Some(body) => entry
                        .collider
                        .world_transform(&predicted_transform(body, dt)),
                    None => start.clone(),
                };
                (start, end)
            })
            .collect();

        // Swept AABB = union of the (rotation-aware) boxes at start and end.
        let swept: Vec<Aabb> = entries
            .iter()
            .zip(&motion)
            .map(|(entry, (start, end))| swept_aabb(entry.collider, start, end))
            .collect();
        let pairs = self.broadphase.find_pairs(&swept);
        let broadphase_time = t_broadphase.elapsed();

        let t_toi = Instant::now();
        let mut events = Vec::new();
        for pair in &pairs {
            if skipped.contains(&(pair.a, pair.b)) {
                continue;
            }
            let (Some(entry_a), Some(entry_b)) = (entries.get(pair.a), entries.get(pair.b)) else {
                continue;
            };
            let (Some((a_start, a_end)), Some((b_start, b_end))) =
                (motion.get(pair.a), motion.get(pair.b))
            else {
                continue;
            };
            if entry_a.collider.is_sensor
                || entry_b.collider.is_sensor
                || shares_body(entry_a.collider, entry_b.collider)
            {
                continue;
            }
            let body_a = entry_a.body(bodies);
            let body_b = entry_b.body(bodies);
            if !pair_needs_resolution(body_a, body_b) {
                continue;
            }
            let relative_velocity = swept_motion(body_a).0 - swept_motion(body_b).0;
            if relative_velocity.norm() * step_dt < CCD_SPEED_THRESHOLD {
                continue;
            }
            if let Some(toi) = time_of_impact(
                entry_a.collider.shape.as_ref(),
                a_start,
                a_end,
                entry_b.collider.shape.as_ref(),
                b_start,
                b_end,
            ) && toi.toi < 1.0
            {
                // Only an impact at which the pair is still closing counts; a
                // pair that is touching but separating (e.g. one that was just
                // resolved) is left to the contact solver.
                let point = (toi.witness_a + toi.witness_b) * 0.5;
                let normal = -toi.normal;
                if approach_speed(body_a, body_b, &point, &normal, toi.toi * dt)
                    < -MIN_APPROACH_SPEED
                {
                    events.push(CcdEvent {
                        toi,
                        a: pair.a,
                        b: pair.b,
                    });
                }
            }
        }

        CcdDetection {
            events,
            broadphase_time,
            toi_time: t_toi.elapsed(),
        }
    }

    /// Resolve the times of impact of fast pairs within a step of length `dt`
    /// (forces are integrated and the contacts are solved already).
    ///
    /// Resolves up to [`max_ccd_substeps`](PhysicsPipeline::max_ccd_substeps)
    /// time-of-impact events, each by advancing the whole world to the time of
    /// impact and applying an impulse there. Returns the length of the part of
    /// the step that is still to be integrated.
    pub(super) fn resolve_ccd_events(
        &mut self,
        dt: Real,
        bodies: &mut RigidBodySet,
        entries: &[ColliderEntry<'_>],
    ) -> Real {
        let mut remaining = dt;
        let mut skipped: HashSet<(usize, usize)> = HashSet::new();
        let mut events_resolved = 0_usize;

        for _ in 0..self.max_ccd_substeps {
            if remaining <= MIN_CCD_INTERVAL {
                break;
            }
            let detection = self.detect_ccd_events(remaining, dt, bodies, entries, &skipped);
            self.record_and_fire("broadphase", detection.broadphase_time);
            self.record_and_fire("narrowphase", detection.toi_time);

            let Some(event) = detection
                .events
                .into_iter()
                .min_by(|x, y| x.toi.toi.partial_cmp(&y.toi.toi).unwrap_or(Ordering::Equal))
            else {
                break;
            };
            let advance = event.toi.toi.clamp(0.0, 1.0) * remaining;

            // Advance the whole world to the time of impact ...
            let t_integration = Instant::now();
            advance_bodies(bodies, advance);
            self.record_and_fire("integration", t_integration.elapsed());

            // ... and resolve the impact there.
            let t_solver = Instant::now();
            let mut applied = false;
            if let (Some(entry_a), Some(entry_b)) = (entries.get(event.a), entries.get(event.b)) {
                let impact = Impact {
                    body_a: entry_a.collider.body_handle,
                    body_b: entry_b.collider.body_handle,
                    point: (event.toi.witness_a + event.toi.witness_b) * 0.5,
                    // The TOI normal points from A to B; contacts point from B to A.
                    normal: -event.toi.normal,
                    restitution: (entry_a.collider.restitution + entry_b.collider.restitution)
                        * 0.5,
                };
                applied = resolve_impact(bodies, &impact);
            }
            self.record_and_fire("solver", t_solver.elapsed());

            if !applied {
                // The pair could not be resolved (nothing to push); do not
                // spend more of the sub-step budget on it.
                skipped.insert((event.a, event.b));
            }
            remaining = (remaining - advance).max(0.0);
            events_resolved += 1;
        }

        self.stats.last_ccd_events = events_resolved;
        remaining
    }
}
