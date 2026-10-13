// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Position-based (XPBD) enforcement of the constraints stored in a
//! [`PhysicsWorld`].
//!
//! The solver follows *Detailed Rigid Body Simulation with Extended Position
//! Based Dynamics* (Müller, Macklin, Chentanez, Jeschke, Kim, 2020). Once the
//! bodies have been integrated for a substep, every supported constraint is
//! reduced to a scalar constraint function `C` that is corrected along a unit
//! direction `n`:
//!
//! * [`ConstraintKind::BallSocket`]: `C = |p_a - p_b|`;
//! * the distance kinds: `C = |p_a - p_b| - target`, where
//!   [`ConstraintKind::DistanceMin`] only ever pushes the anchors apart and
//!   [`ConstraintKind::DistanceMax`] only ever pulls them together;
//!
//! with `n = (p_a - p_b) / |p_a - p_b|`, `p_a = x_a + R_a * anchor_a` and
//! `p_b = x_b + R_b * anchor_b` (or the world point `anchor_b` when the
//! constraint has no second body).
//!
//! For every Gauss-Seidel iteration the multiplier increment is
//!
//! ```text
//! dlambda = (-C - alpha * lambda) / (w_a + w_b + alpha),     alpha = compliance / dt^2
//! w_i     = 1/m_i + |r_i x n|^2 / I_i                        (generalised inverse mass)
//! ```
//!
//! and the positional impulse `p = dlambda * n` translates body `a` by
//! `p / m_a`, body `b` by `-p / m_b` and rotates them by `(r_a x p) / I_a` and
//! `-(r_b x p) / I_b`, so an anchor that is offset from the centre of mass
//! produces the proper coupling between translation and rotation. The
//! multipliers start from zero in every substep, which makes `compliance` a
//! physical, time-step independent softness (the inverse of a spring
//! stiffness). Finally the solver converts the applied displacement and
//! rotation into velocity changes (`v += dx / dt`, `w += dtheta / dt`), so
//! constrained bodies keep a velocity consistent with where the constraint
//! left them.
//!
//! Only [`BodyMode::Dynamic`] bodies are moved; static, kinematic and
//! sleeping bodies act as immovable references. Inertia is the scalar
//! `Body::inertia`, in line with the rest of the lightweight world.

use super::types::{Body, BodyArena, BodyMode, BodySlot, ConstraintKind, PhysicsWorld};
use crate::math::{Quat, Real, Vec3};
use crate::types::BodyHandle;
use nalgebra::Quaternion;

/// Anchor separations shorter than this count as "coincident".
const COINCIDENT_EPSILON: Real = 1e-12;
/// Constraints whose effective inverse mass is below this cannot be solved.
const MIN_EFFECTIVE_INVERSE_MASS: Real = 1e-12;

/// Private working copy of one body taking part in the solve.
struct WorkingBody {
    /// Slot of the body inside the arena.
    slot: usize,
    position: Vec3,
    rotation: Quat,
    /// Inverse mass; zero for immovable bodies.
    inv_mass: Real,
    /// Inverse (scalar) inertia; zero for immovable bodies.
    inv_inertia: Real,
    /// Total translation the solver applied in this substep.
    shift: Vec3,
    /// Total world-frame rotation vector the solver applied in this substep.
    turn: Vec3,
}

impl WorkingBody {
    fn from_body(slot: usize, body: &Body) -> Self {
        let dynamic = body.mode == BodyMode::Dynamic;
        let inv_mass = if dynamic && body.inv_mass.is_finite() && body.inv_mass > 0.0 {
            body.inv_mass
        } else {
            0.0
        };
        let inv_inertia = if dynamic && body.inertia.is_finite() && body.inertia > 0.0 {
            1.0 / body.inertia
        } else {
            0.0
        };
        Self {
            slot,
            position: body.transform.position,
            rotation: body.transform.rotation,
            inv_mass,
            inv_inertia,
            shift: Vec3::zeros(),
            turn: Vec3::zeros(),
        }
    }

    /// Whether any correction can move this body.
    fn is_movable(&self) -> bool {
        self.inv_mass > 0.0 || self.inv_inertia > 0.0
    }

    /// Translate by `shift` and rotate by the world-frame rotation vector `turn`.
    fn apply(&mut self, shift: Vec3, turn: Vec3) {
        self.position += shift;
        self.shift += shift;
        if turn.norm_squared() > 0.0 {
            let q = self.rotation.quaternion();
            let dq = Quaternion::from_parts(0.0, turn * 0.5) * q;
            self.rotation = Quat::from_quaternion(q + dq);
            self.turn += turn;
        }
    }

    /// Whether the working pose is free of NaN and infinities.
    fn is_finite(&self) -> bool {
        self.position.iter().all(|c| c.is_finite())
            && self
                .rotation
                .quaternion()
                .coords
                .iter()
                .all(|c| c.is_finite())
    }
}

/// A constraint resolved against the working bodies of the current solve.
struct WorkingConstraint {
    kind: ConstraintKind,
    /// Index into the working bodies.
    a: usize,
    /// Index into the working bodies; `None` anchors body `a` to the world.
    b: Option<usize>,
    /// Anchor in the local frame of body `a`.
    anchor_a: Vec3,
    /// Anchor in the local frame of body `b`, or a world point without body `b`.
    anchor_b: Vec3,
    /// Target separation of the anchors (zero for ball-sockets).
    target: Real,
    /// Time-step scaled compliance `compliance / dt^2`.
    alpha: Real,
    /// Accumulated XPBD multiplier of this substep.
    lambda: Real,
}

/// Find (or create) the working copy of the body behind `handle`.
///
/// Returns `None` when the handle is stale.
fn gather_body(
    bodies: &BodyArena,
    handle: BodyHandle,
    local_of_slot: &mut [Option<usize>],
    working: &mut Vec<WorkingBody>,
) -> Option<usize> {
    let body = bodies.get(handle)?;
    let slot = handle.index as usize;
    let entry = local_of_slot.get_mut(slot)?;
    if let Some(index) = *entry {
        return Some(index);
    }
    working.push(WorkingBody::from_body(slot, body));
    let index = working.len() - 1;
    *entry = Some(index);
    Some(index)
}

/// Apply one Gauss-Seidel update of `constraint` to the working bodies.
fn solve_constraint(constraint: &mut WorkingConstraint, bodies: &mut [WorkingBody]) {
    let (arm_a, point_a, centre_a, inv_mass_a, inv_inertia_a) = {
        let body = &bodies[constraint.a];
        let arm = body.rotation * constraint.anchor_a;
        (
            arm,
            body.position + arm,
            body.position,
            body.inv_mass,
            body.inv_inertia,
        )
    };
    let (arm_b, point_b, centre_b, inv_mass_b, inv_inertia_b) = match constraint.b {
        Some(index) => {
            let body = &bodies[index];
            let arm = body.rotation * constraint.anchor_b;
            (
                arm,
                body.position + arm,
                body.position,
                body.inv_mass,
                body.inv_inertia,
            )
        }
        None => (
            Vec3::zeros(),
            constraint.anchor_b,
            constraint.anchor_b,
            0.0,
            0.0,
        ),
    };

    let separation = point_a - point_b;
    let distance = separation.norm();
    let direction = if distance > COINCIDENT_EPSILON {
        separation / distance
    } else if constraint.kind == ConstraintKind::BallSocket {
        // The anchors already coincide: there is nothing to correct and no
        // direction to correct it along.
        return;
    } else {
        // Coincident anchors under a distance constraint: separate along the
        // line between the body centres, or an arbitrary axis if those
        // coincide as well.
        (centre_a - centre_b)
            .try_normalize(COINCIDENT_EPSILON)
            .unwrap_or_else(Vec3::x)
    };

    let generalised_inverse_mass_a =
        inv_mass_a + inv_inertia_a * arm_a.cross(&direction).norm_squared();
    let generalised_inverse_mass_b =
        inv_mass_b + inv_inertia_b * arm_b.cross(&direction).norm_squared();
    let denominator = generalised_inverse_mass_a + generalised_inverse_mass_b + constraint.alpha;
    if denominator.is_nan() || denominator <= MIN_EFFECTIVE_INVERSE_MASS {
        return;
    }

    let violation = distance - constraint.target;
    let raw_lambda =
        constraint.lambda + (-violation - constraint.alpha * constraint.lambda) / denominator;
    // One-sided constraints may only push (min) or pull (max); clamping the
    // accumulated multiplier also lets them release an earlier overshoot.
    let new_lambda = match constraint.kind {
        ConstraintKind::DistanceMin => raw_lambda.max(0.0),
        ConstraintKind::DistanceMax => raw_lambda.min(0.0),
        _ => raw_lambda,
    };
    let delta_lambda = new_lambda - constraint.lambda;
    constraint.lambda = new_lambda;

    let impulse = direction * delta_lambda;
    bodies[constraint.a].apply(impulse * inv_mass_a, arm_a.cross(&impulse) * inv_inertia_a);
    if let Some(index) = constraint.b {
        bodies[index].apply(
            -impulse * inv_mass_b,
            -arm_b.cross(&impulse) * inv_inertia_b,
        );
    }
}

impl PhysicsWorld {
    /// Enforce the stored constraints after the bodies were integrated.
    ///
    /// Runs [`PhysicsConfig::solver_iterations`](crate::types::PhysicsConfig)
    /// Gauss-Seidel sweeps of the XPBD update described in the module
    /// documentation and converts the resulting position and orientation
    /// corrections into velocity changes. Inactive constraints, constraints of
    /// a kind for which [`ConstraintKind::is_solved_by_step`] is `false`,
    /// constraints naming a stale body handle (or the same body twice) and
    /// constraints with non-finite data are skipped.
    pub(super) fn solve_constraints(&mut self, dt: Real) {
        if self.constraints.is_empty() || !dt.is_finite() || dt <= 0.0 {
            return;
        }
        let mut local_of_slot: Vec<Option<usize>> = vec![None; self.bodies.slots.len()];
        let mut working: Vec<WorkingBody> = Vec::new();
        let mut resolved: Vec<WorkingConstraint> = Vec::with_capacity(self.constraints.len());
        let inv_dt_sq = 1.0 / (dt * dt);

        for constraint in &self.constraints {
            if !constraint.active
                || !constraint.kind.is_solved_by_step()
                || constraint.body_b == Some(constraint.body_a)
            {
                continue;
            }
            let anchor_a = Vec3::from(constraint.anchor_a);
            let anchor_b = Vec3::from(constraint.anchor_b);
            let data_is_finite = anchor_a.iter().all(|c| c.is_finite())
                && anchor_b.iter().all(|c| c.is_finite())
                && constraint.target.is_finite()
                && constraint.compliance.is_finite();
            if !data_is_finite {
                continue;
            }
            let Some(a) = gather_body(
                &self.bodies,
                constraint.body_a,
                &mut local_of_slot,
                &mut working,
            ) else {
                continue;
            };
            let b = match constraint.body_b {
                Some(handle) => {
                    match gather_body(&self.bodies, handle, &mut local_of_slot, &mut working) {
                        Some(index) => Some(index),
                        None => continue,
                    }
                }
                None => None,
            };
            let b_is_movable = b.is_some_and(|index| working[index].is_movable());
            if !working[a].is_movable() && !b_is_movable {
                continue;
            }
            resolved.push(WorkingConstraint {
                kind: constraint.kind,
                a,
                b,
                anchor_a,
                anchor_b,
                target: match constraint.kind {
                    ConstraintKind::BallSocket => 0.0,
                    _ => constraint.target.max(0.0),
                },
                alpha: constraint.compliance.max(0.0) * inv_dt_sq,
                lambda: 0.0,
            });
        }

        for _ in 0..self.config.solver_iterations.max(1) {
            for constraint in &mut resolved {
                solve_constraint(constraint, &mut working);
            }
        }

        let inv_dt = 1.0 / dt;
        for body_state in &working {
            let moved =
                body_state.shift.norm_squared() > 0.0 || body_state.turn.norm_squared() > 0.0;
            if !moved || !body_state.is_finite() {
                continue;
            }
            if let Some(BodySlot::Occupied { body, .. }) =
                self.bodies.slots.get_mut(body_state.slot)
            {
                body.transform.position = body_state.position;
                body.transform.rotation = body_state.rotation;
                for (v, d) in body.velocity.linear.iter_mut().zip(body_state.shift.iter()) {
                    *v += d * inv_dt;
                }
                for (w, d) in body.velocity.angular.iter_mut().zip(body_state.turn.iter()) {
                    *w += d * inv_dt;
                }
            }
        }
    }
}
