// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conservative advancement of a translating convex shape towards a static one.
//!
//! Let `g(t)` be the gap between the moving shape A, translated by `t * vel`,
//! and the static shape B. For convex shapes `g` is a convex function of `t`
//! (the distance between a translating convex set and a fixed one is convex,
//! and subtracting a constant radius keeps it convex). Writing `n` for the unit
//! vector from the closest point of B to the closest point of A, the slope of
//! `g` at `t` is `dot(vel, n)`, and by convexity the tangent line at `t` lies
//! below the whole curve. The Newton step
//!
//! ```text
//! t' = t + (g(t) - target) / closing,    closing = -dot(vel, n)
//! ```
//!
//! therefore never overshoots the first time `g` drops to `target`: no contact
//! can happen before `t'`. Compared with stepping by `g / |vel|`, which only
//! shrinks the gap by the factor `1 - cos(theta)` per step for an approach at
//! angle `theta` to the contact normal, the Newton step reaches a planar
//! contact in a single iteration and converges quadratically otherwise. If
//! `closing <= 0` the shapes are not approaching, and by convexity they never
//! will: the cast misses.

use super::{ShapeCastConfig, ShapeCastResult, dot3, len3, scale3, sub3};

/// Upper bound on the number of advancement steps of one cast.
///
/// The Newton iteration converges quadratically (transversal contacts) or by
/// halving the remaining distance (tangential contacts), so this is far more
/// than any cast needs; it only guards against rounding-induced stalls.
const MAX_STEPS: usize = 128;

/// Relative size of the smallest separation whose direction is trusted.
///
/// The witness points of a nearly touching pair carry an absolute rounding
/// error of about `f64::EPSILON * scale`; below `DIRECTION_EPS * (1 + scale)`
/// the direction of their difference is noise.
const DIRECTION_EPS: f64 = 1e-10;

/// Closest-feature data of the moving shape A and the static shape B at one
/// instant of the motion.
pub(super) struct Probe {
    /// Gap between the surfaces: positive while separated, non-positive when
    /// the shapes touch or overlap.
    pub gap: f64,
    /// Unit vector from B towards A along the closest features, or `None` when
    /// the closest features coincide (overlap, touching within rounding).
    pub direction: Option<[f64; 3]>,
    /// Direction to report when `direction` is `None` and the advancement has
    /// no earlier direction to fall back on (for example the outward face
    /// normal at a point inside a box).
    pub fallback: Option<[f64; 3]>,
    /// Witness point reported for shape A.
    pub witness_a: [f64; 3],
    /// Witness point reported for shape B.
    pub witness_b: [f64; 3],
}

/// Unit vector of `from_b_to_a = a - b`, or `None` if the two points are closer
/// than the rounding noise of their coordinates.
pub(super) fn separation_direction(a: [f64; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let diff = sub3(a, b);
    let length = len3(diff);
    let scale = a
        .iter()
        .chain(b.iter())
        .fold(0.0_f64, |largest, c| largest.max(c.abs()));
    if length.is_finite() && length > DIRECTION_EPS * (1.0 + scale) {
        Some(scale3(diff, 1.0 / length))
    } else {
        None
    }
}

/// Unit vector along `v`, or `None` for a (numerically) zero vector.
pub(super) fn unit_vector(v: [f64; 3]) -> Option<[f64; 3]> {
    let length = len3(v);
    if length.is_finite() && length > 1e-300 {
        Some(scale3(v, 1.0 / length))
    } else {
        None
    }
}

/// Casts a shape that moves by `t * vel` against a static shape.
///
/// `probe_at(t)` measures the shapes with the moving one displaced by
/// `t * vel`; it returns `None` when no measurement is possible (for example
/// an empty point set), which ends the cast without a hit.
///
/// The result is a hit at the first iterate whose gap is at most
/// `cfg.target_distance + cfg.toi_tolerance`. If the shapes already satisfy that
/// at `t = 0` and actually touch or overlap (`gap <= 0`) the result is flagged
/// as penetrating. The contact normal points from B to A: it is the direction
/// of the closest features, or, when those coincide, the direction at the last
/// iterate that had one, the probe's own fallback, or the direction opposite to
/// the motion.
pub(super) fn advance<F>(vel: [f64; 3], cfg: &ShapeCastConfig, mut probe_at: F) -> ShapeCastResult
where
    F: FnMut(f64) -> Option<Probe>,
{
    let target = cfg.target_distance;
    let hit_gap = target + cfg.toi_tolerance.max(0.0);
    let speed = len3(vel);
    let approach = unit_vector(scale3(vel, -1.0));
    let mut last_direction: Option<[f64; 3]> = None;
    let mut t = 0.0_f64;
    for step in 0..MAX_STEPS {
        let Some(probe) = probe_at(t) else {
            return ShapeCastResult::no_hit(cfg.max_toc);
        };
        if !probe.gap.is_finite() {
            return ShapeCastResult::no_hit(cfg.max_toc);
        }
        let normal = probe
            .direction
            .or(last_direction)
            .or(probe.fallback)
            .or(approach)
            .unwrap_or([0.0, 1.0, 0.0]);
        if probe.gap <= hit_gap {
            return if step == 0 && probe.gap <= 0.0 {
                ShapeCastResult::penetrating(probe.witness_a, probe.witness_b, normal)
            } else {
                ShapeCastResult::hit(t, normal, probe.witness_a, probe.witness_b)
            };
        }
        if probe.direction.is_some() {
            last_direction = probe.direction;
        }
        // An earlier direction can only underestimate the step: by convexity the
        // closing speed never grows as the shapes get closer. Without any
        // direction the speed itself bounds the closing speed.
        let closing = match last_direction {
            Some(direction) => -dot3(vel, direction),
            None => speed,
        };
        if closing.is_nan() || closing <= 0.0 {
            return ShapeCastResult::no_hit(cfg.max_toc);
        }
        t += (probe.gap - target) / closing;
        if t > cfg.max_toc {
            return ShapeCastResult::no_hit(cfg.max_toc);
        }
    }
    ShapeCastResult::no_hit(cfg.max_toc)
}
