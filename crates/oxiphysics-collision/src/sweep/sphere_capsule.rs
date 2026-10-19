// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact swept sphere vs static capsule time of impact.

use super::point_cast::{capsule_entry, closest_point_on_segment};
use super::{SweepResult, add, dot, len, len_sq, normalize, scale, sub};

/// Contact normal (from the capsule towards the sphere) for a sphere centre
/// that lies on the capsule axis, where the offset to the axis has no
/// direction: the direction opposing the motion perpendicular to the axis, or
/// else any direction perpendicular to the axis.
fn perpendicular_fallback(axis: [f64; 3], motion: [f64; 3]) -> [f64; 3] {
    let axis_dir = normalize(axis);
    let perpendicular = |v: [f64; 3]| sub(v, scale(axis_dir, dot(v, axis_dir)));
    let against_motion = perpendicular(scale(motion, -1.0));
    // Motion (anti)parallel to the axis leaves only rounding noise here.
    if len_sq(against_motion) > 1e-24 * len_sq(motion) {
        return normalize(against_motion);
    }
    let mut best = [1.0, 0.0, 0.0];
    let mut best_len_sq = -1.0;
    for unit in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
        let candidate = perpendicular(unit);
        let candidate_len_sq = len_sq(candidate);
        if candidate_len_sq > best_len_sq {
            best_len_sq = candidate_len_sq;
            best = candidate;
        }
    }
    normalize(best)
}

/// Swept sphere vs static capsule time-of-impact.
///
/// The sphere of radius `r_s` moves with its centre from `c0_s` to `c1_s`. The
/// capsule is static with axis segment `(cap0, cap1)` and radius `r_c`. A
/// sphere touches the capsule exactly when its centre comes within
/// `R = r_s + r_c` of the axis segment, so the time of impact is the first time
/// the centre path enters the capsule of radius `R`: the centre starting within
/// `R` of the axis (`toi = 0`), the centre reaching the cylinder of radius `R`
/// around the axis at a point between the end points, or reaching one of the
/// two end-point spheres of radius `R`. Each case is solved exactly and the
/// earliest is returned, so a path that crosses the capsule anywhere along its
/// length is detected.
///
/// The result's `normal` is the unit vector from the capsule towards the sphere
/// (the capsule is shape B, the sphere shape A), `point_a` is the contact point
/// on the sphere surface and `point_b` the contact point on the capsule
/// surface. A zero-length axis makes the capsule a sphere. Returns `None` when
/// the sphere never touches the capsule during the step, and for non-finite
/// input or a negative radius.
pub fn swept_sphere_capsule(
    c0_s: [f64; 3],
    c1_s: [f64; 3],
    r_s: f64,
    cap0: [f64; 3],
    cap1: [f64; 3],
    r_c: f64,
) -> Option<SweepResult> {
    let finite = c0_s
        .iter()
        .chain(c1_s.iter())
        .chain(cap0.iter())
        .chain(cap1.iter())
        .chain([r_s, r_c].iter())
        .all(|x| x.is_finite());
    if !finite || r_s < 0.0 || r_c < 0.0 {
        return None;
    }
    let motion = sub(c1_s, c0_s);
    let toi = capsule_entry(c0_s, motion, cap0, cap1, r_s + r_c)?;

    let centre = add(c0_s, scale(motion, toi));
    let (on_axis, _) = closest_point_on_segment(cap0, cap1, centre);
    let offset = sub(centre, on_axis);
    let dist = len(offset);
    let normal = if dist > 1e-12 * (1.0 + r_s + r_c) {
        scale(offset, 1.0 / dist)
    } else {
        perpendicular_fallback(sub(cap1, cap0), motion)
    };
    Some(SweepResult {
        toi,
        normal,
        point_a: sub(centre, scale(normal, r_s)),
        point_b: add(on_axis, scale(normal, r_c)),
    })
}
