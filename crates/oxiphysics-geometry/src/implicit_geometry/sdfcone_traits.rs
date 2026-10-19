//! # SdfCone - Trait Implementations
//!
//! This module contains trait implementations for `SdfCone`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfCone;

/// Largest half-apex angle the cone SDF works with.  A half-angle of exactly
/// `pi / 2` is a flat disc of infinite radius; the clamp keeps the base radius
/// `height * tan(half_angle)` finite.
const MAX_HALF_ANGLE: f64 = std::f64::consts::FRAC_PI_2 - 1.0e-9;

/// Exact signed distance to a finite right circular cone.
///
/// The apex sits at the origin, the axis is `+Y`, and the base is the disc at
/// `y = height` with radius `height * tan(half_angle)`.  The solid is a body of
/// revolution, so the distance depends only on the meridian-plane coordinates
/// `w = (rho, y)` with `rho = |(x, z)|`, and the boundary seen in that half-plane
/// is the slant segment from the apex to the rim `B = (height * tan(half_angle),
/// height)` followed by the base segment from `(0, height)` to `B` (Quilez'
/// `sdCone`):
///
/// * the magnitude is the smaller of the distances to the two segments;
/// * the sign is negative only strictly inside, i.e. below the slant line
///   (`height * rho - B.rho * y < 0`) and below the base plane (`y < height`).
///
/// A negative `height` is treated as zero (a point at the apex) and a
/// `half_angle` outside `[0, pi / 2)` is clamped to that range (NaN counts as
/// zero).  Such a cone has no interior (the solid is the apex or the segment from
/// the apex to `(0, height, 0)`), so its distance is positive everywhere except on
/// the solid itself.
impl Sdf for SdfCone {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let height = self.height.max(0.0);
        let half_angle = if self.half_angle > 0.0 {
            self.half_angle.min(MAX_HALF_ANGLE)
        } else {
            0.0
        };
        let rim = [height * half_angle.tan(), height];
        let w = [len2([p[0], p[2]]), p[1]];
        let rim_len_sq = dot2(rim, rim);
        let along_slant = if rim_len_sq > 0.0 {
            (dot2(w, rim) / rim_len_sq).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let to_slant = [w[0] - along_slant * rim[0], w[1] - along_slant * rim[1]];
        let to_base = [(w[0] - rim[0]).max(0.0), w[1] - height];
        let magnitude = dot2(to_slant, to_slant).min(dot2(to_base, to_base)).sqrt();
        let beyond = (height * w[0] - rim[0] * w[1]).max(w[1] - height);
        // On the boundary `beyond` is zero and so is the magnitude.  A cone without
        // interior also has `beyond == 0` along its axis, where the magnitude is
        // the positive distance to the solid, so only a strictly negative `beyond`
        // flips the sign.
        if beyond < 0.0 { -magnitude } else { magnitude }
    }
}
