//! # SdfPyramid - Trait Implementations
//!
//! This module contains trait implementations for `SdfPyramid`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfPyramid;

/// Signed distance to the four lateral faces of the pyramid with a unit-side
/// base on `y = 0` and the apex at `(0, h, 0)` (Quilez' `sdPyramid`).
///
/// The point is folded into the octant `x >= z >= 0`, expressed in the frame of
/// the single lateral face that remains (`q.2 / sqrt(m2)` is the signed distance
/// to that face plane) and the in-plane distance to the face triangle is
/// recovered from its two slanted edges.
///
/// Outside the solid, and for `p.y >= 0`, this is the exact Euclidean distance.
/// Inside the solid it is minus the distance to the nearest *lateral* face plane;
/// the base is not considered, so the caller must still apply `max(., -p.y)`.
fn lateral_distance(p: [f64; 3], h: f64) -> f64 {
    let m2 = h * h + 0.25;
    let (ax, az) = (p[0].abs(), p[2].abs());
    let (x, z) = if az > ax { (az, ax) } else { (ax, az) };
    let (x, z) = (x - 0.5, z - 0.5);
    let q = [z, h * p[1] - 0.5 * x, h * x + 0.5 * p[1]];
    let s = (-q[0]).max(0.0);
    let t = ((q[1] - 0.5 * q[0]) / (m2 + 0.25)).clamp(0.0, 1.0);
    let a = m2 * (q[0] + s) * (q[0] + s) + q[1] * q[1];
    let b = m2 * (q[0] + 0.5 * t) * (q[0] + 0.5 * t) + (q[1] - m2 * t) * (q[1] - m2 * t);
    let in_face_region = q[1].min(-q[0] * m2 - q[1] * 0.5) > 0.0;
    let d2 = if in_face_region { 0.0 } else { a.min(b) };
    let magnitude = ((d2 + q[2] * q[2]) / m2).sqrt();
    // `q.2 > 0`: outside the face plane.  `-p.y > 0`: below the base.
    let side = q[2].max(-p[1]);
    if side > 0.0 {
        magnitude
    } else if side < 0.0 {
        -magnitude
    } else {
        0.0
    }
}

/// Exact signed distance to the square-base pyramid.
///
/// The base is the square `|x|, |z| <= half_base` on the plane `y = 0` and the
/// apex is `(0, height, 0)`.
///
/// * Below the base plane (`y <= 0`) the nearest point of the solid lies on the
///   base, so the distance is that to the base square.
/// * Above it, the lateral-face distance `lateral_distance` (evaluated on
///   the pyramid scaled to a unit base) is combined with the base plane through
///   `max(lateral, -y)`, which gives the exact interior distance
///   `-min(distance to the lateral planes, y)`.
///
/// A non-positive extent describes a degenerate pyramid: a zero `half_base` is the
/// vertical segment from the origin to the apex (a point when `height` is zero as
/// well) and a zero `height` is the flat base square.
impl Sdf for SdfPyramid {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let half_base = self.half_base.max(0.0);
        let height = self.height.max(0.0);
        if half_base == 0.0 {
            let along_axis = p[1] - p[1].clamp(0.0, height);
            return (p[0] * p[0] + along_axis * along_axis + p[2] * p[2]).sqrt();
        }
        if p[1] <= 0.0 || height == 0.0 {
            let dx = (p[0].abs() - half_base).max(0.0);
            let dz = (p[2].abs() - half_base).max(0.0);
            return (dx * dx + p[1] * p[1] + dz * dz).sqrt();
        }
        let size = 2.0 * half_base;
        let lateral =
            size * lateral_distance([p[0] / size, p[1] / size, p[2] / size], height / size);
        lateral.max(-p[1])
    }
}
