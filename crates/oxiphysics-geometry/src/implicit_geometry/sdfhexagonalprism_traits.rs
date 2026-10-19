//! # SdfHexagonalPrism - Trait Implementations
//!
//! This module contains trait implementations for `SdfHexagonalPrism`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfHexagonalPrism;

/// Exact signed distance to a regular hexagonal prism (Quilez' `sdHexPrism`).
///
/// The axis is `Y` and the cross-section lies in the XZ plane with its vertices
/// on the X axis at `(+-radius, 0)`; the flat sides are the lines `|z| = a`
/// (and their rotations by 60 degrees) with the apothem `a = radius * sqrt(3) / 2`.
/// `radius` is the circumscribed radius, as documented on [`SdfHexagonalPrism`].
///
/// The point is folded into the first quadrant and reflected across the line
/// through the origin whose unit normal is `k = (-sqrt(3) / 2, 1 / 2)` when it lies
/// on the far side of it, which maps the three pairs of sides of the hexagon onto
/// the single side `z = a`, `|x| <= a / sqrt(3)`.  The in-plane distance is then
/// the distance to that side (negative below it), and the prism distance combines
/// it with the axial distance `|y| - half_height` as for any extrusion.
///
/// A non-positive `radius` or `half_height` is treated as zero.
impl Sdf for SdfHexagonalPrism {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let sqrt3 = 3.0_f64.sqrt();
        let k0 = -0.5 * sqrt3;
        let k1 = 0.5_f64;
        let k2 = 1.0_f64 / sqrt3;
        let apothem = self.radius.max(0.0) * 0.5 * sqrt3;
        let half_height = self.half_height.max(0.0);
        let qa = [p[0].abs(), p[2].abs()];
        let reflect = (k0 * qa[0] + k1 * qa[1]).min(0.0);
        let qx = qa[0] - 2.0 * k0 * reflect;
        let qz = qa[1] - 2.0 * k1 * reflect;
        let half_side = k2 * apothem;
        let ex = (qx.abs() - half_side).max(0.0);
        let ez = qz - apothem;
        let d_xz = len2([ex, ez]) * if qz > apothem { 1.0 } else { -1.0 };
        let d_y = p[1].abs() - half_height;
        let d_xz_pos = d_xz.max(0.0);
        let d_y_pos = d_y.max(0.0);
        len2([d_xz_pos, d_y_pos]) + d_xz.min(0.0).max(d_y.min(0.0))
    }
}
