//! # SdfTriangularPrism - Trait Implementations
//!
//! This module contains trait implementations for `SdfTriangularPrism`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfTriangularPrism;

/// Exact signed distance to an equilateral triangle in the plane `(x, y)`.
///
/// The triangle has its centroid at the origin, its base edge on the line
/// `y = -r / sqrt(3)` between `x = -r` and `x = r`, and its apex at
/// `(0, 2 r / sqrt(3))`, so `r` is half the side length (Quilez'
/// `sdEquilateralTriangle`).  The distance is negative inside the triangle.
fn equilateral_triangle_distance(x: f64, y: f64, r: f64) -> f64 {
    let k = 3.0_f64.sqrt();
    let mut px = x.abs() - r;
    let mut py = y + r / k;
    if px + k * py > 0.0 {
        let reflected = ((px - k * py) / 2.0, (-k * px - py) / 2.0);
        px = reflected.0;
        py = reflected.1;
    }
    px -= px.clamp(-2.0 * r, 0.0);
    let distance = px.hypot(py);
    // `py > 0` is the inside of the folded triangle.  A point exactly on the
    // line through an edge is either on the edge (`distance == 0`) or outside.
    if py > 0.0 { -distance } else { distance }
}

/// Exact signed distance to the prism whose cross-section is an equilateral
/// triangle of side `side` in the `(x, z)` plane, extruded along `y` between
/// `-half_height` and `half_height`.
///
/// The triangle is centred on the axis: its centroid is `(0, 0)`, its apex points
/// towards `+z` and its base edge lies at `z = -side / (2 sqrt(3))`.  The 2-D
/// distance is combined with the slab `|y| <= half_height` by the standard
/// extrusion rule `min(max(d_xz, d_y), 0) + |max((d_xz, d_y), 0)|`, which is
/// negative everywhere inside and exact outside.
///
/// A non-positive `side` collapses the cross-section to the axis and a
/// non-positive `half_height` to the plane `y = 0`.
impl Sdf for SdfTriangularPrism {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let r = 0.5 * self.side.max(0.0);
        let half_height = self.half_height.max(0.0);
        let d_xz = if r > 0.0 {
            equilateral_triangle_distance(p[0], p[2], r)
        } else {
            len2([p[0], p[2]])
        };
        let d_y = p[1].abs() - half_height;
        d_xz.max(d_y).min(0.0) + len2([d_xz.max(0.0), d_y.max(0.0)])
    }
}
