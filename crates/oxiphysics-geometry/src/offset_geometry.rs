#![allow(clippy::if_same_then_else, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Offset / buffer geometry operations.
//!
//! Provides algorithms for:
//!
//! - **Polygon offset** via Minkowski sum approach (inward/outward)
//! - **Curve offset** (equidistant curves for polylines)
//! - **Fillet** (rounding interior corners)
//! - **Chamfer** (beveling corners)
//! - **Rounded rectangles** generation
//! - **Minkowski sum / difference** for convex 2-D shapes
//! - **Morphological operations** (erosion, dilation, opening, closing)
//! - **Tubular neighborhoods** (thickening of curves/surfaces)

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 2-D point type
// ---------------------------------------------------------------------------

/// A 2-D point / vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point2 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}

impl Point2 {
    /// Create a new point.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    /// Euclidean distance to another point.
    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx * dx + dy * dy).sqrt()
    }

    /// Length (norm) of the vector.
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Normalise to unit length. Returns zero if degenerate.
    pub fn normalized(&self) -> Self {
        let l = self.length();
        if l < 1.0e-14 {
            return Self::zero();
        }
        Self {
            x: self.x / l,
            y: self.y / l,
        }
    }

    /// Dot product.
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// 2-D cross product (scalar).
    pub fn cross(&self, other: &Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Perpendicular (rotate 90 degrees CCW).
    pub fn perp(&self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
        }
    }

    /// Add two points.
    pub fn add(&self, other: &Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    /// Subtract.
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    /// Scale.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
        }
    }

    /// Lerp between self and other.
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        Self {
            x: self.x * (1.0 - t) + other.x * t,
            y: self.y * (1.0 - t) + other.y * t,
        }
    }

    /// Angle in radians (atan2).
    pub fn angle(&self) -> f64 {
        self.y.atan2(self.x)
    }

    /// Create from angle and length.
    pub fn from_polar(angle: f64, length: f64) -> Self {
        Self {
            x: length * angle.cos(),
            y: length * angle.sin(),
        }
    }
}

// ---------------------------------------------------------------------------
// Polygon utilities
// ---------------------------------------------------------------------------

/// Signed area of a simple polygon (positive = CCW).
pub fn signed_area(poly: &[Point2]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    area * 0.5
}

/// Absolute area.
pub fn polygon_area(poly: &[Point2]) -> f64 {
    signed_area(poly).abs()
}

/// Perimeter of a polygon.
pub fn polygon_perimeter(poly: &[Point2]) -> f64 {
    let n = poly.len();
    if n < 2 {
        return 0.0;
    }
    let mut perim = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        perim += poly[i].distance(&poly[j]);
    }
    perim
}

/// Centroid of a simple polygon.
pub fn polygon_centroid(poly: &[Point2]) -> Point2 {
    let n = poly.len();
    if n == 0 {
        return Point2::zero();
    }
    let a = signed_area(poly);
    if a.abs() < 1.0e-14 {
        // Degenerate: return average
        let mut cx = 0.0;
        let mut cy = 0.0;
        for p in poly {
            cx += p.x;
            cy += p.y;
        }
        return Point2::new(cx / n as f64, cy / n as f64);
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = poly[i].x * poly[j].y - poly[j].x * poly[i].y;
        cx += (poly[i].x + poly[j].x) * cross;
        cy += (poly[i].y + poly[j].y) * cross;
    }
    let factor = 1.0 / (6.0 * a);
    Point2::new(cx * factor, cy * factor)
}

/// Ensure polygon is in CCW order.
pub fn ensure_ccw(poly: &mut Vec<Point2>) {
    if signed_area(poly) < 0.0 {
        poly.reverse();
    }
}

/// Check if a polygon is convex.
pub fn is_convex(poly: &[Point2]) -> bool {
    let n = poly.len();
    if n < 3 {
        return true;
    }
    let mut sign = 0i32;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let c = poly[(i + 2) % n];
        let cross = b.sub(&a).cross(&c.sub(&b));
        if cross.abs() < 1.0e-12 {
            continue;
        }
        let s = if cross > 0.0 { 1 } else { -1 };
        if sign == 0 {
            sign = s;
        } else if sign != s {
            return false;
        }
    }
    true
}

/// Point-in-polygon test (ray casting).
pub fn point_in_polygon(pt: &Point2, poly: &[Point2]) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let pi = &poly[i];
        let pj = &poly[j];
        if ((pi.y > pt.y) != (pj.y > pt.y))
            && (pt.x < (pj.x - pi.x) * (pt.y - pi.y) / (pj.y - pi.y) + pi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

// ---------------------------------------------------------------------------
// Polygon offset (Minkowski sum with a disk approach)
// ---------------------------------------------------------------------------

/// Offset a simple polygon by `distance`.
///
/// Positive distance = outward (dilation), negative = inward (erosion).
/// Uses edge normal offsetting with miter joins.
pub fn polygon_offset(poly: &[Point2], distance: f64) -> Vec<Point2> {
    let n = poly.len();
    if n < 3 {
        return poly.to_vec();
    }

    // Get edge normals (outward-pointing for CCW polygon)
    let area = signed_area(poly);
    let flip = if area < 0.0 { -1.0 } else { 1.0 };

    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;

        // Edge vectors
        let e_prev = poly[i].sub(&poly[prev]).normalized();
        let e_next = poly[next].sub(&poly[i]).normalized();

        // Outward normals (for CCW: right-hand normal = (dy, -dx))
        let n_prev = Point2::new(e_prev.y * flip, -e_prev.x * flip);
        let n_next = Point2::new(e_next.y * flip, -e_next.x * flip);

        // Bisector
        let bisector = n_prev.add(&n_next);
        let bl = bisector.length();
        if bl < 1.0e-14 {
            // Degenerate: use one normal
            result.push(poly[i].add(&n_prev.scale(distance)));
            continue;
        }
        let bisector = bisector.scale(1.0 / bl);

        // Miter factor: 1 / cos(half-angle)
        let cos_half = bisector.dot(&n_prev).max(0.1); // clamp to avoid extreme miter
        let offset_len = distance / cos_half;

        result.push(poly[i].add(&bisector.scale(offset_len)));
    }

    result
}

/// Polygon offset with arc joins at convex corners.
///
/// `segments_per_corner` controls arc resolution.
pub fn polygon_offset_round(
    poly: &[Point2],
    distance: f64,
    segments_per_corner: usize,
) -> Vec<Point2> {
    let n = poly.len();
    if n < 3 {
        return poly.to_vec();
    }

    let area = signed_area(poly);
    let flip = if area < 0.0 { -1.0 } else { 1.0 };
    let segs = segments_per_corner.max(1);

    let mut result = Vec::new();

    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;

        let e_prev = poly[i].sub(&poly[prev]).normalized();
        let e_next = poly[next].sub(&poly[i]).normalized();

        let n_prev = Point2::new(e_prev.y * flip, -e_prev.x * flip);
        let n_next = Point2::new(e_next.y * flip, -e_next.x * flip);

        let cross = e_prev.cross(&e_next) * flip;

        if (distance > 0.0 && cross < -1.0e-10) || (distance < 0.0 && cross > 1.0e-10) {
            // Convex corner: insert arc
            let angle_start = n_prev.angle();
            let mut angle_end = n_next.angle();
            if distance > 0.0 {
                while angle_end < angle_start {
                    angle_end += 2.0 * PI;
                }
            } else {
                while angle_end > angle_start {
                    angle_end -= 2.0 * PI;
                }
            }
            for k in 0..=segs {
                let t = k as f64 / segs as f64;
                let angle = angle_start + t * (angle_end - angle_start);
                let offset = Point2::from_polar(angle, distance.abs());
                let sign = if distance > 0.0 { 1.0 } else { -1.0 };
                result.push(poly[i].add(&offset.scale(sign)));
            }
        } else {
            // Concave or near-straight: miter
            let bisector = n_prev.add(&n_next);
            let bl = bisector.length();
            if bl < 1.0e-14 {
                result.push(poly[i].add(&n_prev.scale(distance)));
            } else {
                let bisector = bisector.scale(1.0 / bl);
                let cos_half = bisector.dot(&n_prev).max(0.1);
                let offset_len = distance / cos_half;
                result.push(poly[i].add(&bisector.scale(offset_len)));
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Curve offset (polyline offset)
// ---------------------------------------------------------------------------

/// Offset an open polyline by `distance`.
///
/// Returns the offset polyline. Positive distance = left side (for CCW).
pub fn polyline_offset(polyline: &[Point2], distance: f64) -> Vec<Point2> {
    let n = polyline.len();
    if n < 2 {
        return polyline.to_vec();
    }
    if n == 2 {
        let e = polyline[1].sub(&polyline[0]).normalized();
        let normal = Point2::new(-e.y, e.x);
        let offset = normal.scale(distance);
        return vec![polyline[0].add(&offset), polyline[1].add(&offset)];
    }

    let mut result = Vec::with_capacity(n);

    // First point: use first edge normal
    let e0 = polyline[1].sub(&polyline[0]).normalized();
    let n0 = Point2::new(-e0.y, e0.x);
    result.push(polyline[0].add(&n0.scale(distance)));

    // Interior points: bisector
    for i in 1..n - 1 {
        let e_prev = polyline[i].sub(&polyline[i - 1]).normalized();
        let e_next = polyline[i + 1].sub(&polyline[i]).normalized();
        let n_prev = Point2::new(-e_prev.y, e_prev.x);
        let n_next = Point2::new(-e_next.y, e_next.x);
        let bisector = n_prev.add(&n_next);
        let bl = bisector.length();
        if bl < 1.0e-14 {
            result.push(polyline[i].add(&n_prev.scale(distance)));
        } else {
            let bisector = bisector.scale(1.0 / bl);
            let cos_half = bisector.dot(&n_prev).max(0.1);
            let offset_len = distance / cos_half;
            result.push(polyline[i].add(&bisector.scale(offset_len)));
        }
    }

    // Last point: use last edge normal
    let e_last = polyline[n - 1].sub(&polyline[n - 2]).normalized();
    let n_last = Point2::new(-e_last.y, e_last.x);
    result.push(polyline[n - 1].add(&n_last.scale(distance)));

    result
}

// ---------------------------------------------------------------------------
// Fillet and chamfer
// ---------------------------------------------------------------------------

/// Apply a fillet (arc rounding) at each vertex of a polygon.
///
/// `radius` is the fillet radius; `segments` controls arc resolution.
pub fn fillet_polygon(poly: &[Point2], radius: f64, segments: usize) -> Vec<Point2> {
    let n = poly.len();
    if n < 3 || radius <= 0.0 {
        return poly.to_vec();
    }
    let segs = segments.max(2);
    let mut result = Vec::new();

    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;

        let to_prev = poly[prev].sub(&poly[i]).normalized();
        let to_next = poly[next].sub(&poly[i]).normalized();

        let dot = to_prev.dot(&to_next).clamp(-1.0, 1.0);
        let half_angle = ((1.0 + dot) / 2.0).sqrt().acos().max(1.0e-6);
        let tan_half = half_angle.tan().max(1.0e-6);
        let trim = (radius / tan_half).min(
            poly[i]
                .distance(&poly[prev])
                .min(poly[i].distance(&poly[next]))
                * 0.49,
        );
        let actual_radius = trim * tan_half;

        // Trim points
        let p_start = poly[i].add(&to_prev.scale(trim));
        let p_end = poly[i].add(&to_next.scale(trim));

        // Arc center: offset from vertex along bisector
        let bisector = to_prev.add(&to_next).normalized();
        let center_dist = actual_radius / half_angle.sin().max(1.0e-6);
        let center = poly[i].add(&bisector.scale(center_dist));

        // Arc from p_start to p_end
        let angle_start = p_start.sub(&center).angle();
        let angle_end = p_end.sub(&center).angle();
        let mut sweep = angle_end - angle_start;
        // Choose the short arc
        while sweep > PI {
            sweep -= 2.0 * PI;
        }
        while sweep < -PI {
            sweep += 2.0 * PI;
        }

        for k in 0..=segs {
            let t = k as f64 / segs as f64;
            let a = angle_start + t * sweep;
            result.push(Point2::new(
                center.x + actual_radius * a.cos(),
                center.y + actual_radius * a.sin(),
            ));
        }
    }

    result
}

/// Apply a chamfer (bevel) at each vertex of a polygon.
///
/// `distance` is how far from the vertex to cut on each adjacent edge.
pub fn chamfer_polygon(poly: &[Point2], distance: f64) -> Vec<Point2> {
    let n = poly.len();
    if n < 3 || distance <= 0.0 {
        return poly.to_vec();
    }
    let mut result = Vec::new();

    for i in 0..n {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;

        let to_prev = poly[prev].sub(&poly[i]).normalized();
        let to_next = poly[next].sub(&poly[i]).normalized();

        let d_prev = poly[i].distance(&poly[prev]);
        let d_next = poly[i].distance(&poly[next]);
        let trim = distance.min(d_prev * 0.49).min(d_next * 0.49);

        let p1 = poly[i].add(&to_prev.scale(trim));
        let p2 = poly[i].add(&to_next.scale(trim));

        result.push(p1);
        result.push(p2);
    }

    result
}

// ---------------------------------------------------------------------------
// Rounded rectangle
// ---------------------------------------------------------------------------

/// Generate a rounded rectangle polygon.
///
/// `width`, `height` are the full dimensions, `radius` is corner radius.
/// The rectangle is centred at the origin.
pub fn rounded_rectangle(
    width: f64,
    height: f64,
    radius: f64,
    segments_per_corner: usize,
) -> Vec<Point2> {
    let r = radius.min(width / 2.0).min(height / 2.0).max(0.0);
    let hw = width / 2.0;
    let hh = height / 2.0;
    let segs = segments_per_corner.max(1);

    let mut result = Vec::new();

    // Four corners: top-right, top-left, bottom-left, bottom-right
    let corners = [
        (hw - r, hh - r, 0.0),
        (-(hw - r), hh - r, PI / 2.0),
        (-(hw - r), -(hh - r), PI),
        (hw - r, -(hh - r), 3.0 * PI / 2.0),
    ];

    for &(cx, cy, start_angle) in &corners {
        for k in 0..=segs {
            let t = k as f64 / segs as f64;
            let a = start_angle + t * (PI / 2.0);
            result.push(Point2::new(cx + r * a.cos(), cy + r * a.sin()));
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Minkowski sum / difference (convex shapes)
// ---------------------------------------------------------------------------

/// Minkowski sum of two convex polygons (both assumed CCW).
///
/// Result is a convex polygon.
pub fn minkowski_sum_convex(a: &[Point2], b: &[Point2]) -> Vec<Point2> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let na = a.len();
    let nb = b.len();

    // Find bottom-most points
    let mut ia = 0;
    for i in 1..na {
        if a[i].y < a[ia].y || (a[i].y == a[ia].y && a[i].x < a[ia].x) {
            ia = i;
        }
    }
    let mut ib = 0;
    for i in 1..nb {
        if b[i].y < b[ib].y || (b[i].y == b[ib].y && b[i].x < b[ib].x) {
            ib = i;
        }
    }

    let mut result = Vec::with_capacity(na + nb);
    let mut i = 0;
    let mut j = 0;

    while i < na || j < nb {
        let ai = (ia + i) % na;
        let bi = (ib + j) % nb;
        result.push(a[ai].add(&b[bi]));

        let a_next = (ia + i + 1) % na;
        let b_next = (ib + j + 1) % nb;

        let ea = a[a_next].sub(&a[ai]);
        let eb = b[b_next].sub(&b[bi]);

        let cross = ea.cross(&eb);

        if i >= na {
            j += 1;
        } else if j >= nb {
            i += 1;
        } else if cross > 0.0 {
            i += 1;
        } else if cross < 0.0 {
            j += 1;
        } else {
            // Parallel edges: advance both
            i += 1;
            j += 1;
        }
    }

    result
}

/// Minkowski difference of two convex polygons.
///
/// This is the Minkowski sum of A and the reflection of B about the origin.
pub fn minkowski_difference_convex(a: &[Point2], b: &[Point2]) -> Vec<Point2> {
    let b_reflected: Vec<Point2> = b.iter().map(|p| Point2::new(-p.x, -p.y)).collect();
    // Re-order reflected polygon to CCW
    let mut b_ref = b_reflected;
    ensure_ccw(&mut b_ref);
    minkowski_sum_convex(a, &b_ref)
}

// ---------------------------------------------------------------------------
// Morphological operations
// ---------------------------------------------------------------------------

/// Dilate a polygon by `radius` (outward offset).
pub fn morphological_dilation(poly: &[Point2], radius: f64) -> Vec<Point2> {
    polygon_offset(poly, radius.abs())
}

/// Erode a polygon by `radius` (inward offset).
pub fn morphological_erosion(poly: &[Point2], radius: f64) -> Vec<Point2> {
    polygon_offset(poly, -radius.abs())
}

/// Morphological opening: erosion followed by dilation.
///
/// Removes small protrusions.
pub fn morphological_opening(poly: &[Point2], radius: f64) -> Vec<Point2> {
    let eroded = morphological_erosion(poly, radius);
    morphological_dilation(&eroded, radius)
}

/// Morphological closing: dilation followed by erosion.
///
/// Fills small concavities.
pub fn morphological_closing(poly: &[Point2], radius: f64) -> Vec<Point2> {
    let dilated = morphological_dilation(poly, radius);
    morphological_erosion(&dilated, radius)
}

// ---------------------------------------------------------------------------
// Tubular neighborhoods
// ---------------------------------------------------------------------------

/// Generate a tubular neighborhood (thickened polyline).
///
/// Returns a closed polygon that is the union of disks of given `radius`
/// centred along the polyline (approximation via offset + closure).
pub fn tubular_neighborhood(polyline: &[Point2], radius: f64) -> Vec<Point2> {
    if polyline.len() < 2 {
        if polyline.len() == 1 {
            // Single point: return a circle
            let mut circle = Vec::new();
            let segs = 16;
            for k in 0..segs {
                let a = 2.0 * PI * k as f64 / segs as f64;
                circle.push(Point2::new(
                    polyline[0].x + radius * a.cos(),
                    polyline[0].y + radius * a.sin(),
                ));
            }
            return circle;
        }
        return Vec::new();
    }

    let left = polyline_offset(polyline, radius);
    let right = polyline_offset(polyline, -radius);

    // Close into a loop: left forward, right reversed, semicircular caps
    let mut result = Vec::new();

    // Start cap (semicircle around first point)
    let e0 = polyline[1].sub(&polyline[0]).normalized();
    let n0 = Point2::new(-e0.y, e0.x);
    let angle_start = n0.scale(radius).add(&polyline[0]).sub(&polyline[0]).angle();
    let cap_segs = 8;
    // Left side forward
    result.extend_from_slice(&left);

    // End cap
    let n = polyline.len();
    let e_end = polyline[n - 1].sub(&polyline[n - 2]).normalized();
    let n_end = Point2::new(-e_end.y, e_end.x);
    let cap_start = n_end.angle();
    for k in 0..=cap_segs {
        let t = k as f64 / cap_segs as f64;
        let a = cap_start - t * PI;
        result.push(Point2::new(
            polyline[n - 1].x + radius * a.cos(),
            polyline[n - 1].y + radius * a.sin(),
        ));
    }

    // Right side reversed
    for p in right.iter().rev() {
        result.push(*p);
    }

    // Start cap
    let _n0_angle = Point2::new(e0.y, -e0.x).angle();
    let cap_start2 = Point2::new(e0.y, -e0.x).angle();
    for k in 0..=cap_segs {
        let t = k as f64 / cap_segs as f64;
        let a = cap_start2 - t * PI;
        result.push(Point2::new(
            polyline[0].x + radius * a.cos(),
            polyline[0].y + radius * a.sin(),
        ));
    }

    let _ = angle_start; // suppress unused
    result
}

/// Compute the distance from a point to a line segment.
pub fn point_to_segment_distance(pt: &Point2, a: &Point2, b: &Point2) -> f64 {
    let ab = b.sub(a);
    let ap = pt.sub(a);
    let t = (ap.dot(&ab) / ab.dot(&ab)).clamp(0.0, 1.0);
    let proj = a.add(&ab.scale(t));
    pt.distance(&proj)
}

/// Check if a point is within distance `d` of a polyline.
pub fn point_in_tubular(pt: &Point2, polyline: &[Point2], radius: f64) -> bool {
    for i in 0..polyline.len().saturating_sub(1) {
        if point_to_segment_distance(pt, &polyline[i], &polyline[i + 1]) <= radius {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// 3-D offset helpers
// ---------------------------------------------------------------------------

/// Offset a 3-D polygon (planar) by `distance` along edge normals.
///
/// The polygon lives in `[f64; 3]` space and is assumed planar.
/// This projects to 2-D, offsets, and lifts back.
pub fn polygon_offset_3d(
    vertices: &[[f64; 3]],
    distance: f64,
    plane_normal: [f64; 3],
) -> Vec<[f64; 3]> {
    if vertices.len() < 3 {
        return vertices.to_vec();
    }

    // Build a local 2-D frame
    let n = normalize_3d(plane_normal);
    let u = find_perpendicular(n);
    let v = cross_3d(n, u);

    // Project to 2-D
    let origin = vertices[0];
    let poly2d: Vec<Point2> = vertices
        .iter()
        .map(|p| {
            let rel = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
            Point2::new(dot_3d(rel, u), dot_3d(rel, v))
        })
        .collect();

    // Offset in 2-D
    let offset2d = polygon_offset(&poly2d, distance);

    // Lift back to 3-D
    offset2d
        .iter()
        .map(|p| {
            [
                origin[0] + p.x * u[0] + p.y * v[0],
                origin[1] + p.x * u[1] + p.y * v[1],
                origin[2] + p.x * u[2] + p.y * v[2],
            ]
        })
        .collect()
}

/// Normalise a 3-D vector.
fn normalize_3d(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if l < 1.0e-14 {
        return [0.0, 0.0, 1.0];
    }
    [v[0] / l, v[1] / l, v[2] / l]
}

/// Dot product of 3-D vectors.
fn dot_3d(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of 3-D vectors.
fn cross_3d(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Find a perpendicular vector to `n`.
fn find_perpendicular(n: [f64; 3]) -> [f64; 3] {
    let candidate = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let c = cross_3d(n, candidate);
    normalize_3d(c)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    fn square() -> Vec<Point2> {
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]
    }

    fn triangle() -> Vec<Point2> {
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, 1.732),
        ]
    }

    // 1. Signed area of CCW square is positive
    #[test]
    fn test_signed_area_ccw() {
        let sq = square();
        let a = signed_area(&sq);
        assert!((a - 1.0).abs() < 1.0e-6, "Area should be 1.0, got {:.6}", a);
    }

    // 2. Signed area of CW square is negative
    #[test]
    fn test_signed_area_cw() {
        let mut sq = square();
        sq.reverse();
        let a = signed_area(&sq);
        assert!(
            (a - (-1.0)).abs() < 1.0e-6,
            "CW area should be -1.0, got {:.6}",
            a
        );
    }

    // 3. Polygon area is absolute
    #[test]
    fn test_polygon_area() {
        let mut sq = square();
        sq.reverse();
        assert!((polygon_area(&sq) - 1.0).abs() < 1.0e-6);
    }

    // 4. Perimeter of unit square
    #[test]
    fn test_perimeter() {
        let sq = square();
        let p = polygon_perimeter(&sq);
        assert!(
            (p - 4.0).abs() < 1.0e-6,
            "Perimeter should be 4.0, got {:.6}",
            p
        );
    }

    // 5. Centroid of unit square at (0.5, 0.5)
    #[test]
    fn test_centroid() {
        let sq = square();
        let c = polygon_centroid(&sq);
        assert!(
            (c.x - 0.5).abs() < 1.0e-6,
            "cx should be 0.5, got {:.6}",
            c.x
        );
        assert!(
            (c.y - 0.5).abs() < 1.0e-6,
            "cy should be 0.5, got {:.6}",
            c.y
        );
    }

    // 6. Square is convex
    #[test]
    fn test_square_convex() {
        assert!(is_convex(&square()));
    }

    // 7. L-shape is not convex
    #[test]
    fn test_l_shape_not_convex() {
        let l = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 2.0),
            Point2::new(0.0, 2.0),
        ];
        assert!(!is_convex(&l));
    }

    // 8. Point inside polygon
    #[test]
    fn test_point_inside() {
        let sq = square();
        assert!(point_in_polygon(&Point2::new(0.5, 0.5), &sq));
    }

    // 9. Point outside polygon
    #[test]
    fn test_point_outside() {
        let sq = square();
        assert!(!point_in_polygon(&Point2::new(2.0, 2.0), &sq));
    }

    // 10. Outward offset increases area
    #[test]
    fn test_offset_increases_area() {
        let sq = square();
        let offset = polygon_offset(&sq, 0.1);
        let area_orig = polygon_area(&sq);
        let area_offset = polygon_area(&offset);
        assert!(
            area_offset > area_orig,
            "Offset area {:.6} should be > original {:.6}",
            area_offset,
            area_orig
        );
    }

    // 11. Inward offset decreases area
    #[test]
    fn test_offset_decreases_area() {
        let sq = square();
        let offset = polygon_offset(&sq, -0.1);
        let area_orig = polygon_area(&sq);
        let area_offset = polygon_area(&offset);
        assert!(
            area_offset < area_orig,
            "Inward offset area {:.6} should be < original {:.6}",
            area_offset,
            area_orig
        );
    }

    // 12. Zero offset preserves polygon
    #[test]
    fn test_zero_offset() {
        let sq = square();
        let offset = polygon_offset(&sq, 0.0);
        assert_eq!(offset.len(), sq.len());
        for (a, b) in offset.iter().zip(sq.iter()) {
            assert!(a.distance(b) < 1.0e-6);
        }
    }

    // 13. Polyline offset preserves point count
    #[test]
    fn test_polyline_offset_count() {
        let line = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 1.0),
        ];
        let offset = polyline_offset(&line, 0.1);
        assert_eq!(offset.len(), line.len());
    }

    // 14. Fillet produces more points
    #[test]
    fn test_fillet_more_points() {
        let sq = square();
        let filleted = fillet_polygon(&sq, 0.1, 4);
        assert!(
            filleted.len() > sq.len(),
            "Fillet should add points: {} vs {}",
            filleted.len(),
            sq.len()
        );
    }

    // 15. Chamfer doubles vertex count
    #[test]
    fn test_chamfer_doubles() {
        let sq = square();
        let chamfered = chamfer_polygon(&sq, 0.1);
        assert_eq!(
            chamfered.len(),
            2 * sq.len(),
            "Chamfer should produce 2 points per vertex"
        );
    }

    // 16. Rounded rectangle has correct number of points
    #[test]
    fn test_rounded_rect_points() {
        let rect = rounded_rectangle(4.0, 2.0, 0.5, 4);
        // 4 corners * (4+1) = 20 points
        assert_eq!(rect.len(), 20, "Should have 20 points, got {}", rect.len());
    }

    // 17. Rounded rectangle area is between rectangle and full roundedness
    #[test]
    fn test_rounded_rect_area() {
        let rect = rounded_rectangle(4.0, 2.0, 0.5, 16);
        let area = polygon_area(&rect);
        let rect_area = 4.0 * 2.0;
        // Area should be slightly less than full rectangle (corners are rounded)
        assert!(
            area < rect_area + 0.1 && area > rect_area - 2.0,
            "Rounded rect area {:.6} should be near {:.6}",
            area,
            rect_area
        );
    }

    // 18. Minkowski sum of two squares
    #[test]
    fn test_minkowski_sum_squares() {
        let a = square();
        let b = vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 0.0),
            Point2::new(0.5, 0.5),
            Point2::new(0.0, 0.5),
        ];
        let sum = minkowski_sum_convex(&a, &b);
        assert!(!sum.is_empty(), "Minkowski sum should not be empty");
        let area = polygon_area(&sum);
        // Should be ~(1+0.5)*(1+0.5) = 2.25
        assert!(
            (area - 2.25).abs() < 0.5,
            "Minkowski sum area should be ~2.25, got {:.6}",
            area
        );
    }

    // 19. Minkowski difference
    #[test]
    fn test_minkowski_difference() {
        let a = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.0, 2.0),
        ];
        let b = vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 0.0),
            Point2::new(0.5, 0.5),
            Point2::new(0.0, 0.5),
        ];
        let diff = minkowski_difference_convex(&a, &b);
        assert!(!diff.is_empty(), "Minkowski difference should not be empty");
    }

    // 20. Dilation increases area
    #[test]
    fn test_dilation_area() {
        let sq = square();
        let dilated = morphological_dilation(&sq, 0.2);
        assert!(polygon_area(&dilated) > polygon_area(&sq));
    }

    // 21. Erosion decreases area
    #[test]
    fn test_erosion_area() {
        let sq = square();
        let eroded = morphological_erosion(&sq, 0.1);
        assert!(
            polygon_area(&eroded) < polygon_area(&sq),
            "Erosion should reduce area"
        );
    }

    // 22. Opening area <= original area
    #[test]
    fn test_opening_area() {
        let sq = square();
        let opened = morphological_opening(&sq, 0.05);
        assert!(polygon_area(&opened) <= polygon_area(&sq) + 1.0e-3);
    }

    // 23. Closing area >= original area
    #[test]
    fn test_closing_area() {
        let sq = square();
        let closed = morphological_closing(&sq, 0.05);
        // For a convex shape, closing ~ original
        assert!(polygon_area(&closed) >= polygon_area(&sq) - 1.0e-3);
    }

    // 24. Tubular neighborhood produces closed polygon
    #[test]
    fn test_tubular_closed() {
        let line = vec![Point2::new(0.0, 0.0), Point2::new(3.0, 0.0)];
        let tube = tubular_neighborhood(&line, 0.5);
        assert!(tube.len() >= 4, "Tube should have at least 4 points");
    }

    // 25. Point-to-segment distance
    #[test]
    fn test_point_segment_dist() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(2.0, 0.0);
        let pt = Point2::new(1.0, 1.0);
        let d = point_to_segment_distance(&pt, &a, &b);
        assert!(
            (d - 1.0).abs() < 1.0e-6,
            "Distance should be 1.0, got {:.6}",
            d
        );
    }

    // 26. Point in tubular
    #[test]
    fn test_point_in_tubular() {
        let line = vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)];
        assert!(point_in_tubular(&Point2::new(1.0, 0.3), &line, 0.5));
        assert!(!point_in_tubular(&Point2::new(1.0, 1.0), &line, 0.5));
    }

    // 27. Round offset produces more points than miter
    #[test]
    fn test_round_offset_more_points() {
        let sq = square();
        let miter = polygon_offset(&sq, 0.2);
        let round = polygon_offset_round(&sq, 0.2, 4);
        assert!(
            round.len() >= miter.len(),
            "Round offset should have >= miter points: {} vs {}",
            round.len(),
            miter.len()
        );
    }

    // 28. ensure_ccw makes CW polygon CCW
    #[test]
    fn test_ensure_ccw() {
        let mut poly = square();
        poly.reverse();
        assert!(signed_area(&poly) < 0.0);
        ensure_ccw(&mut poly);
        assert!(signed_area(&poly) > 0.0);
    }

    // 29. 3-D polygon offset preserves point count
    #[test]
    fn test_3d_offset_count() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let offset = polygon_offset_3d(&verts, 0.1, [0.0, 0.0, 1.0]);
        assert_eq!(offset.len(), verts.len());
    }

    // 30. Point2 distance
    #[test]
    fn test_point2_distance() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(3.0, 4.0);
        assert!((a.distance(&b) - 5.0).abs() < EPS);
    }

    // 31. Point2 normalized
    #[test]
    fn test_point2_normalized() {
        let p = Point2::new(3.0, 4.0);
        let n = p.normalized();
        assert!((n.length() - 1.0).abs() < EPS);
    }

    // 32. Triangle area
    #[test]
    fn test_triangle_area() {
        let tri = triangle();
        let a = polygon_area(&tri);
        // equilateral-ish triangle with base 2, height ~1.732 => area ~ 1.732
        assert!(
            a > 1.5 && a < 2.0,
            "Triangle area should be ~1.732, got {:.6}",
            a
        );
    }
}
