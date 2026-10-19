// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Trimmed NURBS surfaces.
//!
//! A trimmed surface is a NURBS surface together with trim curves in its
//! `(u, v)` parameter space. The trim curves bound the region of the surface
//! that is kept:
//!
//! - The curves are flattened into polylines and chained end to start into
//!   closed loops; a loop that does not close by itself is closed by a straight
//!   segment, so the order of the curves does not matter but each curve must run
//!   in the direction of its loop.
//! - A parameter point is kept when the **sum of the winding numbers** of all
//!   loops around it is non-zero (non-zero winding rule), or when it lies on a
//!   loop. With the usual orientation, the outer boundary runs counter-clockwise
//!   and the boundary of every hole runs clockwise, so a point inside a hole has
//!   winding number `+1 - 1 = 0` and is trimmed away.

use super::{NurbsBasis, NurbsControlPoint, NurbsKnotVector, NurbsSurface};

/// Largest flattening error of a trim curve relative to the extent of the
/// control polygons of the trim curves.
const FLATTEN_REL_TOL: f64 = 1.0e-6;

/// Two curve ends closer than this (relative to the same extent) are one point.
const JOIN_REL_TOL: f64 = 1.0e-9;

/// A point closer than this (relative to the same extent) to a loop lies on it.
const BOUNDARY_REL_TOL: f64 = 1.0e-12;

/// Every knot span is split into at least `2^MIN_FLATTEN_DEPTH` pieces, so a
/// curve cannot cheat the flatness test at the midpoint of a span.
const MIN_FLATTEN_DEPTH: u32 = 3;

/// No knot span is split into more than `2^MAX_FLATTEN_DEPTH` pieces.
const MAX_FLATTEN_DEPTH: u32 = 14;

/// A trim curve in the (u, v) parameter space of a NURBS surface.
#[derive(Debug, Clone)]
pub struct NurbsTrimCurve {
    /// Parameter-space control points \[u, v, 0\] with weights.
    pub control_points: Vec<NurbsControlPoint>,
    /// Knot vector for the trim curve.
    pub knots: NurbsKnotVector,
    /// Degree.
    pub degree: usize,
}

impl NurbsTrimCurve {
    /// Create a trim curve in parameter space.
    pub fn new(
        control_points: Vec<NurbsControlPoint>,
        knots: NurbsKnotVector,
        degree: usize,
    ) -> Self {
        Self {
            control_points,
            knots,
            degree,
        }
    }

    /// Evaluate trim curve at parameter `t`, returning (u, v).
    ///
    /// A curve with fewer than `degree + 1` control points or a knot vector of
    /// the wrong length cannot be evaluated and yields `[NaN, NaN]`.
    pub fn evaluate(&self, t: f64) -> [f64; 2] {
        let n = self.control_points.len();
        let p = self.degree;
        if n < p + 1 || self.knots.len() != n + p + 1 {
            return [f64::NAN; 2];
        }
        let basis = NurbsBasis::evaluate(&self.knots, t, p, n);
        let span = self.knots.find_span(t, p, n);
        let mut w_sum = 0.0_f64;
        let mut uv = [0.0_f64; 2];
        for (i, b_i) in basis.iter().enumerate() {
            let idx = span - p + i;
            if idx < n {
                let cp = &self.control_points[idx];
                let wn = b_i * cp.weight;
                w_sum += wn;
                uv[0] += wn * cp.point[0];
                uv[1] += wn * cp.point[1];
            }
        }
        if w_sum.abs() > 1.0e-15 {
            uv[0] /= w_sum;
            uv[1] /= w_sum;
        }
        uv
    }

    /// Whether the curve can be used as (part of) a trim loop.
    ///
    /// A usable trim curve has degree of at least 1, at least `degree + 1`
    /// control points, a finite non-decreasing knot vector of length
    /// `control points + degree + 1` with a non-empty parameter domain, and
    /// finite control points with finite positive weights. Trimming ignores
    /// curves that are not usable.
    pub fn is_valid(&self) -> bool {
        let n = self.control_points.len();
        let p = self.degree;
        if p == 0 || n < p + 1 || self.knots.len() != n + p + 1 {
            return false;
        }
        let knots = &self.knots.knots;
        if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[1] < w[0]) {
            return false;
        }
        if knots[n] <= knots[p] {
            return false;
        }
        self.control_points.iter().all(|cp| {
            cp.weight.is_finite()
                && cp.weight > 0.0
                && cp.point[0].is_finite()
                && cp.point[1].is_finite()
        })
    }

    /// The parameter interval `[knots[degree], knots[n]]` on which the curve is
    /// defined, or `None` if the curve is not [valid](Self::is_valid).
    pub fn domain(&self) -> Option<(f64, f64)> {
        if self.is_valid() {
            Some((
                self.knots.knots[self.degree],
                self.knots.knots[self.control_points.len()],
            ))
        } else {
            None
        }
    }

    /// The same curve traversed in the opposite direction.
    ///
    /// Reverses the control points and mirrors the knot vector within its own
    /// range, so the new curve at `t` is the old curve at `t_first + t_last - t`.
    /// Use it to turn a counter-clockwise loop into the clockwise boundary of a
    /// hole.
    pub fn reversed(&self) -> Self {
        let knots = &self.knots.knots;
        let mirrored = match (knots.first(), knots.last()) {
            (Some(first), Some(last)) => knots.iter().rev().map(|k| first + last - k).collect(),
            _ => Vec::new(),
        };
        Self {
            control_points: self.control_points.iter().rev().copied().collect(),
            knots: NurbsKnotVector { knots: mirrored },
            degree: self.degree,
        }
    }

    /// Flatten the curve into a polyline in the (u, v) plane.
    ///
    /// Every knot span is split recursively until the curve at the middle of a
    /// piece is within `tolerance` of the chord of the piece. The polyline starts
    /// at the start point of the curve and ends at its end point. A curve that is
    /// not [valid](Self::is_valid) gives an empty polyline.
    pub fn flatten(&self, tolerance: f64) -> Vec<[f64; 2]> {
        let Some((t_start, t_end)) = self.domain() else {
            return Vec::new();
        };
        let mut breaks = vec![t_start];
        for &knot in &self.knots.knots {
            if breaks.last().is_some_and(|&last| knot > last) && knot < t_end {
                breaks.push(knot);
            }
        }
        breaks.push(t_end);

        let mut flattener = Flattener {
            curve: self,
            tolerance,
            polyline: vec![self.evaluate(t_start)],
        };
        for span in breaks.windows(2) {
            let start = flattener.polyline.last().copied().unwrap_or_default();
            let end = self.evaluate(span[1]);
            flattener.split(span[0], span[1], start, end, 0);
        }
        flattener.polyline
    }
}

/// Recursive flattening of one trim curve.
struct Flattener<'a> {
    curve: &'a NurbsTrimCurve,
    tolerance: f64,
    polyline: Vec<[f64; 2]>,
}

impl Flattener<'_> {
    /// Append the flattened piece of the curve between `t_a` and `t_b` (without
    /// its start point `p_a`, with its end point `p_b`) to the polyline.
    fn split(&mut self, t_a: f64, t_b: f64, p_a: [f64; 2], p_b: [f64; 2], depth: u32) {
        let t_m = 0.5 * (t_a + t_b);
        let p_m = self.curve.evaluate(t_m);
        let flat = distance_to_segment(p_m, p_a, p_b) <= self.tolerance;
        if depth >= MAX_FLATTEN_DEPTH || (depth >= MIN_FLATTEN_DEPTH && flat) {
            self.polyline.push(p_b);
        } else {
            self.split(t_a, t_m, p_a, p_m, depth + 1);
            self.split(t_m, t_b, p_m, p_b, depth + 1);
        }
    }
}

/// Distance between the point `p` and the segment `a`-`b`.
fn distance_to_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let length_sq = ab[0] * ab[0] + ab[1] * ab[1];
    let s = if length_sq > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / length_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let dx = p[0] - (a[0] + s * ab[0]);
    let dy = p[1] - (a[1] + s * ab[1]);
    (dx * dx + dy * dy).sqrt()
}

/// Distance between two points of the plane.
fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

/// The closed polygons that bound the kept region of a trimmed surface.
///
/// Built once by [`NurbsTrimming::trim_loops`] (or from curves with
/// [`TrimLoops::from_curves`]), so that many parameter points can be classified
/// without flattening the trim curves again. It is a snapshot: later changes of
/// the trim curves are not seen by an existing value.
#[derive(Debug, Clone, Default)]
pub struct TrimLoops {
    /// Each loop is a polygon; the edge from the last to the first vertex closes it.
    loops: Vec<Vec<[f64; 2]>>,
    /// Distance below which a point counts as lying on a loop.
    boundary_tolerance: f64,
}

impl TrimLoops {
    /// Flatten `curves` and chain them into closed loops.
    ///
    /// Curves that are not [valid](NurbsTrimCurve::is_valid) are ignored, and so
    /// are chains that enclose no area.
    pub fn from_curves(curves: &[NurbsTrimCurve]) -> Self {
        let mut low = [f64::INFINITY; 2];
        let mut high = [f64::NEG_INFINITY; 2];
        for curve in curves.iter().filter(|c| c.is_valid()) {
            for cp in &curve.control_points {
                for axis in 0..2 {
                    low[axis] = low[axis].min(cp.point[axis]);
                    high[axis] = high[axis].max(cp.point[axis]);
                }
            }
        }
        let extent = (high[0] - low[0]).max(high[1] - low[1]);
        if !extent.is_finite() || extent <= 0.0 {
            return Self::default();
        }

        let polylines: Vec<Vec<[f64; 2]>> = curves
            .iter()
            .map(|curve| curve.flatten(FLATTEN_REL_TOL * extent))
            .filter(|polyline| polyline.len() >= 2)
            .collect();
        Self {
            loops: chain_loops(&polylines, JOIN_REL_TOL * extent),
            boundary_tolerance: BOUNDARY_REL_TOL * extent,
        }
    }

    /// Number of loops.
    pub fn len(&self) -> usize {
        self.loops.len()
    }

    /// `true` when there is no loop, that is, nothing is trimmed.
    pub fn is_empty(&self) -> bool {
        self.loops.is_empty()
    }

    /// Sum of the winding numbers of all loops around `(u, v)`.
    ///
    /// Counter-clockwise loops contribute `+1` to the points they enclose and
    /// clockwise loops `-1`. The value for a point on a loop is not meaningful;
    /// use [`contains`](Self::contains) to classify points.
    pub fn winding_number(&self, u: f64, v: f64) -> i32 {
        self.loops
            .iter()
            .map(|polygon| polygon_winding_number(polygon, [u, v]))
            .sum()
    }

    /// Whether `(u, v)` belongs to the kept region.
    ///
    /// A point is kept when it lies on a loop or when the sum of the winding
    /// numbers of the loops around it is non-zero. Without loops every point is
    /// kept.
    pub fn contains(&self, u: f64, v: f64) -> bool {
        if self.loops.is_empty() {
            return true;
        }
        let p = [u, v];
        let on_loop = self.loops.iter().any(|polygon| {
            (0..polygon.len()).any(|i| {
                let next = polygon[(i + 1) % polygon.len()];
                distance_to_segment(p, polygon[i], next) <= self.boundary_tolerance
            })
        });
        on_loop || self.winding_number(u, v) != 0
    }
}

/// Chain polylines end to start into closed polygons.
///
/// The polylines are taken in the given order: a chain starts at the first
/// polyline not used yet and takes on the first unused polyline that starts where
/// the chain ends, until the chain is closed or cannot be continued. A chain that
/// is not closed is closed by a straight segment. Chains with fewer than three
/// vertices enclose no area and are dropped.
fn chain_loops(polylines: &[Vec<[f64; 2]>], join_tolerance: f64) -> Vec<Vec<[f64; 2]>> {
    let mut used = vec![false; polylines.len()];
    let mut loops = Vec::new();
    for first in 0..polylines.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut chain = polylines[first].clone();
        while let (Some(&head), Some(&tail)) = (chain.first(), chain.last())
            && distance(head, tail) > join_tolerance
        {
            let continuation = (0..polylines.len()).find(|&j| {
                !used[j]
                    && polylines[j]
                        .first()
                        .is_some_and(|&s| distance(s, tail) <= join_tolerance)
            });
            let Some(next) = continuation else {
                break;
            };
            used[next] = true;
            chain.extend(polylines[next].iter().skip(1).copied());
        }
        if let (Some(&head), Some(&tail)) = (chain.first(), chain.last())
            && distance(head, tail) <= join_tolerance
        {
            chain.pop();
        }
        if chain.len() >= 3 {
            loops.push(chain);
        }
    }
    loops
}

/// Winding number of the closed polygon around the point `p`.
///
/// Counts the signed crossings of the ray from `p` in the +u direction (half-open
/// rule, so a vertex on the ray is counted once). The result is only meaningful
/// for points that do not lie on the polygon.
fn polygon_winding_number(polygon: &[[f64; 2]], p: [f64; 2]) -> i32 {
    let is_left = |a: [f64; 2], b: [f64; 2]| -> f64 {
        (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])
    };
    let mut winding = 0;
    for (i, &a) in polygon.iter().enumerate() {
        let b = polygon[(i + 1) % polygon.len()];
        if a[1] <= p[1] {
            if b[1] > p[1] && is_left(a, b) > 0.0 {
                winding += 1;
            }
        } else if b[1] <= p[1] && is_left(a, b) < 0.0 {
            winding -= 1;
        }
    }
    winding
}

/// Trimmed NURBS surface with one or more trim curves.
#[derive(Debug, Clone)]
pub struct NurbsTrimming {
    /// The underlying NURBS surface.
    pub surface: NurbsSurface,
    /// Trim curves defining the trimmed boundary.
    pub trim_curves: Vec<NurbsTrimCurve>,
}

impl NurbsTrimming {
    /// Create a trimmed NURBS surface.
    pub fn new(surface: NurbsSurface) -> Self {
        Self {
            surface,
            trim_curves: Vec::new(),
        }
    }

    /// Add a trim curve.
    pub fn add_trim_curve(&mut self, trim: NurbsTrimCurve) {
        self.trim_curves.push(trim);
    }

    /// The closed loops formed by the trim curves, for classifying many points.
    ///
    /// See [`TrimLoops`]. [`is_inside`](Self::is_inside) builds the loops anew on
    /// every call, which costs the flattening of all trim curves; build the loops
    /// once with this method when many points are to be classified.
    pub fn trim_loops(&self) -> TrimLoops {
        TrimLoops::from_curves(&self.trim_curves)
    }

    /// Sum of the winding numbers of the trim loops around `(u, v)`.
    ///
    /// See [`TrimLoops::winding_number`].
    pub fn winding_number(&self, u: f64, v: f64) -> i32 {
        self.trim_loops().winding_number(u, v)
    }

    /// Check if a parameter-space point (u, v) lies in the kept region.
    ///
    /// The trim curves are chained into closed loops and the point is kept when
    /// the winding numbers of the loops around it sum to a non-zero value, or when
    /// it lies on a loop (non-zero winding rule, see the [module
    /// documentation](self)). Outer boundaries must run counter-clockwise and the
    /// boundaries of holes clockwise. Without trim curves every point is inside;
    /// trim curves that are not [valid](NurbsTrimCurve::is_valid) are ignored.
    pub fn is_inside(&self, u: f64, v: f64) -> bool {
        if self.trim_curves.is_empty() {
            return true;
        }
        self.trim_loops().contains(u, v)
    }

    /// Evaluate the surface at (u, v), returning None if trimmed out.
    pub fn evaluate(&self, u: f64, v: f64) -> Option<[f64; 3]> {
        if self.is_inside(u, v) {
            Some(self.surface.evaluate(u, v))
        } else {
            None
        }
    }
}
