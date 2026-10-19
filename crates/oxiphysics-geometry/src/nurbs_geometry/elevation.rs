// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bézier decomposition and degree elevation of NURBS curves.
//!
//! Both algorithms work on **homogeneous** control points `(w·x, w·y, w·z, w)`.
//! In that space a rational curve is an ordinary polynomial B-spline, so every
//! blending step below is exact for rational curves (conic arcs, circles) too.
//! The control points of the results are never samples of the input curve.
//!
//! - [`BezierDecomposition`] implements Piegl & Tiller, *The NURBS Book*
//!   (2nd ed.), Algorithm A5.6 (`DecomposeCurve`): every interior knot is
//!   inserted until its multiplicity equals the degree, which splits the curve
//!   into its Bézier segments.
//! - [`DegreeElevation`] raises the degree of a Bézier segment with the closed
//!   form `Q_i = i/(p+1)·P_{i-1} + (1 - i/(p+1))·P_i` and the degree of a whole
//!   curve with Algorithm A5.9 (`DegreeElevateCurve`), which keeps the original
//!   continuity: every distinct knot gains exactly one copy.
//!
//! Curves whose knot vector is not clamped (for example those built with
//! [`NurbsKnotVector::uniform`](super::NurbsKnotVector::uniform)) are first
//! clamped by knot insertion at both ends of the parameter domain, so the
//! algorithms act on the part of the curve that
//! [`NurbsCurve::evaluate`](super::NurbsCurve::evaluate) can reach.

use super::{NurbsControlPoint, NurbsCurve, NurbsKnotVector};
use crate::{Error, Result};
use oxiphysics_fem::isogeometric::{
    DegreeElevation as FemDegreeElevation, NurbsCurve as FemNurbsCurve,
};

/// Homogeneous control point `(w·x, w·y, w·z, w)`.
type Homogeneous = [f64; 4];

/// Relative tolerance below which two knot values are treated as one knot.
const KNOT_TOLERANCE: f64 = 1.0e-14;

/// `true` when two knot values coincide within [`KNOT_TOLERANCE`].
fn knots_equal(a: f64, b: f64) -> bool {
    (a - b).abs() <= KNOT_TOLERANCE * a.abs().max(b.abs()).max(1.0)
}

/// Affine combination `alpha·a + (1 - alpha)·b` of two homogeneous points.
fn blend(alpha: f64, a: Homogeneous, b: Homogeneous) -> Homogeneous {
    let beta = 1.0 - alpha;
    [
        alpha * a[0] + beta * b[0],
        alpha * a[1] + beta * b[1],
        alpha * a[2] + beta * b[2],
        alpha * a[3] + beta * b[3],
    ]
}

/// Back from homogeneous coordinates: divide the position by the weight.
fn from_homogeneous(h: Homogeneous) -> NurbsControlPoint {
    NurbsControlPoint {
        point: [h[0] / h[3], h[1] / h[3], h[2] / h[3]],
        weight: h[3],
    }
}

/// One Bézier segment of a decomposed NURBS curve.
#[derive(Debug, Clone)]
pub struct BezierSegment {
    /// Parameter of the curve at the start of the segment.
    pub t_start: f64,
    /// Parameter of the curve at the end of the segment.
    pub t_end: f64,
    /// The `degree + 1` weighted Bézier control points of the segment. The
    /// segment parameter `s ∈ [0, 1]` corresponds to the curve parameter
    /// `t_start + s·(t_end - t_start)`.
    pub control_points: Vec<NurbsControlPoint>,
}

/// A curve prepared for decomposition and elevation: validated, with snapped
/// knots, homogeneous control points and a clamped knot vector.
struct ClampedCurve {
    degree: usize,
    knots: Vec<f64>,
    points: Vec<Homogeneous>,
}

/// Insert the knot `u` once into a homogeneous B-spline (Boehm).
///
/// `u` must lie in the half-open span `[knots[k], knots[k + 1])` of the
/// parameter domain; a knot that already occurs is simply repeated.
fn insert_knot_once(degree: usize, knots: &mut Vec<f64>, points: &mut Vec<Homogeneous>, u: f64) {
    let span = knots.partition_point(|&k| k <= u).saturating_sub(1);
    let mut refined = Vec::with_capacity(points.len() + 1);
    refined.extend_from_slice(&points[..=span - degree]);
    for i in (span - degree + 1)..=span {
        let alpha = (u - knots[i]) / (knots[i + degree] - knots[i]);
        refined.push(blend(alpha, points[i], points[i - 1]));
    }
    refined.extend_from_slice(&points[span..]);
    *points = refined;
    knots.insert(span + 1, u);
}

/// Make the left end of a homogeneous B-spline clamped: the first knot of the
/// parameter domain `knots[degree]` gets multiplicity `degree + 1` and the
/// control points and knots in front of it are dropped.
fn clamp_left(degree: usize, knots: &mut Vec<f64>, points: &mut Vec<Homogeneous>) {
    let start = knots[degree];
    let copies = |knots: &[f64]| knots.iter().filter(|&&k| k == start).count();
    while copies(knots) < degree + 1 {
        insert_knot_once(degree, knots, points, start);
    }
    let last_copy = knots.partition_point(|&k| k <= start) - 1;
    let cut = last_copy - degree;
    knots.drain(..cut);
    points.drain(..cut);
}

/// Mirror a homogeneous B-spline: reverse the control points and the knots and
/// negate the knots, so that the right end becomes the left end.
fn mirror(knots: &mut [f64], points: &mut [Homogeneous]) {
    knots.reverse();
    for knot in knots.iter_mut() {
        *knot = -*knot;
    }
    points.reverse();
}

impl ClampedCurve {
    /// Validate `curve` and bring it into clamped homogeneous form.
    fn new(curve: &NurbsCurve) -> Result<Self> {
        let degree = curve.degree;
        let n_ctrl = curve.control_points.len();
        if n_ctrl < degree + 1 {
            return Err(Error::too_few_points(degree + 1, n_ctrl));
        }
        let expected = n_ctrl + degree + 1;
        let raw_knots = &curve.knots.knots;
        if raw_knots.len() != expected {
            return Err(Error::length_mismatch(expected, raw_knots.len()));
        }
        if raw_knots.iter().any(|k| !k.is_finite()) || raw_knots.windows(2).any(|w| w[1] < w[0]) {
            return Err(Error::degenerate_geometry(
                "knot vector must be finite and non-decreasing",
            ));
        }
        if curve
            .control_points
            .iter()
            .any(|cp| !cp.weight.is_finite() || cp.weight <= 0.0)
        {
            return Err(Error::degenerate_geometry(
                "NURBS weights must be finite and positive",
            ));
        }
        if curve
            .control_points
            .iter()
            .any(|cp| cp.point.iter().any(|x| !x.is_finite()))
        {
            return Err(Error::degenerate_geometry(
                "NURBS control points must be finite",
            ));
        }

        // Knots that differ by rounding noise only are one knot.
        let mut knots = raw_knots.clone();
        for i in 1..knots.len() {
            if knots_equal(knots[i], knots[i - 1]) {
                knots[i] = knots[i - 1];
            }
        }
        if knots[n_ctrl] <= knots[degree] {
            return Err(Error::degenerate_geometry(
                "the parameter domain of the curve is empty",
            ));
        }

        let mut points: Vec<Homogeneous> = curve
            .control_points
            .iter()
            .map(NurbsControlPoint::homogeneous)
            .collect();
        clamp_left(degree, &mut knots, &mut points);
        mirror(&mut knots, &mut points);
        clamp_left(degree, &mut knots, &mut points);
        mirror(&mut knots, &mut points);

        let clamped = Self {
            degree,
            knots,
            points,
        };
        clamped.check_interior_multiplicities()?;
        Ok(clamped)
    }

    /// Interior knots may repeat at most `degree + 1` times (a jump).
    fn check_interior_multiplicities(&self) -> Result<()> {
        let last = self.knots.len() - 1;
        let mut i = self.degree + 1;
        while i < last - self.degree {
            let mut multiplicity = 1;
            while i + multiplicity < self.knots.len()
                && self.knots[i + multiplicity] == self.knots[i]
            {
                multiplicity += 1;
            }
            if multiplicity > self.degree + 1 {
                return Err(Error::degenerate_geometry(format!(
                    "interior knot {} has multiplicity {multiplicity}, more than degree + 1 = {}",
                    self.knots[i],
                    self.degree + 1
                )));
            }
            i += multiplicity;
        }
        Ok(())
    }

    /// Piegl & Tiller A5.6: split the clamped curve into Bézier segments.
    fn decompose(&self) -> Vec<(f64, f64, Vec<Homogeneous>)> {
        let p = self.degree;
        let u = &self.knots;
        let m = u.len() - 1;
        let mut segments = Vec::new();
        let mut current: Vec<Homogeneous> = self.points[..=p].to_vec();
        let mut alphas = vec![0.0_f64; p.max(1)];
        let mut a = p;
        let mut b = p + 1;
        while b < m {
            let first = b;
            while b < m && u[b + 1] == u[b] {
                b += 1;
            }
            let multiplicity = b - first + 1;
            let mut next = vec![[0.0_f64; 4]; p + 1];
            if multiplicity < p {
                // Insert the knot u[b] until its multiplicity equals the degree.
                let numer = u[b] - u[a];
                for j in ((multiplicity + 1)..=p).rev() {
                    alphas[j - multiplicity - 1] = numer / (u[a + j] - u[a]);
                }
                let inserted = p - multiplicity;
                for j in 1..=inserted {
                    let save = inserted - j;
                    let s = multiplicity + j;
                    for k in (s..=p).rev() {
                        current[k] = blend(alphas[k - s], current[k], current[k - 1]);
                    }
                    next[save] = current[p];
                }
            }
            segments.push((u[a], u[b], current));
            if b < m {
                for (i, slot) in next
                    .iter_mut()
                    .enumerate()
                    .skip(p.saturating_sub(multiplicity))
                {
                    *slot = self.points[b - p + i];
                }
                current = next;
                a = b;
                b += 1;
            } else {
                break;
            }
        }
        segments
    }

    /// The clamped curve as a regular [`FemNurbsCurve`].
    fn to_elevatable(&self) -> FemNurbsCurve {
        let points = self
            .points
            .iter()
            .map(|h| [h[0] / h[3], h[1] / h[3], h[2] / h[3]])
            .collect();
        let weights = self.points.iter().map(|h| h[3]).collect();
        FemNurbsCurve::new(points, weights, self.degree, self.knots.clone())
    }
}

/// Decompose a NURBS curve into Bézier segments.
///
/// Implements Piegl & Tiller Algorithm A5.6 on homogeneous control points, so the
/// control points of the segments carry the true weights of the curve.
pub struct BezierDecomposition;

impl BezierDecomposition {
    /// Extract the Bézier control points of every segment of a NURBS curve.
    ///
    /// Returns one `Vec` of `degree + 1` weighted control points per non-empty
    /// knot span of the parameter domain. A malformed curve (see
    /// [`try_decompose`](Self::try_decompose)) yields an empty `Vec`.
    pub fn decompose(curve: &NurbsCurve) -> Vec<Vec<NurbsControlPoint>> {
        Self::try_decompose(curve)
            .map(|segments| {
                segments
                    .into_iter()
                    .map(|segment| segment.control_points)
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Decompose a NURBS curve into Bézier segments, reporting malformed input.
    ///
    /// A knot of multiplicity `degree + 1` (a jump of the curve) separates two
    /// segments that do not share an end point; otherwise neighbouring segments
    /// share their common end point.
    ///
    /// # Errors
    ///
    /// Fails for a curve with fewer than `degree + 1` control points, a knot
    /// vector of the wrong length or one that is decreasing or not finite,
    /// non-finite or non-positive weights, non-finite control points, an empty
    /// parameter domain or an interior knot of multiplicity above `degree + 1`.
    pub fn try_decompose(curve: &NurbsCurve) -> Result<Vec<BezierSegment>> {
        let clamped = ClampedCurve::new(curve)?;
        Ok(clamped
            .decompose()
            .into_iter()
            .map(|(t_start, t_end, points)| BezierSegment {
                t_start,
                t_end,
                control_points: points.into_iter().map(from_homogeneous).collect(),
            })
            .collect())
    }
}

/// Degree elevation for NURBS curves.
///
/// Raising the degree of a curve leaves its geometry and its parametrisation
/// unchanged; only the representation gets richer.
pub struct DegreeElevation;

impl DegreeElevation {
    /// Elevate the degree of one Bézier segment by one.
    ///
    /// For `p + 1` control points `P_0 … P_p` the result has the `p + 2` points
    /// `Q_0 = P_0`, `Q_i = i/(p+1)·P_{i-1} + (1 - i/(p+1))·P_i` for
    /// `1 ≤ i ≤ p` and `Q_{p+1} = P_p`, blended in homogeneous coordinates so
    /// that rational segments keep their shape. An empty slice gives an empty
    /// `Vec`.
    pub fn elevate_bezier(pts: &[NurbsControlPoint]) -> Vec<NurbsControlPoint> {
        let (Some(first), Some(last)) = (pts.first(), pts.last()) else {
            return Vec::new();
        };
        let degree = pts.len() - 1;
        let homogeneous: Vec<Homogeneous> =
            pts.iter().map(NurbsControlPoint::homogeneous).collect();
        let mut elevated = Vec::with_capacity(pts.len() + 1);
        elevated.push(*first);
        for i in 1..=degree {
            let alpha = i as f64 / (degree + 1) as f64;
            elevated.push(from_homogeneous(blend(
                alpha,
                homogeneous[i - 1],
                homogeneous[i],
            )));
        }
        elevated.push(*last);
        elevated
    }

    /// Elevate a NURBS curve's degree by 1 without changing its geometry.
    ///
    /// Every distinct knot gains one copy, so the continuity of the curve at
    /// each knot is unchanged. This is the infallible form of
    /// [`try_elevate_curve`](Self::try_elevate_curve): for a curve that cannot be
    /// elevated, the returned curve keeps the degree and the knot vector of the
    /// input and all of its control points and weights are `NaN`, so a failed
    /// elevation can never be mistaken for a valid curve.
    pub fn elevate_curve(curve: &NurbsCurve) -> NurbsCurve {
        match Self::try_elevate_curve(curve) {
            Ok(elevated) => elevated,
            Err(_) => NurbsCurve {
                knots: curve.knots.clone(),
                control_points: vec![
                    NurbsControlPoint {
                        point: [f64::NAN; 3],
                        weight: f64::NAN,
                    };
                    curve.control_points.len()
                ],
                degree: curve.degree,
            },
        }
    }

    /// Elevate a NURBS curve's degree by 1, reporting malformed input.
    ///
    /// Shorthand for [`try_elevate_curve_times`](Self::try_elevate_curve_times)
    /// with `times = 1`.
    ///
    /// # Errors
    ///
    /// See [`try_elevate_curve_times`](Self::try_elevate_curve_times).
    pub fn try_elevate_curve(curve: &NurbsCurve) -> Result<NurbsCurve> {
        Self::try_elevate_curve_times(curve, 1)
    }

    /// Elevate a NURBS curve's degree by `times` in a single pass of Piegl &
    /// Tiller Algorithm A5.9.
    ///
    /// The result is clamped, has degree `p + times`, represents exactly the
    /// same curve with the same parametrisation, and each distinct knot of the
    /// clamped input gains `times` copies. Curves whose knot vector is not clamped
    /// are clamped first, so the result covers the parameter domain
    /// `[knots[degree], knots[n_ctrl]]` of the input. `times = 0` returns the
    /// clamped curve itself.
    ///
    /// # Errors
    ///
    /// Fails for degree 0, for the malformed curves listed at
    /// [`BezierDecomposition::try_decompose`] and when the elevation itself
    /// breaks down numerically.
    pub fn try_elevate_curve_times(curve: &NurbsCurve, times: usize) -> Result<NurbsCurve> {
        if curve.degree == 0 {
            return Err(Error::unsupported(
                "degree elevation",
                "a curve of degree 0 cannot be elevated",
            ));
        }
        let clamped = ClampedCurve::new(curve)?;
        let elevated = FemDegreeElevation::try_elevate(&clamped.to_elevatable(), times)
            .map_err(|err| Error::general(err.to_string()))?;
        let control_points = elevated
            .control_points
            .iter()
            .zip(elevated.weights.iter())
            .map(|(&point, &weight)| NurbsControlPoint { point, weight })
            .collect();
        Ok(NurbsCurve {
            knots: NurbsKnotVector {
                knots: elevated.basis.knot,
            },
            control_points,
            degree: elevated.basis.degree,
        })
    }
}
