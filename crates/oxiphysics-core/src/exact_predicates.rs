// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact (adaptive-filter) geometric predicates for robust computational geometry.
//!
//! This module implements the orientation and in-sphere/in-circle predicates that
//! underpin Delaunay triangulation, convex-hull construction, and mesh generation.
//! Each predicate first evaluates a cheap `f64` *fast filter* together with a
//! conservative forward-error bound; only when the floating-point determinant is
//! too close to zero to trust the sign does the code fall back to an *exact*
//! evaluation carried out in arbitrary-precision floating-point *expansions*. The
//! exact path is therefore taken rarely, yet the reported sign is always provably
//! correct -- there are no robustness failures, no spurious degeneracies, and no
//! inconsistent answers across permutations of the inputs.
//!
//! The expansion arithmetic (multi-term nonoverlapping `f64` sequences, the
//! `grow`/`scale`/`fast-expansion-sum` operations, and the zero-elimination
//! variants used here) and the adaptive-filtering strategy follow:
//!
//! * Jonathan Richard Shewchuk, "Adaptive Precision Floating-Point Arithmetic and
//!   Fast Robust Geometric Predicates", Discrete & Computational Geometry
//!   18(3):305-363, 1997.
//!
//! All exact arithmetic is built from the error-free transforms in
//! [`crate::extended_precision`] ([`two_sum`], [`two_prod`], [`fast_two_sum`], and
//! a local two-difference transform). **No fused multiply-add (FMA) is used
//! anywhere**, so every result is bit-for-bit deterministic across platforms.
//!
//! The stage-A error-bound constants are those derived in Shewchuk §4; the current
//! implementation pairs each fast filter directly with the fully exact expansion
//! fallback. The intermediate adaptive refinement stages (the B/C/D stages of the
//! reference paper) are a future performance optimization -- their omission costs
//! some speed in near-degenerate cases but never correctness, because the exact
//! expansion is the authority for the sign.

use crate::extended_precision::{fast_two_sum, two_prod, two_sum};

// =====================================================================
// (A) Expansion arithmetic
// =====================================================================

/// Merges two nonoverlapping, increasing-magnitude expansions into their exact
/// sum, eliminating zero components.
///
/// Both `e` and `f` must be *nonoverlapping* expansions whose components are
/// sorted by increasing magnitude (the form produced by every routine in this
/// module). The result is another nonoverlapping increasing-magnitude expansion
/// equal to `sum(e) + sum(f)` exactly.
///
/// # Algorithm
///
/// This is Shewchuk's `O(m + n)` `fast_expansion_sum_zeroelim`. First the two
/// inputs are linearly merged by ascending component magnitude into a single
/// ordered stream `g` (an index-based two-pointer merge followed by draining the
/// remaining tail of whichever input is left). The merged stream is then swept
/// once with [`two_sum`], carrying the running high-order term `q` forward and
/// emitting every nonzero low-order error term. The leading carry is appended
/// last (and an all-zero result collapses to a single `0.0`), so the output is a
/// canonical nonoverlapping expansion with no spurious zero components.
fn fast_expansion_sum_zeroelim(e: &[f64], f: &[f64]) -> Vec<f64> {
    let m = e.len();
    let n = f.len();
    let mut g: Vec<f64> = Vec::with_capacity(m + n);

    // Index-based two-pointer merge by ascending |component|.
    let mut i = 0usize;
    let mut j = 0usize;
    while i < m && j < n {
        if e[i].abs() < f[j].abs() {
            g.push(e[i]);
            i += 1;
        } else {
            g.push(f[j]);
            j += 1;
        }
    }
    while i < m {
        g.push(e[i]);
        i += 1;
    }
    while j < n {
        g.push(f[j]);
        j += 1;
    }

    if g.is_empty() {
        return vec![0.0];
    }

    let mut q = g[0];
    let mut out: Vec<f64> = Vec::with_capacity(g.len());
    for &gk in &g[1..] {
        let (qnew, hh) = two_sum(q, gk);
        if hh != 0.0 {
            out.push(hh);
        }
        q = qnew;
    }
    if q != 0.0 || out.is_empty() {
        out.push(q);
    }
    out
}

/// Multiplies an expansion `e` by a single scalar `b`, eliminating zero
/// components.
///
/// `e` must be a nonoverlapping increasing-magnitude expansion; the result is
/// another such expansion equal to `sum(e) * b` exactly. This is Shewchuk's
/// `scale_expansion_zeroelim`: each component is multiplied by `b` with
/// [`two_prod`], and the doubled-length partial products are summed in order with
/// [`two_sum`]/[`fast_two_sum`], dropping zero error terms.
fn scale_expansion_zeroelim(e: &[f64], b: f64) -> Vec<f64> {
    let mut h = Vec::with_capacity(e.len() * 2);
    if e.is_empty() {
        return vec![0.0];
    }
    let (mut q, hh) = two_prod(e[0], b);
    if hh != 0.0 {
        h.push(hh);
    }
    for &e_i in &e[1..] {
        let (product1, product0) = two_prod(e_i, b);
        let (sum, hh1) = two_sum(q, product0);
        if hh1 != 0.0 {
            h.push(hh1);
        }
        let (qnew, hh2) = fast_two_sum(product1, sum);
        q = qnew;
        if hh2 != 0.0 {
            h.push(hh2);
        }
    }
    if q != 0.0 || h.is_empty() {
        h.push(q);
    }
    h
}

/// Returns the ordinary `f64` approximation of an expansion's value.
///
/// This is a plain left-fold sum of all components. Because an expansion's value
/// is `sum(e[i])` and the components are sorted by increasing magnitude, the fold
/// gives a good (though not correctly-rounded) estimate of the represented real
/// number. It is exposed publicly as a convenience for callers that only need an
/// approximate magnitude rather than an exact sign.
pub fn estimate(e: &[f64]) -> f64 {
    e.iter().sum()
}

/// Returns the exact sign of an expansion: `1`, `-1`, or `0`.
///
/// For a nonoverlapping increasing-magnitude expansion the component of largest
/// magnitude is the last nonzero one, and it dominates the sum, so the sign of
/// the whole expansion equals the sign of that component. This function scans
/// from the high-order end and returns `1` for the first strictly-positive
/// component, `-1` for the first strictly-negative component, and `0` if every
/// component is zero. This is the exact sign authority used by all predicates.
fn expansion_sign(e: &[f64]) -> i32 {
    for &c in e.iter().rev() {
        if c > 0.0 {
            return 1;
        }
        if c < 0.0 {
            return -1;
        }
    }
    0
}

/// Local two-difference error-free transform of `a - b`.
///
/// Returns `(s, err)` where `s = fl(a - b)` is the correctly-rounded difference
/// and `err` is the exact rounding error, so that `s + err == a - b` holds
/// exactly in real arithmetic. This mirrors Knuth's two-sum, adapted to
/// subtraction; it is reimplemented locally because the analogous transform in
/// [`crate::extended_precision`] is private. No FMA is used.
fn two_diff_local(a: f64, b: f64) -> (f64, f64) {
    let s = a - b;
    let bb = s - a;
    let err = (a - (s - bb)) - (b + bb);
    (s, err)
}

/// Represents the exact difference `a - b` as a two-term expansion.
///
/// [`two_diff_local`] yields `(s, err)` with `s` the larger-magnitude rounded
/// difference and `err` the small correction. An expansion stores its components
/// in *increasing* magnitude order, so the correction comes first and the rounded
/// difference second: the returned `[err, s]` is exactly `a - b` and is already a
/// valid nonoverlapping increasing-magnitude two-term expansion.
fn diff2(a: f64, b: f64) -> [f64; 2] {
    let (s, err) = two_diff_local(a, b);
    [err, s]
}

/// Negates every component of an expansion.
///
/// Floating-point negation is exact, so the result is a nonoverlapping
/// increasing-magnitude expansion equal to `-sum(e)`.
fn negate(e: &[f64]) -> Vec<f64> {
    e.iter().map(|&x| -x).collect()
}

/// Returns the exact difference `sum(e) - sum(f)` as an expansion.
///
/// Implemented as the fast expansion sum of `e` with the negation of `f`; both
/// inputs must be nonoverlapping increasing-magnitude expansions and the result
/// is another such expansion.
fn expansion_diff(e: &[f64], f: &[f64]) -> Vec<f64> {
    fast_expansion_sum_zeroelim(e, &negate(f))
}

/// Returns the exact product `sum(e) * sum(f)` as an expansion.
///
/// Both inputs must be nonoverlapping increasing-magnitude expansions. The
/// product is accumulated by scaling `e` by each component of `f` with
/// [`scale_expansion_zeroelim`] and summing the partial products with
/// [`fast_expansion_sum_zeroelim`]. The result is a nonoverlapping
/// increasing-magnitude expansion equal to the exact product.
fn expansion_mul(e: &[f64], f: &[f64]) -> Vec<f64> {
    let mut acc: Vec<f64> = vec![];
    for &f_j in f {
        let part = scale_expansion_zeroelim(e, f_j);
        acc = if acc.is_empty() {
            part
        } else {
            fast_expansion_sum_zeroelim(&acc, &part)
        };
    }
    if acc.is_empty() { vec![0.0] } else { acc }
}

/// Returns the exact `2x2` determinant `a*b - c*d` of four expansions.
///
/// Each argument is a nonoverlapping increasing-magnitude expansion; the two
/// products are formed with [`expansion_mul`] and subtracted with
/// [`expansion_diff`], yielding the exact determinant as an expansion.
fn two_by_two_det(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let ab = expansion_mul(a, b);
    let cd = expansion_mul(c, d);
    expansion_diff(&ab, &cd)
}

// =====================================================================
// (B) Orientation enum
// =====================================================================

/// The qualitative outcome of a geometric orientation or incidence predicate.
///
/// For an orientation test this records on which side of an oriented flat a query
/// point lies; for an in-circle/in-sphere test it records whether the query point
/// is inside, on, or outside the circle/sphere. The mapping from the underlying
/// signed determinant to this enum is fixed by [`from_sign`]/[`from_isign`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// The signed determinant is strictly positive.
    Positive,
    /// The signed determinant is strictly negative.
    Negative,
    /// The signed determinant is exactly zero (collinear/coplanar/cocircular/cospherical).
    Degenerate,
}

/// Maps the sign of an `f64` determinant onto an [`Orientation`].
///
/// Returns [`Orientation::Positive`] for a strictly positive value,
/// [`Orientation::Negative`] for a strictly negative value, and
/// [`Orientation::Degenerate`] otherwise (including for `0.0`, `-0.0`, and NaN).
fn from_sign(s: f64) -> Orientation {
    if s > 0.0 {
        Orientation::Positive
    } else if s < 0.0 {
        Orientation::Negative
    } else {
        Orientation::Degenerate
    }
}

/// Maps the sign of an `i32` determinant onto an [`Orientation`].
///
/// Returns [`Orientation::Positive`] for a positive value,
/// [`Orientation::Negative`] for a negative value, and
/// [`Orientation::Degenerate`] for zero.
fn from_isign(s: i32) -> Orientation {
    if s > 0 {
        Orientation::Positive
    } else if s < 0 {
        Orientation::Negative
    } else {
        Orientation::Degenerate
    }
}

// =====================================================================
// (C) Error-bound constants (Shewchuk §4)
// =====================================================================

/// Machine epsilon for f64 (2^-53), unit roundoff in Shewchuk's filters.
const EPSILON: f64 = 1.1102230246251565e-16; // 2^-53
/// Stage-A relative error bound for the `orient2d` fast filter (Shewchuk §4).
const CCWERRBOUND_A: f64 = (3.0 + 16.0 * EPSILON) * EPSILON; // orient2d stage A
/// Stage-A relative error bound for the `orient3d` fast filter (Shewchuk §4).
const O3DERRBOUND_A: f64 = (7.0 + 56.0 * EPSILON) * EPSILON; // orient3d stage A
/// Stage-A relative error bound for the `incircle` fast filter (Shewchuk §4).
const ICCERRBOUND_A: f64 = (10.0 + 96.0 * EPSILON) * EPSILON; // incircle stage A
/// Stage-A relative error bound for the `insphere` fast filter (Shewchuk §4).
const ISPERRBOUND_A: f64 = (16.0 + 224.0 * EPSILON) * EPSILON; // insphere stage A

// =====================================================================
// (D) Predicates
// =====================================================================

/// Determines the orientation of the ordered triple of 2D points `(pa, pb, pc)`.
///
/// Returns [`Orientation::Positive`] when the points make a counterclockwise
/// (left) turn, [`Orientation::Negative`] when they make a clockwise (right)
/// turn, and [`Orientation::Degenerate`] when they are collinear. The result is
/// the exact sign of the determinant
/// `(ax - cx)*(by - cy) - (ay - cy)*(bx - cx)`.
///
/// A cheap `f64` filter is tried first; only when its forward-error bound cannot
/// certify the sign does the call fall back to the exact expansion path, so the
/// answer is always correct yet usually fast.
pub fn orient2d(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> Orientation {
    let detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    let detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    let det = detleft - detright;
    let detsum = if detleft > 0.0 {
        if detright <= 0.0 {
            return from_sign(det);
        } else {
            detleft + detright
        }
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return from_sign(det);
        } else {
            -detleft - detright
        }
    } else {
        return from_sign(det);
    };
    let errbound = CCWERRBOUND_A * detsum;
    if det.abs() >= errbound {
        return from_sign(det);
    }
    // TODO(perf): intermediate adaptive C-stage refinement (Shewchuk) is a follow-up optimization; exact expansion guarantees correctness now.
    from_isign(orient2d_exact(pa, pb, pc))
}

/// Computes the exact sign of the 2D orientation determinant of `(pa, pb, pc)`.
///
/// Each coordinate difference is represented exactly with [`diff2`], and the
/// determinant `(ax - cx)*(by - cy) - (ay - cy)*(bx - cx)` is assembled exactly
/// via [`two_by_two_det`]. The returned value is `1`, `-1`, or `0`.
fn orient2d_exact(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> i32 {
    let acx = diff2(pa[0], pc[0]);
    let bcy = diff2(pb[1], pc[1]);
    let acy = diff2(pa[1], pc[1]);
    let bcx = diff2(pb[0], pc[0]);
    let det = two_by_two_det(&acx, &bcy, &acy, &bcx);
    expansion_sign(&det)
}

/// Determines the orientation of the ordered tuple of 3D points `(pa, pb, pc, pd)`.
///
/// Returns [`Orientation::Positive`] when `pd` lies below the plane through
/// `pa, pb, pc` (the three viewed counterclockwise from above), and
/// [`Orientation::Negative`] when `pd` lies above it; [`Orientation::Degenerate`]
/// means the four points are coplanar. The result is the exact sign of the `3x3`
/// determinant of the coordinate differences relative to `pd`.
///
/// A cheap `f64` filter with a forward-error bound is tried first, falling back
/// to the exact expansion path only when the sign cannot otherwise be certified.
pub fn orient3d(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> Orientation {
    let adx = pa[0] - pd[0];
    let ady = pa[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdx = pb[0] - pd[0];
    let bdy = pb[1] - pd[1];
    let bdz = pb[2] - pd[2];
    let cdx = pc[0] - pd[0];
    let cdy = pc[1] - pd[1];
    let cdz = pc[2] - pd[2];
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
        + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
        + (adxbdy.abs() + bdxady.abs()) * cdz.abs();
    let errbound = O3DERRBOUND_A * permanent;
    if det.abs() >= errbound {
        return from_sign(det);
    }
    // TODO(perf): intermediate adaptive C-stage refinement (Shewchuk) is a follow-up optimization; exact expansion guarantees correctness now.
    from_isign(orient3d_exact(pa, pb, pc, pd))
}

/// Computes the exact sign of the 3D orientation determinant of `(pa, pb, pc, pd)`.
///
/// The nine coordinate differences relative to `pd` are formed exactly with
/// [`diff2`]; the three `2x2` minors are built with [`two_by_two_det`], scaled by
/// the corresponding z-difference expansions via [`expansion_mul`], and summed
/// with [`fast_expansion_sum_zeroelim`]. The returned value is `1`, `-1`, or `0`.
fn orient3d_exact(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> i32 {
    let adx = diff2(pa[0], pd[0]);
    let ady = diff2(pa[1], pd[1]);
    let adz = diff2(pa[2], pd[2]);
    let bdx = diff2(pb[0], pd[0]);
    let bdy = diff2(pb[1], pd[1]);
    let bdz = diff2(pb[2], pd[2]);
    let cdx = diff2(pc[0], pd[0]);
    let cdy = diff2(pc[1], pd[1]);
    let cdz = diff2(pc[2], pd[2]);

    let m_a = two_by_two_det(&bdx, &cdy, &cdx, &bdy); // bdx*cdy - cdx*bdy
    let m_b = two_by_two_det(&cdx, &ady, &adx, &cdy); // cdx*ady - adx*cdy
    let m_c = two_by_two_det(&adx, &bdy, &bdx, &ady); // adx*bdy - bdx*ady

    let t_a = expansion_mul(&adz, &m_a);
    let t_b = expansion_mul(&bdz, &m_b);
    let t_c = expansion_mul(&cdz, &m_c);

    let s1 = fast_expansion_sum_zeroelim(&t_a, &t_b);
    let det = fast_expansion_sum_zeroelim(&s1, &t_c);
    expansion_sign(&det)
}

/// Tests the position of `pd` relative to the circle through `pa, pb, pc`.
///
/// Assuming `pa, pb, pc` are given in counterclockwise order (Shewchuk's
/// convention), returns [`Orientation::Positive`] when `pd` lies strictly
/// *inside* that circle, [`Orientation::Negative`] when it lies strictly
/// *outside*, and [`Orientation::Degenerate`] when the four points are
/// cocircular. The result is the exact sign of the standard in-circle
/// determinant.
///
/// A cheap `f64` filter with a forward-error bound is tried first; the exact
/// expansion path is taken only when the sign cannot otherwise be certified.
pub fn incircle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> Orientation {
    let adx = pa[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdx = pb[0] - pd[0];
    let bdy = pb[1] - pd[1];
    let cdx = pc[0] - pd[0];
    let cdy = pc[1] - pd[1];
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;
    let errbound = ICCERRBOUND_A * permanent;
    if det.abs() >= errbound {
        return from_sign(det);
    }
    // TODO(perf): incircle adaptive B/C/D stages (Shewchuk) are a follow-up; exact expansion is the correctness authority.
    from_isign(incircle_exact(pa, pb, pc, pd))
}

/// Computes the exact sign of the in-circle determinant of `(pa, pb, pc, pd)`.
///
/// The six coordinate differences relative to `pd` are formed exactly with
/// [`diff2`]. Each "lift" `(x^2 + y^2)` is built as an exact sum of two squared
/// expansions, the three `2x2` minors via [`two_by_two_det`], and the final
/// determinant by multiplying each lift by its minor with [`expansion_mul`] and
/// summing with [`fast_expansion_sum_zeroelim`]. The returned value is `1`, `-1`,
/// or `0`.
fn incircle_exact(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> i32 {
    let adx = diff2(pa[0], pd[0]);
    let ady = diff2(pa[1], pd[1]);
    let bdx = diff2(pb[0], pd[0]);
    let bdy = diff2(pb[1], pd[1]);
    let cdx = diff2(pc[0], pd[0]);
    let cdy = diff2(pc[1], pd[1]);

    let m_bc = two_by_two_det(&bdx, &cdy, &cdx, &bdy);
    let m_ca = two_by_two_det(&cdx, &ady, &adx, &cdy);
    let m_ab = two_by_two_det(&adx, &bdy, &bdx, &ady);

    let alift = fast_expansion_sum_zeroelim(&expansion_mul(&adx, &adx), &expansion_mul(&ady, &ady));
    let blift = fast_expansion_sum_zeroelim(&expansion_mul(&bdx, &bdx), &expansion_mul(&bdy, &bdy));
    let clift = fast_expansion_sum_zeroelim(&expansion_mul(&cdx, &cdx), &expansion_mul(&cdy, &cdy));

    let ta = expansion_mul(&alift, &m_bc);
    let tb = expansion_mul(&blift, &m_ca);
    let tc = expansion_mul(&clift, &m_ab);

    let s1 = fast_expansion_sum_zeroelim(&ta, &tb);
    let det = fast_expansion_sum_zeroelim(&s1, &tc);
    expansion_sign(&det)
}

/// Tests the position of a query point relative to the sphere through four others.
///
/// The argument `points` holds five 3D points: `points[0..4]` define a sphere and
/// are assumed positively oriented (consistent with [`orient3d`]), while
/// `points[4]` is the query point `pe`. Returns [`Orientation::Positive`] when
/// `pe` lies strictly *inside* the sphere, [`Orientation::Negative`] when it lies
/// strictly *outside*, and [`Orientation::Degenerate`] when all five points are
/// cospherical. The reported sign is the exact sign of the in-sphere determinant
/// in the same form as the fast filter, and the double-double oracle used in the
/// tests evaluates that identical formula -- it is the definitive source of truth
/// for the sign convention.
///
/// A cheap `f64` filter with a forward-error bound is tried first; the exact
/// expansion path is taken only when the sign cannot otherwise be certified.
pub fn insphere(points: [[f64; 3]; 5]) -> Orientation {
    let pa = points[0];
    let pb = points[1];
    let pc = points[2];
    let pd = points[3];
    let pe = points[4];
    let aex = pa[0] - pe[0];
    let aey = pa[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bex = pb[0] - pe[0];
    let bey = pb[1] - pe[1];
    let bez = pb[2] - pe[2];
    let cex = pc[0] - pe[0];
    let cey = pc[1] - pe[1];
    let cez = pc[2] - pe[2];
    let dex = pd[0] - pe[0];
    let dey = pd[1] - pe[1];
    let dez = pd[2] - pe[2];
    let ab = aex * bey - bex * aey;
    let bc = bex * cey - cex * bey;
    let cd = cex * dey - dex * cey;
    let da = dex * aey - aex * dey;
    let ac = aex * cey - cex * aey;
    let bd = bex * dey - dex * bey;
    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;
    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;
    let det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
    let abc_p = aez.abs() * bc.abs() + bez.abs() * ac.abs() + cez.abs() * ab.abs();
    let bcd_p = bez.abs() * cd.abs() + cez.abs() * bd.abs() + dez.abs() * bc.abs();
    let cda_p = cez.abs() * da.abs() + dez.abs() * ac.abs() + aez.abs() * cd.abs();
    let dab_p = dez.abs() * ab.abs() + aez.abs() * bd.abs() + bez.abs() * da.abs();
    let permanent = (dlift * abc_p + clift * dab_p) + (blift * cda_p + alift * bcd_p);
    let errbound = ISPERRBOUND_A * permanent;
    if det.abs() >= errbound {
        return from_sign(det);
    }
    // TODO(perf): insphere adaptive stages are a follow-up; exact expansion is the correctness authority.
    from_isign(insphere_exact(points))
}

/// Computes the exact sign of the in-sphere determinant for the five `points`.
///
/// The twelve coordinate differences relative to `pe` (= `points[4]`) are formed
/// exactly with [`diff2`]. The six `2x2` minors are built with
/// [`two_by_two_det`]; the four triple products `abc`, `bcd`, `cda`, `dab` are
/// assembled with the same sign pattern as the fast filter via [`expansion_mul`],
/// [`fast_expansion_sum_zeroelim`], and [`expansion_diff`]; each lift is an exact
/// sum of three squared expansions; and the final determinant
/// `(dlift*abc - clift*dab) + (blift*cda - alift*bcd)` is evaluated exactly. The
/// returned value is `1`, `-1`, or `0` and agrees with the double-double oracle.
fn insphere_exact(points: [[f64; 3]; 5]) -> i32 {
    let pa = points[0];
    let pb = points[1];
    let pc = points[2];
    let pd = points[3];
    let pe = points[4];

    let aex = diff2(pa[0], pe[0]);
    let aey = diff2(pa[1], pe[1]);
    let aez = diff2(pa[2], pe[2]);
    let bex = diff2(pb[0], pe[0]);
    let bey = diff2(pb[1], pe[1]);
    let bez = diff2(pb[2], pe[2]);
    let cex = diff2(pc[0], pe[0]);
    let cey = diff2(pc[1], pe[1]);
    let cez = diff2(pc[2], pe[2]);
    let dex = diff2(pd[0], pe[0]);
    let dey = diff2(pd[1], pe[1]);
    let dez = diff2(pd[2], pe[2]);

    let ab = two_by_two_det(&aex, &bey, &bex, &aey); // aex*bey - bex*aey
    let bc = two_by_two_det(&bex, &cey, &cex, &bey);
    let cd = two_by_two_det(&cex, &dey, &dex, &cey);
    let da = two_by_two_det(&dex, &aey, &aex, &dey);
    let ac = two_by_two_det(&aex, &cey, &cex, &aey);
    let bd = two_by_two_det(&bex, &dey, &dex, &bey);

    // abc = aez*bc - bez*ac + cez*ab
    let abc = {
        let t1 = expansion_mul(&aez, &bc);
        let t2 = expansion_mul(&bez, &ac);
        let t3 = expansion_mul(&cez, &ab);
        let s = expansion_diff(&t1, &t2);
        fast_expansion_sum_zeroelim(&s, &t3)
    };
    // bcd = bez*cd - cez*bd + dez*bc
    let bcd = {
        let t1 = expansion_mul(&bez, &cd);
        let t2 = expansion_mul(&cez, &bd);
        let t3 = expansion_mul(&dez, &bc);
        let s = expansion_diff(&t1, &t2);
        fast_expansion_sum_zeroelim(&s, &t3)
    };
    // cda = cez*da + dez*ac + aez*cd
    let cda = {
        let t1 = expansion_mul(&cez, &da);
        let t2 = expansion_mul(&dez, &ac);
        let t3 = expansion_mul(&aez, &cd);
        let s = fast_expansion_sum_zeroelim(&t1, &t2);
        fast_expansion_sum_zeroelim(&s, &t3)
    };
    // dab = dez*ab + aez*bd + bez*da
    let dab = {
        let t1 = expansion_mul(&dez, &ab);
        let t2 = expansion_mul(&aez, &bd);
        let t3 = expansion_mul(&bez, &da);
        let s = fast_expansion_sum_zeroelim(&t1, &t2);
        fast_expansion_sum_zeroelim(&s, &t3)
    };

    let alift = fast_expansion_sum_zeroelim(
        &fast_expansion_sum_zeroelim(&expansion_mul(&aex, &aex), &expansion_mul(&aey, &aey)),
        &expansion_mul(&aez, &aez),
    );
    let blift = fast_expansion_sum_zeroelim(
        &fast_expansion_sum_zeroelim(&expansion_mul(&bex, &bex), &expansion_mul(&bey, &bey)),
        &expansion_mul(&bez, &bez),
    );
    let clift = fast_expansion_sum_zeroelim(
        &fast_expansion_sum_zeroelim(&expansion_mul(&cex, &cex), &expansion_mul(&cey, &cey)),
        &expansion_mul(&cez, &cez),
    );
    let dlift = fast_expansion_sum_zeroelim(
        &fast_expansion_sum_zeroelim(&expansion_mul(&dex, &dex), &expansion_mul(&dey, &dey)),
        &expansion_mul(&dez, &dez),
    );

    let p1 = expansion_diff(&expansion_mul(&dlift, &abc), &expansion_mul(&clift, &dab));
    let p2 = expansion_diff(&expansion_mul(&blift, &cda), &expansion_mul(&alift, &bcd));
    let det = fast_expansion_sum_zeroelim(&p1, &p2);
    expansion_sign(&det)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extended_precision::Dd;

    /// Maps a double-double value to its sign as `1`, `-1`, or `0`.
    fn dd_sign(x: Dd) -> i32 {
        match x.to_f64().partial_cmp(&0.0) {
            Some(core::cmp::Ordering::Greater) => 1,
            Some(core::cmp::Ordering::Less) => -1,
            _ => 0,
        }
    }

    /// Double-double oracle for the 2D orientation determinant sign.
    fn orient2d_oracle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> i32 {
        let acx = Dd::from_f64(pa[0]) - Dd::from_f64(pc[0]);
        let bcy = Dd::from_f64(pb[1]) - Dd::from_f64(pc[1]);
        let acy = Dd::from_f64(pa[1]) - Dd::from_f64(pc[1]);
        let bcx = Dd::from_f64(pb[0]) - Dd::from_f64(pc[0]);
        let det = acx * bcy - acy * bcx;
        dd_sign(det)
    }

    /// Double-double oracle for the 3D orientation determinant sign.
    fn orient3d_oracle(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> i32 {
        let adx = Dd::from_f64(pa[0]) - Dd::from_f64(pd[0]);
        let ady = Dd::from_f64(pa[1]) - Dd::from_f64(pd[1]);
        let adz = Dd::from_f64(pa[2]) - Dd::from_f64(pd[2]);
        let bdx = Dd::from_f64(pb[0]) - Dd::from_f64(pd[0]);
        let bdy = Dd::from_f64(pb[1]) - Dd::from_f64(pd[1]);
        let bdz = Dd::from_f64(pb[2]) - Dd::from_f64(pd[2]);
        let cdx = Dd::from_f64(pc[0]) - Dd::from_f64(pd[0]);
        let cdy = Dd::from_f64(pc[1]) - Dd::from_f64(pd[1]);
        let cdz = Dd::from_f64(pc[2]) - Dd::from_f64(pd[2]);
        let det = adz * (bdx * cdy - cdx * bdy)
            + bdz * (cdx * ady - adx * cdy)
            + cdz * (adx * bdy - bdx * ady);
        dd_sign(det)
    }

    /// Double-double oracle for the in-circle determinant sign.
    fn incircle_oracle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> i32 {
        let adx = Dd::from_f64(pa[0]) - Dd::from_f64(pd[0]);
        let ady = Dd::from_f64(pa[1]) - Dd::from_f64(pd[1]);
        let bdx = Dd::from_f64(pb[0]) - Dd::from_f64(pd[0]);
        let bdy = Dd::from_f64(pb[1]) - Dd::from_f64(pd[1]);
        let cdx = Dd::from_f64(pc[0]) - Dd::from_f64(pd[0]);
        let cdy = Dd::from_f64(pc[1]) - Dd::from_f64(pd[1]);
        let alift = adx * adx + ady * ady;
        let blift = bdx * bdx + bdy * bdy;
        let clift = cdx * cdx + cdy * cdy;
        let det = alift * (bdx * cdy - cdx * bdy)
            + blift * (cdx * ady - adx * cdy)
            + clift * (adx * bdy - bdx * ady);
        dd_sign(det)
    }

    /// Double-double oracle for the in-sphere determinant sign.
    ///
    /// Mirrors the exact formula used by [`insphere`] stage A so that the fast
    /// filter, the exact path, and this oracle all share one sign convention.
    fn insphere_oracle(points: [[f64; 3]; 5]) -> i32 {
        let pa = points[0];
        let pb = points[1];
        let pc = points[2];
        let pd = points[3];
        let pe = points[4];
        let aex = Dd::from_f64(pa[0]) - Dd::from_f64(pe[0]);
        let aey = Dd::from_f64(pa[1]) - Dd::from_f64(pe[1]);
        let aez = Dd::from_f64(pa[2]) - Dd::from_f64(pe[2]);
        let bex = Dd::from_f64(pb[0]) - Dd::from_f64(pe[0]);
        let bey = Dd::from_f64(pb[1]) - Dd::from_f64(pe[1]);
        let bez = Dd::from_f64(pb[2]) - Dd::from_f64(pe[2]);
        let cex = Dd::from_f64(pc[0]) - Dd::from_f64(pe[0]);
        let cey = Dd::from_f64(pc[1]) - Dd::from_f64(pe[1]);
        let cez = Dd::from_f64(pc[2]) - Dd::from_f64(pe[2]);
        let dex = Dd::from_f64(pd[0]) - Dd::from_f64(pe[0]);
        let dey = Dd::from_f64(pd[1]) - Dd::from_f64(pe[1]);
        let dez = Dd::from_f64(pd[2]) - Dd::from_f64(pe[2]);
        let ab = aex * bey - bex * aey;
        let bc = bex * cey - cex * bey;
        let cd = cex * dey - dex * cey;
        let da = dex * aey - aex * dey;
        let ac = aex * cey - cex * aey;
        let bd = bex * dey - dex * bey;
        let abc = aez * bc - bez * ac + cez * ab;
        let bcd = bez * cd - cez * bd + dez * bc;
        let cda = cez * da + dez * ac + aez * cd;
        let dab = dez * ab + aez * bd + bez * da;
        let alift = aex * aex + aey * aey + aez * aez;
        let blift = bex * bex + bey * bey + bez * bez;
        let clift = cex * cex + cey * cey + cez * cez;
        let dlift = dex * dex + dey * dey + dez * dez;
        let det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
        dd_sign(det)
    }

    /// Returns the opposite orientation, leaving `Degenerate` unchanged.
    fn opp(o: Orientation) -> Orientation {
        match o {
            Orientation::Positive => Orientation::Negative,
            Orientation::Negative => Orientation::Positive,
            Orientation::Degenerate => Orientation::Degenerate,
        }
    }

    /// Minimal LCG (Knuth MMIX constants) producing f64 in [0,1) for deterministic torture tests.
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    // ---- (1) Permutation consistency ----

    #[test]
    fn orient2d_permutation_consistency() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        // CCW triangle is positive.
        assert_eq!(orient2d(a, b, c), Orientation::Positive);
        // Swapping two vertices flips the sign.
        assert_eq!(orient2d(b, a, c), opp(orient2d(a, b, c)));
        assert_eq!(orient2d(a, c, b), opp(orient2d(a, b, c)));
        assert_eq!(orient2d(c, b, a), opp(orient2d(a, b, c)));
        // A cyclic (even) permutation preserves the sign.
        assert_eq!(orient2d(b, c, a), orient2d(a, b, c));
        assert_eq!(orient2d(c, a, b), orient2d(a, b, c));
    }

    #[test]
    fn orient3d_permutation_consistency() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        let base = orient3d(a, b, c, d);
        assert_ne!(base, Orientation::Degenerate);
        // Swapping two vertices flips the sign.
        assert_eq!(orient3d(b, a, c, d), opp(base));
        assert_eq!(orient3d(a, c, b, d), opp(base));
        assert_eq!(orient3d(a, b, d, c), opp(base));
        // An even permutation (rotate b->c->d->b) preserves the sign.
        assert_eq!(orient3d(a, c, d, b), base);
    }

    #[test]
    fn incircle_permutation_consistency() {
        // Non-cocircular configuration: pa,pb,pc CCW, pd off the circle.
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        let d = [0.3, 0.3];
        let base = incircle(a, b, c, d);
        assert_ne!(base, Orientation::Degenerate);
        // Swapping two of pa,pb,pc flips the sign.
        assert_eq!(incircle(b, a, c, d), opp(base));
        assert_eq!(incircle(a, c, b, d), opp(base));
        assert_eq!(incircle(c, b, a, d), opp(base));
    }

    #[test]
    fn insphere_permutation_consistency() {
        // Non-cospherical configuration of five 3D points.
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.2, 0.2, 0.2],
        ];
        let base = insphere(points);
        assert_ne!(base, Orientation::Degenerate);
        // An odd permutation (single swap) of the first four points flips the sign.
        let mut swapped = points;
        swapped.swap(0, 1);
        assert_eq!(insphere(swapped), opp(base));
        let mut swapped2 = points;
        swapped2.swap(1, 2);
        assert_eq!(insphere(swapped2), opp(base));
    }

    // ---- (2) Torture tests vs Dd oracle ----

    #[test]
    fn torture_orient2d() {
        let mut rng = Lcg(0x1234_5678_9abc_def0);
        for i in 0..100_000u64 {
            // Power-of-two scale.
            let exp = (rng.next_f64() * 40.0) as i32 - 20;
            let s = 2f64.powi(exp);
            let k = ((rng.next_f64() * 8.0) as i64 + 1) as f64;
            let mut pa = [s, s];
            let mut pb = [s * 2.0, s * 2.0];
            let mut pc = [s * k, s * k];
            // Nudge ONE coordinate of one point by +/- 1 ulp.
            let which = (rng.next_f64() * 6.0) as usize % 6;
            let up = rng.next_f64() < 0.5;
            let target: &mut f64 = match which {
                0 => &mut pa[0],
                1 => &mut pa[1],
                2 => &mut pb[0],
                3 => &mut pb[1],
                4 => &mut pc[0],
                _ => &mut pc[1],
            };
            *target = if up {
                target.next_up()
            } else {
                target.next_down()
            };
            let got = orient2d(pa, pb, pc);
            let oracle = from_isign(orient2d_oracle(pa, pb, pc));
            assert_eq!(
                got, oracle,
                "mismatch at iter {i}: pa={pa:?} pb={pb:?} pc={pc:?}"
            );
        }
    }

    #[test]
    fn torture_incircle() {
        let mut rng = Lcg(0x0bad_c0de_dead_beef);
        for i in 0..100_000u64 {
            let exp = (rng.next_f64() * 30.0) as i32 - 10;
            let s = 2f64.powi(exp);
            // Axis-aligned square (cocircular).
            let mut pa = [0.0, 0.0];
            let mut pb = [s, 0.0];
            let mut pc = [s, s];
            let mut pd = [0.0, s];
            let which = (rng.next_f64() * 8.0) as usize % 8;
            let up = rng.next_f64() < 0.5;
            let target: &mut f64 = match which {
                0 => &mut pa[0],
                1 => &mut pa[1],
                2 => &mut pb[0],
                3 => &mut pb[1],
                4 => &mut pc[0],
                5 => &mut pc[1],
                6 => &mut pd[0],
                _ => &mut pd[1],
            };
            *target = if up {
                target.next_up()
            } else {
                target.next_down()
            };
            let got = incircle(pa, pb, pc, pd);
            let oracle = from_isign(incircle_oracle(pa, pb, pc, pd));
            assert_eq!(
                got, oracle,
                "mismatch at iter {i}: pa={pa:?} pb={pb:?} pc={pc:?} pd={pd:?}"
            );
        }
    }

    #[test]
    fn torture_orient3d() {
        let mut rng = Lcg(0xfeed_face_cafe_b00b);
        for i in 0..20_000u64 {
            let exp = (rng.next_f64() * 20.0) as i32 - 5;
            let s = 2f64.powi(exp);
            // Coplanar points (all z equal).
            let z = s * 3.0;
            let mut pa = [0.0, 0.0, z];
            let mut pb = [s, 0.0, z];
            let mut pc = [0.0, s, z];
            let mut pd = [s, s, z];
            let which = (rng.next_f64() * 12.0) as usize % 12;
            let up = rng.next_f64() < 0.5;
            let target: &mut f64 = match which {
                0 => &mut pa[0],
                1 => &mut pa[1],
                2 => &mut pa[2],
                3 => &mut pb[0],
                4 => &mut pb[1],
                5 => &mut pb[2],
                6 => &mut pc[0],
                7 => &mut pc[1],
                8 => &mut pc[2],
                9 => &mut pd[0],
                10 => &mut pd[1],
                _ => &mut pd[2],
            };
            *target = if up {
                target.next_up()
            } else {
                target.next_down()
            };
            let got = orient3d(pa, pb, pc, pd);
            let oracle = from_isign(orient3d_oracle(pa, pb, pc, pd));
            assert_eq!(
                got, oracle,
                "mismatch at iter {i}: pa={pa:?} pb={pb:?} pc={pc:?} pd={pd:?}"
            );
        }
    }

    #[test]
    fn torture_insphere() {
        let mut rng = Lcg(0x5151_5151_2727_2727);
        for i in 0..20_000u64 {
            let exp = (rng.next_f64() * 12.0) as i32 - 2;
            let s = 2f64.powi(exp);
            // Axis-aligned cospherical config on a sphere of radius s about origin.
            let mut points = [
                [s, 0.0, 0.0],
                [-s, 0.0, 0.0],
                [0.0, s, 0.0],
                [0.0, 0.0, s],
                [0.0, -s, 0.0],
            ];
            let pt = (rng.next_f64() * 5.0) as usize % 5;
            let co = (rng.next_f64() * 3.0) as usize % 3;
            let up = rng.next_f64() < 0.5;
            let target = &mut points[pt][co];
            *target = if up {
                target.next_up()
            } else {
                target.next_down()
            };
            let got = insphere(points);
            let oracle = from_isign(insphere_oracle(points));
            assert_eq!(got, oracle, "mismatch at iter {i}: points={points:?}");
        }
    }

    /// Heavy 1,000,000-iteration torture run for `orient2d` and `incircle`.
    ///
    /// Skipped unless the `OXIPHYSICS_TORTURE_1M` environment variable is set, so
    /// ordinary test runs stay fast while the exhaustive sweep is still available
    /// (and always compiled). Every iteration asserts the predicate matches the
    /// double-double oracle exactly.
    #[test]
    fn torture_1m_orient2d_incircle() {
        if std::env::var("OXIPHYSICS_TORTURE_1M").is_err() {
            return;
        }
        let mut rng = Lcg(0xa5a5_5a5a_3c3c_c3c3);
        for i in 0..1_000_000u64 {
            // orient2d sweep.
            {
                let exp = (rng.next_f64() * 40.0) as i32 - 20;
                let s = 2f64.powi(exp);
                let k = ((rng.next_f64() * 8.0) as i64 + 1) as f64;
                let mut pa = [s, s];
                let mut pb = [s * 2.0, s * 2.0];
                let mut pc = [s * k, s * k];
                let which = (rng.next_f64() * 6.0) as usize % 6;
                let up = rng.next_f64() < 0.5;
                let target: &mut f64 = match which {
                    0 => &mut pa[0],
                    1 => &mut pa[1],
                    2 => &mut pb[0],
                    3 => &mut pb[1],
                    4 => &mut pc[0],
                    _ => &mut pc[1],
                };
                *target = if up {
                    target.next_up()
                } else {
                    target.next_down()
                };
                assert_eq!(
                    orient2d(pa, pb, pc),
                    from_isign(orient2d_oracle(pa, pb, pc)),
                    "orient2d mismatch at iter {i}"
                );
            }
            // incircle sweep.
            {
                let exp = (rng.next_f64() * 30.0) as i32 - 10;
                let s = 2f64.powi(exp);
                let mut pa = [0.0, 0.0];
                let mut pb = [s, 0.0];
                let mut pc = [s, s];
                let mut pd = [0.0, s];
                let which = (rng.next_f64() * 8.0) as usize % 8;
                let up = rng.next_f64() < 0.5;
                let target: &mut f64 = match which {
                    0 => &mut pa[0],
                    1 => &mut pa[1],
                    2 => &mut pb[0],
                    3 => &mut pb[1],
                    4 => &mut pc[0],
                    5 => &mut pc[1],
                    6 => &mut pd[0],
                    _ => &mut pd[1],
                };
                *target = if up {
                    target.next_up()
                } else {
                    target.next_down()
                };
                assert_eq!(
                    incircle(pa, pb, pc, pd),
                    from_isign(incircle_oracle(pa, pb, pc, pd)),
                    "incircle mismatch at iter {i}"
                );
            }
        }
    }

    // ---- (3) Classic naive-f64 failure ----

    #[test]
    fn classic_naive_f64_failure() {
        // Truly collinear points on the line y = x.
        let a = [0.5, 0.5];
        let b = [12.0, 12.0];
        let c = [24.0, 24.0];
        // The exact predicate must report exact collinearity.
        assert_eq!(orient2d(a, b, c), Orientation::Degenerate);

        // A near-collinear point off the line y = x: perturbing the x-coordinate
        // of a point near `b` yields a tiny-but-nonzero naive determinant. The
        // exact predicate must agree with the Dd oracle on this off-line point.
        // (Note: nudging the *far* corner c = [24, 24] by one ulp in x leaves the
        // naive determinant rounding back to exactly 0.0 -- precisely the
        // catastrophic misrounding the exact predicate defends against -- so the
        // demonstration point is chosen near `b` where the naive value survives.)
        let c_near = [(12.0f64).next_up(), 12.0];
        let naive_c_near =
            (a[0] - c_near[0]) * (b[1] - c_near[1]) - (a[1] - c_near[1]) * (b[0] - c_near[0]);
        // Demonstrate the naive determinant is nonzero for this off-line point.
        assert!(
            naive_c_near != 0.0,
            "near-collinear point should have nonzero naive determinant"
        );
        assert_eq!(
            orient2d(a, b, c_near),
            from_isign(orient2d_oracle(a, b, c_near))
        );
        // c_near is genuinely NOT collinear, so the exact predicate is decisive.
        assert_ne!(orient2d(a, b, c_near), Orientation::Degenerate);

        // Build a batch of triples (collinear + near-collinear) and verify the
        // exact predicate agrees with the Dd oracle across the whole batch, and
        // that at least one naive determinant in the batch is nonzero while the
        // genuinely collinear triple stays Degenerate.
        let batch: [[[f64; 2]; 3]; 4] = [
            [a, b, c],
            [a, b, c_near],
            [a, b, [12.0, (12.0f64).next_down()]],
            [a, [12.0, (12.0f64).next_down()], c],
        ];
        let mut any_naive_nonzero = false;
        for triple in &batch {
            let [pa, pb, pc] = *triple;
            assert_eq!(
                orient2d(pa, pb, pc),
                from_isign(orient2d_oracle(pa, pb, pc))
            );
            let naive = (pa[0] - pc[0]) * (pb[1] - pc[1]) - (pa[1] - pc[1]) * (pb[0] - pc[0]);
            if naive != 0.0 {
                any_naive_nonzero = true;
            }
        }
        assert!(
            any_naive_nonzero,
            "expected at least one nonzero naive determinant"
        );
        // The genuinely collinear triple stays exactly Degenerate.
        assert_eq!(orient2d(a, b, c), Orientation::Degenerate);
    }

    // ---- (4) Degenerate exact inputs ----

    #[test]
    fn degenerate_orient2d() {
        assert_eq!(
            orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]),
            Orientation::Degenerate
        );
    }

    #[test]
    fn degenerate_orient3d() {
        assert_eq!(
            orient3d(
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0]
            ),
            Orientation::Degenerate
        );
    }

    #[test]
    fn degenerate_incircle() {
        // Unit square corners (cocircular), pa,pb,pc CCW.
        assert_eq!(
            incircle([0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]),
            Orientation::Degenerate
        );
    }

    #[test]
    fn degenerate_insphere() {
        // Five cospherical points (radius 1 about origin).
        let points = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        assert_eq!(insphere(points), Orientation::Degenerate);
        assert_eq!(insphere_oracle(points), 0);
    }

    // ---- Expansion arithmetic sanity (exercises estimate + EFT helpers) ----

    #[test]
    fn expansion_estimate_and_sign() {
        // (3 - 1) exactly as an expansion.
        let e = diff2(3.0, 1.0);
        assert_eq!(estimate(&e), 2.0);
        assert_eq!(expansion_sign(&e), 1);
        let neg = negate(&e);
        assert_eq!(expansion_sign(&neg), -1);
        // Exact product 2 * 3 = 6.
        let prod = expansion_mul(&diff2(3.0, 1.0), &diff2(5.0, 2.0));
        assert_eq!(estimate(&prod), 6.0);
        assert_eq!(expansion_sign(&prod), 1);
        // Zero expansion.
        let zero = expansion_diff(&e, &e);
        assert_eq!(expansion_sign(&zero), 0);
        assert_eq!(estimate(&zero), 0.0);
    }
}
