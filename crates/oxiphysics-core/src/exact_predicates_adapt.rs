// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive intermediate stages (B/C) for the robust geometric predicates.
//!
//! This file is a private sibling of [`super`] (`exact_predicates`), wired in
//! with `#[path = "exact_predicates_adapt.rs"] mod adapt;`. It implements the
//! incremental adaptive-precision refinement of Shewchuk's predicates -- the
//! stages that sit *between* the cheap stage-A fast filter and the fully exact
//! expansion fallback:
//!
//! * Jonathan Richard Shewchuk, "Adaptive Precision Floating-Point Arithmetic and
//!   Fast Robust Geometric Predicates", Discrete & Computational Geometry
//!   18(3):305-363, 1997, §4-5.
//!
//! Each `*_adapt` routine first evaluates the determinant of the *rounded*
//! coordinate differences as an exact `f64` expansion (the **stage-B** value),
//! tested against the stage-B error bound. If that does not certify the sign, it
//! adds the first-order tail correction in plain `f64` and tests the
//! **result** bound (`RESULTERRBOUND`) together with the stage-C bound. Only when
//! *both* incremental stages remain uncertain does it fall through to the
//! parent module's fully exact expansion (`super::*_exact`), which is the final
//! authority for the sign. Consequently the adaptive path always returns exactly
//! the same sign as the exact expansion, never a looser or wrong answer -- it is
//! only *faster* on the easy near-degenerate cases, which it resolves without
//! paying for the full exact computation.
//!
//! All arithmetic is built from the error-free transforms of
//! [`crate::extended_precision`] and the parent module's expansion routines; **no
//! fused multiply-add (FMA) is used**, so results are bit-for-bit deterministic.

use super::{
    EPSILON, estimate, expansion_sign, fast_expansion_sum_zeroelim, incircle_exact, insphere_exact,
    orient2d_exact, orient3d_exact, scale_expansion_zeroelim, two_diff_local,
};
use crate::extended_precision::{two_prod, two_sum};

// =====================================================================
// Stage-B / stage-C / result error-bound constants (Shewchuk §4-5).
//
// Each is the canonical `(b + c * EPSILON) * EPSILON`-form (or
// `(b + c * EPSILON) * EPSILON * EPSILON` for the second-order C bounds)
// constant from Shewchuk's `exactinit`, with `EPSILON = 2^-53` the f64 unit
// roundoff. They are NOT loosened: they are exactly the reference constants.
// =====================================================================

/// Stage-B relative error bound for `orient2d` (Shewchuk `ccwerrboundB`).
const CCWERRBOUND_B: f64 = (2.0 + 12.0 * EPSILON) * EPSILON;
/// Stage-C relative error bound for `orient2d` (Shewchuk `ccwerrboundC`).
const CCWERRBOUND_C: f64 = (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
/// Result error bound shared by all predicates (Shewchuk `resulterrbound`).
const RESULTERRBOUND: f64 = (3.0 + 8.0 * EPSILON) * EPSILON;
/// Stage-B relative error bound for `orient3d` (Shewchuk `o3derrboundB`).
const O3DERRBOUND_B: f64 = (3.0 + 28.0 * EPSILON) * EPSILON;
/// Stage-C relative error bound for `orient3d` (Shewchuk `o3derrboundC`).
const O3DERRBOUND_C: f64 = (26.0 + 288.0 * EPSILON) * EPSILON * EPSILON;
/// Stage-B relative error bound for `incircle` (Shewchuk `iccerrboundB`).
const ICCERRBOUND_B: f64 = (4.0 + 48.0 * EPSILON) * EPSILON;
/// Stage-C relative error bound for `incircle` (Shewchuk `iccerrboundC`).
const ICCERRBOUND_C: f64 = (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;
/// Stage-B relative error bound for `insphere` (Shewchuk `isperrboundB`).
const ISPERRBOUND_B: f64 = (5.0 + 72.0 * EPSILON) * EPSILON;
/// Stage-C relative error bound for `insphere` (Shewchuk `isperrboundC`).
const ISPERRBOUND_C: f64 = (71.0 + 1408.0 * EPSILON) * EPSILON * EPSILON;

// =====================================================================
// Low-level error-free transforms used by the stage-B expansions.
// =====================================================================

/// Returns the sign of an `f64` as `1`, `-1`, or `0`.
///
/// Used to map a certified stage-B/stage-C determinant estimate onto an integer
/// sign in the same convention as [`super::expansion_sign`].
fn sign_of(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// Shewchuk's `Two_One_Diff`: exact three-term expansion of `(a1 + a0) - b`.
///
/// `(a1, a0)` is a two-term expansion (high then low) and `b` a single `f64`.
/// Returns `(x2, x1, x0)` in decreasing magnitude, a nonoverlapping expansion
/// whose components sum exactly to `a1 + a0 - b`. No FMA is used.
fn two_one_diff(a1: f64, a0: f64, b: f64) -> (f64, f64, f64) {
    let (i, x0) = two_diff_local(a0, b);
    let (x2, x1) = two_sum(a1, i);
    (x2, x1, x0)
}

/// Shewchuk's `Two_Two_Diff`: exact four-term expansion of `(a1+a0) - (b1+b0)`.
///
/// Both `(a1, a0)` and `(b1, b0)` are two-term expansions (high then low). The
/// returned `(x3, x2, x1, x0)` is a nonoverlapping expansion in decreasing
/// magnitude summing exactly to `(a1 + a0) - (b1 + b0)`. No FMA is used.
fn two_two_diff(a1: f64, a0: f64, b1: f64, b0: f64) -> (f64, f64, f64, f64) {
    let (j, j0, x0) = two_one_diff(a1, a0, b0);
    let (x3, x2, x1) = two_one_diff(j, j0, b1);
    (x3, x2, x1, x0)
}

/// Exact four-term expansion of the `2x2` cross-difference `a*b - c*d`.
///
/// Each factor is a single `f64`. The two products are formed with Dekker's
/// [`two_prod`] error-free transform and subtracted with [`two_two_diff`],
/// yielding `[x0, x1, x2, x3]` -- a nonoverlapping increasing-magnitude
/// expansion equal to `a*b - c*d` exactly. Its leading (largest) component is
/// `[3]`, which the in-sphere stage-C correction reuses.
fn cross_diff_expansion(a: f64, b: f64, c: f64, d: f64) -> [f64; 4] {
    let (hi1, lo1) = two_prod(a, b);
    let (hi2, lo2) = two_prod(c, d);
    let (x3, x2, x1, x0) = two_two_diff(hi1, lo1, hi2, lo2);
    [x0, x1, x2, x3]
}

/// Exact expansion `k0*e0 + k1*e1 + k2*e2` of three scaled expansions.
///
/// Each `e_i` is a nonoverlapping increasing-magnitude expansion and `k_i` a
/// scalar; the result is another such expansion. Used to assemble the in-sphere
/// triple products (`abc`, `bcd`, `cda`, `dab`) from the six `2x2` minors.
fn triple_sum(e0: &[f64], k0: f64, e1: &[f64], k1: f64, e2: &[f64], k2: f64) -> Vec<f64> {
    let t0 = scale_expansion_zeroelim(e0, k0);
    let t1 = scale_expansion_zeroelim(e1, k1);
    let t2 = scale_expansion_zeroelim(e2, k2);
    let s = fast_expansion_sum_zeroelim(&t0, &t1);
    fast_expansion_sum_zeroelim(&s, &t2)
}

/// Exact expansion `(dx^2 + dy^2) * e` (the planar "lift" scaling).
///
/// Each coordinate-squared scaling is performed by scaling `e` twice with
/// [`scale_expansion_zeroelim`], so the squared factor is represented exactly;
/// the two contributions are summed with [`fast_expansion_sum_zeroelim`].
fn lift_scale_2d(e: &[f64], dx: f64, dy: f64) -> Vec<f64> {
    let xd = scale_expansion_zeroelim(&scale_expansion_zeroelim(e, dx), dx);
    let yd = scale_expansion_zeroelim(&scale_expansion_zeroelim(e, dy), dy);
    fast_expansion_sum_zeroelim(&xd, &yd)
}

/// Exact expansion `sign * (dx^2 + dy^2 + dz^2) * e` (the spatial "lift").
///
/// `sign` must be `+1.0` or `-1.0`; it is folded into the *second* scaling of
/// each axis so that, e.g., `sign = -1.0` yields `-dx^2 * e` exactly. The three
/// axis contributions are summed with [`fast_expansion_sum_zeroelim`].
fn lift_scale_3d(e: &[f64], dx: f64, dy: f64, dz: f64, sign: f64) -> Vec<f64> {
    let xd = scale_expansion_zeroelim(&scale_expansion_zeroelim(e, dx), sign * dx);
    let yd = scale_expansion_zeroelim(&scale_expansion_zeroelim(e, dy), sign * dy);
    let zd = scale_expansion_zeroelim(&scale_expansion_zeroelim(e, dz), sign * dz);
    let xy = fast_expansion_sum_zeroelim(&xd, &yd);
    fast_expansion_sum_zeroelim(&xy, &zd)
}

// =====================================================================
// Adaptive predicates (stages B and C).
// =====================================================================

/// Adaptive `orient2d` (Shewchuk's `orient2dadapt`), returning a `1/-1/0` sign.
///
/// `detsum` is the magnitude estimate from the stage-A filter. The routine
/// computes the determinant of the rounded differences as an exact four-term
/// expansion (stage B), tests it against [`CCWERRBOUND_B`] `* detsum`, then -- if
/// the coordinate differences carry nonzero tails -- adds the first-order tail
/// correction and tests the result bound. It defers to [`orient2d_exact`] only
/// when both incremental stages remain uncertain, so the returned sign always
/// equals the exact sign.
pub(super) fn orient2d_adapt(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], detsum: f64) -> i32 {
    let (acx, acxtail) = two_diff_local(pa[0], pc[0]);
    let (bcx, bcxtail) = two_diff_local(pb[0], pc[0]);
    let (acy, acytail) = two_diff_local(pa[1], pc[1]);
    let (bcy, bcytail) = two_diff_local(pb[1], pc[1]);

    // Stage B: exact 4-term expansion of acx*bcy - acy*bcx (rounded diffs).
    let (dl_hi, dl_lo) = two_prod(acx, bcy);
    let (dr_hi, dr_lo) = two_prod(acy, bcx);
    let (b3, b2, b1, b0) = two_two_diff(dl_hi, dl_lo, dr_hi, dr_lo);
    let b = [b0, b1, b2, b3];

    let mut det = estimate(&b);
    if det.abs() >= CCWERRBOUND_B * detsum {
        return sign_of(det);
    }

    // No tails => the rounded differences are exact => `b` is the exact result.
    if acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0 {
        return expansion_sign(&b);
    }

    // Stage C: add the first-order tail correction, test the result bound.
    let errbound = CCWERRBOUND_C * detsum + RESULTERRBOUND * det.abs();
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if det.abs() >= errbound {
        return sign_of(det);
    }

    // Still uncertain: the exact expansion is the authority.
    orient2d_exact(pa, pb, pc)
}

/// Adaptive `orient3d` (Shewchuk's `orient3dadapt`), returning a `1/-1/0` sign.
///
/// `permanent` is the stage-A magnitude estimate. Stage B forms the three signed
/// `2x2` minors as exact four-term expansions, scales each by the corresponding
/// `z`-difference, and sums them; the estimate is tested against
/// [`O3DERRBOUND_B`] `* permanent`. Stage C adds the first-order tail correction
/// and tests the result bound. It defers to [`orient3d_exact`] only when both
/// stages remain uncertain.
pub(super) fn orient3d_adapt(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
    permanent: f64,
) -> i32 {
    let (adx, adxtail) = two_diff_local(pa[0], pd[0]);
    let (ady, adytail) = two_diff_local(pa[1], pd[1]);
    let (adz, adztail) = two_diff_local(pa[2], pd[2]);
    let (bdx, bdxtail) = two_diff_local(pb[0], pd[0]);
    let (bdy, bdytail) = two_diff_local(pb[1], pd[1]);
    let (bdz, bdztail) = two_diff_local(pb[2], pd[2]);
    let (cdx, cdxtail) = two_diff_local(pc[0], pd[0]);
    let (cdy, cdytail) = two_diff_local(pc[1], pd[1]);
    let (cdz, cdztail) = two_diff_local(pc[2], pd[2]);

    // Stage B: minors of the rounded differences, scaled and summed.
    let bc = cross_diff_expansion(bdx, cdy, cdx, bdy); // bdx*cdy - cdx*bdy
    let ca = cross_diff_expansion(cdx, ady, adx, cdy); // cdx*ady - adx*cdy
    let ab = cross_diff_expansion(adx, bdy, bdx, ady); // adx*bdy - bdx*ady
    let adet = scale_expansion_zeroelim(&bc, adz);
    let bdet = scale_expansion_zeroelim(&ca, bdz);
    let cdet = scale_expansion_zeroelim(&ab, cdz);
    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet);
    let fin1 = fast_expansion_sum_zeroelim(&abdet, &cdet);

    let mut det = estimate(&fin1);
    if det.abs() >= O3DERRBOUND_B * permanent {
        return sign_of(det);
    }

    if adxtail == 0.0
        && bdxtail == 0.0
        && cdxtail == 0.0
        && adytail == 0.0
        && bdytail == 0.0
        && cdytail == 0.0
        && adztail == 0.0
        && bdztail == 0.0
        && cdztail == 0.0
    {
        return expansion_sign(&fin1);
    }

    let errbound = O3DERRBOUND_C * permanent + RESULTERRBOUND * det.abs();
    det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + adztail * (bdx * cdy - bdy * cdx))
        + (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + bdztail * (cdx * ady - cdy * adx))
        + (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + cdztail * (adx * bdy - ady * bdx));
    if det.abs() >= errbound {
        return sign_of(det);
    }

    orient3d_exact(pa, pb, pc, pd)
}

/// Adaptive `incircle` (Shewchuk's `incircleadapt`), returning a `1/-1/0` sign.
///
/// `permanent` is the stage-A magnitude estimate. Stage B forms the three signed
/// `2x2` minors as exact expansions, multiplies each by the exact squared
/// "lift" `(dx^2 + dy^2)` of the opposite vertex, and sums them; the estimate is
/// tested against [`ICCERRBOUND_B`] `* permanent`. Stage C adds the first-order
/// tail correction (including the `2 * (d . dtail)` lift derivative) and tests
/// the result bound. It defers to [`incircle_exact`] only when both stages
/// remain uncertain.
pub(super) fn incircle_adapt(
    pa: [f64; 2],
    pb: [f64; 2],
    pc: [f64; 2],
    pd: [f64; 2],
    permanent: f64,
) -> i32 {
    let (adx, adxtail) = two_diff_local(pa[0], pd[0]);
    let (ady, adytail) = two_diff_local(pa[1], pd[1]);
    let (bdx, bdxtail) = two_diff_local(pb[0], pd[0]);
    let (bdy, bdytail) = two_diff_local(pb[1], pd[1]);
    let (cdx, cdxtail) = two_diff_local(pc[0], pd[0]);
    let (cdy, cdytail) = two_diff_local(pc[1], pd[1]);

    // Stage B: lift-weighted minors of the rounded differences.
    let bc = cross_diff_expansion(bdx, cdy, cdx, bdy); // bdx*cdy - cdx*bdy
    let ca = cross_diff_expansion(cdx, ady, adx, cdy); // cdx*ady - adx*cdy
    let ab = cross_diff_expansion(adx, bdy, bdx, ady); // adx*bdy - bdx*ady
    let adet = lift_scale_2d(&bc, adx, ady); // (adx^2 + ady^2) * bc
    let bdet = lift_scale_2d(&ca, bdx, bdy); // (bdx^2 + bdy^2) * ca
    let cdet = lift_scale_2d(&ab, cdx, cdy); // (cdx^2 + cdy^2) * ab
    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet);
    let fin1 = fast_expansion_sum_zeroelim(&abdet, &cdet);

    let mut det = estimate(&fin1);
    if det.abs() >= ICCERRBOUND_B * permanent {
        return sign_of(det);
    }

    if adxtail == 0.0
        && bdxtail == 0.0
        && cdxtail == 0.0
        && adytail == 0.0
        && bdytail == 0.0
        && cdytail == 0.0
    {
        return expansion_sign(&fin1);
    }

    let errbound = ICCERRBOUND_C * permanent + RESULTERRBOUND * det.abs();
    det += ((adx * adx + ady * ady)
        * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
        + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
        + ((bdx * bdx + bdy * bdy)
            * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
            + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
        + ((cdx * cdx + cdy * cdy)
            * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
            + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
    if det.abs() >= errbound {
        return sign_of(det);
    }

    incircle_exact(pa, pb, pc, pd)
}

/// Adaptive `insphere` (Shewchuk's `insphereadapt`), returning a `1/-1/0` sign.
///
/// `points[0..4]` define the sphere and `points[4]` is the query point;
/// `permanent` is the stage-A magnitude estimate. Stage B forms the six `2x2`
/// minors as exact four-term expansions, assembles the four triple products
/// (`abc`, `bcd`, `cda`, `dab`), multiplies each by the exact spatial lift of
/// the opposite vertex with the correct sign, and sums them into the
/// determinant expansion `(dlift*abc - clift*dab) + (blift*cda - alift*bcd)`; the
/// estimate is tested against [`ISPERRBOUND_B`] `* permanent`. Stage C adds the
/// first-order tail correction (reusing the minors' leading components) and
/// tests the result bound. It defers to [`insphere_exact`] only when both stages
/// remain uncertain.
pub(super) fn insphere_adapt(points: [[f64; 3]; 5], permanent: f64) -> i32 {
    let pa = points[0];
    let pb = points[1];
    let pc = points[2];
    let pd = points[3];
    let pe = points[4];

    let (aex, aextail) = two_diff_local(pa[0], pe[0]);
    let (aey, aeytail) = two_diff_local(pa[1], pe[1]);
    let (aez, aeztail) = two_diff_local(pa[2], pe[2]);
    let (bex, bextail) = two_diff_local(pb[0], pe[0]);
    let (bey, beytail) = two_diff_local(pb[1], pe[1]);
    let (bez, beztail) = two_diff_local(pb[2], pe[2]);
    let (cex, cextail) = two_diff_local(pc[0], pe[0]);
    let (cey, ceytail) = two_diff_local(pc[1], pe[1]);
    let (cez, ceztail) = two_diff_local(pc[2], pe[2]);
    let (dex, dextail) = two_diff_local(pd[0], pe[0]);
    let (dey, deytail) = two_diff_local(pd[1], pe[1]);
    let (dez, deztail) = two_diff_local(pd[2], pe[2]);

    // Stage B: six minors of the rounded differences.
    let ab = cross_diff_expansion(aex, bey, bex, aey); // aex*bey - bex*aey
    let bc = cross_diff_expansion(bex, cey, cex, bey); // bex*cey - cex*bey
    let cd = cross_diff_expansion(cex, dey, dex, cey); // cex*dey - dex*cey
    let da = cross_diff_expansion(dex, aey, aex, dey); // dex*aey - aex*dey
    let ac = cross_diff_expansion(aex, cey, cex, aey); // aex*cey - cex*aey
    let bd = cross_diff_expansion(bex, dey, dex, bey); // bex*dey - dex*bey

    // Triple products with the same sign pattern as the fast filter.
    let bcd = triple_sum(&cd, bez, &bd, -cez, &bc, dez); // bez*cd - cez*bd + dez*bc
    let cda = triple_sum(&da, cez, &ac, dez, &cd, aez); // cez*da + dez*ac + aez*cd
    let dab = triple_sum(&ab, dez, &bd, aez, &da, bez); // dez*ab + aez*bd + bez*da
    let abc = triple_sum(&bc, aez, &ac, -bez, &ab, cez); // aez*bc - bez*ac + cez*ab

    // Lift-weighted, signed so the sum is the in-sphere determinant.
    let adet = lift_scale_3d(&bcd, aex, aey, aez, -1.0); // -alift * bcd
    let bdet = lift_scale_3d(&cda, bex, bey, bez, 1.0); //  blift * cda
    let cdet = lift_scale_3d(&dab, cex, cey, cez, -1.0); // -clift * dab
    let ddet = lift_scale_3d(&abc, dex, dey, dez, 1.0); //  dlift * abc
    let abdet = fast_expansion_sum_zeroelim(&adet, &bdet);
    let cddet = fast_expansion_sum_zeroelim(&cdet, &ddet);
    let fin1 = fast_expansion_sum_zeroelim(&abdet, &cddet);

    let mut det = estimate(&fin1);
    if det.abs() >= ISPERRBOUND_B * permanent {
        return sign_of(det);
    }

    if aextail == 0.0
        && aeytail == 0.0
        && aeztail == 0.0
        && bextail == 0.0
        && beytail == 0.0
        && beztail == 0.0
        && cextail == 0.0
        && ceytail == 0.0
        && ceztail == 0.0
        && dextail == 0.0
        && deytail == 0.0
        && deztail == 0.0
    {
        return expansion_sign(&fin1);
    }

    // Leading (largest) components of the minors, reused by the correction.
    let ab3 = ab[3];
    let bc3 = bc[3];
    let cd3 = cd[3];
    let da3 = da[3];
    let ac3 = ac[3];
    let bd3 = bd[3];

    let errbound = ISPERRBOUND_C * permanent + RESULTERRBOUND * det.abs();
    let abeps = (aex * beytail + bey * aextail) - (aey * bextail + bex * aeytail);
    let bceps = (bex * ceytail + cey * bextail) - (bey * cextail + cex * beytail);
    let cdeps = (cex * deytail + dey * cextail) - (cey * dextail + dex * ceytail);
    let daeps = (dex * aeytail + aey * dextail) - (dey * aextail + aex * deytail);
    let aceps = (aex * ceytail + cey * aextail) - (aey * cextail + cex * aeytail);
    let bdeps = (bex * deytail + dey * bextail) - (bey * dextail + dex * beytail);
    det += (((bex * bex + bey * bey + bez * bez)
        * ((cez * daeps + dez * aceps + aez * cdeps)
            + (ceztail * da3 + deztail * ac3 + aeztail * cd3))
        + (dex * dex + dey * dey + dez * dez)
            * ((aez * bceps - bez * aceps + cez * abeps)
                + (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
        - ((aex * aex + aey * aey + aez * aez)
            * ((bez * cdeps - cez * bdeps + dez * bceps)
                + (beztail * cd3 - ceztail * bd3 + deztail * bc3))
            + (cex * cex + cey * cey + cez * cez)
                * ((dez * abeps + aez * bdeps + bez * daeps)
                    + (deztail * ab3 + aeztail * bd3 + beztail * da3))))
        + 2.0
            * (((bex * bextail + bey * beytail + bez * beztail)
                * (cez * da3 + dez * ac3 + aez * cd3)
                + (dex * dextail + dey * deytail + dez * deztail)
                    * (aez * bc3 - bez * ac3 + cez * ab3))
                - ((aex * aextail + aey * aeytail + aez * aeztail)
                    * (bez * cd3 - cez * bd3 + dez * bc3)
                    + (cex * cextail + cey * ceytail + cez * ceztail)
                        * (dez * ab3 + aez * bd3 + bez * da3)));
    if det.abs() >= errbound {
        return sign_of(det);
    }

    insphere_exact(points)
}
