//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ContactResult, EnhBox, EnhSimplex, EnhSphere, EpaEnhancedResult, GjkMetrics, GjkResult,
    SupportCache, SupportPoint, ToiResult, TranslatedShape, WarmStartCache,
};
use crate::gjk_epa::{SupportPoint as EpaPoint, complete_tetrahedron, expand_polytope};

/// Add two 3-vectors.
#[inline]
pub fn vadd(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors: a − b.
#[inline]
pub fn vsub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
#[inline]
pub fn vscale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product.
#[inline]
pub fn vdot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product.
#[inline]
pub fn vcross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Squared length of a 3-vector.
#[inline]
pub fn vlen_sq(a: [f64; 3]) -> f64 {
    vdot(a, a)
}
/// Length of a 3-vector.
#[inline]
pub fn vlen(a: [f64; 3]) -> f64 {
    vlen_sq(a).sqrt()
}
/// Normalise a 3-vector; returns `[0,0,0]` if near-zero.
#[inline]
pub fn vnorm(a: [f64; 3]) -> [f64; 3] {
    let l = vlen(a);
    if l < 1e-300 {
        return [0.0; 3];
    }
    vscale(a, 1.0 / l)
}
/// Negate a 3-vector.
#[inline]
pub fn vneg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
/// Linear interpolation between two 3-vectors.
#[inline]
pub fn vlerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    vadd(vscale(a, 1.0 - t), vscale(b, t))
}
/// Trait for convex shapes that can return a support point in a given direction.
pub trait ConvexShape: Send + Sync {
    /// Return the point on the shape that is furthest in direction `dir`.
    fn support(&self, dir: [f64; 3]) -> [f64; 3];
    /// Return the centre (approximate) of the shape.
    fn centre(&self) -> [f64; 3];
}
/// Johnson's distance: closest point on a line segment to the origin.
///
/// Returns the barycentric coordinates `(λ_A, λ_B)` and the closest point `w`.
/// The segment is `A + λ*(B-A)` with `λ ∈ [0,1]`.
pub fn johnson_segment(a: [f64; 3], b: [f64; 3]) -> (f64, f64, [f64; 3]) {
    let ab = vsub(b, a);
    let len_sq = vlen_sq(ab);
    if len_sq < 1e-24 {
        return (1.0, 0.0, a);
    }
    let t = (-vdot(a, ab) / len_sq).clamp(0.0, 1.0);
    let pt = vadd(a, vscale(ab, t));
    (1.0 - t, t, pt)
}
/// Johnson's distance: closest point on a triangle to the origin.
///
/// Returns barycentric coordinates `(λ_A, λ_B, λ_C)` and the closest point.
pub fn johnson_triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> (f64, f64, f64, [f64; 3]) {
    let ab = vsub(b, a);
    let ac = vsub(c, a);
    let neg_a = vneg(a);
    let d1 = vdot(ab, neg_a);
    let d2 = vdot(ac, neg_a);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (1.0, 0.0, 0.0, a);
    }
    let neg_b = vneg(b);
    let d3 = vdot(ab, neg_b);
    let d4 = vdot(ac, neg_b);
    if d3 >= 0.0 && d4 <= d3 {
        return (0.0, 1.0, 0.0, b);
    }
    let neg_c = vneg(c);
    let d5 = vdot(ab, neg_c);
    let d6 = vdot(ac, neg_c);
    if d6 >= 0.0 && d5 <= d6 {
        return (0.0, 0.0, 1.0, c);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let pt = vadd(a, vscale(ab, v));
        return (1.0 - v, v, 0.0, pt);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let pt = vadd(a, vscale(ac, w));
        return (1.0 - w, 0.0, w, pt);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let pt = vadd(b, vscale(vsub(c, b), w));
        return (0.0, 1.0 - w, w, pt);
    }
    let denom = 1.0 / (va + vb + vc);
    let u = va * denom;
    let v = vb * denom;
    let w = 1.0 - u - v;
    let pt = vadd(vadd(vscale(a, u), vscale(b, v)), vscale(c, w));
    (u, v, w, pt)
}
/// Johnson's distance: closest point on a tetrahedron to the origin.
///
/// Returns the closest point and whether the origin is inside the tetrahedron.
pub fn johnson_tetrahedron(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> ([f64; 3], bool) {
    let inside = |p: [f64; 3], q: [f64; 3], r: [f64; 3], apex: [f64; 3]| -> bool {
        let n = vcross(vsub(q, p), vsub(r, p));
        let to_apex = vdot(n, vsub(apex, p));
        let to_origin = vdot(n, vneg(p));
        to_apex * to_origin >= 0.0
    };
    // A flat tetrahedron has no interior: every apex lies in the plane of the
    // opposite face, so the sign tests would pass for any origin.
    let edge = vlen(vsub(b, a)).max(vlen(vsub(c, a))).max(vlen(vsub(d, a)));
    let volume = vdot(vcross(vsub(b, a), vsub(c, a)), vsub(d, a)).abs();
    let flat = volume.is_nan() || volume <= 1e-12 * edge * edge * edge;
    if !flat && inside(a, b, c, d) && inside(a, b, d, c) && inside(a, c, d, b) && inside(b, c, d, a)
    {
        return ([0.0; 3], true);
    }
    let faces = [(a, b, c), (a, b, d), (a, c, d), (b, c, d)];
    let mut best_dist = f64::INFINITY;
    let mut best_pt = a;
    for &(p, q, r) in &faces {
        let (_, _, _, pt) = johnson_triangle(p, q, r);
        let d = vlen_sq(pt);
        if d < best_dist {
            best_dist = d;
            best_pt = pt;
        }
    }
    (best_pt, false)
}
/// Relative accuracy at which the GJK iteration stops: it has converged once the
/// newest support point cannot bring the simplex closer to the origin by more
/// than `GJK_CONVERGENCE_TOLERANCE * (1 + distance)`.
const GJK_CONVERGENCE_TOLERANCE: f64 = 1e-9;

/// Relative decrease of the squared simplex distance below which an iteration
/// counts as having made no progress (rounding noise or a cycling simplex).
///
/// It must be far below the improvement a still-unconverged iteration makes: a
/// step that is `d` short of the exact distance improves the squared distance
/// only by about `(d / edge)^2`, so a coarse threshold such as `1e-8` stops
/// the iteration while the distance is still off by several `1e-6`.
const GJK_STAGNATION_TOLERANCE: f64 = 1e-12;

/// A [`GjkResult`] together with a separation certificate.
///
/// The certificate is a unit `axis` pointing from shape A toward shape B and a
/// `gap`: every point of B lies at least `gap` further along `axis` than every
/// point of A, so the shapes are at least `gap` apart and stay separated when A
/// is moved by any displacement `d` with `d . axis < gap`. It is a *lower*
/// bound of the distance (the closest points of the result give the upper
/// bound) and is what a conservative advancement has to step on.
pub(super) struct GjkCore {
    /// The public result.
    pub(super) result: GjkResult,
    /// Unit separating axis pointing from A toward B (zero if `gap == 0`).
    pub(super) axis: [f64; 3],
    /// Certified lower bound of the distance, `0.0` when no separating axis was
    /// found (for example because the shapes intersect).
    pub(super) gap: f64,
}

/// Enhanced GJK algorithm with warm starting.
///
/// Accepts an optional `WarmStartCache` from the previous frame. A valid cache
/// seeds the first search direction with its cached separating axis; the
/// search always starts from an empty simplex, because the cached vertices
/// belong to the previous configuration of the shapes. The cache itself is not
/// modified: call [`WarmStartCache::update`] with the returned simplex.
///
/// For separated shapes the iteration runs until the distance has converged,
/// so [`GjkResult::distance`] is the true distance (to within roughly
/// `1e-9 * (1 + distance)`), not merely an upper bound. If `max_iter`
/// iterations are not enough, the best simplex found so far is reported and its
/// distance is an upper bound of the true distance (for shapes that overlap it
/// can therefore be positive). At least one iteration is always performed, so
/// the result never describes an empty simplex.
pub fn gjk_enhanced(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: Option<&WarmStartCache>,
    max_iter: u32,
) -> GjkResult {
    gjk_enhanced_core(shape_a, shape_b, cache, max_iter).result
}
/// GJK distance iteration that also records the best separation certificate.
///
/// See [`gjk_enhanced`]; the extra output is described on [`GjkCore`].
pub(super) fn gjk_enhanced_core(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: Option<&WarmStartCache>,
    max_iter: u32,
) -> GjkCore {
    let mut metrics = GjkMetrics::new();
    // A warm start seeds the search direction only. The cached vertices belong to
    // the previous configuration of the shapes: starting from them would measure
    // a stale Minkowski difference, and a cached intersection (four vertices)
    // plus the first new support point is no simplex at all.
    let warm_direction = cache
        .filter(|c| c.valid)
        .map(|c| c.direction)
        .filter(|d| d.iter().all(|x| x.is_finite()) && vlen_sq(*d) > 1e-18);
    metrics.warm_started = warm_direction.is_some();
    let dir = warm_direction.unwrap_or_else(|| vsub(shape_b.centre(), shape_a.centre()));
    let mut support = |counter: &mut GjkMetrics, d: [f64; 3]| {
        counter.support_calls += 2;
        let sa = shape_a.support(d);
        let sb = shape_b.support(vneg(d));
        SupportPoint::new(vsub(sa, sb), sa, sb)
    };
    gjk_iterate(dir, max_iter, metrics, &mut support)
}
/// The GJK distance iteration over a support oracle.
///
/// `dir` is the first search direction (a zero or non-finite vector falls back
/// to `+x`); `support(metrics, d)` returns the support point of the Minkowski
/// difference `A - B` along the unit direction `d` and counts its own support
/// calls in `metrics`. See [`gjk_enhanced`] for the termination rules.
fn gjk_iterate(
    mut dir: [f64; 3],
    max_iter: u32,
    mut metrics: GjkMetrics,
    support: &mut dyn FnMut(&mut GjkMetrics, [f64; 3]) -> SupportPoint,
) -> GjkCore {
    let mut simplex = EnhSimplex::new();
    if vlen_sq(dir) < 1e-18 {
        dir = [1.0, 0.0, 0.0];
    }
    dir = vnorm(dir);
    if vlen_sq(dir) < 1e-18 {
        dir = [1.0, 0.0, 0.0];
    }
    let mut prev_dist_sq = f64::INFINITY;
    let mut axis = [0.0; 3];
    let mut gap = 0.0_f64;
    for _iter in 0..max_iter.max(1) {
        metrics.iterations += 1;
        let point = support(&mut metrics, dir);
        let w = point.w;
        // `dir` is a unit vector: `reach` is how far the Minkowski difference
        // extends along it, measured from the origin. A negative value proves
        // that A and B are separated along `dir` by at least `-reach`.
        let reach = vdot(w, dir);
        if reach < 0.0 && -reach > gap {
            gap = -reach;
            axis = dir;
        }
        if prev_dist_sq.is_finite() {
            // The simplex is `closest` away from the origin and `dir` points
            // from its closest point to the origin, so `reach >= -closest`.
            // When equality holds up to the tolerance no support point can
            // bring the simplex closer: the distance has converged.
            let closest = prev_dist_sq.sqrt();
            if reach + closest <= GJK_CONVERGENCE_TOLERANCE * (1.0 + closest) {
                metrics.converged_early = true;
                return GjkCore {
                    result: build_separated_result(simplex, metrics),
                    axis,
                    gap,
                };
            }
        }
        simplex.add(point);
        let (new_dir, contains_origin) = simplex.reduce();
        metrics.simplex_reductions += 1;
        if contains_origin {
            metrics.final_dist_sq = 0.0;
            return GjkCore {
                result: GjkResult {
                    intersecting: true,
                    closest_a: [0.0; 3],
                    closest_b: [0.0; 3],
                    dist_sq: 0.0,
                    simplex: simplex.verts,
                    metrics,
                },
                axis: [0.0; 3],
                gap: 0.0,
            };
        }
        let new_dist_sq = vlen_sq(new_dir);
        if new_dist_sq < 1e-18 {
            metrics.converged_early = true;
            return GjkCore {
                result: build_separated_result(simplex, metrics),
                axis,
                gap,
            };
        }
        if new_dist_sq >= prev_dist_sq * (1.0 - GJK_STAGNATION_TOLERANCE)
            && !simplex.verts.is_empty()
        {
            metrics.converged_early = true;
            return GjkCore {
                result: build_separated_result(simplex, metrics),
                axis,
                gap,
            };
        }
        prev_dist_sq = new_dist_sq;
        dir = vnorm(new_dir);
        if vlen_sq(dir) < 1e-18 {
            metrics.converged_early = true;
            return GjkCore {
                result: build_separated_result(simplex, metrics),
                axis,
                gap,
            };
        }
    }
    GjkCore {
        result: build_separated_result(simplex, metrics),
        axis,
        gap,
    }
}
pub(super) fn build_separated_result(simplex: EnhSimplex, mut metrics: GjkMetrics) -> GjkResult {
    let (closest_a, closest_b) = extract_closest_points(&simplex);
    let diff = vsub(closest_a, closest_b);
    let dist_sq = vlen_sq(diff);
    metrics.final_dist_sq = dist_sq;
    GjkResult {
        intersecting: false,
        closest_a,
        closest_b,
        dist_sq,
        simplex: simplex.verts,
        metrics,
    }
}
pub(super) fn extract_closest_points(simplex: &EnhSimplex) -> ([f64; 3], [f64; 3]) {
    match simplex.verts.len() {
        0 => ([0.0; 3], [0.0; 3]),
        1 => (simplex.verts[0].a, simplex.verts[0].b),
        2 => {
            let pa = simplex.verts[1].w;
            let pb = simplex.verts[0].w;
            let (la, lb, _) = johnson_segment(pa, pb);
            let ca = vadd(
                vscale(simplex.verts[1].a, la),
                vscale(simplex.verts[0].a, lb),
            );
            let cb = vadd(
                vscale(simplex.verts[1].b, la),
                vscale(simplex.verts[0].b, lb),
            );
            (ca, cb)
        }
        3 => {
            let wa = simplex.verts[2].w;
            let wb = simplex.verts[1].w;
            let wc = simplex.verts[0].w;
            let (la, lb, lc, _) = johnson_triangle(wa, wb, wc);
            let ca = vadd(
                vadd(
                    vscale(simplex.verts[2].a, la),
                    vscale(simplex.verts[1].a, lb),
                ),
                vscale(simplex.verts[0].a, lc),
            );
            let cb = vadd(
                vadd(
                    vscale(simplex.verts[2].b, la),
                    vscale(simplex.verts[1].b, lb),
                ),
                vscale(simplex.verts[0].b, lc),
            );
            (ca, cb)
        }
        _ => (simplex.verts[0].a, simplex.verts[0].b),
    }
}
/// Compute the distance between two convex shapes using GJK.
///
/// Returns `0.0` if the shapes are intersecting.
pub fn gjk_distance_enhanced(shape_a: &dyn ConvexShape, shape_b: &dyn ConvexShape) -> f64 {
    let result = gjk_enhanced(shape_a, shape_b, None, 64);
    result.distance()
}
/// Compute the closest points between two separated convex shapes.
///
/// Returns `None` if the shapes are intersecting (distance = 0).
pub fn gjk_closest_points(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
) -> Option<([f64; 3], [f64; 3])> {
    let result = gjk_enhanced(shape_a, shape_b, None, 64);
    if result.intersecting {
        None
    } else {
        Some((result.closest_a, result.closest_b))
    }
}
/// Run the enhanced EPA after GJK confirmed intersection.
///
/// `initial_simplex` is the simplex GJK ended with. It may hold any number of
/// points: the distinct, affinely independent ones are kept and completed to a
/// tetrahedron of support points of the Minkowski difference `A - B`. The
/// tetrahedron does not have to enclose the origin; the polytope first grows
/// toward it.
///
/// Each step takes the face nearest to the origin, asks for the support point
/// along its normal, removes **every** face that can see that point and joins
/// the point to the horizon (the edges of exactly one removed face), so the
/// polytope stays convex. The iteration stops when the support point lies
/// within `tolerance` beyond the nearest face, or after `max_iter` steps (the
/// depth is then a lower bound).
///
/// The result is that of the nearest facet: `normal` points from B into A (the
/// direction in which A has to move to separate), `depth` is the distance of
/// the facet from the origin and `contact_point` is the midpoint of the two
/// deepest points, the barycentric blends of the support points of the facet
/// that reproduce the foot of the perpendicular from the origin. The deepest
/// point of A is `contact_point - normal * depth / 2` and that of B
/// `contact_point + normal * depth / 2`.
///
/// If no penetration can be derived (the Minkowski difference is flat, or the
/// shapes do not overlap) the depth is zero, `normal` points along the line of
/// centres from B to A (or +y when they coincide), `contact_point` is the
/// midpoint of the centres and `iterations` is 0.
pub fn epa_enhanced(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    initial_simplex: Vec<SupportPoint>,
    max_iter: u32,
    tolerance: f64,
) -> EpaEnhancedResult {
    let oracle = |direction: [f64; 3]| {
        let a = shape_a.support(direction);
        let b = shape_b.support(vneg(direction));
        EpaPoint {
            w: vsub(a, b),
            a,
            b,
        }
    };
    let known: Vec<EpaPoint> = initial_simplex
        .iter()
        .map(|p| EpaPoint {
            w: p.w,
            a: p.a,
            b: p.b,
        })
        .collect();
    let penetration = <[EpaPoint; 4]>::try_from(complete_tetrahedron(&known, &oracle).as_slice())
        .ok()
        .and_then(|seed| {
            expand_polytope(
                seed,
                &oracle,
                tolerance,
                usize::try_from(max_iter).unwrap_or(usize::MAX),
            )
        });
    match penetration {
        Some(found) => EpaEnhancedResult {
            depth: found.depth,
            normal: vneg(found.normal),
            contact_point: vscale(vadd(found.contact_a, found.contact_b), 0.5),
            iterations: u32::try_from(found.iterations).unwrap_or(u32::MAX),
        },
        None => {
            let toward_a = vnorm(vsub(shape_a.centre(), shape_b.centre()));
            EpaEnhancedResult {
                depth: 0.0,
                normal: if vlen_sq(toward_a) > 0.0 {
                    toward_a
                } else {
                    [0.0, 1.0, 0.0]
                },
                contact_point: vscale(vadd(shape_a.centre(), shape_b.centre()), 0.5),
                iterations: 0,
            }
        }
    }
}
/// Run the combined GJK-EPA pipeline.
///
/// If GJK reports intersection, EPA is run to compute penetration depth.
/// Otherwise the GJK closest-point result is returned.
pub fn gjk_epa_pipeline(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: Option<&WarmStartCache>,
    max_gjk_iter: u32,
    max_epa_iter: u32,
    epa_tol: f64,
) -> ContactResult {
    let gjk_res = gjk_enhanced(shape_a, shape_b, cache, max_gjk_iter);
    if !gjk_res.intersecting {
        let sep = gjk_res.distance();
        let diff = vsub(gjk_res.closest_a, gjk_res.closest_b);
        let normal = vnorm(diff);
        return ContactResult {
            in_contact: false,
            signed_distance: sep,
            point_a: gjk_res.closest_a,
            point_b: gjk_res.closest_b,
            normal,
            gjk_metrics: gjk_res.metrics,
            epa_iterations: 0,
        };
    }
    let epa_res = epa_enhanced(
        shape_a,
        shape_b,
        gjk_res.simplex.clone(),
        max_epa_iter,
        epa_tol,
    );
    ContactResult {
        in_contact: true,
        signed_distance: -epa_res.depth,
        // The deepest point of A lies half a depth against the B-to-A normal
        // from the midpoint, that of B half a depth along it.
        point_a: vsub(
            epa_res.contact_point,
            vscale(epa_res.normal, epa_res.depth * 0.5),
        ),
        point_b: vadd(
            epa_res.contact_point,
            vscale(epa_res.normal, epa_res.depth * 0.5),
        ),
        normal: epa_res.normal,
        gjk_metrics: gjk_res.metrics,
        epa_iterations: epa_res.iterations,
    }
}
/// Number of GJK iterations allowed for each distance query of the time-of-impact
/// search.
const TOI_GJK_ITERATIONS: u32 = 64;

/// Closing speeds (along the separating axis) at or below this value never close
/// the gap: the shapes are treated as separating.
const MIN_CLOSING_SPEED: f64 = 1e-12;

/// GJK-based time-of-impact query.
///
/// Sweeps `shape_a` along the relative velocity `vel_rel` over the time
/// interval \[0, 1\] while `shape_b` stays in place, and determines when the two
/// first touch.
///
/// Uses conservative advancement: at time `t` the shapes are separated by a
/// certified gap along a separating axis `n` (A toward B), so they cannot touch
/// before `t + gap / (vel_rel . n)`; time is advanced by exactly that amount
/// until the shapes are within `tolerance` of each other (a hit at that time)
/// or the advance passes `t = 1` (no hit). If the separating axis is not
/// approached (`vel_rel . n <= 0`) the shapes are known to stay apart and the
/// sweep misses. Unlike a test of the end positions this cannot tunnel: a shape
/// that starts and ends the interval separated from the other, but passes
/// through it in between, is reported as a hit.
///
/// Shapes that already intersect, or are within `tolerance`, hit at `toi = 0`.
/// The reported `toi` is the first sampled time at which the gap is within
/// `tolerance`; it never exceeds the true time of impact. `normal` points from
/// B toward A. If the iteration budget `max_iter` is exhausted before the sweep
/// is decided, a hit at the last conservative time is reported (the true impact
/// cannot be earlier).
///
/// `vel_rel` is the displacement of A relative to B over the interval, so
/// `vel_a - vel_b` for two moving shapes.
pub fn gjk_time_of_impact(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    vel_rel: [f64; 3],
    max_iter: u32,
    tolerance: f64,
) -> ToiResult {
    let touching_dist_sq = tolerance.max(0.0) * tolerance.max(0.0);
    let centre_normal = |moved: &dyn ConvexShape| vnorm(vsub(moved.centre(), shape_b.centre()));
    if !vel_rel.iter().all(|c| c.is_finite()) {
        return ToiResult {
            hit: false,
            toi: 1.0,
            normal: centre_normal(shape_a),
            iterations: 0,
        };
    }

    let mut t = 0.0_f64;
    let mut normal = centre_normal(shape_a);
    let mut iterations = 0u32;
    for _ in 0..max_iter.max(1) {
        iterations += 1;
        let moved = TranslatedShape {
            shape: shape_a,
            offset: vscale(vel_rel, t),
        };
        let query = gjk_enhanced_core(&moved, shape_b, None, TOI_GJK_ITERATIONS);
        let result = &query.result;
        if result.intersecting {
            return ToiResult {
                hit: true,
                toi: t,
                normal: centre_normal(&moved),
                iterations,
            };
        }
        if vlen_sq(vsub(result.closest_a, result.closest_b)) > 1e-24 {
            normal = vnorm(vsub(result.closest_a, result.closest_b));
        } else if query.gap > 0.0 {
            normal = vneg(query.axis);
        }
        if result.dist_sq <= touching_dist_sq {
            return ToiResult {
                hit: true,
                toi: t,
                normal,
                iterations,
            };
        }
        if query.gap <= 0.0 {
            // No separation could be certified, so a contact before the next
            // sample cannot be ruled out: stay conservative.
            return ToiResult {
                hit: true,
                toi: t,
                normal,
                iterations,
            };
        }
        let closing_speed = vdot(vel_rel, query.axis);
        if closing_speed <= MIN_CLOSING_SPEED {
            return ToiResult {
                hit: false,
                toi: 1.0,
                normal: vneg(query.axis),
                iterations,
            };
        }
        t += query.gap / closing_speed;
        if t > 1.0 {
            return ToiResult {
                hit: false,
                toi: 1.0,
                normal: vneg(query.axis),
                iterations,
            };
        }
    }
    ToiResult {
        hit: true,
        toi: t,
        normal,
        iterations,
    }
}
/// Compute the closest point on a sphere to a given external point.
///
/// Returns the point on the sphere surface closest to `p`.
pub fn closest_point_on_sphere(sphere: &EnhSphere, p: [f64; 3]) -> [f64; 3] {
    let dir = vsub(p, sphere.centre);
    let l = vlen(dir);
    if l < 1e-12 {
        return vadd(sphere.centre, [sphere.radius, 0.0, 0.0]);
    }
    vadd(sphere.centre, vscale(dir, sphere.radius / l))
}
/// Compute the closest point on an AABB (box) to a given point.
pub fn closest_point_on_box(b: &EnhBox, p: [f64; 3]) -> [f64; 3] {
    [
        p[0].clamp(b.centre[0] - b.half[0], b.centre[0] + b.half[0]),
        p[1].clamp(b.centre[1] - b.half[1], b.centre[1] + b.half[1]),
        p[2].clamp(b.centre[2] - b.half[2], b.centre[2] + b.half[2]),
    ]
}
/// Compute the closest point on a line segment to a given point.
///
/// Returns the closest point and the parameter `t ∈ [0,1]`.
pub fn closest_point_on_segment(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> ([f64; 3], f64) {
    let ab = vsub(b, a);
    let ap = vsub(p, a);
    let len_sq = vlen_sq(ab);
    if len_sq < 1e-24 {
        return (a, 0.0);
    }
    let t = (vdot(ap, ab) / len_sq).clamp(0.0, 1.0);
    (vadd(a, vscale(ab, t)), t)
}
/// Compute the distance between two line segments.
///
/// Returns `(pa, pb, dist)` where `pa` is the closest point on segment 1,
/// `pb` on segment 2, and `dist` is their distance.
pub fn segment_segment_distance(
    a0: [f64; 3],
    a1: [f64; 3],
    b0: [f64; 3],
    b1: [f64; 3],
) -> ([f64; 3], [f64; 3], f64) {
    let d1 = vsub(a1, a0);
    let d2 = vsub(b1, b0);
    let r = vsub(a0, b0);
    let a = vdot(d1, d1);
    let e = vdot(d2, d2);
    let f = vdot(d2, r);
    let (s, t) = if a < 1e-12 && e < 1e-12 {
        (0.0, 0.0)
    } else if a < 1e-12 {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = vdot(d1, r);
        if e < 1e-12 {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b_dot = vdot(d1, d2);
            let denom = a * e - b_dot * b_dot;
            let s0 = if denom.abs() > 1e-12 {
                ((b_dot * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t0 = (b_dot * s0 + f) / e;
            if t0 < 0.0 {
                ((-c / a).clamp(0.0, 1.0), 0.0)
            } else if t0 > 1.0 {
                (((b_dot - c) / a).clamp(0.0, 1.0), 1.0)
            } else {
                (s0, t0)
            }
        }
    };
    let pa = vadd(a0, vscale(d1, s));
    let pb = vadd(b0, vscale(d2, t));
    let dist = vlen(vsub(pa, pb));
    (pa, pb, dist)
}
/// GJK variant that maintains a `SupportCache` to skip redundant support calls.
///
/// It runs the same distance iteration as [`gjk_enhanced`] (so the distance is
/// the true distance, not the first separating-axis estimate); a support query
/// whose direction matches the cached one within `1e-6` reuses the cached
/// support points instead of calling the shapes.
pub fn gjk_with_support_cache(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: &mut SupportCache,
    max_iter: u32,
) -> GjkResult {
    let dir = vsub(shape_b.centre(), shape_a.centre());
    let mut support = |metrics: &mut GjkMetrics, d: [f64; 3]| {
        let (sa, sb) = match cache.get(d, 1e-6) {
            Some(cached) => cached,
            None => {
                let sa = shape_a.support(d);
                let sb = shape_b.support(vneg(d));
                cache.put(d, sa, sb);
                metrics.support_calls += 2;
                (sa, sb)
            }
        };
        SupportPoint::new(vsub(sa, sb), sa, sb)
    };
    gjk_iterate(dir, max_iter, GjkMetrics::new(), &mut support).result
}
#[cfg(test)]
mod tests_vector_math {
    use super::*;

    #[test]
    fn test_vadd() {
        assert_eq!(vadd([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_vsub() {
        assert_eq!(vsub([5.0, 7.0, 9.0], [4.0, 5.0, 6.0]), [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_vdot() {
        assert!((vdot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])).abs() < 1e-12);
        assert!((vdot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vcross_orthogonal() {
        let c = vcross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vnorm_unit() {
        let v = vnorm([3.0, 4.0, 0.0]);
        assert!((vlen(v) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_vnorm_zero() {
        let v = vnorm([0.0, 0.0, 0.0]);
        assert_eq!(v, [0.0; 3]);
    }
    #[test]
    fn test_vlerp_endpoints() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let l0 = vlerp(a, b, 0.0);
        let l1 = vlerp(a, b, 1.0);
        assert!((l0[0] - 1.0).abs() < 1e-12);
        assert!((l1[1] - 1.0).abs() < 1e-12);
    }
}
#[cfg(test)]
mod tests_support_functions {
    use super::*;
    use crate::gjk_enhanced::ConvexHull;
    use crate::gjk_enhanced::EnhBox;
    use crate::gjk_enhanced::EnhCapsule;
    use crate::gjk_enhanced::EnhSphere;
    use crate::gjk_enhanced::TorusApprox;

    #[test]
    fn test_sphere_support_max_x() {
        let s = EnhSphere::new([0.0; 3], 2.0);
        let sup = s.support([1.0, 0.0, 0.0]);
        assert!((sup[0] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_box_support_positive() {
        let b = EnhBox::new([0.0; 3], [1.0, 2.0, 3.0]);
        let sup = b.support([1.0, 1.0, 1.0]);
        assert_eq!(sup, [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_box_support_negative() {
        let b = EnhBox::new([0.0; 3], [1.0, 2.0, 3.0]);
        let sup = b.support([-1.0, -1.0, -1.0]);
        assert_eq!(sup, [-1.0, -2.0, -3.0]);
    }
    #[test]
    fn test_capsule_support_along_axis() {
        let cap = EnhCapsule::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let sup = cap.support([0.0, 1.0, 0.0]);
        assert!((sup[1] - 1.5).abs() < 1e-9, "sup[1]={}", sup[1]);
    }
    #[test]
    fn test_convex_hull_support() {
        let hull = ConvexHull::new(vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]);
        let sup = hull.support([1.0, 0.0, 0.0]);
        assert!((sup[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_torus_support_in_plane() {
        let t = TorusApprox::new([0.0; 3], 2.0, 0.5);
        let sup = t.support([1.0, 0.0, 0.0]);
        assert!((sup[0] - 2.5).abs() < 1e-9, "got {}", sup[0]);
    }
}
#[cfg(test)]
mod tests_johnson {

    use crate::gjk_enhanced::johnson_segment;
    use crate::gjk_enhanced::johnson_tetrahedron;
    use crate::gjk_enhanced::johnson_triangle;
    use crate::gjk_enhanced::vlen_sq;
    #[test]
    fn test_johnson_segment_midpoint() {
        let (la, lb, pt) = johnson_segment([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((pt[0]).abs() < 1e-9, "pt={pt:?}");
        assert!((la - 0.5).abs() < 1e-9);
        assert!((lb - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_johnson_segment_endpoint_a() {
        let (la, _lb, pt) = johnson_segment([2.0, 0.0, 0.0], [4.0, 0.0, 0.0]);
        assert!((la - 1.0).abs() < 1e-9);
        assert!((pt[0] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_johnson_segment_endpoint_b() {
        let (_la, lb, pt) = johnson_segment([-4.0, 0.0, 0.0], [-2.0, 0.0, 0.0]);
        assert!((lb - 1.0).abs() < 1e-9);
        assert!((pt[0] - (-2.0)).abs() < 1e-9);
    }
    #[test]
    fn test_johnson_triangle_vertex_a() {
        let a = [1.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [1.5, 1.0, 0.0];
        let (la, lb, lc, pt) = johnson_triangle(a, b, c);
        assert!(
            (la - 1.0).abs() < 0.1,
            "la={la}, lb={lb}, lc={lc}, pt={pt:?}"
        );
    }
    #[test]
    fn test_johnson_triangle_interior() {
        let a = [1.0, 0.0, 0.0];
        let b = [-1.0, 1.0, 0.0];
        let c = [-1.0, -1.0, 0.0];
        let (_la, _lb, _lc, pt) = johnson_triangle(a, b, c);
        assert!(vlen_sq(pt) < 0.2, "pt not near origin: {pt:?}");
    }
    #[test]
    fn test_johnson_tetrahedron_contains_origin() {
        let a = [1.0, 0.0, 0.0];
        let b = [-1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        let (_pt, inside) = johnson_tetrahedron(a, b, c, d);
        assert!(inside, "origin should be inside tetrahedron");
    }
    #[test]
    fn test_johnson_tetrahedron_outside() {
        let a = [2.0, 0.0, 0.0];
        let b = [3.0, 0.0, 0.0];
        let c = [2.5, 1.0, 0.0];
        let d = [2.5, 0.5, 1.0];
        let (_pt, inside) = johnson_tetrahedron(a, b, c, d);
        assert!(!inside);
    }
}
#[cfg(test)]
mod tests_gjk_enhanced {

    use crate::gjk_enhanced::EnhBox;

    use crate::gjk_enhanced::EnhSphere;

    use crate::gjk_enhanced::WarmStartCache;
    use crate::gjk_enhanced::gjk_enhanced;

    #[test]
    fn test_gjk_separated_spheres() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 1.0);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
        let res = gjk_enhanced(&s1, &s2, None, 64);
        assert!(!res.intersecting, "separated spheres should not intersect");
        assert!(
            (res.distance() - 3.0).abs() < 0.01,
            "dist={}",
            res.distance()
        );
    }
    #[test]
    fn test_gjk_overlapping_spheres() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 2.0);
        let s2 = EnhSphere::new([1.0, 0.0, 0.0], 2.0);
        let res = gjk_enhanced(&s1, &s2, None, 64);
        assert!(res.intersecting, "overlapping spheres should intersect");
    }
    #[test]
    fn test_gjk_touching_spheres() {
        let s1 = EnhSphere::new([0.0, 0.0, 0.0], 1.0);
        let s2 = EnhSphere::new([2.0, 0.0, 0.0], 1.0);
        let res = gjk_enhanced(&s1, &s2, None, 64);
        assert!(
            res.dist_sq < 0.01 || res.intersecting,
            "touching spheres dist={}",
            res.distance()
        );
    }
    #[test]
    fn test_gjk_box_vs_box_separated() {
        let b1 = EnhBox::new([0.0; 3], [1.0; 3]);
        let b2 = EnhBox::new([5.0, 0.0, 0.0], [1.0; 3]);
        let res = gjk_enhanced(&b1, &b2, None, 64);
        assert!(!res.intersecting);
        assert!(
            (res.distance() - 3.0).abs() < 0.01,
            "dist={}",
            res.distance()
        );
    }
    #[test]
    fn test_gjk_box_vs_box_overlapping() {
        let b1 = EnhBox::new([0.0; 3], [1.0; 3]);
        let b2 = EnhBox::new([1.0, 0.0, 0.0], [1.0; 3]);
        let res = gjk_enhanced(&b1, &b2, None, 64);
        assert!(res.intersecting);
    }
    #[test]
    fn test_gjk_metrics_iteration_count() {
        let s1 = EnhSphere::new([0.0; 3], 1.0);
        let s2 = EnhSphere::new([10.0, 0.0, 0.0], 1.0);
        let res = gjk_enhanced(&s1, &s2, None, 64);
        assert!(res.metrics.iterations >= 1);
    }
    #[test]
    fn test_gjk_warm_start_reduces_iterations() {
        let s1 = EnhSphere::new([0.0; 3], 1.0);
        let s2 = EnhSphere::new([5.0, 0.0, 0.0], 1.0);
        let res1 = gjk_enhanced(&s1, &s2, None, 64);
        let mut cache = WarmStartCache::new();
        cache.update(&res1.simplex, [1.0, 0.0, 0.0]);
        let s2b = EnhSphere::new([5.1, 0.0, 0.0], 1.0);
        let res2 = gjk_enhanced(&s1, &s2b, Some(&cache), 64);
        assert!(res2.metrics.warm_started);
        // Radii 1 + 1 against a centre distance of 5.1.
        assert!(!res2.intersecting);
        assert!(
            (res2.distance() - 3.1).abs() < 1e-9,
            "dist={}",
            res2.distance()
        );
    }
}
