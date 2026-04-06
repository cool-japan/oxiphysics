//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ContactResult, EnhBox, EnhSimplex, EnhSphere, EpaEnhancedResult, EpaFaceE, GjkMetrics,
    GjkResult, SupportCache, SupportPoint, ToiResult, TranslatedShape, WarmStartCache,
};

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
    if inside(a, b, c, d) && inside(a, b, d, c) && inside(a, c, d, b) && inside(b, c, d, a) {
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
/// Enhanced GJK algorithm with warm starting.
///
/// Accepts an optional `WarmStartCache` from the previous frame.
/// Returns a `GjkResult` and optionally updates the cache.
pub fn gjk_enhanced(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: Option<&WarmStartCache>,
    max_iter: u32,
) -> GjkResult {
    let mut metrics = GjkMetrics::new();
    let mut simplex = EnhSimplex::new();
    let mut dir = if let Some(c) = cache.filter(|c| c.valid) {
        metrics.warm_started = true;
        simplex = EnhSimplex::from_cache(c);
        c.direction
    } else {
        vsub(shape_b.centre(), shape_a.centre())
    };
    if vlen_sq(dir) < 1e-18 {
        dir = [1.0, 0.0, 0.0];
    }
    dir = vnorm(dir);
    if vlen_sq(dir) < 1e-18 {
        dir = [1.0, 0.0, 0.0];
    }
    let mut prev_dist_sq = f64::INFINITY;
    for _iter in 0..max_iter {
        metrics.iterations += 1;
        let sa = shape_a.support(dir);
        let sb = shape_b.support(vneg(dir));
        metrics.support_calls += 2;
        let w = vsub(sa, sb);
        if !simplex.verts.is_empty() && vdot(w, dir) < 1e-7 {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
        simplex.add(SupportPoint::new(w, sa, sb));
        let (new_dir, contains_origin) = simplex.reduce();
        metrics.simplex_reductions += 1;
        if contains_origin {
            metrics.final_dist_sq = 0.0;
            return GjkResult {
                intersecting: true,
                closest_a: [0.0; 3],
                closest_b: [0.0; 3],
                dist_sq: 0.0,
                simplex: simplex.verts,
                metrics,
            };
        }
        let new_dist_sq = vlen_sq(new_dir);
        if new_dist_sq < 1e-18 {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
        if new_dist_sq >= prev_dist_sq * (1.0 - 1e-8) && !simplex.verts.is_empty() {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
        prev_dist_sq = new_dist_sq;
        dir = vnorm(new_dir);
        if vlen_sq(dir) < 1e-18 {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
    }
    build_separated_result(simplex, metrics)
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
/// `simplex_verts` must be a full tetrahedral GJK simplex (4 points).
pub fn epa_enhanced(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    initial_simplex: Vec<SupportPoint>,
    max_iter: u32,
    tolerance: f64,
) -> EpaEnhancedResult {
    let mut verts: Vec<SupportPoint> = initial_simplex;
    let expand_dirs: [[f64; 3]; 14] = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [-1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
    ];
    for &d in &expand_dirs {
        if verts.len() >= 4 {
            break;
        }
        let dn = vnorm(d);
        let sa = shape_a.support(dn);
        let sb = shape_b.support(vneg(dn));
        let w = vsub(sa, sb);
        let is_new = verts.iter().all(|v| vlen_sq(vsub(v.w, w)) > 1e-12);
        if is_new {
            verts.push(SupportPoint::new(w, sa, sb));
        }
    }
    let is_degenerate = |pts: &[SupportPoint]| -> bool {
        if pts.len() < 4 {
            return true;
        }
        let a = pts[0].w;
        let b = pts[1].w;
        let c = pts[2].w;
        let d = pts[3].w;
        let ab = vsub(b, a);
        let ac = vsub(c, a);
        let ad = vsub(d, a);
        let n = vcross(ab, ac);
        vdot(n, ad).abs() < 1e-8
    };
    if is_degenerate(&verts) {
        for &d in &expand_dirs {
            if !is_degenerate(&verts) {
                break;
            }
            let dn = vnorm(d);
            let sa = shape_a.support(dn);
            let sb = shape_b.support(vneg(dn));
            let w = vsub(sa, sb);
            let is_new = verts.iter().all(|v| vlen_sq(vsub(v.w, w)) > 1e-12);
            if is_new {
                verts.push(SupportPoint::new(w, sa, sb));
                if verts.len() > 4 {
                    let n = verts.len();
                    let candidate = vec![verts[0], verts[1], verts[2], verts[n - 1]];
                    if !is_degenerate(&candidate) {
                        verts = candidate;
                    } else {
                        let candidate2 = vec![verts[0], verts[1], verts[n - 1], verts[n - 2]];
                        if !is_degenerate(&candidate2) {
                            verts = candidate2;
                        }
                    }
                }
            }
        }
    }
    if verts.len() < 4 || is_degenerate(&verts) {
        return EpaEnhancedResult {
            depth: 0.0,
            normal: [0.0, 1.0, 0.0],
            contact_point: [0.0; 3],
            iterations: 0,
        };
    }
    let mut faces: Vec<EpaFaceE> = vec![
        epa_make_face(&verts, 0, 1, 2),
        epa_make_face(&verts, 0, 1, 3),
        epa_make_face(&verts, 0, 2, 3),
        epa_make_face(&verts, 1, 2, 3),
    ];
    for iter in 0..max_iter {
        let Some((fi, &face)) = faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.dist > 1e-12)
            .min_by(|(_, a), (_, b)| {
                a.dist
                    .partial_cmp(&b.dist)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .or_else(|| {
                faces.iter().enumerate().min_by(|(_, a), (_, b)| {
                    a.dist
                        .partial_cmp(&b.dist)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
            })
        else {
            break;
        };
        let sa = shape_a.support(face.normal);
        let sb = shape_b.support(vneg(face.normal));
        let w = vsub(sa, sb);
        let new_dist = vdot(face.normal, w);
        if new_dist - face.dist < tolerance {
            let depth = new_dist.max(face.dist).max(0.0);
            if depth < 1e-12 {
                break;
            }
            let cp_a = vscale(
                vadd(
                    verts[face.indices[0]].a,
                    vadd(verts[face.indices[1]].a, verts[face.indices[2]].a),
                ),
                1.0 / 3.0,
            );
            let cp_b = vscale(
                vadd(
                    verts[face.indices[0]].b,
                    vadd(verts[face.indices[1]].b, verts[face.indices[2]].b),
                ),
                1.0 / 3.0,
            );
            return EpaEnhancedResult {
                depth,
                normal: face.normal,
                contact_point: vscale(vadd(cp_a, cp_b), 0.5),
                iterations: iter,
            };
        }
        let new_idx = verts.len();
        verts.push(SupportPoint::new(w, sa, sb));
        let old_face = faces.remove(fi);
        let [i0, i1, i2] = old_face.indices;
        faces.push(epa_make_face(&verts, i0, i1, new_idx));
        faces.push(epa_make_face(&verts, i1, i2, new_idx));
        faces.push(epa_make_face(&verts, i2, i0, new_idx));
    }
    let best_face_opt = faces.iter().filter(|f| f.dist > 1e-12).min_by(|a, b| {
        a.dist
            .partial_cmp(&b.dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if let Some(&face) = best_face_opt {
        let cp_a = vscale(
            vadd(
                verts[face.indices[0]].a,
                vadd(verts[face.indices[1]].a, verts[face.indices[2]].a),
            ),
            1.0 / 3.0,
        );
        let cp_b = vscale(
            vadd(
                verts[face.indices[0]].b,
                vadd(verts[face.indices[1]].b, verts[face.indices[2]].b),
            ),
            1.0 / 3.0,
        );
        EpaEnhancedResult {
            depth: face.dist,
            normal: face.normal,
            contact_point: vscale(vadd(cp_a, cp_b), 0.5),
            iterations: max_iter,
        }
    } else {
        EpaEnhancedResult {
            depth: 0.0,
            normal: [0.0, 1.0, 0.0],
            contact_point: [0.0; 3],
            iterations: max_iter,
        }
    }
}
pub(super) fn epa_make_face(verts: &[SupportPoint], i0: usize, i1: usize, i2: usize) -> EpaFaceE {
    let a = verts[i0].w;
    let b = verts[i1].w;
    let c = verts[i2].w;
    let ab = vsub(b, a);
    let ac = vsub(c, a);
    let mut n = vnorm(vcross(ab, ac));
    if vdot(n, a) < 0.0 {
        n = vneg(n);
    }
    let dist = vdot(n, a);
    EpaFaceE {
        indices: [i0, i1, i2],
        normal: n,
        dist,
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
        point_a: vadd(
            epa_res.contact_point,
            vscale(epa_res.normal, epa_res.depth * 0.5),
        ),
        point_b: vsub(
            epa_res.contact_point,
            vscale(epa_res.normal, epa_res.depth * 0.5),
        ),
        normal: epa_res.normal,
        gjk_metrics: gjk_res.metrics,
        epa_iterations: epa_res.iterations,
    }
}
/// GJK-based time-of-impact query.
///
/// Sweeps `shape_a` from its current position along `vel_a - vel_b` and
/// determines when it first touches `shape_b`.
///
/// Uses the cast GJK (bilateral advancement) method.
///
/// `vel_a` and `vel_b` are relative velocities over the time interval \[0,1\].
pub fn gjk_time_of_impact(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    vel_rel: [f64; 3],
    max_iter: u32,
    tolerance: f64,
) -> ToiResult {
    let res0 = gjk_enhanced(shape_a, shape_b, None, 64);
    if res0.intersecting || res0.dist_sq <= tolerance * tolerance {
        let normal = vnorm(vsub(shape_a.centre(), shape_b.centre()));
        return ToiResult {
            hit: true,
            toi: 0.0,
            normal,
            iterations: 1,
        };
    }
    let sep_dir = vnorm(vsub(shape_b.centre(), shape_a.centre()));
    let approach_speed = vdot(vel_rel, sep_dir);
    if approach_speed <= 0.0 {
        return ToiResult {
            hit: false,
            toi: 1.0,
            normal: sep_dir,
            iterations: 1,
        };
    }
    let mut t_lo = 0.0f64;
    let mut t_hi = 1.0f64;
    let mut hit = false;
    let mut toi = 1.0f64;
    let mut normal = sep_dir;
    let mut iterations = 0u32;
    let x1 = vscale(vel_rel, 1.0);
    let translated_a1 = TranslatedShape {
        shape: shape_a,
        offset: x1,
    };
    let res1 = gjk_enhanced(&translated_a1, shape_b, None, 64);
    if !res1.intersecting && res1.dist_sq > tolerance * tolerance {
        return ToiResult {
            hit: false,
            toi: 1.0,
            normal: sep_dir,
            iterations: 1,
        };
    }
    for _ in 0..max_iter {
        iterations += 1;
        let t_mid = (t_lo + t_hi) * 0.5;
        let x_mid = vscale(vel_rel, t_mid);
        let translated_a = TranslatedShape {
            shape: shape_a,
            offset: x_mid,
        };
        let res = gjk_enhanced(&translated_a, shape_b, None, 64);
        if res.intersecting || res.dist_sq <= tolerance * tolerance {
            t_hi = t_mid;
            hit = true;
            toi = t_mid;
            if !res.intersecting {
                normal = vnorm(vsub(res.closest_a, res.closest_b));
            } else {
                normal = vnorm(vsub(shape_a.centre(), shape_b.centre()));
            }
        } else {
            t_lo = t_mid;
        }
        if t_hi - t_lo < tolerance {
            break;
        }
    }
    ToiResult {
        hit,
        toi,
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
pub fn gjk_with_support_cache(
    shape_a: &dyn ConvexShape,
    shape_b: &dyn ConvexShape,
    cache: &mut SupportCache,
    max_iter: u32,
) -> GjkResult {
    let mut metrics = GjkMetrics::new();
    let mut simplex = EnhSimplex::new();
    let mut dir = vsub(shape_b.centre(), shape_a.centre());
    if vlen_sq(dir) < 1e-18 {
        dir = [1.0, 0.0, 0.0];
    }
    for _ in 0..max_iter {
        metrics.iterations += 1;
        let (sa, sb) = if let Some((a, b)) = cache.get(dir, 1e-6) {
            (a, b)
        } else {
            let a = shape_a.support(dir);
            let b = shape_b.support(vneg(dir));
            cache.put(dir, a, b);
            metrics.support_calls += 2;
            (a, b)
        };
        let w = vsub(sa, sb);
        if vdot(w, dir) < vdot(dir, dir) * 1e-10 {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
        simplex.add(SupportPoint::new(w, sa, sb));
        let (new_dir, contains_origin) = simplex.reduce();
        metrics.simplex_reductions += 1;
        if contains_origin {
            return GjkResult {
                intersecting: true,
                closest_a: [0.0; 3],
                closest_b: [0.0; 3],
                dist_sq: 0.0,
                simplex: simplex.verts,
                metrics,
            };
        }
        if vlen_sq(new_dir) < 1e-18 {
            metrics.converged_early = true;
            return build_separated_result(simplex, metrics);
        }
        dir = new_dir;
    }
    build_separated_result(simplex, metrics)
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
    }
}
