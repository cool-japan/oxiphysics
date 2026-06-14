//! Minkowski Portal Refinement (MPR / XenoCollide).
//!
//! Full Snethen XenoCollide portal-refinement algorithm (3-vertex triangular
//! portal). Used as an EPA fallback when the GJK simplex is degenerate.
//!
//! Reference: Gary Snethen, "XenoCollide: Complex Collision Made Simple",
//! Game Programming Gems 7, 2008.

use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::Shape;

use super::functions::{MAX_ITERATIONS, do_simplex, support};
use super::types::{MprResult, Simplex, SupportPoint};
use crate::types::Contact;

/// Maximum portal-refinement iterations.
const MPR_MAX_ITERS: usize = 64;
/// Portal-convergence tolerance.
const MPR_TOL: f64 = 1e-7;
/// Threshold below which a direction vector is treated as degenerate.
const MPR_EPS: f64 = 1e-10;

/// Normalize `v`, returning `None` when it is shorter than [`MPR_EPS`].
///
/// Used everywhere in place of `Vec3::normalize` (which panics on a zero
/// vector) so the portal search can never panic on degenerate input.
fn safe_normalize(v: Vec3) -> Option<Vec3> {
    let n = v.norm();
    if n < MPR_EPS { None } else { Some(v / n) }
}

/// Internal result of the shared MPR core.
enum MprOutcome {
    /// Shapes overlap; carries the contact normal (B->A), penetration depth and
    /// the witness points on shapes A and B.
    Hit {
        normal: Vec3,
        depth: f64,
        point_a: Vec3,
        point_b: Vec3,
    },
    /// Shapes are separated.
    Miss,
}

/// Project the origin onto the portal triangle and recover the contact points
/// on shapes A and B via barycentric interpolation of the per-vertex supports.
fn barycentric_contact(v1: &SupportPoint, v2: &SupportPoint, v3: &SupportPoint) -> (Vec3, Vec3) {
    let a = v1.point;
    let b = v2.point;
    let c = v3.point;
    let v0v = b - a;
    let v1v = c - a;
    let v2v = -a; // origin - a
    let d00 = v0v.dot(&v0v);
    let d01 = v0v.dot(&v1v);
    let d11 = v1v.dot(&v1v);
    let d20 = v2v.dot(&v0v);
    let d21 = v2v.dot(&v1v);
    let denom = d00 * d11 - d01 * d01;
    let (u, v, w) = if denom.abs() > 1e-12 {
        let v = (d11 * d20 - d01 * d21) / denom;
        let w = (d00 * d21 - d01 * d20) / denom;
        let u = 1.0 - v - w;
        (u.max(0.0), v.max(0.0), w.max(0.0))
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    };
    let total = u + v + w;
    let (u, v, w) = if total > 1e-12 {
        (u / total, v / total, w / total)
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    };
    let pa = v1.support_a * u + v2.support_a * v + v3.support_a * w;
    let pb = v1.support_b * u + v2.support_b * v + v3.support_b * w;
    (pa, pb)
}

/// Squared distance from point `p` to segment `[a, b]`, plus the closest point
/// on the segment. Mirror of libccd `__ccdVec3PointSegmentDist2`.
fn point_segment_dist2(p: Vec3, a: Vec3, b: Vec3) -> (f64, Vec3) {
    let d = b - a;
    let len2 = d.dot(&d);
    if len2 < MPR_EPS {
        return ((p - a).dot(&(p - a)), a);
    }
    let t = ((p - a).dot(&d) / len2).clamp(0.0, 1.0);
    let w = a + d * t;
    ((p - w).dot(&(p - w)), w)
}

/// Squared distance from point `p` to triangle `(a, b, c)`, plus the closest
/// point on the triangle. Faithful port of libccd `ccdVec3PointTriDist2`:
/// solves the unconstrained 2-D minimization in the triangle's parametric
/// plane, and when the solution lies outside the triangle (or the triangle is
/// degenerate) falls back to the minimum over the three edges. This edge
/// clamping is what makes the recovered penetration the TRUE minimum distance
/// from the origin to the portal facet.
fn point_tri_dist2(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> (f64, Vec3) {
    // Parametrize T(s, t) = a + s*d1 + t*d2 and minimize |T - p|^2.
    let d1 = b - a;
    let d2 = c - a;
    let ad = a - p;

    let u = d1.dot(&d1);
    let v = d2.dot(&d2);
    let w = d1.dot(&d2);
    let p_coef = d1.dot(&ad);
    let q_coef = d2.dot(&ad);

    let denom = w * w - u * v;

    if denom.abs() >= MPR_EPS {
        let s = (q_coef * w - p_coef * v) / denom;
        let t = (p_coef * w - q_coef * u) / denom;
        if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
            let witness = a + d1 * s + d2 * t;
            return ((p - witness).dot(&(p - witness)), witness);
        }
    }

    // Closest point lies outside the triangle (or it is degenerate): take the
    // minimum over the three edges.
    let (d_ab, w_ab) = point_segment_dist2(p, a, b);
    let (d_ac, w_ac) = point_segment_dist2(p, a, c);
    let (d_bc, w_bc) = point_segment_dist2(p, b, c);
    let mut best = d_ab;
    let mut witness = w_ab;
    if d_ac < best {
        best = d_ac;
        witness = w_ac;
    }
    if d_bc < best {
        best = d_bc;
        witness = w_bc;
    }
    (best, witness)
}

/// Recover the penetration depth, contact normal and witness pair from a
/// converged portal triangle, following libccd `findPenetr`.
///
/// Depth is the distance from the origin to the closest point on the portal
/// triangle `(v1, v2, v3)`; the normal is that origin->closest-point direction
/// (clamped to the triangle, hence to the true minimum-penetration face). When
/// the closest point coincides with the origin (a grazing/touching contact) the
/// direction is recovered from the portal face normal `face_dir` instead.
fn extract_penetration(
    v1: &SupportPoint,
    v2: &SupportPoint,
    v3: &SupportPoint,
    face_dir: Vec3,
) -> MprOutcome {
    let (dist2, witness) = point_tri_dist2(Vec3::zeros(), v1.point, v2.point, v3.point);
    let depth = dist2.max(0.0).sqrt();
    let normal = match safe_normalize(witness) {
        Some(x) => x,
        None => face_dir,
    };
    let (pa, pb) = barycentric_contact(v1, v2, v3);
    MprOutcome::Hit {
        normal,
        depth,
        point_a: pa,
        point_b: pb,
    }
}

/// Shared XenoCollide core driving both [`mpr_full`] and [`mpr_contact`].
///
/// `v0_point` is the Minkowski-space interior point of A-B, taken as
/// `center_A - center_B` (same A-B space as [`support`], whose `point` field is
/// `world_a - world_b`). The origin lies at 0, so the ray we follow is
/// `v0 -> origin`, i.e. `-v0_point`. Portal sign rules follow Snethen's
/// formulation, validated against EPA on sphere/sphere (see module tests).
fn mpr_core(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> MprOutcome {
    // --- Phase 1: portal discovery -------------------------------------------------
    // Interior point v0 of the Minkowski difference A-B, taken as the difference
    // of the two shape centers of mass (guaranteed to lie inside A-B). Note this
    // is center_A - center_B to live in the same A-B space as `support().point`.
    // Only its position `v0_point` is needed (the portal-discovery ray reference);
    // it is never an EPA seed vertex, so no witness pair is kept.
    let center_a = transform_a.transform_point(&shape_a.center_of_mass());
    let center_b = transform_b.transform_point(&shape_b.center_of_mass());
    let mut v0_point = center_a - center_b;
    if v0_point.norm() < MPR_EPS {
        // Centers coincide — perturb along +x so the search has a direction.
        v0_point = Vec3::new(1e-5, 0.0, 0.0);
    }

    // v1: support toward the origin from the interior point.
    let n = -v0_point;
    let mut v1 = support(shape_a, transform_a, shape_b, transform_b, &n);
    if v1.point.dot(&n) < 0.0 {
        return MprOutcome::Miss;
    }

    // v2: support perpendicular to the plane (origin, v0, v1).
    //
    // libccd `discoverPortal` uses `v0 x v1` (the normal of the origin-v0-v1
    // plane). When v1 is parallel to v0 the cross product collapses to zero; this
    // is libccd's return-code 1/2 case (origin lies on v1, or on the v0-v1
    // segment). The shapes still penetrate straight along the v0->v1 axis, and
    // the exact answer is available without building a full portal:
    // depth = |v1.point|, normal = v1.point / |v1.point| (libccd
    // `findPenetrSegment`); or, when v1 itself reaches the origin, a depth-0
    // touching contact (libccd `findPenetrTouch`).
    let n = v0_point.cross(&v1.point);
    if n.norm() < MPR_EPS {
        match safe_normalize(v1.point) {
            Some(normal) => {
                return MprOutcome::Hit {
                    normal,
                    depth: v1.point.norm(),
                    point_a: v1.support_a,
                    point_b: v1.support_b,
                };
            }
            None => {
                let normal = match safe_normalize(-v0_point) {
                    Some(x) => x,
                    None => Vec3::new(0.0, 1.0, 0.0),
                };
                return MprOutcome::Hit {
                    normal,
                    depth: 0.0,
                    point_a: v1.support_a,
                    point_b: v1.support_b,
                };
            }
        }
    }
    let mut v2 = support(shape_a, transform_a, shape_b, transform_b, &n);
    if v2.point.dot(&n) < 0.0 {
        return MprOutcome::Miss;
    }

    // Orient the portal so face (v1, v2) winds "outside" the origin relative to
    // v0 (libccd `discoverPortal` swap step): if `(v1-v0) x (v2-v0)` aligns with
    // v0, swap v1<->v2 so the maintained winding stays outward.
    let mut n = (v1.point - v0_point).cross(&(v2.point - v0_point));
    if n.dot(&v0_point) > 0.0 {
        std::mem::swap(&mut v1, &mut v2);
        n = -n;
    }

    // Third portal vertex; refined inside the discovery loop. Seeded so it is in
    // scope for Phase 2 even if the loop body never assigns it.
    let mut v3 = support(shape_a, transform_a, shape_b, transform_b, &n);

    // Discovery loop. The portal normal `n` is re-derived each pass; once it can
    // no longer be normalized the search has degenerated and we stop (this is the
    // `while let` exit, equivalent to the classic `match { None => break }`).
    let mut iter = 0usize;
    while let Some(nn) = safe_normalize(n) {
        v3 = support(shape_a, transform_a, shape_b, transform_b, &nn);
        if v3.point.dot(&nn) < 0.0 {
            return MprOutcome::Miss;
        }
        // Does the origin ray (v0 -> origin) pass through triangle (v1, v2, v3)?
        if v1.point.cross(&v3.point).dot(&v0_point) < 0.0 {
            v2 = v3;
            n = (v1.point - v0_point).cross(&(v2.point - v0_point));
        } else if v3.point.cross(&v2.point).dot(&v0_point) < 0.0 {
            v1 = v3;
            n = (v1.point - v0_point).cross(&(v2.point - v0_point));
        } else {
            break; // portal found: (v1, v2, v3) brackets the origin ray.
        }
        iter += 1;
        if iter > MPR_MAX_ITERS {
            break;
        }
    }

    // --- Phase 2: minimum-penetration refinement (EPA) ----------------------------
    // The discovery portal (v1, v2, v3) confirms the overlap and brackets the
    // origin ray from the interior point v0. Plain XenoCollide portal-walking would
    // then report the penetration along whatever face the v0->origin ray pierces,
    // which for deeply overlapping / rotated / curved shapes is frequently NOT the
    // minimum-translation face. We therefore finish with an expanding-polytope
    // (EPA) refinement, which converges to the true global minimum-penetration face.
    //
    // EPA must START from an origin-enclosing tetrahedron of Minkowski-boundary
    // support points whose closest face is already no farther from the origin than
    // the true penetration — EPA only ever pushes the closest face OUTWARD, so a
    // seed hull that is too large would over-report shallow contacts. The natural
    // such seed is the GJK termination simplex (its closest face hugs the contact),
    // so we build it here with a focused GJK pass warm-started from the portal
    // points. The portal already proved overlap, so this GJK pass terminates inside.
    let portal_dir = safe_normalize(-v0_point).unwrap_or_else(|| Vec3::new(1.0, 0.0, 0.0));
    let tetra = match gjk_origin_tetra(shape_a, transform_a, shape_b, transform_b, v1, portal_dir) {
        Some(t) => t,
        // GJK could not refine an origin-enclosing simplex (extremely shallow /
        // degenerate). Recover the penetration from the portal triangle directly.
        None => {
            let face_dir = safe_normalize((v2.point - v1.point).cross(&(v3.point - v1.point)))
                .unwrap_or(portal_dir);
            return extract_penetration(&v1, &v2, &v3, face_dir);
        }
    };

    epa_min_penetration(shape_a, transform_a, shape_b, transform_b, tetra)
}

/// Build a tetrahedral GJK simplex that strictly contains the origin, in
/// Minkowski-difference space, using the standard GJK descent.
///
/// Warm-started from the boundary support `warm` (a portal vertex, near the
/// contact) it runs the shared [`do_simplex`] machinery until the simplex encloses
/// the origin (`do_simplex` returns `None`), then returns those four boundary
/// support points — the tightest origin-enclosing seed for EPA (its closest face
/// hugs the actual contact, so EPA does not over-report shallow penetrations).
/// `initial_dir` seeds the search when `warm` already sits at the origin. Returns
/// `None` if the descent fails to enclose the origin within the iteration budget
/// or collapses to a sub-tetrahedral simplex.
fn gjk_origin_tetra(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    warm: SupportPoint,
    initial_dir: Vec3,
) -> Option<[SupportPoint; 4]> {
    let mut simplex = Simplex::new();
    simplex.push(warm);
    let mut direction = safe_normalize(-warm.point).unwrap_or(initial_dir);

    for _ in 0..MAX_ITERATIONS {
        let d = match safe_normalize(direction) {
            Some(x) => x,
            None => break,
        };
        let w = support(shape_a, transform_a, shape_b, transform_b, &d);
        if w.point.dot(&d) < 0.0 {
            // Should not happen for a confirmed overlap, but guard anyway.
            return None;
        }
        simplex.push(w);
        match do_simplex(&mut simplex) {
            Some(new_dir) => direction = new_dir,
            None => {
                // Origin enclosed: a 4-point simplex is the EPA seed.
                if simplex.len() == 4 {
                    return Some([
                        simplex.points[0],
                        simplex.points[1],
                        simplex.points[2],
                        simplex.points[3],
                    ]);
                }
                return None;
            }
        }
    }
    None
}

/// A polytope face: indices into the EPA vertex list, plus the cached outward
/// unit normal (pointing away from the origin) and origin-plane distance.
struct EpaFace {
    verts: [usize; 3],
    normal: Vec3,
    dist: f64,
}

/// Build an [`EpaFace`] for triangle `(a, b, c)` (vertex indices into `pts`),
/// orienting the normal to point AWAY from the origin (the origin is interior to
/// the polytope, so the outward normal satisfies `normal.dot(a) >= 0`). Returns
/// `None` when the triangle is degenerate (zero area).
fn make_face(pts: &[SupportPoint], a: usize, b: usize, c: usize) -> Option<EpaFace> {
    let pa = pts[a].point;
    let pb = pts[b].point;
    let pc = pts[c].point;
    let normal = safe_normalize((pb - pa).cross(&(pc - pa)))?;
    let dist = normal.dot(&pa);
    let (normal, dist) = if dist < 0.0 {
        (-normal, -dist)
    } else {
        (normal, dist)
    };
    Some(EpaFace {
        verts: [a, b, c],
        normal,
        dist,
    })
}

/// Insert a new vertex (index `wi` into `pts`) into a polytope `faces` by the
/// quickhull horizon step: remove every face the vertex can see and re-stitch the
/// horizon edges to it, preserving convexity. Returns `false` (leaving `faces`
/// untouched) when the vertex sees no face — i.e. it is already inside the hull,
/// so the polytope is unchanged and the EPA expansion has converged.
fn insert_vertex(pts: &[SupportPoint], faces: &mut Vec<EpaFace>, wi: usize) -> bool {
    let wp = pts[wi].point;
    let mut horizon: Vec<(usize, usize)> = Vec::with_capacity(faces.len());
    let mut kept: Vec<EpaFace> = Vec::with_capacity(faces.len());
    let mut saw_any = false;
    for f in faces.drain(..) {
        let fp = pts[f.verts[0]].point;
        if f.normal.dot(&(wp - fp)) > MPR_TOL {
            saw_any = true;
            let [i0, i1, i2] = f.verts;
            for &(s, t) in &[(i0, i1), (i1, i2), (i2, i0)] {
                // Cancel an edge already present (shared by two visible faces);
                // otherwise add it. Both orderings are checked so the cancellation
                // is robust to per-face winding chosen by `make_face`.
                if let Some(pos) = horizon
                    .iter()
                    .position(|&(a, b)| (a == t && b == s) || (a == s && b == t))
                {
                    horizon.swap_remove(pos);
                } else {
                    horizon.push((s, t));
                }
            }
        } else {
            kept.push(f);
        }
    }
    *faces = kept;
    if !saw_any {
        return false;
    }
    for (s, t) in horizon {
        if let Some(f) = make_face(pts, s, t, wi) {
            faces.push(f);
        }
    }
    true
}

/// Expanding Polytope Algorithm seeded from an origin-enclosing tetrahedron of
/// Minkowski-boundary support points.
///
/// Iteratively pushes the closest polytope face toward the boundary until it can
/// no longer expand, yielding the minimum penetration depth, the contact normal
/// (B->A) and the witness pair. Self-contained and bounded.
///
/// The closest face's origin-distance increases monotonically as the polytope is
/// refined toward the boundary, so the answer is the closest face at convergence;
/// the running best is updated each iteration *before* expanding and kept if a
/// numerically degenerate re-stitch ends the loop early.
fn epa_min_penetration(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    tetra: [SupportPoint; 4],
) -> MprOutcome {
    let mut pts: Vec<SupportPoint> = tetra.to_vec();

    // Seed the four tetrahedron faces (each omitting one vertex), oriented outward.
    let mut faces: Vec<EpaFace> = Vec::with_capacity(16);
    for &(a, b, c) in &[(0usize, 1usize, 2usize), (0, 1, 3), (0, 2, 3), (1, 2, 3)] {
        if let Some(f) = make_face(&pts, a, b, c) {
            faces.push(f);
        }
    }
    if faces.len() < 4 {
        // Degenerate tetra: recover a contact from three of its vertices.
        let face_dir =
            safe_normalize((pts[1].point - pts[0].point).cross(&(pts[2].point - pts[0].point)))
                .unwrap_or_else(|| Vec3::new(0.0, 1.0, 0.0));
        return extract_penetration(&pts[0], &pts[1], &pts[2], face_dir);
    }

    let mut best_face = faces[0].verts;
    let mut best_normal = faces[0].normal;
    let mut best_dist = faces[0].dist;

    for _ in 0..MPR_MAX_ITERS {
        // Closest face to the origin in the current polytope becomes the running
        // answer (EPA's closest-face distance is monotone non-decreasing).
        let mut ci = 0usize;
        let mut cdist = f64::INFINITY;
        for (i, f) in faces.iter().enumerate() {
            if f.dist < cdist {
                cdist = f.dist;
                ci = i;
            }
        }
        let close = &faces[ci];
        best_face = close.verts;
        best_normal = close.normal;
        best_dist = close.dist;

        // Support in the outward face-normal direction.
        let w = support(shape_a, transform_a, shape_b, transform_b, &best_normal);
        let w_dist = w.point.dot(&best_normal);

        // Converged: the support adds no measurable extent beyond the closest face.
        if w_dist - best_dist < MPR_TOL {
            break;
        }

        // A duplicate support means the closest face is already on the boundary in
        // its normal direction and cannot expand further; the loop has converged.
        let wp = w.point;
        if pts
            .iter()
            .any(|p| (p.point - wp).norm_squared() < MPR_EPS * MPR_EPS)
        {
            break;
        }

        // Insert the new boundary vertex; on a degenerate re-stitch (sees no face
        // or empties the polytope) stop and keep the closest face found so far.
        let wi = pts.len();
        pts.push(w);
        if !insert_vertex(&pts, &mut faces, wi) || faces.is_empty() {
            break;
        }
    }

    // Recover witnesses from the closest face via barycentric interpolation.
    let [ia, ib, ic] = best_face;
    let (pa, pb) = barycentric_contact(&pts[ia], &pts[ib], &pts[ic]);
    let normal = match safe_normalize(best_normal) {
        Some(x) => x,
        None => Vec3::new(0.0, 1.0, 0.0),
    };
    MprOutcome::Hit {
        normal,
        depth: best_dist.max(0.0),
        point_a: pa,
        point_b: pb,
    }
}

/// Full Snethen XenoCollide MPR (3-vertex portal version).
///
/// Returns penetration depth, the separating/contact normal (pointing from B
/// towards A), and a contact point when the shapes intersect; otherwise
/// returns `MprResult::Separated`.
pub fn mpr_full(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> MprResult {
    match mpr_core(shape_a, transform_a, shape_b, transform_b) {
        MprOutcome::Hit {
            normal,
            depth,
            point_a,
            point_b,
        } => MprResult::Intersecting {
            normal,
            depth,
            point: (point_a + point_b) * 0.5,
        },
        MprOutcome::Miss => MprResult::Separated,
    }
}

/// Like [`mpr_full`] but returns the full witness pair as an EPA-compatible
/// [`Contact`] when intersecting, for use as an EPA fallback.
pub fn mpr_contact(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> Option<Contact> {
    match mpr_core(shape_a, transform_a, shape_b, transform_b) {
        MprOutcome::Hit {
            normal,
            depth,
            point_a,
            point_b,
        } => Some(Contact::new(point_a, point_b, normal, depth)),
        MprOutcome::Miss => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_geometry::Sphere;

    #[test]
    fn mpr_full_unit_spheres_overlap() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        match mpr_full(&s1, &t1, &s2, &t2) {
            MprResult::Intersecting { normal, depth, .. } => {
                assert!(
                    (depth - 1.0).abs() < 1e-3,
                    "expected depth approx 1.0, got {depth}"
                );
                assert!(
                    normal.x.abs() > 0.99,
                    "expected normal approx +/-x, got {normal:?}"
                );
            }
            MprResult::Separated => panic!("unit spheres at distance 1 must intersect"),
        }
    }

    #[test]
    fn mpr_full_far_spheres_separated() {
        let s1 = Sphere::new(0.5);
        let s2 = Sphere::new(0.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(matches!(mpr_full(&s1, &t1, &s2, &t2), MprResult::Separated));
    }

    #[test]
    fn mpr_contact_unit_spheres_witnesses() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let contact = mpr_contact(&s1, &t1, &s2, &t2).expect("spheres overlap");
        assert!(
            (contact.depth - 1.0).abs() < 1e-3,
            "depth {}",
            contact.depth
        );
        assert!(contact.normal.x.abs() > 0.99, "normal {:?}", contact.normal);
    }
}
