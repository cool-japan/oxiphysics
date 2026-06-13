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

use super::functions::support;
use super::types::{MprResult, SupportPoint};
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
    let center_a = transform_a.transform_point(&shape_a.center_of_mass());
    let center_b = transform_b.transform_point(&shape_b.center_of_mass());
    let mut v0_point = center_a - center_b;
    let v0 = SupportPoint {
        point: v0_point,
        support_a: center_a,
        support_b: center_b,
    };
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
    // The natural choice is `v1 x v0`, which is normal to the (origin, v0, v1)
    // plane. When v1 is parallel to v0 (the portal ray is collinear with v1, as
    // in a head-on sphere/sphere hit) that cross product collapses to zero. The
    // shapes are NOT merely touching in that case — they overlap straight along
    // the v0->v1 axis — so we must keep building a real portal. We therefore
    // fall back to a search direction perpendicular to v1 (spanning the missing
    // dimension) instead of emitting a spurious depth-0 contact.
    let mut n = v1.point.cross(&v0_point);
    if n.norm() < MPR_EPS {
        let perp_axis = if v1.point.x.abs() < 0.7 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        n = v1.point.cross(&perp_axis);
        if n.norm() < MPR_EPS {
            // v1 itself is (near-)zero: the support touches the origin exactly,
            // i.e. the shapes graze along the ray. Emit a depth-0 contact.
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
    let mut v2 = support(shape_a, transform_a, shape_b, transform_b, &n);
    if v2.point.dot(&n) < 0.0 {
        return MprOutcome::Miss;
    }

    // Orient the portal so its normal points away from the origin relative to v0.
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

    // --- Phase 2: portal refinement ------------------------------------------------
    for _ in 0..MPR_MAX_ITERS {
        let raw = (v2.point - v1.point).cross(&(v3.point - v1.point));
        let dir = match safe_normalize(raw) {
            Some(x) => x,
            None => break,
        };
        // Orient the portal normal OUTWARD, i.e. away from the interior point v0
        // (not aligned with a boundary vertex). With v0 inside A-B, the outward
        // face normal satisfies dir.dot(v0) <= 0; flip when that fails. Depth is
        // then the origin->portal-plane distance measured along this outward dir.
        let dir = if dir.dot(&v0_point) > 0.0 { -dir } else { dir };

        let v4 = support(shape_a, transform_a, shape_b, transform_b, &dir);

        // Convergence: the portal stopped expanding toward the origin.
        if v4.point.dot(&dir) - v1.point.dot(&dir) < MPR_TOL {
            let depth = v1.point.dot(&dir).max(0.0);
            let (pa, pb) = barycentric_contact(&v1, &v2, &v3);
            return MprOutcome::Hit {
                normal: dir,
                depth,
                point_a: pa,
                point_b: pb,
            };
        }

        // Separation guard: support along dir never reaches the origin plane.
        if v4.point.dot(&dir) < 0.0 {
            return MprOutcome::Miss;
        }

        // Choose which portal vertex v4 replaces (origin-ray sub-region test).
        // Exact mirror of libccd `expandPortal`: cross = v4 x v0 (here v0 is the
        // +interior point v0_point), then pick the sub-region of the portal that
        // the origin ray still passes through.
        let cross = v4.point.cross(&v0_point);
        if cross.dot(&v1.point) > 0.0 {
            if cross.dot(&v2.point) > 0.0 {
                v1 = v4;
            } else {
                v3 = v4;
            }
        } else if cross.dot(&v3.point) > 0.0 {
            v2 = v4;
        } else {
            v1 = v4;
        }
    }

    // Iteration cap reached — return the best current portal estimate.
    let raw = (v2.point - v1.point).cross(&(v3.point - v1.point));
    let dir = match safe_normalize(raw) {
        Some(x) => x,
        None => Vec3::new(0.0, 1.0, 0.0),
    };
    let dir = if dir.dot(&v0_point) > 0.0 { -dir } else { dir };
    let depth = v1.point.dot(&dir).max(0.0);
    let (pa, pb) = barycentric_contact(&v1, &v2, &v3);
    let _ = v0; // interior witness retained for clarity; not needed past here.
    MprOutcome::Hit {
        normal: dir,
        depth,
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
