//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

use super::types::{
    CcdResult, Gjk, GjkDistanceResult, GjkResult, MprResult, RayCastResult, Simplex,
    SpeculativeContact, SupportPoint,
};
use oxiphysics_geometry::Sphere;

/// Maximum iterations for GJK convergence.
pub(crate) const MAX_ITERATIONS: usize = 64;
/// Tolerance for GJK termination.
pub(crate) const TOLERANCE: Real = 1e-8;
/// Compute support point in Minkowski difference.
pub fn support(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    direction: &Vec3,
) -> SupportPoint {
    let local_dir_a = transform_a.rotation.inverse() * direction;
    let local_dir_b = transform_b.rotation.inverse() * (-direction);
    let local_a = shape_a.support_point(&local_dir_a);
    let local_b = shape_b.support_point(&local_dir_b);
    let world_a = transform_a.transform_point(&local_a);
    let world_b = transform_b.transform_point(&local_b);
    SupportPoint {
        point: world_a - world_b,
        support_a: world_a,
        support_b: world_b,
    }
}
/// Process the simplex and return a new search direction, or None if origin is contained.
pub fn do_simplex(simplex: &mut Simplex) -> Option<Vec3> {
    match simplex.len() {
        2 => do_simplex_line(simplex),
        3 => do_simplex_triangle(simplex),
        4 => do_simplex_tetrahedron(simplex),
        _ => Some(Vec3::new(1.0, 0.0, 0.0)),
    }
}
pub fn do_simplex_line(simplex: &mut Simplex) -> Option<Vec3> {
    let a = simplex.points[1].point;
    let b = simplex.points[0].point;
    let ab = b - a;
    let ao = -a;
    if ab.dot(&ao) > 0.0 {
        let dir = ab.cross(&ao).cross(&ab);
        if dir.norm_squared() < TOLERANCE {
            let perp = if ab.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            Some(ab.cross(&perp))
        } else {
            Some(dir)
        }
    } else {
        simplex.points = vec![simplex.points[1]];
        Some(ao)
    }
}
pub fn do_simplex_triangle(simplex: &mut Simplex) -> Option<Vec3> {
    let a = simplex.points[2].point;
    let b = simplex.points[1].point;
    let c = simplex.points[0].point;
    let ab = b - a;
    let ac = c - a;
    let ao = -a;
    let abc = ab.cross(&ac);
    if abc.cross(&ac).dot(&ao) > 0.0 {
        if ac.dot(&ao) > 0.0 {
            simplex.points = vec![simplex.points[0], simplex.points[2]];
            Some(ac.cross(&ao).cross(&ac))
        } else {
            simplex.points = vec![simplex.points[1], simplex.points[2]];
            do_simplex_line(simplex)
        }
    } else if ab.cross(&abc).dot(&ao) > 0.0 {
        simplex.points = vec![simplex.points[1], simplex.points[2]];
        do_simplex_line(simplex)
    } else if abc.dot(&ao) > 0.0 {
        Some(abc)
    } else {
        simplex.points = vec![simplex.points[1], simplex.points[0], simplex.points[2]];
        Some(-abc)
    }
}
pub fn do_simplex_tetrahedron(simplex: &mut Simplex) -> Option<Vec3> {
    let a = simplex.points[3].point;
    let b = simplex.points[2].point;
    let c = simplex.points[1].point;
    let d = simplex.points[0].point;
    let ao = -a;
    let ab = b - a;
    let ac = c - a;
    let ad = d - a;
    let abc = ab.cross(&ac);
    let acd = ac.cross(&ad);
    let adb = ad.cross(&ab);
    if abc.dot(&ao) > 0.0 {
        simplex.points = vec![simplex.points[1], simplex.points[2], simplex.points[3]];
        do_simplex_triangle(simplex)
    } else if acd.dot(&ao) > 0.0 {
        simplex.points = vec![simplex.points[0], simplex.points[1], simplex.points[3]];
        do_simplex_triangle(simplex)
    } else if adb.dot(&ao) > 0.0 {
        simplex.points = vec![simplex.points[2], simplex.points[0], simplex.points[3]];
        do_simplex_triangle(simplex)
    } else {
        None
    }
}
pub fn closest_point_line(a: &Vec3, b: &Vec3) -> (Vec3, Vec<Real>) {
    let ab = b - a;
    let t = (-a).dot(&ab) / ab.dot(&ab);
    let t = t.clamp(0.0, 1.0);
    let point = a + ab * t;
    (point, vec![1.0 - t, t])
}
pub fn closest_point_triangle(a: &Vec3, b: &Vec3, c: &Vec3) -> (Vec3, Vec<Real>) {
    let ab = b - a;
    let ac = c - a;
    let ao = -a;
    let d1 = ab.dot(&ao);
    let d2 = ac.dot(&ao);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (*a, vec![1.0, 0.0, 0.0]);
    }
    let bo = -b;
    let d3 = ab.dot(&bo);
    let d4 = ac.dot(&bo);
    if d3 >= 0.0 && d4 <= d3 {
        return (*b, vec![0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (a + ab * v, vec![1.0 - v, v, 0.0]);
    }
    let co = -c;
    let d5 = ab.dot(&co);
    let d6 = ac.dot(&co);
    if d6 >= 0.0 && d5 <= d6 {
        return (*c, vec![0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (a + ac * w, vec![1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (b + (c - b) * w, vec![0.0, 1.0 - w, w]);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    (a + ab * v + ac * w, vec![1.0 - v - w, v, w])
}
/// Check if a simplex is degenerate (has near-zero volume/area).
#[allow(dead_code)]
pub fn is_simplex_degenerate(simplex: &Simplex) -> bool {
    match simplex.len() {
        0 | 1 => false,
        2 => {
            let edge = simplex.points[1].point - simplex.points[0].point;
            edge.norm_squared() < TOLERANCE * TOLERANCE
        }
        3 => {
            let ab = simplex.points[1].point - simplex.points[0].point;
            let ac = simplex.points[2].point - simplex.points[0].point;
            let cross = ab.cross(&ac);
            cross.norm_squared() < TOLERANCE * TOLERANCE
        }
        4 => {
            let ab = simplex.points[1].point - simplex.points[0].point;
            let ac = simplex.points[2].point - simplex.points[0].point;
            let ad = simplex.points[3].point - simplex.points[0].point;
            let vol = ab.dot(&ac.cross(&ad));
            vol.abs() < TOLERANCE * TOLERANCE * TOLERANCE
        }
        _ => true,
    }
}
/// Attempt to fix a degenerate simplex by removing problematic points.
///
/// Returns `true` if the simplex was modified.
#[allow(dead_code)]
pub fn fix_degenerate_simplex(simplex: &mut Simplex) -> bool {
    if simplex.len() < 2 {
        return false;
    }
    let mut i = 0;
    let mut modified = false;
    while i < simplex.points.len() {
        let mut j = i + 1;
        while j < simplex.points.len() {
            let diff = simplex.points[j].point - simplex.points[i].point;
            if diff.norm_squared() < TOLERANCE * TOLERANCE {
                simplex.points.remove(j);
                modified = true;
            } else {
                j += 1;
            }
        }
        i += 1;
    }
    modified
}
/// Compute the signed distance between two convex shapes using GJK.
///
/// Returns the Minkowski difference distance (positive if separated,
/// negative estimate if overlapping).
#[allow(dead_code)]
pub fn gjk_distance(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> f64 {
    match Gjk::query(shape_a, transform_a, shape_b, transform_b) {
        GjkResult::Separated { distance, .. } => distance,
        GjkResult::Intersecting(_) => -1.0,
    }
}
/// Compute the Minkowski difference support in a given direction
/// and return both the support point and the witness points.
#[allow(dead_code)]
pub fn minkowski_support(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    direction: &Vec3,
) -> SupportPoint {
    support(shape_a, transform_a, shape_b, transform_b, direction)
}
/// Compute support point for a *uniformly scaled* shape.
///
/// Scaling `s > 1` inflates the shape; `s < 1` shrinks it.
#[allow(dead_code)]
pub fn scaled_support(
    shape: &dyn Shape,
    transform: &Transform,
    direction: &Vec3,
    scale: f64,
) -> SupportPoint {
    let local_dir = transform.rotation.inverse() * direction;
    let local_pt = shape.support_point(&local_dir);
    let world_pt = transform.transform_point(&(local_pt * scale));
    SupportPoint {
        point: world_pt,
        support_a: world_pt,
        support_b: Vec3::zeros(),
    }
}
/// Compute Minkowski difference support for scaled shapes.
pub fn scaled_support_minkowski(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    scale_a: f64,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    scale_b: f64,
    direction: &Vec3,
) -> SupportPoint {
    let local_dir_a = transform_a.rotation.inverse() * direction;
    let local_dir_b = transform_b.rotation.inverse() * (-direction);
    let local_a = shape_a.support_point(&local_dir_a) * scale_a;
    let local_b = shape_b.support_point(&local_dir_b) * scale_b;
    let world_a = transform_a.transform_point(&local_a);
    let world_b = transform_b.transform_point(&local_b);
    SupportPoint {
        point: world_a - world_b,
        support_a: world_a,
        support_b: world_b,
    }
}
/// Compute the time-of-impact (TOI) between two moving convex shapes.
///
/// Uses a conservative advancement (linear cast) approach:
/// iteratively advances `t` until the shapes touch or `max_toi` is reached.
///
/// Returns `Some(t)` where `t ∈ [0, max_toi]` is the first time the shapes
/// touch, or `None` if they do not collide within `[0, max_toi]`.
#[allow(dead_code)]
pub fn gjk_time_of_impact(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    velocity_a: Vec3,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    velocity_b: Vec3,
    max_toi: f64,
) -> Option<f64> {
    if Gjk::intersect(shape_a, transform_a, shape_b, transform_b) {
        return None;
    }
    let rel_vel = velocity_a - velocity_b;
    let rel_speed = rel_vel.norm();
    if rel_speed < 1e-14 {
        return None;
    }
    let mut t = 0.0_f64;
    let max_iter = 64;
    for _ in 0..max_iter {
        let ta = Transform {
            position: transform_a.position + velocity_a * t,
            rotation: transform_a.rotation,
        };
        let tb = Transform {
            position: transform_b.position + velocity_b * t,
            rotation: transform_b.rotation,
        };
        match Gjk::query(shape_a, &ta, shape_b, &tb) {
            GjkResult::Intersecting(_) => {
                return Some(t);
            }
            GjkResult::Separated { distance, .. } => {
                if distance < 1e-6 {
                    return Some(t);
                }
                let step = distance / rel_speed;
                t += step;
                if t > max_toi {
                    return None;
                }
            }
        }
    }
    None
}
/// Refine the penetration direction from a GJK simplex.
///
/// When GJK terminates with an intersecting simplex, this function returns
/// an approximate penetration direction by finding the face of the simplex
/// that is closest to the origin and returning its outward normal.
///
/// For a full penetration vector, use the EPA algorithm after GJK.
#[allow(dead_code)]
pub fn refine_penetration_direction(simplex: &Simplex) -> Vec3 {
    match simplex.len() {
        0 | 1 => Vec3::new(0.0, 0.0, 1.0),
        2 => {
            let edge = simplex.points[1].point - simplex.points[0].point;
            let perp = if edge.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let normal = edge.cross(&perp);
            if normal.norm_squared() < TOLERANCE {
                perp
            } else {
                normal.normalize()
            }
        }
        3 => {
            let a = simplex.points[0].point;
            let b = simplex.points[1].point;
            let c = simplex.points[2].point;
            let normal = (b - a).cross(&(c - a));
            if normal.norm_squared() < TOLERANCE {
                Vec3::new(0.0, 0.0, 1.0)
            } else {
                let n = normal.normalize();
                if n.dot(&(-a)) < 0.0 { -n } else { n }
            }
        }
        _ => {
            let pts = &simplex.points;
            let faces: [(usize, usize, usize); 4] = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];
            let mut best_normal = Vec3::new(0.0, 0.0, 1.0);
            let mut best_dist = f64::MAX;
            for &(i, j, k) in &faces {
                let a = pts[i].point;
                let b = pts[j].point;
                let c = pts[k].point;
                let raw = (b - a).cross(&(c - a));
                if raw.norm_squared() < TOLERANCE {
                    continue;
                }
                let n = raw.normalize();
                let n = if n.dot(&(-a)) < 0.0 { -n } else { n };
                let dist = n.dot(&a).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best_normal = n;
                }
            }
            best_normal
        }
    }
}
/// Estimate the penetration depth between two overlapping shapes.
///
/// Uses the GJK simplex to approximate the minimum translation distance.
/// For a precise depth, use EPA.
#[allow(dead_code)]
pub fn estimate_penetration_depth(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> f64 {
    match Gjk::query(shape_a, transform_a, shape_b, transform_b) {
        GjkResult::Separated { distance, .. } => -distance,
        GjkResult::Intersecting(simplex) => {
            let dirs: [Vec3; 6] = [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(-1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, -1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -1.0),
            ];
            let _ = simplex;
            let mut min_depth = f64::MAX;
            for dir in &dirs {
                let sp = support(shape_a, transform_a, shape_b, transform_b, dir);
                let proj = sp.point.dot(dir);
                if proj > 0.0 && proj < min_depth {
                    min_depth = proj;
                }
            }
            if min_depth == f64::MAX {
                0.0
            } else {
                min_depth
            }
        }
    }
}
/// Perform a precise GJK distance query tracking closest-point witnesses.
///
/// Unlike the basic `Gjk::query`, this variant accumulates barycentric
/// weights at each step to reconstruct world-space witness points accurately.
#[allow(dead_code)]
pub fn gjk_distance_query(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> GjkDistanceResult {
    let mut simplex = Simplex::new();
    let mut direction = transform_b.position - transform_a.position;
    if direction.norm_squared() < TOLERANCE {
        direction = Vec3::new(1.0, 0.0, 0.0);
    }
    let first = support(shape_a, transform_a, shape_b, transform_b, &direction);
    simplex.push(first);
    direction = -first.point;
    let mut iterations = 0usize;
    for _ in 0..MAX_ITERATIONS {
        iterations += 1;
        let new_pt = support(shape_a, transform_a, shape_b, transform_b, &direction);
        if new_pt.point.dot(&direction) < -TOLERANCE {
            let (closest_md, bary) = simplex.closest_point_to_origin();
            let dist = closest_md.norm();
            let (closest_a, closest_b) = barycentric_witnesses(&simplex, &bary);
            return GjkDistanceResult {
                distance: dist,
                closest_a,
                closest_b,
                iterations,
            };
        }
        simplex.push(new_pt);
        if let Some(new_dir) = do_simplex(&mut simplex) {
            direction = new_dir;
        } else {
            return GjkDistanceResult {
                distance: -1.0,
                closest_a: transform_a.position,
                closest_b: transform_b.position,
                iterations,
            };
        }
    }
    GjkDistanceResult {
        distance: -1.0,
        closest_a: transform_a.position,
        closest_b: transform_b.position,
        iterations,
    }
}
/// Reconstruct world-space witness points from a simplex and barycentric coordinates.
pub fn barycentric_witnesses(simplex: &Simplex, bary: &[Real]) -> (Vec3, Vec3) {
    let n = simplex.len().min(bary.len());
    if n == 0 {
        return (Vec3::zeros(), Vec3::zeros());
    }
    let mut wa = Vec3::zeros();
    let mut wb = Vec3::zeros();
    for i in 0..n {
        wa += simplex.points[i].support_a * bary[i];
        wb += simplex.points[i].support_b * bary[i];
    }
    (wa, wb)
}
/// Cast a ray against a convex shape using the GJK algorithm.
///
/// `ray_origin` and `ray_dir` define the ray in world space.
/// `max_t` limits the ray length.
///
/// Returns `Some(RayCastResult)` if the ray hits the shape within `[0, max_t]`.
#[allow(dead_code)]
pub fn gjk_ray_cast(
    shape: &dyn Shape,
    transform: &Transform,
    ray_origin: Vec3,
    ray_dir: Vec3,
    max_t: f64,
) -> Option<RayCastResult> {
    let dir_len = ray_dir.norm();
    if dir_len < 1e-14 {
        return None;
    }
    let unit_dir = ray_dir / dir_len;
    let mut t = 0.0_f64;
    let mut normal = Vec3::new(0.0, 1.0, 0.0);
    let start_origin = ray_origin;
    for _ in 0..MAX_ITERATIONS {
        let query_pos = start_origin + unit_dir * t;
        let query_t = Transform {
            position: query_pos,
            rotation: transform.rotation,
        };
        let sp_fwd = {
            let local = transform.rotation.inverse() * unit_dir;
            let lp = shape.support_point(&local);
            transform.transform_point(&lp)
        };
        let sp_bwd = {
            let local = transform.rotation.inverse() * (-unit_dir);
            let lp = shape.support_point(&local);
            transform.transform_point(&lp)
        };
        let to_fwd = sp_fwd - query_t.position;
        let to_bwd = sp_bwd - query_t.position;
        if to_fwd.dot(&unit_dir) < 0.0 && (-to_bwd).dot(&unit_dir) < 0.0 {
            return Some(RayCastResult {
                t,
                point: query_pos,
                normal,
            });
        }
        let diff = sp_bwd - query_pos;
        let proj = diff.dot(&unit_dir);
        if proj.abs() < 1e-6 {
            let surf_diff = query_pos - sp_bwd;
            let surf_len = surf_diff.norm();
            normal = if surf_len > 1e-10 {
                surf_diff / surf_len
            } else {
                unit_dir
            };
            return Some(RayCastResult {
                t,
                point: query_pos,
                normal,
            });
        }
        if proj > 0.0 {
            t += proj;
            let n = diff;
            let nl = n.norm();
            normal = if nl > 1e-10 { n / nl } else { normal };
        } else {
            break;
        }
        if t > max_t {
            break;
        }
    }
    let local_origin = transform.inverse().transform_point(&start_origin);
    let local_dir = transform.rotation.inverse() * unit_dir;
    if let Some(hit) = shape.ray_cast(&local_origin, &local_dir, max_t) {
        let world_point = transform.transform_point(&hit.point);
        let world_normal = transform.transform_vector(&hit.normal);
        let nl = world_normal.norm();
        let world_normal = if nl > 1e-10 {
            world_normal / nl
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        return Some(RayCastResult {
            t: hit.toi,
            point: world_point,
            normal: world_normal,
        });
    }
    None
}
/// Cast a ray against the Minkowski sum of two convex shapes (GJK portal ray-cast).
///
/// Returns the `t` value at which the cast ray first touches the Minkowski difference.
#[allow(dead_code)]
pub fn gjk_cast_ray_minkowski(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    _ray_origin: Vec3,
    ray_dir: Vec3,
    max_t: f64,
) -> Option<f64> {
    let dir_len = ray_dir.norm();
    if dir_len < 1e-14 {
        return None;
    }
    let unit_dir = ray_dir / dir_len;
    let mut t = 0.0_f64;
    for _ in 0..MAX_ITERATIONS {
        let ta_adv = Transform {
            position: transform_a.position - unit_dir * t,
            rotation: transform_a.rotation,
        };
        match Gjk::query(shape_a, &ta_adv, shape_b, transform_b) {
            GjkResult::Intersecting(_) => return Some(t),
            GjkResult::Separated { distance, .. } => {
                if distance < 1e-6 {
                    return Some(t);
                }
                t += distance;
                if t > max_t {
                    return None;
                }
            }
        }
    }
    None
}
/// Minkowski Portal Refinement algorithm for convex-convex intersection.
///
/// MPR (XenoCollide) finds a "portal" through the Minkowski difference that
/// contains the origin, confirming intersection and providing a contact normal.
///
/// Reference: Gary Snethen, "XenoCollide: Complex Collision Made Simple",
/// Game Programming Gems 7, 2008.
#[allow(dead_code)]
pub fn mpr_query(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> MprResult {
    let v0_vec = transform_b.position - transform_a.position;
    if v0_vec.norm_squared() < TOLERANCE * TOLERANCE {
        return MprResult::Intersecting {
            normal: Vec3::new(0.0, 1.0, 0.0),
            depth: 0.0,
            point: transform_a.position,
        };
    }
    let v0 = SupportPoint {
        point: v0_vec,
        support_a: transform_a.position,
        support_b: transform_b.position,
    };
    let dir0 = -v0_vec;
    let v1 = support(shape_a, transform_a, shape_b, transform_b, &dir0);
    if v1.point.dot(&dir0) < 0.0 {
        return MprResult::Separated;
    }
    let dir1 = dir0.cross(&v1.point).cross(&v1.point);
    let dir1 = if dir1.norm_squared() < TOLERANCE {
        let perp = if dir0.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        dir0.cross(&perp)
    } else {
        dir1
    };
    let mut v1 = v1;
    let mut v2 = support(shape_a, transform_a, shape_b, transform_b, &dir1);
    for _ in 0..MAX_ITERATIONS {
        let portal_normal = (v1.point - v0.point).cross(&(v2.point - v0.point));
        let pn_len = portal_normal.norm();
        if pn_len < TOLERANCE {
            break;
        }
        let pn = portal_normal / pn_len;
        if pn.dot(&v1.point) >= 0.0 {
            let depth = pn.dot(&v1.point).abs();
            let point = (v1.point + v2.point) * 0.5;
            return MprResult::Intersecting {
                normal: pn,
                depth,
                point,
            };
        }
        let v3 = support(shape_a, transform_a, shape_b, transform_b, &pn);
        if v3.point.dot(&pn) < 0.0 {
            return MprResult::Separated;
        }
        let d1 = (v3.point - v0.point).cross(&(v1.point - v0.point));
        if d1.dot(&(-v0.point)) > 0.0 {
            v2 = v3;
        } else {
            v1 = v3;
        }
    }
    MprResult::Separated
}
/// Internal recursive portal refinement step for MPR.
pub fn mpr_refine_portal(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    v0: &SupportPoint,
    v1: &SupportPoint,
    v2: &SupportPoint,
    depth: usize,
) -> MprResult {
    if depth > MAX_ITERATIONS {
        return MprResult::Separated;
    }
    let portal_normal = (v1.point - v0.point).cross(&(v2.point - v0.point));
    let pn_len = portal_normal.norm();
    if pn_len < TOLERANCE {
        return MprResult::Separated;
    }
    let pn = portal_normal / pn_len;
    if pn.dot(&v1.point) >= 0.0 {
        let penetration_depth = pn.dot(&v1.point).abs();
        let contact_point = (v1.point + v2.point) * 0.5;
        return MprResult::Intersecting {
            normal: pn,
            depth: penetration_depth,
            point: contact_point,
        };
    }
    let v3 = support(shape_a, transform_a, shape_b, transform_b, &pn);
    if v3.point.dot(&pn) < 0.0 {
        return MprResult::Separated;
    }
    let d1 = (v3.point - v0.point).cross(&(v1.point - v0.point));
    if d1.dot(&(-v0.point)) > 0.0 {
        mpr_refine_portal(
            shape_a,
            transform_a,
            shape_b,
            transform_b,
            v0,
            v1,
            &v3,
            depth + 1,
        )
    } else {
        mpr_refine_portal(
            shape_a,
            transform_a,
            shape_b,
            transform_b,
            v0,
            &v3,
            v2,
            depth + 1,
        )
    }
}
/// Compute time-of-impact using a bisection refinement after conservative advancement.
///
/// Unlike `gjk_time_of_impact` (which uses only conservative advancement), this
/// function performs a bisection step once an interval `[t_lo, t_hi]` bracketing
/// the TOI is found, giving sub-millisecond accuracy.
#[allow(dead_code)]
pub fn gjk_toi_bisection(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    velocity_a: Vec3,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    velocity_b: Vec3,
    max_toi: f64,
    tolerance: f64,
) -> Option<f64> {
    if Gjk::intersect(shape_a, transform_a, shape_b, transform_b) {
        return None;
    }
    let rel_vel = velocity_a - velocity_b;
    let rel_speed = rel_vel.norm();
    if rel_speed < 1e-14 {
        return None;
    }
    let mut t_lo = 0.0_f64;
    let mut t_hi = max_toi;
    let mut found = false;
    let mut t = 0.0_f64;
    for _ in 0..MAX_ITERATIONS {
        let ta = Transform {
            position: transform_a.position + velocity_a * t,
            rotation: transform_a.rotation,
        };
        let tb = Transform {
            position: transform_b.position + velocity_b * t,
            rotation: transform_b.rotation,
        };
        match Gjk::query(shape_a, &ta, shape_b, &tb) {
            GjkResult::Intersecting(_) => {
                t_hi = t;
                found = true;
                break;
            }
            GjkResult::Separated { distance, .. } => {
                if distance < 1e-5 {
                    t_hi = t;
                    found = true;
                    break;
                }
                t_lo = t;
                t += distance / rel_speed;
                if t > max_toi {
                    break;
                }
            }
        }
    }
    if !found {
        return None;
    }
    for _ in 0..32 {
        if t_hi - t_lo < tolerance {
            break;
        }
        let mid = (t_lo + t_hi) * 0.5;
        let ta = Transform {
            position: transform_a.position + velocity_a * mid,
            rotation: transform_a.rotation,
        };
        let tb = Transform {
            position: transform_b.position + velocity_b * mid,
            rotation: transform_b.rotation,
        };
        if Gjk::intersect(shape_a, &ta, shape_b, &tb) {
            t_hi = mid;
        } else {
            t_lo = mid;
        }
    }
    Some(t_hi)
}
/// Compute the time-of-impact for an angular sweep (rotation + translation).
///
/// Both bodies may rotate and translate.  The function samples `n_substeps`
/// intermediate states and returns the smallest TOI found.
#[allow(dead_code)]
pub fn gjk_toi_angular_sweep(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    velocity_a: Vec3,
    angular_a: Vec3,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    velocity_b: Vec3,
    angular_b: Vec3,
    max_toi: f64,
    n_substeps: usize,
) -> Option<f64> {
    let n = n_substeps.max(2);
    let dt = max_toi / (n as f64);
    let mut prev_ta = transform_a.clone();
    let mut prev_tb = transform_b.clone();
    for i in 1..=n {
        let t = i as f64 * dt;
        let ang_a_vec = angular_a * dt;
        let ang_a_len = ang_a_vec.norm();
        let rot_a_delta = if ang_a_len > 1e-14 {
            oxiphysics_core::math::quat_from_axis_angle(&(ang_a_vec / ang_a_len), ang_a_len)
        } else {
            oxiphysics_core::math::Quat::identity()
        };
        let ang_b_vec = angular_b * dt;
        let ang_b_len = ang_b_vec.norm();
        let rot_b_delta = if ang_b_len > 1e-14 {
            oxiphysics_core::math::quat_from_axis_angle(&(ang_b_vec / ang_b_len), ang_b_len)
        } else {
            oxiphysics_core::math::Quat::identity()
        };
        let next_ta = Transform {
            position: transform_a.position + velocity_a * t,
            rotation: rot_a_delta * transform_a.rotation,
        };
        let next_tb = Transform {
            position: transform_b.position + velocity_b * t,
            rotation: rot_b_delta * transform_b.rotation,
        };
        if Gjk::intersect(shape_a, &next_ta, shape_b, &next_tb) {
            let toi_refined = gjk_toi_bisection(
                shape_a, &prev_ta, velocity_a, shape_b, &prev_tb, velocity_b, dt, 1e-6,
            );
            return Some((i - 1) as f64 * dt + toi_refined.unwrap_or(dt));
        }
        prev_ta = next_ta;
        prev_tb = next_tb;
    }
    None
}
/// Perform a linear CCD sweep (conservative advancement) returning contact info.
///
/// This version computes a surface normal at the TOI in addition to the time value.
#[allow(dead_code)]
pub fn gjk_ccd_linear(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    velocity_a: Vec3,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    velocity_b: Vec3,
) -> Option<CcdResult> {
    let toi = gjk_toi_bisection(
        shape_a,
        transform_a,
        velocity_a,
        shape_b,
        transform_b,
        velocity_b,
        1.0,
        1e-6,
    )?;
    let ta_toi = Transform {
        position: transform_a.position + velocity_a * toi,
        rotation: transform_a.rotation,
    };
    let tb_toi = Transform {
        position: transform_b.position + velocity_b * toi,
        rotation: transform_b.rotation,
    };
    let rel_dir = (tb_toi.position - ta_toi.position).normalize();
    let sp = support(shape_a, &ta_toi, shape_b, &tb_toi, &rel_dir);
    let normal = if sp.point.norm() > 1e-10 {
        sp.point.normalize()
    } else {
        rel_dir
    };
    let point = (sp.support_a + sp.support_b) * 0.5;
    Some(CcdResult { toi, normal, point })
}
/// Estimate penetration depth and axis by sampling many support directions.
///
/// More accurate than `estimate_penetration_depth` because it samples along
/// 26 directions (faces + edges + vertices of an octahedron) rather than 6.
#[allow(dead_code)]
pub fn estimate_penetration_depth_26(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
) -> (f64, Vec3) {
    let inv_sqrt2 = 1.0_f64 / 2.0_f64.sqrt();
    let inv_sqrt3 = 1.0_f64 / 3.0_f64.sqrt();
    let dirs: [Vec3; 26] = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(inv_sqrt2, inv_sqrt2, 0.0),
        Vec3::new(-inv_sqrt2, inv_sqrt2, 0.0),
        Vec3::new(inv_sqrt2, -inv_sqrt2, 0.0),
        Vec3::new(-inv_sqrt2, -inv_sqrt2, 0.0),
        Vec3::new(inv_sqrt2, 0.0, inv_sqrt2),
        Vec3::new(-inv_sqrt2, 0.0, inv_sqrt2),
        Vec3::new(inv_sqrt2, 0.0, -inv_sqrt2),
        Vec3::new(-inv_sqrt2, 0.0, -inv_sqrt2),
        Vec3::new(0.0, inv_sqrt2, inv_sqrt2),
        Vec3::new(0.0, -inv_sqrt2, inv_sqrt2),
        Vec3::new(0.0, inv_sqrt2, -inv_sqrt2),
        Vec3::new(0.0, -inv_sqrt2, -inv_sqrt2),
        Vec3::new(inv_sqrt3, inv_sqrt3, inv_sqrt3),
        Vec3::new(-inv_sqrt3, inv_sqrt3, inv_sqrt3),
        Vec3::new(inv_sqrt3, -inv_sqrt3, inv_sqrt3),
        Vec3::new(-inv_sqrt3, -inv_sqrt3, inv_sqrt3),
        Vec3::new(inv_sqrt3, inv_sqrt3, -inv_sqrt3),
        Vec3::new(-inv_sqrt3, inv_sqrt3, -inv_sqrt3),
        Vec3::new(inv_sqrt3, -inv_sqrt3, -inv_sqrt3),
        Vec3::new(-inv_sqrt3, -inv_sqrt3, -inv_sqrt3),
    ];
    let mut min_depth = f64::MAX;
    let mut min_axis = Vec3::new(0.0, 1.0, 0.0);
    for dir in &dirs {
        let sp = support(shape_a, transform_a, shape_b, transform_b, dir);
        let proj = sp.point.dot(dir);
        if proj > 0.0 && proj < min_depth {
            min_depth = proj;
            min_axis = *dir;
        }
    }
    if min_depth == f64::MAX {
        (0.0, Vec3::new(0.0, 1.0, 0.0))
    } else {
        (min_depth, min_axis)
    }
}
/// Generate a speculative contact if two shapes are within `threshold` of each other.
#[allow(dead_code)]
pub fn gjk_speculative_contact(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    threshold: f64,
) -> Option<SpeculativeContact> {
    let result = gjk_distance_query(shape_a, transform_a, shape_b, transform_b);
    if result.distance > threshold {
        return None;
    }
    let diff = result.closest_b - result.closest_a;
    let dist = diff.norm();
    let normal = if dist > 1e-10 {
        diff / dist
    } else {
        let d = transform_b.position - transform_a.position;
        if d.norm_squared() > 1e-10 {
            d.normalize()
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        }
    };
    Some(SpeculativeContact {
        point_a: result.closest_a,
        point_b: result.closest_b,
        normal,
        gap: result.distance,
    })
}
/// Compute the GJK hull distance (minimum distance between convex hull vertices).
///
/// This is a simplified version that works directly on vertex arrays without
/// needing `Shape` trait objects.
#[allow(dead_code)]
pub fn gjk_hull_distance(verts_a: &[[f64; 3]], verts_b: &[[f64; 3]]) -> f64 {
    if verts_a.is_empty() || verts_b.is_empty() {
        return f64::MAX;
    }
    let support_arr = |verts: &[[f64; 3]], dir: &Vec3| -> Vec3 {
        let mut best = Vec3::new(verts[0][0], verts[0][1], verts[0][2]);
        let mut best_dot = best.dot(dir);
        for v in &verts[1..] {
            let vv = Vec3::new(v[0], v[1], v[2]);
            let d = vv.dot(dir);
            if d > best_dot {
                best_dot = d;
                best = vv;
            }
        }
        best
    };
    let mut dir = Vec3::new(1.0, 0.0, 0.0);
    let mut simplex: Vec<Vec3> = Vec::new();
    let first = support_arr(verts_a, &dir) - support_arr(verts_b, &(-dir));
    simplex.push(first);
    dir = -first;
    for _ in 0..MAX_ITERATIONS {
        if dir.norm_squared() < TOLERANCE * TOLERANCE {
            return 0.0;
        }
        let new_pt = support_arr(verts_a, &dir) - support_arr(verts_b, &(-dir));
        if new_pt.dot(&dir) < TOLERANCE {
            let (closest, _bary) = match simplex.len() {
                1 => (simplex[0], vec![1.0]),
                2 => {
                    let ab = simplex[1] - simplex[0];
                    let t = (-simplex[0]).dot(&ab) / ab.dot(&ab);
                    let t = t.clamp(0.0, 1.0);
                    (simplex[0] + ab * t, vec![1.0 - t, t])
                }
                _ => (dir, vec![1.0]),
            };
            return closest.norm();
        }
        simplex.push(new_pt);
        if simplex.len() > 3 {
            simplex = simplex[simplex.len() - 3..].to_vec();
        }
        let a = simplex[simplex.len() - 1];
        let ao = -a;
        if simplex.len() >= 2 {
            let b = simplex[simplex.len() - 2];
            let ab = b - a;
            if ab.dot(&ao) > 0.0 {
                dir = ab.cross(&ao).cross(&ab);
                if dir.norm_squared() < TOLERANCE {
                    return 0.0;
                }
            } else {
                simplex = vec![a];
                dir = ao;
            }
        } else {
            dir = ao;
        }
    }
    0.0
}
/// Test whether a point lies inside a convex shape.
#[allow(dead_code)]
pub fn gjk_point_inside_convex(point: Vec3, shape: &dyn Shape, transform: &Transform) -> bool {
    use oxiphysics_geometry::Sphere;
    let point_sphere = Sphere::new(1e-6);
    let point_transform = Transform {
        position: point,
        rotation: transform.rotation,
    };
    Gjk::intersect(&point_sphere, &point_transform, shape, transform)
}
/// Compute the distance from a point to the surface of a convex shape.
#[allow(dead_code)]
pub fn gjk_point_to_convex_distance(point: Vec3, shape: &dyn Shape, transform: &Transform) -> f64 {
    let point_sphere = Sphere::new(1e-6);
    let point_transform = Transform {
        position: point,
        rotation: transform.rotation,
    };
    gjk_distance(shape, transform, &point_sphere, &point_transform)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::narrowphase::ClosestFeature;
    use crate::narrowphase::GjkDebugTrace;
    use crate::narrowphase::VoronoiSimplexSolver;
    use crate::narrowphase::gjk::types::{GjkSolver, SupportCache};
    use oxiphysics_geometry::{BoxShape, Sphere};
    #[test]
    fn test_gjk_sphere_sphere_intersect() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        assert!(Gjk::intersect(&s1, &t1, &s2, &t2));
    }
    #[test]
    fn test_gjk_sphere_sphere_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(!Gjk::intersect(&s1, &t1, &s2, &t2));
    }
    #[test]
    fn test_gjk_box_box_intersect() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        assert!(Gjk::intersect(&b1, &t1, &b2, &t2));
    }
    #[test]
    fn test_gjk_box_box_separated() {
        let b1 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let b2 = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(!Gjk::intersect(&b1, &t1, &b2, &t2));
    }
    #[test]
    fn test_gjk_sphere_box_intersect() {
        let s = Sphere::new(1.0);
        let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        assert!(Gjk::intersect(&s, &t1, &b, &t2));
    }
    #[test]
    fn test_gjk_closest_sphere_sphere_touching() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(2.0, 0.0, 0.0));
        match Gjk::query(&s1, &t1, &s2, &t2) {
            GjkResult::Separated {
                point_a,
                point_b,
                distance,
            } => {
                assert!(point_a.x.is_finite() && point_b.x.is_finite());
                assert!(distance >= 0.0, "distance = {}", distance);
            }
            GjkResult::Intersecting(_) => {
                assert!(Gjk::closest_points(&s1, &t1, &s2, &t2).is_none());
            }
        }
    }
    #[test]
    fn test_gjk_closest_sphere_sphere_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        let result = Gjk::closest_points(&s1, &t1, &s2, &t2);
        assert!(
            result.is_some(),
            "separated spheres must yield closest points"
        );
        match Gjk::query(&s1, &t1, &s2, &t2) {
            GjkResult::Separated { distance, .. } => {
                assert!(distance > 0.0, "distance = {}", distance);
            }
            GjkResult::Intersecting(_) => {
                panic!("spheres 3 m apart should not intersect");
            }
        }
    }
    /// Single point: closest feature is the point itself.
    #[test]
    fn test_voronoi_single_point() {
        let mut solver = VoronoiSimplexSolver::new();
        solver.add_point(Vec3::new(1.0, 0.0, 0.0));
        let (pt, feature) = solver.solve();
        assert!((pt - Vec3::new(1.0, 0.0, 0.0)).norm() < 1e-12);
        assert_eq!(feature, ClosestFeature::Vertex(0));
    }
    /// Two points: closest feature should be on the edge or a vertex.
    #[test]
    fn test_voronoi_line_segment() {
        let mut solver = VoronoiSimplexSolver::new();
        solver.add_point(Vec3::new(-1.0, 1.0, 0.0));
        solver.add_point(Vec3::new(1.0, 1.0, 0.0));
        let (pt, feature) = solver.solve();
        assert!(
            (pt - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-10,
            "closest = {:?}",
            pt
        );
        assert!(
            matches!(feature, ClosestFeature::Edge(_, _)),
            "expected edge feature, got {feature:?}"
        );
    }
    /// Closest point to origin on a triangle.
    #[test]
    fn test_voronoi_triangle() {
        let mut solver = VoronoiSimplexSolver::new();
        solver.add_point(Vec3::new(0.0, 2.0, 0.0));
        solver.add_point(Vec3::new(-1.0, 1.0, 0.0));
        solver.add_point(Vec3::new(1.0, 1.0, 0.0));
        let (pt, _feature) = solver.solve();
        assert!(
            pt.norm() < 2.0,
            "closest point should be reasonably close to origin, got {:?}",
            pt
        );
    }
    /// Cache should store and retrieve support points.
    #[test]
    fn test_support_cache_basic() {
        let mut cache = SupportCache::new(10);
        assert!(cache.is_empty());
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let sp = SupportPoint {
            point: Vec3::new(1.0, 0.0, 0.0),
            support_a: Vec3::new(1.0, 0.0, 0.0),
            support_b: Vec3::zeros(),
        };
        cache.insert(dir, sp);
        assert_eq!(cache.len(), 1);
        let result = cache.lookup(&Vec3::new(1.0, 0.0, 0.0));
        assert!(result.is_some());
        let result = cache.lookup(&Vec3::new(0.0, 1.0, 0.0));
        assert!(result.is_none());
    }
    /// Cache should evict oldest entry when full.
    #[test]
    fn test_support_cache_eviction() {
        let mut cache = SupportCache::new(2);
        let sp = SupportPoint {
            point: Vec3::zeros(),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        };
        cache.insert(Vec3::new(1.0, 0.0, 0.0), sp);
        cache.insert(Vec3::new(0.0, 1.0, 0.0), sp);
        cache.insert(Vec3::new(0.0, 0.0, 1.0), sp);
        assert_eq!(cache.len(), 2);
        assert!(cache.lookup(&Vec3::new(1.0, 0.0, 0.0)).is_none());
    }
    /// Cache clear should empty the cache.
    #[test]
    fn test_support_cache_clear() {
        let mut cache = SupportCache::new(10);
        let sp = SupportPoint {
            point: Vec3::zeros(),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        };
        cache.insert(Vec3::new(1.0, 0.0, 0.0), sp);
        cache.clear();
        assert!(cache.is_empty());
    }
    /// A single-point simplex is not degenerate.
    #[test]
    fn test_non_degenerate_single_point() {
        let mut s = Simplex::new();
        s.push(SupportPoint {
            point: Vec3::new(1.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        });
        assert!(!is_simplex_degenerate(&s));
    }
    /// Two identical points form a degenerate simplex.
    #[test]
    fn test_degenerate_duplicate_points() {
        let mut s = Simplex::new();
        let p = SupportPoint {
            point: Vec3::new(1.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        };
        s.push(p);
        s.push(p);
        assert!(is_simplex_degenerate(&s));
    }
    /// Collinear triangle is degenerate.
    #[test]
    fn test_degenerate_collinear_triangle() {
        let mut s = Simplex::new();
        s.push(SupportPoint {
            point: Vec3::new(0.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        });
        s.push(SupportPoint {
            point: Vec3::new(1.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        });
        s.push(SupportPoint {
            point: Vec3::new(2.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        });
        assert!(is_simplex_degenerate(&s));
    }
    /// Fix degenerate simplex removes duplicate points.
    #[test]
    fn test_fix_degenerate_removes_duplicates() {
        let mut s = Simplex::new();
        let p = SupportPoint {
            point: Vec3::new(1.0, 0.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        };
        s.push(p);
        s.push(p);
        s.push(SupportPoint {
            point: Vec3::new(0.0, 1.0, 0.0),
            support_a: Vec3::zeros(),
            support_b: Vec3::zeros(),
        });
        let modified = fix_degenerate_simplex(&mut s);
        assert!(modified);
        assert_eq!(s.len(), 2);
    }
    /// GJK distance for separated spheres should be positive.
    #[test]
    fn test_gjk_distance_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let d = gjk_distance(&s1, &t1, &s2, &t2);
        assert!(d > 0.0, "distance should be positive, got {d}");
    }
    /// GJK distance for overlapping spheres should be negative.
    #[test]
    fn test_gjk_distance_overlapping() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let d = gjk_distance(&s1, &t1, &s2, &t2);
        assert!(d < 0.0, "distance should be negative for overlap, got {d}");
    }
    /// Minkowski support should return a valid support point.
    #[test]
    fn test_minkowski_support() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let sp = minkowski_support(&s, &t, &s, &t, &dir);
        assert!(sp.point.x.is_finite());
    }
    /// Simplex Default trait should create an empty simplex.
    #[test]
    fn test_simplex_default() {
        let s = Simplex::default();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }
    /// Closest point on a line segment.
    #[test]
    fn test_closest_point_line_midpoint() {
        let a = Vec3::new(-1.0, 1.0, 0.0);
        let b = Vec3::new(1.0, 1.0, 0.0);
        let (pt, bary) = closest_point_line(&a, &b);
        assert!((pt.x).abs() < 1e-10, "x = {}", pt.x);
        assert!((pt.y - 1.0).abs() < 1e-10, "y = {}", pt.y);
        assert!((bary[0] - 0.5).abs() < 1e-10);
        assert!((bary[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_scaled_support_sphere() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let sp = scaled_support(&s, &t, &dir, 2.0);
        assert!(
            (sp.support_a.x - 2.0).abs() < 1e-10,
            "scaled support x = {}",
            sp.support_a.x
        );
    }
    #[test]
    fn test_scaled_support_different_scales() {
        let s = Sphere::new(1.0);
        let t = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let sp1 = scaled_support(&s, &t, &dir, 1.0);
        let sp3 = scaled_support(&s, &t, &dir, 3.0);
        assert!(
            sp3.support_a.y > sp1.support_a.y,
            "larger scale gives larger support"
        );
    }
    #[test]
    fn test_gjk_with_scale_intersect() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(3.0, 0.0, 0.0));
        assert!(Gjk::intersect_scaled(&s1, &t1, 2.0, &s2, &t2, 2.0));
    }
    #[test]
    fn test_gjk_with_scale_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        assert!(!Gjk::intersect_scaled(&s1, &t1, 1.0, &s2, &t2, 1.0));
    }
    #[test]
    fn test_gjk_toi_approaching() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let vel1 = Vec3::new(1.0, 0.0, 0.0);
        let vel2 = Vec3::new(0.0, 0.0, 0.0);
        let toi = gjk_time_of_impact(&s1, &t1, vel1, &s2, &t2, vel2, 10.0);
        assert!(toi.is_some(), "should find toi for approaching objects");
        let t = toi.unwrap();
        assert!(t > 0.0 && t < 10.0, "toi={t} should be in (0,10)");
    }
    #[test]
    fn test_gjk_toi_receding() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let vel1 = Vec3::new(-1.0, 0.0, 0.0);
        let vel2 = Vec3::new(1.0, 0.0, 0.0);
        let toi = gjk_time_of_impact(&s1, &t1, vel1, &s2, &t2, vel2, 10.0);
        assert!(toi.is_none(), "receding objects should not have toi");
    }
    #[test]
    fn test_gjk_toi_already_overlapping() {
        let s1 = Sphere::new(2.0);
        let s2 = Sphere::new(2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let vel1 = Vec3::new(0.0, 0.0, 0.0);
        let vel2 = Vec3::new(0.0, 0.0, 0.0);
        let toi = gjk_time_of_impact(&s1, &t1, vel1, &s2, &t2, vel2, 10.0);
        assert!(toi.is_none(), "already overlapping → toi = None (use EPA)");
    }
    #[test]
    fn test_gjk_debug_trace_separated() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let trace = GjkDebugTrace::collect(&s1, &t1, &s2, &t2);
        assert!(
            trace.n_iterations() > 0,
            "should have at least one iteration"
        );
        assert!(
            !trace.terminated_as_intersecting(),
            "separated → no intersection"
        );
    }
    #[test]
    fn test_gjk_debug_trace_intersecting() {
        let s1 = Sphere::new(2.0);
        let s2 = Sphere::new(2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let trace = GjkDebugTrace::collect(&s1, &t1, &s2, &t2);
        assert!(
            trace.terminated_as_intersecting(),
            "overlapping → intersection"
        );
    }
    #[test]
    fn test_gjk_debug_directions_nonempty() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let trace = GjkDebugTrace::collect(&s1, &t1, &s2, &t2);
        assert!(!trace.directions().is_empty());
    }
    #[test]
    fn test_penetration_direction_refinement_basic() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        if let GjkResult::Intersecting(simplex) = result {
            let dir = refine_penetration_direction(&simplex);
            assert!(
                dir.norm() > 1e-10,
                "penetration direction should not be zero"
            );
        }
    }
    #[test]
    fn test_penetration_direction_axis_aligned() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        if let GjkResult::Intersecting(simplex) = result {
            let dir = refine_penetration_direction(&simplex);
            assert!(
                dir.norm() > 1e-10,
                "penetration direction should not be zero"
            );
        }
    }
    #[test]
    fn test_penetration_estimate_depth() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let depth = estimate_penetration_depth(&s1, &t1, &s2, &t2);
        assert!(
            depth > 0.0,
            "overlapping spheres should have positive depth"
        );
    }
    #[test]
    fn test_penetration_estimate_depth_no_overlap() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let depth = estimate_penetration_depth(&s1, &t1, &s2, &t2);
        assert!(
            depth <= 0.0,
            "separated shapes should have non-positive depth"
        );
    }
    #[test]
    fn test_gjk_solver_signed_distance_separated_positive() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let d = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        assert!(
            d > 0.0,
            "separated shapes should have positive signed distance, got {d}"
        );
    }
    #[test]
    fn test_gjk_solver_signed_distance_overlapping_negative() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let d = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        assert!(
            d <= 0.0,
            "overlapping shapes should have non-positive signed distance, got {d}"
        );
    }
    #[test]
    fn test_gjk_solver_signed_distance_stores_simplex_on_overlap() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let d = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        assert!(d <= 0.0);
        assert!(
            solver.has_cached_simplex(),
            "solver should cache the simplex on overlap"
        );
    }
    #[test]
    fn test_gjk_solver_signed_distance_clears_simplex_on_separation() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        let _d = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        assert!(
            !solver.has_cached_simplex(),
            "no cached simplex for separated case"
        );
    }
    #[test]
    fn test_gjk_solver_witness_points_separated() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        let result = solver.compute_witness_points(&s1, &t1, &s2, &t2);
        assert!(
            result.is_some(),
            "separated spheres must yield witness points"
        );
        let (pa, pb) = result.unwrap();
        assert!(
            pa.x.is_finite() && pb.x.is_finite(),
            "witness points must be finite"
        );
    }
    #[test]
    fn test_gjk_solver_witness_points_intersecting_returns_none() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(2.0);
        let s2 = Sphere::new(2.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = solver.compute_witness_points(&s1, &t1, &s2, &t2);
        assert!(
            result.is_none(),
            "overlapping shapes must return None for witness points"
        );
    }
    #[test]
    fn test_gjk_solver_witness_points_sensible_distance() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(4.0, 0.0, 0.0));
        if let Some((pa, pb)) = solver.compute_witness_points(&s1, &t1, &s2, &t2) {
            let gap = (pb - pa).norm();
            assert!(gap > 0.0 && gap < 10.0, "gap = {gap}");
        }
    }
    #[test]
    fn test_gjk_solver_warm_start_stores_direction() {
        let mut solver = GjkSolver::new();
        assert!(
            solver.warm_direction().is_none(),
            "no warm direction initially"
        );
        solver.warm_start(Vec3::new(1.0, 0.0, 0.0));
        let d = solver
            .warm_direction()
            .expect("warm direction should be stored");
        assert!(
            (d.norm() - 1.0).abs() < 1e-10,
            "warm direction should be unit, norm={}",
            d.norm()
        );
    }
    #[test]
    fn test_gjk_solver_warm_start_zero_ignored() {
        let mut solver = GjkSolver::new();
        solver.warm_start(Vec3::new(0.0, 0.0, 0.0));
        assert!(
            solver.warm_direction().is_none(),
            "zero direction should not be stored"
        );
    }
    #[test]
    fn test_gjk_solver_warm_start_normalises_input() {
        let mut solver = GjkSolver::new();
        solver.warm_start(Vec3::new(3.0, 4.0, 0.0));
        let d = solver.warm_direction().unwrap();
        assert!(
            (d.norm() - 1.0).abs() < 1e-10,
            "stored direction must be unit, norm={}",
            d.norm()
        );
        assert!(
            (d.x - 0.6).abs() < 1e-10 && (d.y - 0.8).abs() < 1e-10,
            "direction not normalised correctly: {:?}",
            d
        );
    }
    #[test]
    fn test_gjk_solver_reset_clears_state() {
        let mut solver = GjkSolver::new();
        solver.warm_start(Vec3::new(1.0, 0.0, 0.0));
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let _ = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        solver.reset();
        assert!(
            !solver.has_cached_simplex(),
            "reset should clear simplex cache"
        );
        assert!(
            solver.warm_direction().is_none(),
            "reset should clear warm direction"
        );
    }
    #[test]
    fn test_gjk_solver_take_simplex_consumes_cache() {
        let mut solver = GjkSolver::new();
        let s1 = Sphere::new(1.5);
        let s2 = Sphere::new(1.5);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(0.5, 0.0, 0.0));
        let _ = solver.compute_signed_distance(&s1, &t1, &s2, &t2);
        if solver.has_cached_simplex() {
            let s = solver.take_simplex();
            assert!(s.is_some(), "take_simplex should return the cached simplex");
            assert!(
                !solver.has_cached_simplex(),
                "cache should be empty after take"
            );
        }
    }
    #[test]
    fn test_gjk_solver_default_equivalent_to_new() {
        let solver: GjkSolver = GjkSolver::default();
        assert!(!solver.has_cached_simplex());
        assert!(solver.warm_direction().is_none());
    }
}
