//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    GjkResult, RaySphereHit, RestitutionModel, RotatingSweptBody, SweptBody, SweptCapsule,
    TranslatedSupport,
};

#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n > 1e-12 {
        scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}
/// Iterative conservative advancement for two moving spheres.
///
/// Returns an approximate TOI in `[0, 1]` (normalised), or `None` if the
/// spheres never come within the sum of their radii.
pub fn conservative_advancement_rigid(
    pos_a: [f64; 3],
    vel_a: [f64; 3],
    radius_a: f64,
    pos_b: [f64; 3],
    vel_b: [f64; 3],
    radius_b: f64,
    max_iter: usize,
) -> Option<f64> {
    let sum_r = radius_a + radius_b;
    let rel_vel = sub(vel_a, vel_b);
    let max_speed = norm(rel_vel);
    let mut t = 0.0_f64;
    for _ in 0..max_iter {
        let pa = add(pos_a, scale(vel_a, t));
        let pb = add(pos_b, scale(vel_b, t));
        let d = norm(sub(pa, pb));
        let gap = d - sum_r;
        if gap <= 1e-6 {
            return Some(t);
        }
        if max_speed < 1e-12 {
            return None;
        }
        let step = gap / max_speed;
        t += step;
        if t > 1.0 {
            return None;
        }
    }
    None
}
/// Compute the swept AABB that bounds a box moving linearly over `[0, dt]`.
///
/// Returns `(min, max)` corners of the bounding AABB.
pub fn linear_sweep_aabb(
    pos: [f64; 3],
    vel: [f64; 3],
    half_extents: [f64; 3],
    dt: f64,
) -> ([f64; 3], [f64; 3]) {
    let end_pos = add(pos, scale(vel, dt));
    let mut mn = [0.0_f64; 3];
    let mut mx = [0.0_f64; 3];
    for i in 0..3 {
        let lo_start = pos[i] - half_extents[i];
        let hi_start = pos[i] + half_extents[i];
        let lo_end = end_pos[i] - half_extents[i];
        let hi_end = end_pos[i] + half_extents[i];
        mn[i] = lo_start.min(lo_end);
        mx[i] = hi_start.max(hi_end);
    }
    (mn, mx)
}
/// Return the time of impact given an initial separation and closing speed.
///
/// Returns `Some(separation / closing_speed)` when `closing_speed > 0` and
/// the bodies are currently separated, `Some(0.0)` when already overlapping,
/// and `None` when diverging.
pub fn time_of_impact_linear(separation: f64, closing_speed: f64) -> Option<f64> {
    if separation <= 0.0 {
        return Some(0.0);
    }
    if closing_speed <= 0.0 {
        return None;
    }
    Some(separation / closing_speed)
}
/// Speculative contact: expand the AABB of each body by its velocity * dt,
/// then detect potential contacts before they happen.
///
/// Returns `true` if the swept sphere of `a` (expanded by `dt * |vel_a|`) and
/// the swept sphere of `b` might contact.
pub fn speculative_contact_candidate(a: &SweptBody, b: &SweptBody, dt: f64) -> bool {
    let disp_a = norm(a.vel) * dt;
    let disp_b = norm(b.vel) * dt;
    let combined_reach = a.radius + b.radius + disp_a + disp_b;
    let dist = norm(sub(a.pos, b.pos));
    dist <= combined_reach
}
/// Refine a time-of-impact estimate by bisecting the interval `[t_lo, t_hi]`.
///
/// At each step, the two spheres are evaluated at the midpoint; the sub-interval
/// where they approach is kept.  Returns a refined TOI in `[t_lo, t_hi]` or
/// `None` if the spheres are separated throughout.
pub fn refine_toi_bisection(
    a: &SweptBody,
    b: &SweptBody,
    t_lo: f64,
    t_hi: f64,
    tolerance: f64,
    max_iter: usize,
) -> Option<f64> {
    let sum_r = a.radius + b.radius;
    let dist_lo = norm(sub(a.position_at(t_lo), b.position_at(t_lo)));
    let dist_hi = norm(sub(a.position_at(t_hi), b.position_at(t_hi)));
    if dist_lo <= sum_r {
        return Some(t_lo);
    }
    if dist_hi > sum_r {
        return None;
    }
    let mut lo = t_lo;
    let mut hi = t_hi;
    for _ in 0..max_iter {
        if hi - lo < tolerance {
            break;
        }
        let mid = (lo + hi) * 0.5;
        let d = norm(sub(a.position_at(mid), b.position_at(mid)));
        if d <= sum_r {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some((lo + hi) * 0.5)
}
/// CCD broad-phase candidate test for rotating bodies.
///
/// Returns `true` if the swept bounding volumes of the two bodies might overlap
/// within `[0, dt]`.
pub fn rotating_body_ccd_candidate(a: &RotatingSweptBody, b: &RotatingSweptBody, dt: f64) -> bool {
    let r_a = a.conservative_radius(dt);
    let r_b = b.conservative_radius(dt);
    let mid_a = add(a.pos, scale(a.vel, dt * 0.5));
    let mid_b = add(b.pos, scale(b.vel, dt * 0.5));
    let half_disp_a = norm(a.vel) * dt * 0.5;
    let half_disp_b = norm(b.vel) * dt * 0.5;
    let dist = norm(sub(mid_a, mid_b));
    dist <= r_a + r_b + half_disp_a + half_disp_b
}
#[cfg(test)]
#[inline]
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Closest-point-on-segment computation.
///
/// Given a line segment from `a` to `b` and a point `p`, returns the
/// parameter `t ∈ [0, 1]` such that `a + t*(b-a)` is closest to `p`, and
/// the closest point itself.
pub fn closest_point_on_segment(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> ([f64; 3], f64) {
    let ab = sub(b, a);
    let len_sq = dot(ab, ab);
    if len_sq < 1e-24 {
        return (a, 0.0);
    }
    let t = (dot(sub(p, a), ab) / len_sq).clamp(0.0, 1.0);
    (add(a, scale(ab, t)), t)
}
/// Squared distance between the closest points of two line segments.
///
/// Uses the Shoemake / Ericson algorithm.
///
/// Returns `(sq_dist, s, t)` where `s` and `t` are parameters on the first
/// and second segments respectively.
pub fn segment_segment_sq_dist(
    p0: [f64; 3],
    p1: [f64; 3],
    q0: [f64; 3],
    q1: [f64; 3],
) -> (f64, f64, f64) {
    let d1 = sub(p1, p0);
    let d2 = sub(q1, q0);
    let r = sub(p0, q0);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (mut s, mut t);
    pub(super) const EPS: f64 = 1e-10;
    if a <= EPS && e <= EPS {
        return (dot(r, r), 0.0, 0.0);
    }
    if a <= EPS {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= EPS {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            if denom.abs() > EPS {
                s = ((b * f - c * e) / denom).clamp(0.0, 1.0);
            } else {
                s = 0.0;
            }
            t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
        }
    }
    let closest_p = add(p0, scale(d1, s));
    let closest_q = add(q0, scale(d2, t));
    let diff = sub(closest_p, closest_q);
    (dot(diff, diff), s, t)
}
/// Capsule–capsule CCD test.
///
/// Returns the time of impact `t ∈ [0, dt]` when the two capsules first
/// touch, or `None` if no contact occurs within `[0, dt]`.
///
/// Uses a conservative advancement loop based on the closest-point-on-axis
/// geometry.
pub fn check_capsule_capsule(
    a: &SweptCapsule,
    b: &SweptCapsule,
    dt: f64,
    max_iter: usize,
) -> Option<f64> {
    let sum_r = a.radius + b.radius;
    let sum_r2 = sum_r * sum_r;
    let mut t = 0.0f64;
    for _ in 0..max_iter {
        let ca = a.advance(t);
        let cb = b.advance(t);
        let (sq_d, _, _) = segment_segment_sq_dist(ca.p0, ca.p1, cb.p0, cb.p1);
        if sq_d <= sum_r2 {
            return Some(t);
        }
        let gap = sq_d.sqrt() - sum_r;
        let rel_speed = norm(sub(a.vel, b.vel));
        if rel_speed < 1e-12 {
            return None;
        }
        let step = gap / rel_speed;
        t += step;
        if t > dt {
            return None;
        }
    }
    None
}
/// Cast a ray `origin + t * dir` against a sphere.
///
/// Returns the closest intersection with `t >= t_min` and `t <= t_max`, or
/// `None` if there is no such intersection.
pub fn ray_cast_sphere(
    origin: [f64; 3],
    dir: [f64; 3],
    sphere_center: [f64; 3],
    sphere_radius: f64,
    t_min: f64,
    t_max: f64,
) -> Option<RaySphereHit> {
    let oc = sub(origin, sphere_center);
    let a = dot(dir, dir);
    let half_b = dot(oc, dir);
    let c = dot(oc, oc) - sphere_radius * sphere_radius;
    let disc = half_b * half_b - a * c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t = (-half_b - sqrt_disc) / a;
    let t = if t >= t_min && t <= t_max {
        t
    } else {
        let t2 = (-half_b + sqrt_disc) / a;
        if t2 >= t_min && t2 <= t_max {
            t2
        } else {
            return None;
        }
    };
    let point = add(origin, scale(dir, t));
    let normal = normalize(sub(point, sphere_center));
    Some(RaySphereHit { t, point, normal })
}
/// Cast a ray against an axis-aligned bounding box using the slab method.
///
/// Returns `Some((t_enter, t_exit))` when the ray hits the box within
/// `[t_min, t_max]`, or `None` if it misses.
pub fn ray_cast_aabb(
    origin: [f64; 3],
    dir: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    t_min: f64,
    t_max: f64,
) -> Option<(f64, f64)> {
    let mut tmin = t_min;
    let mut tmax = t_max;
    for i in 0..3 {
        if dir[i].abs() < 1e-14 {
            if origin[i] < aabb_min[i] || origin[i] > aabb_max[i] {
                return None;
            }
        } else {
            let inv = 1.0 / dir[i];
            let t0 = (aabb_min[i] - origin[i]) * inv;
            let t1 = (aabb_max[i] - origin[i]) * inv;
            let (t0, t1) = if inv >= 0.0 { (t0, t1) } else { (t1, t0) };
            tmin = tmin.max(t0);
            tmax = tmax.min(t1);
            if tmin > tmax {
                return None;
            }
        }
    }
    Some((tmin, tmax))
}
/// Cast a ray against an infinite plane `n·x = d`.
///
/// Returns `Some(t)` such that `origin + t * dir` lies on the plane,
/// with `t ∈ [t_min, t_max]`, or `None` if the ray is parallel or misses.
pub fn ray_cast_plane(
    origin: [f64; 3],
    dir: [f64; 3],
    plane_normal: [f64; 3],
    plane_d: f64,
    t_min: f64,
    t_max: f64,
) -> Option<f64> {
    let n = normalize(plane_normal);
    let denom = dot(n, dir);
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = (plane_d - dot(n, origin)) / denom;
    if t >= t_min && t <= t_max {
        Some(t)
    } else {
        None
    }
}
/// A convex shape support function abstraction for GJK.
///
/// Each shape provides a `support(d)` method that returns the point on its
/// surface furthest in direction `d`.
pub trait ConvexSupport {
    /// Return the support point in direction `d`.
    fn support(&self, d: [f64; 3]) -> [f64; 3];
}
/// Minkowski-difference support: for two convex shapes A and B,
/// `support_A(d) - support_B(-d)`.
pub fn minkowski_difference_support<A: ConvexSupport, B: ConvexSupport>(
    a: &A,
    b: &B,
    d: [f64; 3],
) -> [f64; 3] {
    let neg_d = scale(d, -1.0);
    sub(a.support(d), b.support(neg_d))
}
/// Run the GJK distance algorithm between two convex shapes.
///
/// This is a simplified iterative GJK that converges when the improvement
/// between iterations falls below `tolerance`.  Returns the squared
/// minimum distance and closest points.
///
/// # Complexity
/// O(`max_iter`) iterations, each O(1) for sphere/AABB support.
pub fn gjk_distance<A: ConvexSupport, B: ConvexSupport>(
    a: &A,
    b: &B,
    tolerance: f64,
    max_iter: usize,
) -> GjkResult {
    let mut dir = [1.0, 0.0, 0.0_f64];
    let mut closest_a = a.support(dir);
    let mut closest_b = b.support(scale(dir, -1.0));
    let mut sq_dist = f64::MAX;
    for _ in 0..max_iter {
        let w = sub(closest_a, closest_b);
        let new_sq_dist = dot(w, w);
        if new_sq_dist < 1e-28 {
            return GjkResult {
                sq_dist: 0.0,
                closest_a,
                closest_b,
                overlapping: true,
            };
        }
        if (sq_dist - new_sq_dist) / sq_dist.max(1.0) < tolerance {
            break;
        }
        sq_dist = new_sq_dist;
        dir = scale(w, -1.0);
        let nd = normalize(dir);
        closest_a = a.support(nd);
        closest_b = b.support(scale(nd, -1.0));
    }
    let w = sub(closest_a, closest_b);
    let final_sq = dot(w, w);
    GjkResult {
        sq_dist: final_sq,
        closest_a,
        closest_b,
        overlapping: final_sq < tolerance * tolerance,
    }
}
/// Conservative-advancement TOI solver using GJK distance queries.
///
/// Advances from `t = 0` towards `t = dt` by computing the GJK distance at
/// each step and estimating a safe advance as `dist / max_relative_speed`.
///
/// Returns the time `t ∈ [0, dt]` when the two shapes first touch, or `None`
/// if they remain separated throughout.
pub fn gjk_conservative_advancement<A: ConvexSupport + Clone, B: ConvexSupport + Clone>(
    a_shape: &A,
    a_pos: [f64; 3],
    a_vel: [f64; 3],
    b_shape: &B,
    b_pos: [f64; 3],
    b_vel: [f64; 3],
    dt: f64,
    max_iter: usize,
    tolerance: f64,
) -> Option<f64> {
    let rel_vel = sub(a_vel, b_vel);
    let max_speed = norm(rel_vel);
    let mut t = 0.0f64;
    for _ in 0..max_iter {
        let a_at_t = TranslatedSupport {
            shape: a_shape,
            offset: add(a_pos, scale(a_vel, t)),
        };
        let b_at_t = TranslatedSupport {
            shape: b_shape,
            offset: add(b_pos, scale(b_vel, t)),
        };
        let result = gjk_distance(&a_at_t, &b_at_t, tolerance, 32);
        if result.overlapping || result.sq_dist < tolerance * tolerance {
            return Some(t);
        }
        if max_speed < 1e-12 {
            return None;
        }
        let gap = result.sq_dist.sqrt();
        let step = gap / max_speed;
        t += step;
        if t > dt {
            return None;
        }
    }
    None
}
/// Compute a spring–damper penalty force for an overlapping pair of spheres.
///
/// Uses the Hertz-inspired formula:
///   `F = k * overlap + c * v_close`
///
/// where `overlap = sum_radii - distance` (positive when penetrating),
/// `v_close` is the closing speed along the normal, and the force is directed
/// along the contact normal.
///
/// # Returns
/// The force vector to apply to body A (equal and opposite for B).
pub fn penalty_force_sphere_sphere(
    pos_a: [f64; 3],
    vel_a: [f64; 3],
    pos_b: [f64; 3],
    vel_b: [f64; 3],
    radius_a: f64,
    radius_b: f64,
    stiffness: f64,
    damping: f64,
) -> [f64; 3] {
    let diff = sub(pos_a, pos_b);
    let dist = norm(diff);
    let sum_r = radius_a + radius_b;
    if dist >= sum_r || dist < 1e-12 {
        return [0.0, 0.0, 0.0];
    }
    let normal = if dist > 1e-12 {
        scale(diff, 1.0 / dist)
    } else {
        [1.0, 0.0, 0.0]
    };
    let overlap = sum_r - dist;
    let rel_vel = sub(vel_a, vel_b);
    let v_close = dot(rel_vel, normal);
    let magnitude = stiffness * overlap - damping * v_close;
    scale(normal, magnitude.max(0.0))
}
/// Compute the impulse magnitude for a sphere–sphere collision using a given
/// restitution model.
///
/// Returns the scalar impulse `j` (applied along the contact normal):
/// `j = -(1 + e) * v_n / (inv_a + inv_b)`.
pub fn compute_collision_impulse(
    rel_vel: [f64; 3],
    contact_normal: [f64; 3],
    inv_mass_a: f64,
    inv_mass_b: f64,
    model: &RestitutionModel,
) -> f64 {
    let v_n = dot(rel_vel, contact_normal);
    if v_n >= 0.0 {
        return 0.0;
    }
    let impact_speed = (-v_n).abs();
    let e = model.evaluate(impact_speed);
    let denom = inv_mass_a + inv_mass_b;
    if denom < 1e-30 {
        return 0.0;
    }
    -(1.0 + e) * v_n / denom
}
/// Apply a Baumgarte stabilisation position correction to two bodies.
///
/// `correction = beta * overlap / dt` is added to the relative velocity
/// along the normal so that penetration is corrected over a few steps.
///
/// Returns `(delta_vel_a, delta_vel_b)` velocity corrections.
pub fn baumgarte_correction(
    normal: [f64; 3],
    overlap: f64,
    inv_mass_a: f64,
    inv_mass_b: f64,
    beta: f64,
    dt: f64,
) -> ([f64; 3], [f64; 3]) {
    if dt < 1e-12 {
        return ([0.0; 3], [0.0; 3]);
    }
    let correction_speed = beta * overlap / dt;
    let denom = inv_mass_a + inv_mass_b;
    if denom < 1e-30 {
        return ([0.0; 3], [0.0; 3]);
    }
    let j = correction_speed / denom;
    let delta_a = scale(normal, j * inv_mass_a);
    let delta_b = scale(normal, -j * inv_mass_b);
    (delta_a, delta_b)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ccd::AabbSupport;
    use crate::ccd::CcdBodySlot;
    use crate::ccd::CcdBroadphase;
    use crate::ccd::CcdEvent;
    use crate::ccd::CcdFilter;
    use crate::ccd::CcdPipeline;
    use crate::ccd::CcdStats;
    use crate::ccd::CcdSweep;
    use crate::ccd::SphereSupport;
    use crate::ccd::SubstepIntegrator;
    use crate::ccd::TunnelTest;
    fn make_body(pos: [f64; 3], vel: [f64; 3], radius: f64, inv_mass: f64) -> SweptBody {
        SweptBody::new(pos, vel, radius, inv_mass)
    }
    #[test]
    fn swept_body_position_at() {
        let b = make_body([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let p = b.position_at(2.0);
        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
    }
    #[test]
    fn swept_body_position_at_zero() {
        let b = make_body([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], 1.0, 1.0);
        assert_eq!(b.position_at(0.0), [1.0, 2.0, 3.0]);
    }
    #[test]
    fn swept_body_clone_eq() {
        let a = make_body([0.0; 3], [0.0; 3], 1.0, 1.0);
        let b = a.clone();
        assert_eq!(a, b);
    }
    #[test]
    fn sphere_sphere_head_on_collision() {
        let a = make_body([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_sphere(&a, &b, 10.0);
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!(t > 0.0 && t <= 10.0);
        let pa = a.position_at(t);
        let pb = b.position_at(t);
        let dist = norm(sub(pa, pb));
        assert!((dist - 1.0).abs() < 1e-6, "dist={dist}");
    }
    #[test]
    fn sphere_sphere_no_collision_diverging() {
        let a = make_body([-5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_sphere(&a, &b, 10.0);
        assert!(toi.is_none());
    }
    #[test]
    fn sphere_sphere_already_overlapping() {
        let a = make_body([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1.0);
        let b = make_body([0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1.0);
        let toi = TunnelTest::check_sphere_sphere(&a, &b, 1.0);
        assert_eq!(toi, Some(0.0));
    }
    #[test]
    fn sphere_sphere_miss_different_axis() {
        let a = make_body([-5.0, 3.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, -3.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_sphere(&a, &b, 10.0);
        assert!(toi.is_none());
    }
    #[test]
    fn sphere_sphere_collision_outside_dt() {
        let a = make_body([-5.0, 0.0, 0.0], [0.5, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, 0.0, 0.0], [-0.5, 0.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_sphere(&a, &b, 5.0);
        assert!(toi.is_none());
    }
    #[test]
    fn sphere_plane_approaching() {
        let sphere = make_body([0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_plane(&sphere, [0.0, 1.0, 0.0], 0.0, 10.0);
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((t - 4.5).abs() < 1e-6, "t={t}");
    }
    #[test]
    fn sphere_plane_moving_away() {
        let sphere = make_body([0.0, 5.0, 0.0], [0.0, 1.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_plane(&sphere, [0.0, 1.0, 0.0], 0.0, 10.0);
        assert!(toi.is_none());
    }
    #[test]
    fn sphere_plane_already_touching() {
        let sphere = make_body([0.0, 0.5, 0.0], [0.0, -1.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_plane(&sphere, [0.0, 1.0, 0.0], 0.0, 10.0);
        assert_eq!(toi, Some(0.0));
    }
    #[test]
    fn sphere_plane_outside_dt() {
        let sphere = make_body([0.0, 100.0, 0.0], [0.0, -1.0, 0.0], 0.5, 1.0);
        let toi = TunnelTest::check_sphere_plane(&sphere, [0.0, 1.0, 0.0], 0.0, 1.0);
        assert!(toi.is_none());
    }
    #[test]
    fn pipeline_detects_collision() {
        let bodies = vec![
            make_body([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
            make_body([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
        ];
        let ids = vec![0usize, 1usize];
        let pipeline = CcdPipeline::new(10.0);
        let events = pipeline.detect_tunneling(&bodies, &ids);
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].body_a, 0);
        assert_eq!(events[0].body_b, 1);
    }
    #[test]
    fn pipeline_sorted_by_toi() {
        let bodies = vec![
            make_body([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0),
            make_body([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.1, 1.0),
            make_body([-10.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0),
            make_body([10.0, 5.0, 0.0], [-1.0, 0.0, 0.0], 0.1, 1.0),
        ];
        let ids = vec![0, 1, 2, 3];
        let pipeline = CcdPipeline::new(20.0);
        let events = pipeline.detect_tunneling(&bodies, &ids);
        for w in events.windows(2) {
            assert!(w[0].toi <= w[1].toi);
        }
    }
    #[test]
    fn pipeline_earliest_event() {
        let bodies = vec![
            make_body([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
            make_body([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
        ];
        let ids = vec![0, 1];
        let pipeline = CcdPipeline::new(10.0);
        let ev = pipeline.earliest_event(&bodies, &ids);
        assert!(ev.is_some());
    }
    #[test]
    fn pipeline_no_events_empty() {
        let bodies: Vec<SweptBody> = vec![];
        let ids: Vec<usize> = vec![];
        let pipeline = CcdPipeline::new(1.0);
        let events = pipeline.detect_tunneling(&bodies, &ids);
        assert!(events.is_empty());
    }
    #[test]
    fn substep_step_to_toi() {
        let mut bodies = vec![make_body([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0)];
        let integrator = SubstepIntegrator::new(4);
        integrator.step_to_toi(&mut bodies, 2.0);
        assert!((bodies[0].pos[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn substep_resolve_event_bounces() {
        let mut bodies = vec![
            make_body([-0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
            make_body([0.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
        ];
        let event = CcdEvent {
            toi: 0.0,
            body_a: 0,
            body_b: 1,
            contact_normal: [-1.0, 0.0, 0.0],
            contact_point: [0.0, 0.0, 0.0],
        };
        let integrator = SubstepIntegrator::new(1);
        integrator.resolve_event(&mut bodies, &event, 1.0);
        assert!(bodies[0].vel[0] < 0.0);
        assert!(bodies[1].vel[0] > 0.0);
    }
    #[test]
    fn substep_resolve_static_body() {
        let mut bodies = vec![
            make_body([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
            make_body([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.5, 0.0),
        ];
        let event = CcdEvent {
            toi: 0.0,
            body_a: 0,
            body_b: 1,
            contact_normal: [-1.0, 0.0, 0.0],
            contact_point: [0.5, 0.0, 0.0],
        };
        let integrator = SubstepIntegrator::new(1);
        integrator.resolve_event(&mut bodies, &event, 1.0);
        assert!(bodies[0].vel[0] < 0.0);
        assert_eq!(bodies[1].vel, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn conservative_advancement_finds_toi() {
        let toi = conservative_advancement_rigid(
            [-0.6, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            [0.6, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            64,
        );
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "t={t}");
    }
    #[test]
    fn conservative_advancement_no_collision() {
        let toi = conservative_advancement_rigid(
            [-5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            64,
        );
        assert!(toi.is_none());
    }
    #[test]
    fn conservative_advancement_stationary() {
        let toi = conservative_advancement_rigid(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            2.0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            2.0,
            64,
        );
        assert_eq!(toi, Some(0.0));
    }
    #[test]
    fn swept_aabb_basic() {
        let (mn, mx) = linear_sweep_aabb([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5; 3], 2.0);
        assert!((mn[0] - (-0.5)).abs() < 1e-10);
        assert!((mx[0] - 2.5).abs() < 1e-10);
        assert!((mn[1] - (-0.5)).abs() < 1e-10);
        assert!((mx[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn swept_aabb_stationary() {
        let (mn, mx) = linear_sweep_aabb([1.0, 2.0, 3.0], [0.0; 3], [0.5; 3], 5.0);
        assert!((mn[0] - 0.5).abs() < 1e-10);
        assert!((mx[0] - 1.5).abs() < 1e-10);
        assert!((mn[2] - 2.5).abs() < 1e-10);
        assert!((mx[2] - 3.5).abs() < 1e-10);
    }
    #[test]
    fn swept_aabb_negative_velocity() {
        let (mn, mx) = linear_sweep_aabb([0.0; 3], [-2.0, 0.0, 0.0], [1.0; 3], 1.0);
        assert!((mn[0] - (-3.0)).abs() < 1e-10);
        assert!((mx[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn toi_linear_basic() {
        let t = time_of_impact_linear(10.0, 2.0);
        assert_eq!(t, Some(5.0));
    }
    #[test]
    fn toi_linear_diverging() {
        assert_eq!(time_of_impact_linear(5.0, -1.0), None);
    }
    #[test]
    fn toi_linear_already_overlapping() {
        assert_eq!(time_of_impact_linear(-1.0, 5.0), Some(0.0));
    }
    #[test]
    fn toi_linear_zero_closing_speed() {
        assert_eq!(time_of_impact_linear(3.0, 0.0), None);
    }
    #[test]
    fn speculative_candidate_approaching_bodies() {
        let a = make_body([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        assert!(speculative_contact_candidate(&a, &b, 10.0));
    }
    #[test]
    fn speculative_candidate_far_apart() {
        let a = make_body([0.0; 3], [0.0; 3], 0.1, 1.0);
        let b = make_body([100.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        assert!(!speculative_contact_candidate(&a, &b, 0.01));
    }
    #[test]
    fn bisection_refines_toi_head_on() {
        let a = make_body([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        let toi = refine_toi_bisection(&a, &b, 0.0, 1.0, 1e-6, 64);
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((t - 0.5).abs() < 0.01, "t={t}, expected ≈0.5");
    }
    #[test]
    fn bisection_returns_none_when_no_contact() {
        let a = make_body([-5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0);
        let b = make_body([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0);
        let toi = refine_toi_bisection(&a, &b, 0.0, 10.0, 1e-6, 64);
        assert!(toi.is_none());
    }
    #[test]
    fn bisection_already_overlapping() {
        let a = make_body([0.0; 3], [0.0; 3], 2.0, 1.0);
        let b = make_body([0.5, 0.0, 0.0], [0.0; 3], 2.0, 1.0);
        let toi = refine_toi_bisection(&a, &b, 0.0, 1.0, 1e-6, 64);
        assert_eq!(toi, Some(0.0));
    }
    #[test]
    fn rotating_body_conservative_radius_grows() {
        let b = RotatingSweptBody::new([0.0; 3], [0.0; 3], [0.0, 10.0, 0.0], 1.0, 1.0);
        let r = b.conservative_radius(0.1);
        assert!(
            r > 1.0,
            "conservative radius should be larger than base radius"
        );
    }
    #[test]
    fn rotating_body_no_rotation_same_radius() {
        let b = RotatingSweptBody::new([0.0; 3], [0.0; 3], [0.0; 3], 1.5, 1.0);
        assert!((b.conservative_radius(1.0) - 1.5).abs() < 1e-12);
    }
    #[test]
    fn rotating_body_ccd_candidate_close_bodies() {
        let a = RotatingSweptBody::new([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 0.5, 1.0);
        let b = RotatingSweptBody::new([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3], 0.5, 1.0);
        assert!(rotating_body_ccd_candidate(&a, &b, 1.0));
    }
    #[test]
    fn ccd_filter_should_check_matching() {
        let fa = CcdFilter::new(0b0001, 0b0010);
        let fb = CcdFilter::new(0b0010, 0b0001);
        assert!(CcdFilter::should_check(&fa, &fb));
    }
    #[test]
    fn ccd_filter_should_not_check_non_matching() {
        let fa = CcdFilter::new(0b0001, 0b0100);
        let fb = CcdFilter::new(0b0010, 0b1000);
        assert!(!CcdFilter::should_check(&fa, &fb));
    }
    #[test]
    fn ccd_broadphase_produces_candidate_pairs() {
        let bp = CcdBroadphase::new(1.0);
        let filter = CcdFilter::new(0xFFFF, 0xFFFF);
        let slots = vec![
            CcdBodySlot::new(
                make_body([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
                0,
                filter,
            ),
            CcdBodySlot::new(
                make_body([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
                1,
                filter,
            ),
        ];
        let pairs = bp.candidate_pairs(&slots);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 1));
    }
    #[test]
    fn ccd_broadphase_no_pairs_when_ccd_disabled() {
        let bp = CcdBroadphase::new(1.0);
        let filter = CcdFilter::new(0xFFFF, 0xFFFF);
        let slots = vec![
            CcdBodySlot::new(
                make_body([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
                0,
                filter,
            )
            .disable_ccd(),
            CcdBodySlot::new(
                make_body([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
                1,
                filter,
            )
            .disable_ccd(),
        ];
        let pairs = bp.candidate_pairs(&slots);
        assert!(pairs.is_empty());
    }
    #[test]
    fn ccd_stats_initial_state() {
        let s = CcdStats::new();
        assert_eq!(s.events_found, 0);
        assert_eq!(s.pairs_tested, 0);
        assert!(s.earliest_toi.is_infinite());
    }
    #[test]
    fn ccd_stats_record_event() {
        let mut s = CcdStats::new();
        s.record_event(2.5);
        s.record_event(1.0);
        assert_eq!(s.events_found, 2);
        assert!((s.earliest_toi - 1.0).abs() < 1e-12);
    }
    #[test]
    fn ccd_stats_reset() {
        let mut s = CcdStats::new();
        s.record_event(1.0);
        s.record_pair_test();
        s.reset();
        assert_eq!(s.events_found, 0);
        assert_eq!(s.pairs_tested, 0);
    }
    #[test]
    fn run_with_stats_resolves_collision() {
        let mut bodies = vec![
            make_body([-4.5, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1.0),
            make_body([4.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.5, 1.0),
        ];
        let ids = vec![0usize, 1usize];
        let pipeline = CcdPipeline::new(10.0);
        let stats = pipeline.run_with_stats(&mut bodies, &ids, 1.0);
        assert!(stats.events_found > 0, "should find at least one event");
        assert!(stats.earliest_toi.is_finite());
    }
    #[test]
    fn run_with_stats_no_collision() {
        let mut bodies = vec![make_body([0.0; 3], [0.0; 3], 0.1, 1.0)];
        let ids = vec![0usize];
        let pipeline = CcdPipeline::new(1.0);
        let stats = pipeline.run_with_stats(&mut bodies, &ids, 1.0);
        assert_eq!(stats.events_found, 0);
    }
    #[test]
    fn swept_capsule_center() {
        let c = SweptCapsule::new([0.0, 0.0, 0.0], [0.0, 4.0, 0.0], 0.5, [0.0; 3], 1.0);
        let cen = c.center();
        assert!((cen[1] - 2.0).abs() < 1e-10, "center.y={}", cen[1]);
    }
    #[test]
    fn swept_capsule_length() {
        let c = SweptCapsule::new([0.0; 3], [0.0, 3.0, 4.0], 0.5, [0.0; 3], 1.0);
        assert!((c.length() - 5.0).abs() < 1e-10, "length={}", c.length());
    }
    #[test]
    fn swept_capsule_bounding_sphere() {
        let c = SweptCapsule::new([0.0; 3], [2.0, 0.0, 0.0], 0.5, [0.0; 3], 1.0);
        assert!((c.bounding_sphere_radius() - 1.5).abs() < 1e-10);
    }
    #[test]
    fn swept_capsule_advance() {
        let c = SweptCapsule::new([0.0; 3], [1.0, 0.0, 0.0], 0.5, [2.0, 0.0, 0.0], 1.0);
        let moved = c.advance(1.5);
        assert!((moved.p0[0] - 3.0).abs() < 1e-10, "p0.x={}", moved.p0[0]);
        assert!((moved.p1[0] - 4.0).abs() < 1e-10, "p1.x={}", moved.p1[0]);
    }
    #[test]
    fn closest_point_segment_middle() {
        let (pt, t) = closest_point_on_segment([0.0; 3], [4.0, 0.0, 0.0], [2.0, 1.0, 0.0]);
        assert!((pt[0] - 2.0).abs() < 1e-10);
        assert!((t - 0.5).abs() < 1e-10);
    }
    #[test]
    fn closest_point_segment_clamp_start() {
        let (pt, t) = closest_point_on_segment([1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0; 3]);
        assert_eq!(pt, [1.0, 0.0, 0.0]);
        assert_eq!(t, 0.0);
    }
    #[test]
    fn closest_point_segment_clamp_end() {
        let (pt, t) = closest_point_on_segment([0.0; 3], [2.0, 0.0, 0.0], [5.0, 0.0, 0.0]);
        assert_eq!(pt, [2.0, 0.0, 0.0]);
        assert_eq!(t, 1.0);
    }
    #[test]
    fn segment_segment_parallel_same_height() {
        let (sq, _s, _t) = segment_segment_sq_dist(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!((sq - 1.0).abs() < 1e-6, "sq={sq}");
    }
    #[test]
    fn segment_segment_perpendicular_closest() {
        let (sq, _s, _t) = segment_segment_sq_dist(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!(sq < 1e-10, "sq={sq}, expected ~0");
    }
    #[test]
    fn capsule_capsule_head_on_collision() {
        let a = SweptCapsule::new(
            [-5.0, 0.0, 0.0],
            [-4.0, 0.0, 0.0],
            0.5,
            [1.0, 0.0, 0.0],
            1.0,
        );
        let b = SweptCapsule::new([4.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5, [-1.0, 0.0, 0.0], 1.0);
        let toi = check_capsule_capsule(&a, &b, 10.0, 64);
        assert!(toi.is_some(), "should collide");
    }
    #[test]
    fn capsule_capsule_no_collision() {
        let a = SweptCapsule::new(
            [-5.0, 0.0, 0.0],
            [-4.0, 0.0, 0.0],
            0.2,
            [-1.0, 0.0, 0.0],
            1.0,
        );
        let b = SweptCapsule::new([4.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.2, [1.0, 0.0, 0.0], 1.0);
        let toi = check_capsule_capsule(&a, &b, 1.0, 64);
        assert!(toi.is_none(), "should not collide");
    }
    #[test]
    fn ray_cast_sphere_hit() {
        let hit = ray_cast_sphere(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_some());
        let h = hit.unwrap();
        assert!((h.t - 4.0).abs() < 1e-6, "t={}", h.t);
        assert!((h.point[0] - (-1.0)).abs() < 1e-6);
    }
    #[test]
    fn ray_cast_sphere_miss() {
        let hit = ray_cast_sphere(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            100.0,
        );
        assert!(hit.is_none());
    }
    #[test]
    fn ray_cast_sphere_inside() {
        let hit = ray_cast_sphere([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], 2.0, 0.0, 100.0);
        assert!(hit.is_some());
    }
    #[test]
    fn ray_cast_aabb_hit() {
        let hit = ray_cast_aabb(
            [-5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_some());
        let (t_enter, _) = hit.unwrap();
        assert!((t_enter - 4.0).abs() < 1e-6, "t_enter={t_enter}");
    }
    #[test]
    fn ray_cast_aabb_miss() {
        let hit = ray_cast_aabb(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            0.0,
            100.0,
        );
        assert!(hit.is_none());
    }
    #[test]
    fn ray_cast_aabb_parallel_inside() {
        let hit = ray_cast_aabb([0.0; 3], [1.0, 0.0, 0.0], [-1.0; 3], [1.0; 3], 0.0, 100.0);
        assert!(hit.is_some());
    }
    #[test]
    fn ray_cast_plane_hit() {
        let t = ray_cast_plane(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.0,
            0.0,
            100.0,
        );
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-6, "t={}", t.unwrap());
    }
    #[test]
    fn ray_cast_plane_parallel_miss() {
        let t = ray_cast_plane(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            0.0,
            0.0,
            100.0,
        );
        assert!(t.is_none());
    }
    #[test]
    fn sphere_support_basic() {
        let s = SphereSupport {
            center: [0.0; 3],
            radius: 2.0,
        };
        let p = s.support([1.0, 0.0, 0.0]);
        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
    }
    #[test]
    fn aabb_support_basic() {
        let aabb = AabbSupport {
            min: [-1.0; 3],
            max: [1.0; 3],
        };
        let p = aabb.support([1.0, 0.0, 0.0]);
        assert!((p[0] - 1.0).abs() < 1e-10, "p.x={}", p[0]);
        let q = aabb.support([-1.0, 0.0, 0.0]);
        assert!((q[0] - (-1.0)).abs() < 1e-10, "q.x={}", q[0]);
    }
    #[test]
    fn gjk_distance_separated_spheres() {
        let a = SphereSupport {
            center: [-3.0, 0.0, 0.0],
            radius: 0.5,
        };
        let b = SphereSupport {
            center: [3.0, 0.0, 0.0],
            radius: 0.5,
        };
        let result = gjk_distance(&a, &b, 1e-8, 64);
        assert!(!result.overlapping);
        assert!(result.sq_dist > 10.0, "sq_dist={}", result.sq_dist);
    }
    #[test]
    fn gjk_distance_converges_for_touching_spheres() {
        let a = SphereSupport {
            center: [-0.5, 0.0, 0.0],
            radius: 0.5,
        };
        let b = SphereSupport {
            center: [0.5, 0.0, 0.0],
            radius: 0.5,
        };
        let result = gjk_distance(&a, &b, 1e-4, 64);
        assert!(
            result.overlapping || result.sq_dist < 0.05,
            "sq_dist={}",
            result.sq_dist
        );
    }
    #[test]
    fn minkowski_difference_support_test() {
        let a = SphereSupport {
            center: [0.0; 3],
            radius: 1.0,
        };
        let b = SphereSupport {
            center: [3.0, 0.0, 0.0],
            radius: 1.0,
        };
        let d = [1.0, 0.0, 0.0_f64];
        let md = minkowski_difference_support(&a, &b, d);
        assert!((md[0] - (1.0 - 2.0)).abs() < 1e-8, "md.x={}", md[0]);
    }
    #[test]
    fn gjk_ca_spheres_approaching() {
        let a = SphereSupport {
            center: [0.0; 3],
            radius: 0.5,
        };
        let b = SphereSupport {
            center: [0.0; 3],
            radius: 0.5,
        };
        let toi = gjk_conservative_advancement(
            &a,
            [-4.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            &b,
            [4.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            10.0,
            64,
            1e-6,
        );
        assert!(toi.is_some(), "spheres should collide");
        let t = toi.unwrap();
        assert!((0.0..=10.0).contains(&t), "t={t}");
    }
    #[test]
    fn gjk_ca_spheres_diverging() {
        let a = SphereSupport {
            center: [0.0; 3],
            radius: 0.5,
        };
        let b = SphereSupport {
            center: [0.0; 3],
            radius: 0.5,
        };
        let toi = gjk_conservative_advancement(
            &a,
            [-4.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            &b,
            [4.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            10.0,
            64,
            1e-6,
        );
        assert!(toi.is_none(), "spheres diverge — should not collide");
    }
    #[test]
    fn penalty_force_overlapping() {
        let f = penalty_force_sphere_sphere(
            [0.0; 3],
            [0.0; 3],
            [0.5, 0.0, 0.0],
            [0.0; 3],
            0.5,
            0.5,
            1000.0,
            0.0,
        );
        assert!(
            f[0] < 0.0,
            "force should push A in -x away from B: f={:?}",
            f
        );
        assert!((f[0] + 500.0).abs() < 1.0, "f.x={}", f[0]);
    }
    #[test]
    fn penalty_force_separated_zero() {
        let f = penalty_force_sphere_sphere(
            [0.0; 3],
            [0.0; 3],
            [10.0, 0.0, 0.0],
            [0.0; 3],
            0.5,
            0.5,
            1000.0,
            0.0,
        );
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn restitution_constant() {
        let m = RestitutionModel::Constant(0.7);
        assert!((m.evaluate(5.0) - 0.7).abs() < 1e-10);
    }
    #[test]
    fn restitution_speed_dependent_slow() {
        let m = RestitutionModel::SpeedDependent {
            e_max: 1.0,
            v_max: 10.0,
        };
        assert!((m.evaluate(0.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn restitution_speed_dependent_fast() {
        let m = RestitutionModel::SpeedDependent {
            e_max: 1.0,
            v_max: 10.0,
        };
        assert!(m.evaluate(10.0) <= 1e-10);
        assert!(m.evaluate(20.0) <= 1e-10);
    }
    #[test]
    fn compute_collision_impulse_elastic() {
        let j = compute_collision_impulse(
            [2.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            1.0,
            1.0,
            &RestitutionModel::Constant(1.0),
        );
        assert!((j - 2.0).abs() < 1e-10, "j={j}");
    }
    #[test]
    fn compute_collision_impulse_separating() {
        let j = compute_collision_impulse(
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            1.0,
            1.0,
            &RestitutionModel::Constant(1.0),
        );
        assert_eq!(j, 0.0);
    }
    #[test]
    fn baumgarte_pushes_apart() {
        let (da, db) = baumgarte_correction([0.0, 1.0, 0.0], 0.1, 1.0, 1.0, 0.2, 0.01);
        assert!(da[1] > 0.0, "da.y={}", da[1]);
        assert!(db[1] < 0.0, "db.y={}", db[1]);
    }
    #[test]
    fn baumgarte_zero_dt() {
        let (da, db) = baumgarte_correction([0.0, 1.0, 0.0], 0.1, 1.0, 1.0, 0.2, 0.0);
        assert_eq!(da, [0.0; 3]);
        assert_eq!(db, [0.0; 3]);
    }
    #[test]
    fn cross_basic() {
        let c = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-12);
        assert!(c[0].abs() < 1e-12);
        assert!(c[1].abs() < 1e-12);
    }
    #[test]
    fn ccd_sweep_sphere_head_on_collision() {
        let toi = CcdSweep::compute_collision_time_sphere(
            [-5.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            1.0,
            [5.0, 0.0, 0.0],
            [-10.0, 0.0, 0.0],
            1.0,
            1.0,
        );
        assert!(toi.is_some(), "head-on spheres should collide");
        let t = toi.unwrap();
        assert!((t - 0.4).abs() < 1e-6, "toi = {t}");
    }
    #[test]
    fn ccd_sweep_sphere_no_collision_moving_apart() {
        let toi = CcdSweep::compute_collision_time_sphere(
            [-5.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            0.5,
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
            1.0,
        );
        assert!(toi.is_none(), "diverging spheres should not collide");
    }
    #[test]
    fn ccd_sweep_sphere_already_overlapping_returns_zero() {
        let toi = CcdSweep::compute_collision_time_sphere(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            1.0,
            1.0,
        );
        assert_eq!(toi, Some(0.0));
    }
    #[test]
    fn ccd_sweep_sphere_miss_lateral() {
        let toi = CcdSweep::compute_collision_time_sphere(
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.4,
            [0.0, 100.0, 0.0],
            [10.0, 0.0, 0.0],
            0.4,
            1.0,
        );
        assert!(toi.is_none(), "lateral miss should produce no collision");
    }
    #[test]
    fn ccd_substep_count_static_scene() {
        let n = CcdSweep::compute_sub_step_count(0.0, 1.0, 0.016, 0.5, 64);
        assert_eq!(n, 1, "static scene → 1 substep");
    }
    #[test]
    fn ccd_substep_count_fast_small_body() {
        let n = CcdSweep::compute_sub_step_count(100.0, 0.1, 0.016, 0.5, 128);
        assert!(n >= 32, "fast small body needs many substeps: {n}");
    }
    #[test]
    fn ccd_substep_count_clamped_at_max() {
        let n = CcdSweep::compute_sub_step_count(1e9, 1e-6, 1.0, 0.5, 16);
        assert_eq!(n, 16, "substep count clamped at max");
    }
    #[test]
    fn ccd_detect_tunneling_none_when_slow() {
        let bodies = vec![
            SweptBody::new([0.0, 0.0, 0.0], [0.1, 0.0, 0.0], 1.0, 1.0),
            SweptBody::new([5.0, 0.0, 0.0], [0.0, 0.1, 0.0], 1.0, 1.0),
        ];
        let indices = CcdSweep::detect_tunneling(&bodies, 0.016, 1.0);
        assert!(
            indices.is_empty(),
            "slow bodies should not be detected as tunneling"
        );
    }
    #[test]
    fn ccd_detect_tunneling_fast_body() {
        let bodies = vec![SweptBody::new(
            [0.0, 0.0, 0.0],
            [1000.0, 0.0, 0.0],
            0.1,
            1.0,
        )];
        let indices = CcdSweep::detect_tunneling(&bodies, 0.016, 1.0);
        assert!(
            !indices.is_empty(),
            "fast body should be flagged for tunneling"
        );
        assert_eq!(indices[0], 0);
    }
    #[test]
    fn ccd_sweep_sphere_toi_within_dt() {
        let toi = CcdSweep::compute_collision_time_sphere(
            [-3.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            0.5,
            [3.0, 0.0, 0.0],
            [-5.0, 0.0, 0.0],
            0.5,
            1.0,
        );
        assert!(toi.is_some(), "should collide within dt=1");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "toi should be in [0, dt]: {t}");
        assert!((t - 0.5).abs() < 1e-6, "expected toi≈0.5, got {t}");
    }
}
