//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CcdResult, LinearCastResult, SweptSphere};

#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
#[inline]
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(a, 1.0 / l)
    }
}
/// Apply a 4x4 homogeneous transform (row-major) to a 3-D point.
pub(super) fn transform_point(m: [[f64; 4]; 4], p: [f64; 3]) -> [f64; 3] {
    let x = m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3];
    let y = m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3];
    let z = m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3];
    [x, y, z]
}
/// Apply a 4x4 homogeneous transform (row-major) to a 3-D direction (no translation).
#[cfg(test)]
pub(super) fn transform_dir(m: [[f64; 4]; 4], d: [f64; 3]) -> [f64; 3] {
    let x = m[0][0] * d[0] + m[0][1] * d[1] + m[0][2] * d[2];
    let y = m[1][0] * d[0] + m[1][1] * d[1] + m[1][2] * d[2];
    let z = m[2][0] * d[0] + m[2][1] * d[1] + m[2][2] * d[2];
    [x, y, z]
}
/// Interpolate two 4x4 matrices element-wise (linear blend).
pub(super) fn lerp_matrix(a: [[f64; 4]; 4], b: [[f64; 4]; 4], t: f64) -> [[f64; 4]; 4] {
    let mut result = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            result[i][j] = a[i][j] + t * (b[i][j] - a[i][j]);
        }
    }
    result
}
/// Minkowski sum of two spheres: simply the sum of their radii.
pub fn minkowski_sum_spheres(r1: f64, r2: f64) -> f64 {
    r1 + r2
}
/// Minkowski sum of two AABBs.
///
/// For boxes A = \[min1, max1\] and B = \[min2, max2\] the Minkowski sum is
/// `result_min[i] = min1[i] + min2[i]`,
/// `result_max[i] = max1[i] + max2[i]`.
pub fn minkowski_sum_aabbs(
    min1: [f64; 3],
    max1: [f64; 3],
    min2: [f64; 3],
    max2: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let result_min = add3(min1, min2);
    let result_max = add3(max1, max2);
    (result_min, result_max)
}
/// Minkowski difference AABB.
///
/// For boxes A = \[min1, max1\] and B = \[min2, max2\] the Minkowski difference
/// (A minus B) has `result_min[i] = min1[i] - max2[i]`,
/// `result_max[i] = max1[i] - min2[i]`.
pub fn minkowski_diff_aabbs(
    min1: [f64; 3],
    max1: [f64; 3],
    min2: [f64; 3],
    max2: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let result_min = sub3(min1, max2);
    let result_max = sub3(max1, min2);
    (result_min, result_max)
}
/// Support point of the Minkowski difference A + (-B) in direction `dir`.
///
/// Used as a building block for GJK: the support of A-B is
/// `support_A(dir) - support_B(-dir)`.
pub fn minkowski_support(
    support_a: impl Fn([f64; 3]) -> [f64; 3],
    support_b: impl Fn([f64; 3]) -> [f64; 3],
    dir: [f64; 3],
) -> [f64; 3] {
    let sa = support_a(dir);
    let sb = support_b(neg3(dir));
    sub3(sa, sb)
}
/// Time-of-impact for a moving sphere (the `SweptSphere`) against a static
/// sphere.
///
/// The swept sphere travels linearly from `center_start` to `center_end`.
/// Returns the normalised time `t` in `[0, 1]` of first contact, or `None` if
/// there is no collision during the sweep.
pub fn swept_sphere_vs_sphere(
    swept: &SweptSphere,
    sphere_center: [f64; 3],
    sphere_radius: f64,
) -> Option<f64> {
    let combined_radius = swept.radius + sphere_radius;
    let vel = sub3(swept.center_end, swept.center_start);
    let rel = sub3(swept.center_start, sphere_center);
    let a = dot3(vel, vel);
    let b = 2.0 * dot3(rel, vel);
    let c = dot3(rel, rel) - combined_radius * combined_radius;
    if a < 1e-14 {
        return if c <= 0.0 { Some(0.0) } else { None };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let t1 = (-b - sq) / (2.0 * a);
    let t2 = (-b + sq) / (2.0 * a);
    let t = if t1 >= 0.0 { t1 } else { t2 };
    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}
/// Time-of-impact for a moving sphere against an infinite plane.
///
/// The plane is defined by `plane_normal` (unit vector) and scalar `plane_d`
/// such that `dot(normal, x) = plane_d` for points `x` on the plane.
///
/// Returns `t` in `[0, 1]` of first contact, or `None` if no contact during the
/// sweep.
pub fn swept_sphere_vs_plane(
    swept: &SweptSphere,
    plane_normal: [f64; 3],
    plane_d: f64,
) -> Option<f64> {
    let vel = sub3(swept.center_end, swept.center_start);
    let d_start = dot3(swept.center_start, plane_normal) - plane_d;
    let d_vel = dot3(vel, plane_normal);
    let mut t_min = f64::INFINITY;
    for &sign in &[1.0_f64, -1.0_f64] {
        let target = sign * swept.radius;
        if d_vel.abs() > 1e-14 {
            let t = (target - d_start) / d_vel;
            if (0.0..=1.0).contains(&t) && t < t_min {
                t_min = t;
            }
        } else {
            if (d_start - target).abs() < 1e-12 {
                t_min = 0.0;
            }
        }
    }
    if t_min.is_finite() { Some(t_min) } else { None }
}
/// Time-of-impact for two moving spheres.
///
/// Sphere A moves from `a_start` to `a_end` with radius `ra`.
/// Sphere B moves from `b_start` to `b_end` with radius `rb`.
///
/// Returns `t` in `[0, 1]` of first contact, or `None`.
pub fn swept_sphere_vs_swept_sphere(
    a_start: [f64; 3],
    a_end: [f64; 3],
    ra: f64,
    b_start: [f64; 3],
    b_end: [f64; 3],
    rb: f64,
) -> Option<f64> {
    let combined = ra + rb;
    let rel_start = sub3(a_start, b_start);
    let rel_end = sub3(a_end, b_end);
    let vel = sub3(rel_end, rel_start);
    let a = dot3(vel, vel);
    let b = 2.0 * dot3(rel_start, vel);
    let c = dot3(rel_start, rel_start) - combined * combined;
    if a < 1e-14 {
        return if c <= 0.0 { Some(0.0) } else { None };
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let t1 = (-b - sq) / (2.0 * a);
    let t2 = (-b + sq) / (2.0 * a);
    let t = if t1 >= 0.0 { t1 } else { t2 };
    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}
/// Linear cast: move sphere A from `start` along `displacement` and find
/// the time of first contact with a static sphere B.
///
/// Returns `None` if no contact occurs within the displacement.
pub fn linear_cast_sphere_vs_sphere(
    a_center: [f64; 3],
    a_radius: f64,
    displacement: [f64; 3],
    b_center: [f64; 3],
    b_radius: f64,
) -> Option<LinearCastResult> {
    let a_end = add3(a_center, displacement);
    let swept = SweptSphere::new(a_center, a_end, a_radius);
    let toi = swept_sphere_vs_sphere(&swept, b_center, b_radius)?;
    let contact_a = lerp3(a_center, a_end, toi);
    let dir = sub3(contact_a, b_center);
    let dist = len3(dir);
    let normal = if dist > 1e-14 {
        scale3(dir, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    let contact_point = add3(b_center, scale3(normal, b_radius));
    Some(LinearCastResult {
        toi,
        contact_point,
        normal,
    })
}
/// Linear cast: move sphere from `start` along `displacement` and find
/// the time of first contact with a static AABB.
///
/// Uses slab method on the AABB inflated by the sphere radius.
pub fn linear_cast_sphere_vs_aabb(
    center: [f64; 3],
    radius: f64,
    displacement: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> Option<f64> {
    let inf_min = [
        aabb_min[0] - radius,
        aabb_min[1] - radius,
        aabb_min[2] - radius,
    ];
    let inf_max = [
        aabb_max[0] + radius,
        aabb_max[1] + radius,
        aabb_max[2] + radius,
    ];
    ray_vs_aabb(center, displacement, inf_min, inf_max)
}
/// Ray-AABB intersection returning `t` in `[0, 1]` or `None`.
pub(super) fn ray_vs_aabb(
    origin: [f64; 3],
    dir: [f64; 3],
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
) -> Option<f64> {
    let mut t_near = f64::NEG_INFINITY;
    let mut t_far = f64::INFINITY;
    for k in 0..3 {
        if dir[k].abs() < 1e-14 {
            if origin[k] < aabb_min[k] || origin[k] > aabb_max[k] {
                return None;
            }
        } else {
            let inv_d = 1.0 / dir[k];
            let mut t1 = (aabb_min[k] - origin[k]) * inv_d;
            let mut t2 = (aabb_max[k] - origin[k]) * inv_d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            t_near = t_near.max(t1);
            t_far = t_far.min(t2);
            if t_near > t_far {
                return None;
            }
        }
    }
    if (0.0..=1.0).contains(&t_near) {
        Some(t_near)
    } else if t_near < 0.0 && t_far >= 0.0 {
        Some(0.0)
    } else {
        None
    }
}
/// Compute a conservative AABB for an arbitrary shape swept linearly.
///
/// The shape is defined by its local AABB `(local_min, local_max)` at the
/// start position and `displacement` is the total translation.
pub fn swept_volume_aabb(
    local_min: [f64; 3],
    local_max: [f64; 3],
    displacement: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let mut result_min = local_min;
    let mut result_max = local_max;
    for k in 0..3 {
        if displacement[k] > 0.0 {
            result_max[k] += displacement[k];
        } else {
            result_min[k] += displacement[k];
        }
    }
    (result_min, result_max)
}
/// Compute conservative AABB for a moving sphere over a time interval.
///
/// `positions` is a list of sampled center positions; the AABB is expanded
/// by `radius` in all directions.
pub fn swept_sphere_aabb_sampled(positions: &[[f64; 3]], radius: f64) -> ([f64; 3], [f64; 3]) {
    if positions.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = positions[0];
    let mut mx = positions[0];
    for p in &positions[1..] {
        for k in 0..3 {
            if p[k] < mn[k] {
                mn[k] = p[k];
            }
            if p[k] > mx[k] {
                mx[k] = p[k];
            }
        }
    }
    for k in 0..3 {
        mn[k] -= radius;
        mx[k] += radius;
    }
    (mn, mx)
}
/// Closest point on a line segment `[a, b]` to point `p`.
pub fn closest_point_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = sub3(b, a);
    let ap = sub3(p, a);
    let ab_sq = dot3(ab, ab);
    if ab_sq < 1e-14 {
        return a;
    }
    let t = (dot3(ap, ab) / ab_sq).clamp(0.0, 1.0);
    add3(a, scale3(ab, t))
}
/// Distance from point `p` to segment `[a, b]`.
pub fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(p, closest_point_on_segment(p, a, b)))
}
/// Closest distance between two line segments `[a1, b1]` and `[a2, b2]`.
pub fn segment_segment_distance(a1: [f64; 3], b1: [f64; 3], a2: [f64; 3], b2: [f64; 3]) -> f64 {
    let d1 = sub3(b1, a1);
    let d2 = sub3(b2, a2);
    let r = sub3(a1, a2);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    if a < 1e-14 && e < 1e-14 {
        return len3(r);
    }
    let (s, t);
    if a < 1e-14 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot3(d1, r);
        if e < 1e-14 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b_val = dot3(d1, d2);
            let denom = a * e - b_val * b_val;
            s = if denom.abs() > 1e-14 {
                ((b_val * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            t = ((b_val * s + f) / e).clamp(0.0, 1.0);
        }
    }
    let closest1 = add3(a1, scale3(d1, s));
    let closest2 = add3(a2, scale3(d2, t));
    len3(sub3(closest1, closest2))
}
/// Time-of-impact between two swept spheres represented as capsules.
///
/// This uses segment-segment closest-distance to check if two capsules
/// (swept spheres) intersect.
pub fn capsule_capsule_toi(a: &SweptSphere, b: &SweptSphere) -> Option<f64> {
    let combined_radius = a.radius + b.radius;
    let n_steps = 32;
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    let dist0 = len3(sub3(a.center_start, b.center_start));
    if dist0 <= combined_radius {
        return Some(0.0);
    }
    let dist1 = len3(sub3(a.center_end, b.center_end));
    let mid_a = a.center_at(0.5);
    let mid_b = b.center_at(0.5);
    let dist_mid = len3(sub3(mid_a, mid_b));
    let min_dist = dist0.min(dist1).min(dist_mid);
    if min_dist
        > combined_radius
            + len3(sub3(a.center_end, a.center_start))
            + len3(sub3(b.center_end, b.center_start))
    {
        return None;
    }
    for _ in 0..n_steps {
        let t = (lo + hi) * 0.5;
        let pa = a.center_at(t);
        let pb = b.center_at(t);
        let d = len3(sub3(pa, pb));
        if d <= combined_radius {
            hi = t;
        } else {
            lo = t;
        }
    }
    if hi <= 1.0 {
        let pa = a.center_at(hi);
        let pb = b.center_at(hi);
        let d = len3(sub3(pa, pb));
        if d <= combined_radius + 1e-6 {
            return Some(hi);
        }
    }
    None
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::swept::SweptBox;
    use crate::swept::SweptCapsule;

    use crate::swept::SweptSphere;

    fn identity() -> [[f64; 4]; 4] {
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
    fn translate(tx: f64, ty: f64, tz: f64) -> [[f64; 4]; 4] {
        [
            [1.0, 0.0, 0.0, tx],
            [0.0, 1.0, 0.0, ty],
            [0.0, 0.0, 1.0, tz],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
    #[test]
    fn test_minkowski_sum_aabbs_unit_boxes() {
        let (rmin, rmax) = minkowski_sum_aabbs(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert_eq!(rmin, [-2.0, -2.0, -2.0]);
        assert_eq!(rmax, [2.0, 2.0, 2.0]);
    }
    #[test]
    fn test_swept_sphere_vs_sphere_toi() {
        let swept = SweptSphere::new([-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5);
        let toi = swept_sphere_vs_sphere(&swept, [0.0, 0.0, 0.0], 0.5);
        let t = toi.expect("expected a TOI");
        assert!((t - 0.4).abs() < 1e-9, "expected TOI ~ 0.4, got {}", t);
    }
    #[test]
    fn test_swept_sphere_vs_sphere_no_hit() {
        let swept = SweptSphere::new([-5.0, 5.0, 0.0], [5.0, 5.0, 0.0], 0.5);
        assert!(swept_sphere_vs_sphere(&swept, [0.0, 0.0, 0.0], 0.5).is_none());
    }
    #[test]
    fn test_swept_sphere_vs_plane() {
        let swept = SweptSphere::new([0.0, 0.0, 5.0], [0.0, 0.0, -5.0], 1.0);
        let toi = swept_sphere_vs_plane(&swept, [0.0, 0.0, 1.0], 0.0);
        let t = toi.expect("expected a TOI");
        assert!((t - 0.4).abs() < 1e-9, "expected TOI ~ 0.4, got {}", t);
    }
    #[test]
    fn test_swept_sphere_vs_plane_no_hit() {
        let swept = SweptSphere::new([0.0, 0.0, 5.0], [0.0, 0.0, 10.0], 1.0);
        let toi = swept_sphere_vs_plane(&swept, [0.0, 0.0, 1.0], 0.0);
        assert!(toi.is_none(), "expected no TOI, got {:?}", toi);
    }
    #[test]
    fn test_swept_sphere_aabb_contains_endpoints() {
        let start = [-3.0, 1.0, 0.0];
        let end = [3.0, -1.0, 2.0];
        let r = 0.5;
        let swept = SweptSphere::new(start, end, r);
        let (mn, mx) = swept.aabb();
        for k in 0..3 {
            assert!(mn[k] <= start[k] - r + 1e-12, "min[{}] too large", k);
            assert!(
                mx[k] >= start[k] + r - 1e-12,
                "max[{}] too small (start)",
                k
            );
            assert!(mn[k] <= end[k] - r + 1e-12, "min[{}] too large (end)", k);
            assert!(mx[k] >= end[k] + r - 1e-12, "max[{}] too small (end)", k);
        }
    }
    #[test]
    fn test_swept_sphere_ray_hit() {
        let swept = SweptSphere::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        let hit = swept.ray_intersect([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        let t = hit.expect("expected a hit");
        assert!((t - 4.5).abs() < 1e-9, "expected t=4.5, got {}", t);
    }
    #[test]
    fn test_swept_sphere_ray_miss() {
        let swept = SweptSphere::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
        assert!(
            swept
                .ray_intersect([0.0, 5.0, 0.0], [1.0, 0.0, 0.0])
                .is_none()
        );
    }
    #[test]
    fn test_swept_box_aabb_identity() {
        let sb = SweptBox::new(identity(), identity(), [1.0, 1.0, 1.0]);
        let (mn, mx) = sb.aabb();
        for k in 0..3 {
            assert!((mn[k] + 1.0).abs() < 1e-12, "min[{}]={} != -1", k, mn[k]);
            assert!((mx[k] - 1.0).abs() < 1e-12, "max[{}]={} != 1", k, mx[k]);
        }
    }
    #[test]
    fn test_swept_box_aabb_translated() {
        let sb = SweptBox::new(identity(), translate(10.0, 0.0, 0.0), [1.0, 1.0, 1.0]);
        let (mn, mx) = sb.aabb();
        assert!((mn[0] + 1.0).abs() < 1e-12);
        assert!((mx[0] - 11.0).abs() < 1e-12);
    }
    #[test]
    fn test_minkowski_support() {
        let support_a = |d: [f64; 3]| {
            let l = len3(d);
            if l < 1e-14 {
                [1.0, 0.0, 0.0]
            } else {
                scale3(d, 1.0 / l)
            }
        };
        let support_b = |d: [f64; 3]| {
            let l = len3(d);
            if l < 1e-14 {
                [2.0, 0.0, 0.0]
            } else {
                scale3(d, 2.0 / l)
            }
        };
        let s = minkowski_support(support_a, support_b, [1.0, 0.0, 0.0]);
        assert!((s[0] - 3.0).abs() < 1e-12, "support[0]={}", s[0]);
        assert!(s[1].abs() < 1e-12);
        assert!(s[2].abs() < 1e-12);
    }
    #[test]
    fn test_minkowski_sum_spheres() {
        assert!((minkowski_sum_spheres(1.0, 2.0) - 3.0).abs() < 1e-12);
        assert!((minkowski_sum_spheres(0.5, 0.5) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_transform_point_and_dir() {
        let m = translate(1.0, 2.0, 3.0);
        let p = transform_point(m, [0.0, 0.0, 0.0]);
        assert_eq!(p, [1.0, 2.0, 3.0]);
        let d = transform_dir(m, [1.0, 0.0, 0.0]);
        assert_eq!(d, [1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_swept_sphere_volume() {
        let swept = SweptSphere::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0);
        let vol = swept.volume();
        let expected = std::f64::consts::PI * 2.0 + (4.0 / 3.0) * std::f64::consts::PI;
        assert!(
            (vol - expected).abs() < 1e-9,
            "vol={} expected={}",
            vol,
            expected
        );
    }
    #[test]
    fn test_swept_sphere_sweep_length() {
        let swept = SweptSphere::new([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 1.0);
        assert!((swept.sweep_length() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_sphere_center_at() {
        let swept = SweptSphere::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 1.0);
        let mid = swept.center_at(0.5);
        assert!((mid[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_capsule_aabb() {
        let sc = SweptCapsule::new([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5, 1.0);
        let (mn, mx) = sc.aabb();
        assert!(mn[0] <= -0.5 + 1e-12);
        assert!(mx[0] >= 10.5 - 1e-12);
        assert!(mn[1] <= -1.5 + 1e-12);
        assert!(mx[1] >= 1.5 - 1e-12);
    }
    #[test]
    fn test_swept_capsule_toi_vs_sphere() {
        let sc = SweptCapsule::new([-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5, 0.5);
        let toi = sc.toi_vs_sphere([0.0, 0.0, 0.0], 0.5);
        assert!(toi.is_some(), "expected a TOI");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t));
    }
    #[test]
    fn test_swept_capsule_volume() {
        let sc = SweptCapsule::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 1.0);
        let vol = sc.capsule_volume();
        let expected = std::f64::consts::PI * 1.0 * 2.0 + (4.0 / 3.0) * std::f64::consts::PI;
        assert!((vol - expected).abs() < 1e-9);
    }
    #[test]
    fn test_swept_sphere_vs_swept_sphere() {
        let toi = swept_sphere_vs_swept_sphere(
            [-5.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.5,
        );
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&t));
    }
    #[test]
    fn test_swept_sphere_vs_swept_sphere_no_hit() {
        let toi = swept_sphere_vs_swept_sphere(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
            [10.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            0.5,
        );
        assert!(toi.is_none());
    }
    #[test]
    fn test_linear_cast_sphere_vs_sphere_hit() {
        let result = linear_cast_sphere_vs_sphere(
            [-5.0, 0.0, 0.0],
            0.5,
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_some());
        let r = result.unwrap();
        assert!((r.toi - 0.4).abs() < 1e-9);
    }
    #[test]
    fn test_linear_cast_sphere_vs_sphere_miss() {
        let result = linear_cast_sphere_vs_sphere(
            [-5.0, 5.0, 0.0],
            0.5,
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.5,
        );
        assert!(result.is_none());
    }
    #[test]
    fn test_linear_cast_sphere_vs_aabb() {
        let toi = linear_cast_sphere_vs_aabb(
            [-5.0, 0.0, 0.0],
            0.5,
            [10.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(toi.is_some());
        let t = toi.unwrap();
        assert!((t - 0.35).abs() < 1e-9, "t={}", t);
    }
    #[test]
    fn test_swept_volume_aabb() {
        let (mn, mx) = swept_volume_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0], [5.0, -3.0, 0.0]);
        assert!((mn[0] - (-1.0)).abs() < 1e-12);
        assert!((mx[0] - 6.0).abs() < 1e-12);
        assert!((mn[1] - (-4.0)).abs() < 1e-12);
        assert!((mx[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_sphere_aabb_sampled() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [3.0, -1.0, 1.0]];
        let (mn, mx) = swept_sphere_aabb_sampled(&positions, 0.5);
        assert!((mn[0] - (-0.5)).abs() < 1e-12);
        assert!((mx[0] - 3.5).abs() < 1e-12);
        assert!((mn[1] - (-1.5)).abs() < 1e-12);
        assert!((mx[1] - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_closest_point_on_segment() {
        let cp = closest_point_on_segment([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((cp[0]).abs() < 1e-12);
        assert!((cp[1]).abs() < 1e-12);
    }
    #[test]
    fn test_point_segment_distance() {
        let d = point_segment_distance([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_segment_segment_distance_parallel() {
        let d = segment_segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        );
        assert!((d - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_segment_segment_distance_crossing() {
        let d = segment_segment_distance(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!((d - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_minkowski_diff_aabbs() {
        let (mn, mx) = minkowski_diff_aabbs(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [-0.5, -0.5, -0.5],
            [0.5, 0.5, 0.5],
        );
        assert!((mn[0] - (-1.5)).abs() < 1e-12);
        assert!((mx[0] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_swept_box_volume() {
        let sb = SweptBox::new(identity(), identity(), [1.0, 2.0, 3.0]);
        assert!((sb.box_volume() - 48.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_box_displacement() {
        let sb = SweptBox::new(identity(), translate(5.0, 3.0, 1.0), [1.0, 1.0, 1.0]);
        let d = sb.displacement();
        assert!((d[0] - 5.0).abs() < 1e-12);
        assert!((d[1] - 3.0).abs() < 1e-12);
        assert!((d[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_box_aabb_sampled_matches_endpoint() {
        let sb = SweptBox::new(identity(), translate(10.0, 0.0, 0.0), [1.0, 1.0, 1.0]);
        let (mn2, mx2) = sb.aabb();
        let (mn_s, mx_s) = sb.aabb_sampled(10);
        for k in 0..3 {
            assert!(mn_s[k] <= mn2[k] + 1e-12);
            assert!(mx_s[k] >= mx2[k] - 1e-12);
        }
    }
    #[test]
    fn test_capsule_capsule_toi_overlap() {
        let a = SweptSphere::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        let b = SweptSphere::new([0.5, 0.0, 0.0], [0.5, 0.0, 0.0], 1.0);
        let toi = capsule_capsule_toi(&a, &b);
        assert!(toi.is_some());
        assert!((toi.unwrap()).abs() < 1e-6, "overlapping => toi=0");
    }
    #[test]
    fn test_capsule_capsule_toi_no_hit() {
        let a = SweptSphere::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1);
        let b = SweptSphere::new([10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.1);
        assert!(capsule_capsule_toi(&a, &b).is_none());
    }
    #[test]
    fn test_lerp_matrix() {
        let a = identity();
        let b = translate(10.0, 20.0, 30.0);
        let mid = lerp_matrix(a, b, 0.5);
        assert!((mid[0][3] - 5.0).abs() < 1e-12);
        assert!((mid[1][3] - 10.0).abs() < 1e-12);
        assert!((mid[2][3] - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize3() {
        let n = normalize3([3.0, 4.0, 0.0]);
        assert!((n[0] - 0.6).abs() < 1e-12);
        assert!((n[1] - 0.8).abs() < 1e-12);
    }
    #[test]
    fn test_ray_vs_aabb_miss() {
        let result = ray_vs_aabb(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_none());
    }
    #[test]
    fn test_swept_sphere_surface_area() {
        let swept = SweptSphere::new([0.0; 3], [0.0; 3], 1.0);
        let sa = swept.surface_area();
        assert!((sa - 4.0 * std::f64::consts::PI).abs() < 1e-9);
    }
}
/// Conservative CCD: time-of-impact for a moving sphere vs a static capsule.
///
/// The capsule axis runs from `cap_a` to `cap_b` with `cap_radius`.
/// The sphere moves from `sphere_start` to `sphere_end` with `sphere_radius`.
///
/// Returns `Some(CcdResult)` at the earliest time of contact, or `None`.
pub fn ccd_sphere_vs_capsule(
    sphere_start: [f64; 3],
    sphere_end: [f64; 3],
    sphere_radius: f64,
    cap_a: [f64; 3],
    cap_b: [f64; 3],
    cap_radius: f64,
) -> Option<CcdResult> {
    let combined = sphere_radius + cap_radius;
    let n_steps = 32_usize;
    {
        let d = point_segment_distance(sphere_start, cap_a, cap_b);
        if d <= combined {
            let closest = closest_point_on_segment(sphere_start, cap_a, cap_b);
            return Some(CcdResult {
                toi: 0.0,
                contact_point: closest,
            });
        }
    }
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    let d_end = point_segment_distance(sphere_end, cap_a, cap_b);
    if d_end > combined {
        let mid = lerp3(sphere_start, sphere_end, 0.5);
        let d_mid = point_segment_distance(mid, cap_a, cap_b);
        if d_mid > combined {
            return None;
        }
    }
    for _ in 0..n_steps {
        let t = (lo + hi) * 0.5;
        let pos = lerp3(sphere_start, sphere_end, t);
        let d = point_segment_distance(pos, cap_a, cap_b);
        if d <= combined {
            hi = t;
        } else {
            lo = t;
        }
    }
    let contact_pos = lerp3(sphere_start, sphere_end, hi);
    let contact_point = closest_point_on_segment(contact_pos, cap_a, cap_b);
    Some(CcdResult {
        toi: hi,
        contact_point,
    })
}
/// CCD time-of-impact for two moving AABBs.
///
/// Uses the separating axis theorem on the swept union: returns the earliest
/// time `t` in `[0, 1]` at which the AABBs first overlap, or `None`.
pub fn ccd_aabb_vs_aabb(
    a_min_start: [f64; 3],
    a_max_start: [f64; 3],
    a_displacement: [f64; 3],
    b_min: [f64; 3],
    b_max: [f64; 3],
) -> Option<f64> {
    let mut t_enter = 0.0f64;
    let mut t_exit = 1.0f64;
    for k in 0..3 {
        let v = a_displacement[k];
        let a_lo = a_min_start[k];
        let a_hi = a_max_start[k];
        let b_lo = b_min[k];
        let b_hi = b_max[k];
        if v.abs() < 1e-14 {
            if a_hi < b_lo || a_lo > b_hi {
                return None;
            }
        } else {
            let inv_v = 1.0 / v;
            let t0 = (b_lo - a_hi) * inv_v;
            let t1 = (b_hi - a_lo) * inv_v;
            let (t_in, t_out) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
            t_enter = t_enter.max(t_in);
            t_exit = t_exit.min(t_out);
            if t_enter > t_exit {
                return None;
            }
        }
    }
    if (0.0..=1.0).contains(&t_enter) {
        Some(t_enter)
    } else if t_enter < 0.0 && t_exit >= 0.0 {
        Some(0.0)
    } else {
        None
    }
}
/// Estimate the surface area of a swept sphere (capsule) numerically.
///
/// Uses the formula: lateral = 2πr·L, end caps = 4πr².
/// This is a fast closed-form estimate.
pub fn swept_sphere_surface_area(radius: f64, sweep_length: f64) -> f64 {
    2.0 * std::f64::consts::PI * radius * sweep_length
        + 4.0 * std::f64::consts::PI * radius * radius
}
/// Estimate the surface area of a swept box via the start and end AABBs.
///
/// Conservative upper bound: sum of the two individual box surface areas
/// plus the lateral contribution from the sweep path.
pub fn swept_box_surface_area_estimate(half_extents: [f64; 3], displacement: [f64; 3]) -> f64 {
    let [hx, hy, hz] = half_extents;
    let box_sa = 8.0 * (hx * hy + hy * hz + hz * hx);
    let lateral_patch = {
        let [dx, dy, dz] = displacement;
        let nx = dy.abs() * 2.0 * hz + dz.abs() * 2.0 * hy;
        let ny = dz.abs() * 2.0 * hx + dx.abs() * 2.0 * hz;
        let nz = dx.abs() * 2.0 * hy + dy.abs() * 2.0 * hx;
        nx + ny + nz
    };
    box_sa + lateral_patch
}
/// Detect whether a swept sphere (capsule) self-intersects by checking
/// whether any two non-adjacent time steps overlap.
///
/// Samples the sweep path at `n_samples` discrete time steps and checks
/// whether the spheres at any two steps (with index difference > 1) overlap.
///
/// Returns `true` if a self-intersection is detected.
pub fn swept_sphere_detect_self_intersection(ss: &SweptSphere, n_samples: usize) -> bool {
    let n = n_samples.max(3);
    let centers: Vec<[f64; 3]> = (0..n)
        .map(|k| ss.center_at(k as f64 / (n - 1) as f64))
        .collect();
    let r = ss.radius;
    let threshold_sq = (2.0 * r) * (2.0 * r);
    for i in 0..centers.len() {
        for j in (i + 2)..centers.len() {
            let d = sub3(centers[i], centers[j]);
            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if d2 < threshold_sq {
                return true;
            }
        }
    }
    false
}
/// Compute an approximate envelope surface for a `SweptSphere` by sampling
/// the sweep and collecting the outermost surface points.
///
/// Returns a `Vec<[f64; 3]>` of points on the swept-sphere envelope sampled
/// at `n_time` time steps and `n_angle` angular positions around the sweep
/// axis.
pub fn swept_sphere_compute_envelope(
    ss: &SweptSphere,
    n_time: usize,
    n_angle: usize,
) -> Vec<[f64; 3]> {
    let nt = n_time.max(2);
    let na = n_angle.max(3);
    let dir = normalize3(ss.direction());
    let perp1 = {
        let t = if dir[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        normalize3(cross3(dir, t))
    };
    let perp2 = cross3(dir, perp1);
    let mut pts = Vec::with_capacity(nt * na);
    for ti in 0..nt {
        let t = ti as f64 / (nt - 1) as f64;
        let center = ss.center_at(t);
        for ai in 0..na {
            let angle = 2.0 * std::f64::consts::PI * ai as f64 / na as f64;
            let (s, c) = angle.sin_cos();
            let offset = add3(scale3(perp1, ss.radius * c), scale3(perp2, ss.radius * s));
            pts.push(add3(center, offset));
        }
    }
    pts
}
/// Compute the minimum clearance distance between a swept sphere and a
/// set of obstacle points.
///
/// Samples the sweep at `n_samples` time steps and computes the minimum
/// distance from any obstacle point to the nearest swept-sphere surface.
///
/// Returns `f64::INFINITY` if `obstacles` is empty.
pub fn swept_sphere_compute_clearance(
    ss: &SweptSphere,
    obstacles: &[[f64; 3]],
    n_samples: usize,
) -> f64 {
    if obstacles.is_empty() {
        return f64::INFINITY;
    }
    let n = n_samples.max(2);
    let mut min_clearance = f64::INFINITY;
    for &obs in obstacles {
        let pa = ss.center_start;
        let pb = ss.center_end;
        let ab = sub3(pb, pa);
        let ab_len2 = dot3(ab, ab);
        let pt_to_pa = sub3(obs, pa);
        let t_clamped = if ab_len2 < 1e-14 {
            0.0
        } else {
            (dot3(pt_to_pa, ab) / ab_len2).clamp(0.0, 1.0)
        };
        let closest = add3(pa, scale3(ab, t_clamped));
        let dist_to_axis = len3(sub3(obs, closest));
        let clearance = (dist_to_axis - ss.radius).max(0.0);
        if clearance < min_clearance {
            min_clearance = clearance;
        }
        let _ = n;
    }
    min_clearance
}
#[cfg(test)]
mod tests_extended {

    use crate::swept::LinearExtrusion;
    use crate::swept::RotationalSweep;
    use crate::swept::SweptAabb;

    use crate::swept::SweptObb;
    use crate::swept::SweptSphere;
    use crate::swept::ccd_aabb_vs_aabb;
    use crate::swept::ccd_sphere_vs_capsule;
    use crate::swept::swept_box_surface_area_estimate;
    use crate::swept::swept_sphere_compute_clearance;
    use crate::swept::swept_sphere_compute_envelope;
    use crate::swept::swept_sphere_detect_self_intersection;
    use crate::swept::swept_sphere_surface_area;
    fn square_profile() -> Vec<[f64; 2]> {
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    }
    #[test]
    fn test_linear_extrusion_profile_area() {
        let ext = LinearExtrusion::new(square_profile(), [0.0, 0.0, 2.0]);
        assert!(
            (ext.profile_area() - 1.0).abs() < 1e-12,
            "square area = 1.0"
        );
    }
    #[test]
    fn test_linear_extrusion_volume() {
        let ext = LinearExtrusion::new(square_profile(), [0.0, 0.0, 2.0]);
        assert!(
            (ext.volume() - 2.0).abs() < 1e-9,
            "volume = {}",
            ext.volume()
        );
    }
    #[test]
    fn test_linear_extrusion_surface_area() {
        let ext = LinearExtrusion::new(square_profile(), [0.0, 0.0, 2.0]);
        let sa = ext.surface_area();
        assert!((sa - 10.0).abs() < 1e-9, "surface area = {sa}");
    }
    #[test]
    fn test_linear_extrusion_aabb_contains_profile() {
        let ext = LinearExtrusion::new(square_profile(), [3.0, 0.0, 2.0]);
        let (mn, mx) = ext.aabb();
        assert!(mn[0] <= 0.0 + 1e-12);
        assert!(mx[0] >= 4.0 - 1e-12);
        assert!(mn[2] <= 0.0 + 1e-12);
        assert!(mx[2] >= 2.0 - 1e-12);
    }
    #[test]
    fn test_linear_extrusion_empty_profile_area_zero() {
        let ext = LinearExtrusion::new(vec![], [0.0, 0.0, 1.0]);
        assert!((ext.profile_area()).abs() < 1e-12);
    }
    #[test]
    fn test_rotational_sweep_aabb_cylinder() {
        let profile = vec![[1.0, 0.0], [1.0, 2.0]];
        let rs = RotationalSweep::new(profile, 16);
        let (mn, mx) = rs.aabb();
        assert!((mn[0] + 1.0).abs() < 1e-9);
        assert!((mx[0] - 1.0).abs() < 1e-9);
        assert!((mn[1]).abs() < 1e-9);
        assert!((mx[1] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_rotational_sweep_volume_cylinder() {
        let profile = vec![[1.0, 0.0], [1.0, 3.0]];
        let rs = RotationalSweep::new(profile, 64);
        let vol = rs.volume();
        let expected = std::f64::consts::PI * 3.0;
        assert!(
            (vol - expected).abs() < 1e-9,
            "vol={vol} expected={expected}"
        );
    }
    #[test]
    fn test_rotational_sweep_lateral_surface_area_cylinder() {
        let profile = vec![[1.0, 0.0], [1.0, 3.0]];
        let rs = RotationalSweep::new(profile, 64);
        let sa = rs.lateral_surface_area();
        let expected = 2.0 * std::f64::consts::PI * 3.0;
        assert!((sa - expected).abs() < 1e-9, "sa={sa} expected={expected}");
    }
    #[test]
    fn test_rotational_sweep_vertices_count() {
        let profile = vec![[1.0, 0.0], [1.0, 1.0], [0.5, 2.0]];
        let segs = 8;
        let rs = RotationalSweep::new(profile.clone(), segs);
        let verts = rs.vertices();
        assert_eq!(verts.len(), profile.len() * segs);
    }
    #[test]
    fn test_rotational_sweep_min_segments() {
        let rs = RotationalSweep::new(vec![[1.0, 0.0], [1.0, 1.0]], 1);
        assert!(rs.segments >= 3);
    }
    fn identity_axes() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }
    #[test]
    fn test_swept_obb_aabb_stationary() {
        let obb = SweptObb::new([0.0; 3], identity_axes(), [1.0, 1.0, 1.0], [0.0; 3]);
        let (mn, mx) = obb.aabb();
        for k in 0..3 {
            assert!((mn[k] + 1.0).abs() < 1e-12, "mn[{k}]={}", mn[k]);
            assert!((mx[k] - 1.0).abs() < 1e-12, "mx[{k}]={}", mx[k]);
        }
    }
    #[test]
    fn test_swept_obb_aabb_translated() {
        let obb = SweptObb::new([0.0; 3], identity_axes(), [1.0, 1.0, 1.0], [5.0, 0.0, 0.0]);
        let (mn, mx) = obb.aabb();
        assert!((mn[0] + 1.0).abs() < 1e-12, "mn.x should be -1");
        assert!((mx[0] - 6.0).abs() < 1e-12, "mx.x should be 6");
    }
    #[test]
    fn test_swept_obb_volume() {
        let obb = SweptObb::new([0.0; 3], identity_axes(), [2.0, 3.0, 4.0], [0.0; 3]);
        assert!((obb.volume() - 192.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_obb_support_plus_x() {
        let obb = SweptObb::new([0.0; 3], identity_axes(), [1.5, 1.0, 1.0], [0.0; 3]);
        let sp = obb.support([1.0, 0.0, 0.0]);
        assert!((sp[0] - 1.5).abs() < 1e-12, "support in +X should be 1.5");
    }
    #[test]
    fn test_swept_aabb_union_aabb() {
        let sa = SweptAabb::new(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [4.0, -1.0, -1.0],
            [6.0, 1.0, 1.0],
        );
        let (mn, mx) = sa.aabb();
        assert!((mn[0] + 1.0).abs() < 1e-12);
        assert!((mx[0] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_swept_aabb_contains_point() {
        let sa = SweptAabb::new([-1.0; 3], [1.0; 3], [4.0, -1.0, -1.0], [6.0, 1.0, 1.0]);
        assert!(sa.contains_point([0.0, 0.0, 0.0]));
        assert!(sa.contains_point([5.0, 0.0, 0.0]));
        assert!(!sa.contains_point([10.0, 0.0, 0.0]));
    }
    #[test]
    fn test_swept_aabb_displacement() {
        let sa = SweptAabb::new([-1.0; 3], [1.0; 3], [3.0, -1.0, -1.0], [5.0, 1.0, 1.0]);
        let d = sa.displacement();
        assert!((d[0] - 4.0).abs() < 1e-12, "displacement x = {}", d[0]);
    }
    #[test]
    fn test_ccd_sphere_vs_capsule_hit() {
        let result = ccd_sphere_vs_capsule(
            [-5.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            0.5,
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(result.is_some(), "expected CCD hit");
        let r = result.unwrap();
        assert!(r.toi >= 0.0 && r.toi <= 1.0, "TOI in [0,1]");
    }
    #[test]
    fn test_ccd_sphere_vs_capsule_miss() {
        let result = ccd_sphere_vs_capsule(
            [5.0, -10.0, 0.0],
            [5.0, 10.0, 0.0],
            0.5,
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.5,
        );
        assert!(result.is_none(), "expected miss");
    }
    #[test]
    fn test_ccd_sphere_vs_capsule_already_overlapping() {
        let result = ccd_sphere_vs_capsule(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            [0.0, -0.5, 0.0],
            [0.0, 0.5, 0.0],
            1.0,
        );
        assert!(result.is_some());
        assert!(result.unwrap().toi < 1e-12, "already overlapping => TOI=0");
    }
    #[test]
    fn test_ccd_aabb_vs_aabb_hit() {
        let toi = ccd_aabb_vs_aabb(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [10.0, 0.0, 0.0],
            [4.0, -1.0, -1.0],
            [6.0, 1.0, 1.0],
        );
        assert!(toi.is_some(), "expected CCD AABB hit");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "TOI in [0,1]");
    }
    #[test]
    fn test_ccd_aabb_vs_aabb_miss() {
        let toi = ccd_aabb_vs_aabb(
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.0, 0.0],
            [10.0, -1.0, -1.0],
            [12.0, 1.0, 1.0],
        );
        assert!(toi.is_none(), "expected miss");
    }
    #[test]
    fn test_ccd_aabb_vs_aabb_already_overlapping() {
        let toi = ccd_aabb_vs_aabb(
            [-1.0; 3],
            [1.0; 3],
            [0.5, 0.0, 0.0],
            [-0.5, -1.0, -1.0],
            [0.5, 1.0, 1.0],
        );
        assert!(toi.is_some());
        assert!(toi.unwrap() < 1e-9, "already overlapping => TOI=0");
    }
    #[test]
    fn test_swept_sphere_surface_area_sphere() {
        let sa = swept_sphere_surface_area(1.0, 0.0);
        assert!((sa - 4.0 * std::f64::consts::PI).abs() < 1e-9, "sa={sa}");
    }
    #[test]
    fn test_swept_sphere_surface_area_capsule() {
        let sa = swept_sphere_surface_area(1.0, 2.0);
        let expected = 8.0 * std::f64::consts::PI;
        assert!((sa - expected).abs() < 1e-9, "sa={sa}");
    }
    #[test]
    fn test_swept_box_surface_area_estimate_static() {
        let sa = swept_box_surface_area_estimate([1.0, 1.0, 1.0], [0.0; 3]);
        assert!((sa - 24.0).abs() < 1e-9, "sa={sa}");
    }
    #[test]
    fn test_swept_box_surface_area_estimate_positive() {
        let sa = swept_box_surface_area_estimate([1.0, 1.0, 1.0], [5.0, 0.0, 0.0]);
        assert!(sa > 24.0, "swept box SA should be > static SA");
    }
    #[test]
    fn test_swept_self_intersection_long_capsule_no_self_intersect() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [100.0, 0.0, 0.0], 1.0);
        let has_si = swept_sphere_detect_self_intersection(&ss, 4);
        assert!(
            !has_si,
            "long capsule with few samples should not self-intersect"
        );
    }
    #[test]
    fn test_swept_self_intersection_stationary_sphere() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0);
        let has_si = swept_sphere_detect_self_intersection(&ss, 5);
        assert!(
            has_si,
            "zero-length sweep (sphere) should detect self-intersection"
        );
    }
    #[test]
    fn test_swept_self_intersection_very_short_long_radius() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], 2.0);
        let has_si = swept_sphere_detect_self_intersection(&ss, 10);
        assert!(
            has_si,
            "short sweep with large radius should self-intersect"
        );
    }
    #[test]
    fn test_swept_envelope_point_count() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0);
        let env = swept_sphere_compute_envelope(&ss, 8, 12);
        assert_eq!(
            env.len(),
            8 * 12,
            "envelope should have n_time * n_angle points"
        );
    }
    #[test]
    fn test_swept_envelope_points_finite() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], 1.0);
        let env = swept_sphere_compute_envelope(&ss, 6, 8);
        for (i, p) in env.iter().enumerate() {
            for (k, &pk) in p.iter().enumerate() {
                assert!(pk.is_finite(), "envelope[{i}][{k}] not finite: {}", pk);
            }
        }
    }
    #[test]
    fn test_swept_envelope_radius_correct() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0);
        let env = swept_sphere_compute_envelope(&ss, 4, 16);
        for p in &env {
            let r = (p[1].powi(2) + p[2].powi(2)).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-6,
                "envelope point radial distance should be ~1.0, got {r}"
            );
        }
    }
    #[test]
    fn test_swept_clearance_empty_obstacles_is_infinity() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0);
        let cl = swept_sphere_compute_clearance(&ss, &[], 10);
        assert!(
            cl.is_infinite(),
            "clearance with no obstacles should be +inf"
        );
    }
    #[test]
    fn test_swept_clearance_far_obstacle() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0);
        let obstacles = [[2.5, 100.0, 0.0]];
        let cl = swept_sphere_compute_clearance(&ss, &obstacles, 10);
        assert!(
            cl > 0.0,
            "clearance to far obstacle should be positive: {cl}"
        );
        assert!(
            cl < f64::INFINITY,
            "clearance to finite obstacle should be finite: {cl}"
        );
    }
    #[test]
    fn test_swept_clearance_obstacle_on_surface_near_zero() {
        let ss = SweptSphere::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0);
        let obstacles = [[2.5, 1.0, 0.0]];
        let cl = swept_sphere_compute_clearance(&ss, &obstacles, 10);
        assert!(
            cl.abs() < 1e-6,
            "obstacle on surface should give clearance ~0, got {cl}"
        );
    }
}
