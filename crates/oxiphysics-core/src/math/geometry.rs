//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::linear_algebra::*;
#[cfg(test)]
use super::types::*;

/// Extract ZYX Euler angles (yaw, pitch, roll) from a rotation matrix.
///
/// Returns `(roll, pitch, yaw)` in radians, same convention as
/// [`quat_to_euler`].
///
/// Assumes the matrix is a proper rotation (no scaling / reflection).
pub fn mat3_to_euler_zyx(m: &Mat3) -> (Real, Real, Real) {
    let pitch = (-m[(2, 0)]).clamp(-1.0, 1.0).asin();
    let cos_pitch = pitch.cos();
    if cos_pitch.abs() < 1e-6 {
        let roll = 0.0;
        let yaw = m[(0, 1)].atan2(m[(1, 1)]);
        (roll, pitch, yaw)
    } else {
        let roll = m[(2, 1)].atan2(m[(2, 2)]);
        let yaw = m[(1, 0)].atan2(m[(0, 0)]);
        (roll, pitch, yaw)
    }
}
/// Build a rotation matrix from ZYX Euler angles (roll, pitch, yaw) in radians.
///
/// Equivalent to `Rz(yaw) * Ry(pitch) * Rx(roll)`.
pub fn mat3_from_euler_zyx(roll: Real, pitch: Real, yaw: Real) -> Mat3 {
    let (cr, sr) = (roll.cos(), roll.sin());
    let (cp, sp) = (pitch.cos(), pitch.sin());
    let (cy, sy) = (yaw.cos(), yaw.sin());
    Mat3::new(
        cy * cp,
        cy * sp * sr - sy * cr,
        cy * sp * cr + sy * sr,
        sy * cp,
        sy * sp * sr + cy * cr,
        sy * sp * cr - cy * sr,
        -sp,
        cp * sr,
        cp * cr,
    )
}
/// Cofactor matrix of a 3×3 matrix.
///
/// `C[i][j] = (-1)^{i+j} * M_{ij}` where `M_{ij}` is the (i,j) minor.
pub fn mat3_cofactor(m: &Mat3) -> Mat3 {
    let c = |r: usize, c_idx: usize| -> Real {
        let rows: [usize; 2] = match r {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let cols: [usize; 2] = match c_idx {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let det2 = m[(rows[0], cols[0])] * m[(rows[1], cols[1])]
            - m[(rows[0], cols[1])] * m[(rows[1], cols[0])];
        let sign = if (r + c_idx).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign * det2
    };
    Mat3::new(
        c(0, 0),
        c(0, 1),
        c(0, 2),
        c(1, 0),
        c(1, 1),
        c(1, 2),
        c(2, 0),
        c(2, 1),
        c(2, 2),
    )
}
/// Adjugate of a nalgebra `Mat3` (transpose of the cofactor matrix).
///
/// `adj(M) = cofactor(M)^T = det(M) * M^{-1}` when M is invertible.
pub fn mat3_adjugate_na(m: &Mat3) -> Mat3 {
    mat3_cofactor(m).transpose()
}
/// Build a unit quaternion that rotates `from` onto `to`.
///
/// Both `from` and `to` must be unit vectors.  Returns the identity for
/// anti-parallel inputs (180° rotation around a perpendicular axis).
pub fn quat_from_two_vectors(from: &Vec3, to: &Vec3) -> Quat {
    let dot = from.dot(to).clamp(-1.0, 1.0);
    if dot > 1.0 - 1e-10 {
        return Quat::identity();
    }
    if dot < -1.0 + 1e-10 {
        let perp = if from.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0).cross(from)
        } else {
            Vec3::new(0.0, 1.0, 0.0).cross(from)
        };
        let axis = Unit::new_normalize(perp);
        return UnitQuaternion::from_axis_angle(&axis, std::f64::consts::PI);
    }
    let axis = from.cross(to);
    let w = 1.0 + dot;
    let len = (w * w + axis.norm_squared()).sqrt();
    let [x, y, z] = [axis.x / len, axis.y / len, axis.z / len];
    UnitQuaternion::new_normalize(nalgebra::Quaternion::new(w / len, x, y, z))
}
/// Quaternion logarithm of a `UnitQuaternion`.
///
/// Returns a pure quaternion `[θ*nx, θ*ny, θ*nz, 0]` where `θ*n̂` is the
/// rotation axis-angle.
pub fn quat_log(q: &Quat) -> nalgebra::Quaternion<Real> {
    let w = q.w.clamp(-1.0, 1.0);
    let theta = w.acos();
    let sin_t = theta.sin();
    if sin_t.abs() < 1e-10 {
        return nalgebra::Quaternion::new(0.0, 0.0, 0.0, 0.0);
    }
    let s = theta / sin_t;
    nalgebra::Quaternion::new(0.0, q.i * s, q.j * s, q.k * s)
}
/// Quaternion exponential of a pure quaternion `v` (w component ignored).
///
/// Returns a unit quaternion `exp(v) = [sin(|v|)*v̂, cos(|v|)]`.
pub fn quat_exp(v: &nalgebra::Quaternion<Real>) -> Quat {
    let norm = (v.i * v.i + v.j * v.j + v.k * v.k).sqrt();
    if norm < 1e-15 {
        return Quat::identity();
    }
    let s = norm.sin() / norm;
    UnitQuaternion::new_normalize(nalgebra::Quaternion::new(
        norm.cos(),
        v.i * s,
        v.j * s,
        v.k * s,
    ))
}
/// Ensure double-cover consistency: negate `q` if its w < 0 so that the
/// canonical representative always has w ≥ 0.
pub fn quat_double_cover_fix(q: Quat) -> Quat {
    if q.w < 0.0 {
        UnitQuaternion::new_normalize(-q.into_inner())
    } else {
        q
    }
}
/// Spherical cubic interpolation (SQUAD) between two `UnitQuaternion`s.
///
/// Interpolates between `q1` and `q2` at `t ∈ [0,1]` using the
/// inner control quaternions `s1` and `s2` for C1 continuity.
pub fn quat_squad_na(q1: &Quat, q2: &Quat, s1: &Quat, s2: &Quat, t: Real) -> Quat {
    let slerp_q = q1.slerp(q2, t);
    let slerp_s = s1.slerp(s2, t);
    slerp_q.slerp(&slerp_s, 2.0 * t * (1.0 - t))
}
#[cfg(test)]
mod extended_math_tests {
    use super::*;
    use crate::Mat3;
    use crate::Vec3;

    use crate::math::mat3_from_axis_angle;

    use crate::math::quat_from_axis_angle;

    use crate::math::vec3_refract;

    #[test]
    fn test_vec3_project_onto_x_axis() {
        let v = Vec3::new(3.0, 4.0, 5.0);
        let x = Vec3::new(1.0, 0.0, 0.0);
        let p = vec3_project_onto(&v, &x);
        assert!((p.x - 3.0).abs() < 1e-12, "x={}", p.x);
        assert!(p.y.abs() < 1e-12);
        assert!(p.z.abs() < 1e-12);
    }
    #[test]
    fn test_vec3_reject_from_perpendicular() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let x = Vec3::new(1.0, 0.0, 0.0);
        let r = vec3_reject_from(&v, &x);
        assert!(r.x.abs() < 1e-12, "x={}", r.x);
        assert!((r.y - 4.0).abs() < 1e-12, "y={}", r.y);
    }
    #[test]
    fn test_vec3_reflect_about_normal() {
        let v = Vec3::new(1.0, -1.0, 0.0);
        let n = Vec3::new(0.0, 1.0, 0.0);
        let r = vec3_reflect_about(&v, &n);
        assert!((r.x - 1.0).abs() < 1e-12);
        assert!((r.y - 1.0).abs() < 1e-12);
        assert!(r.z.abs() < 1e-12);
    }
    #[test]
    fn test_vec3_refract_no_total_internal_reflection() {
        let v = Vec3::new(0.0, -1.0, 0.0);
        let n = Vec3::new(0.0, 1.0, 0.0);
        let refr = vec3_refract(&v, &n, 1.0).expect("should not TIR");
        assert!((refr.x).abs() < 1e-10);
        assert!((refr.y + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec3_refract_total_internal_reflection() {
        let v = Vec3::new(0.9999, -0.0141, 0.0).normalize();
        let n = Vec3::new(0.0, 1.0, 0.0);
        let result = vec3_refract(&v, &n, 1.5);
        if let Some(r) = result {
            assert!(
                (r.norm() - 1.0).abs() < 1e-8,
                "refracted not unit: {}",
                r.norm()
            );
        }
    }
    #[test]
    fn test_vec3_angle_between_perpendicular() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        let angle = vec3_angle_between(&a, &b);
        assert!(
            (angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "angle={}",
            angle
        );
    }
    #[test]
    fn test_vec3_angle_between_parallel() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let angle = vec3_angle_between(&a, &a);
        assert!(angle.abs() < 1e-10, "angle={}", angle);
    }
    #[test]
    fn test_vec3_rotate_by_quat() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let q = quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
        let r = vec3_rotate_by_quat(&v, &q);
        assert!(r.x.abs() < 1e-10, "x={}", r.x);
        assert!((r.y - 1.0).abs() < 1e-10, "y={}", r.y);
        assert!(r.z.abs() < 1e-10, "z={}", r.z);
    }
    #[test]
    fn test_mat3_from_axis_angle_identity_at_zero() {
        let m = mat3_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), 0.0);
        let diff = (m - Mat3::identity()).norm();
        assert!(diff < 1e-10, "diff={}", diff);
    }
    #[test]
    fn test_mat3_from_axis_angle_90z() {
        let m = mat3_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
        let v = m * Vec3::new(1.0, 0.0, 0.0);
        assert!(v.x.abs() < 1e-10, "x={}", v.x);
        assert!((v.y - 1.0).abs() < 1e-10, "y={}", v.y);
        assert!(v.z.abs() < 1e-10, "z={}", v.z);
    }
    #[test]
    fn test_mat3_euler_roundtrip() {
        let roll = 0.4_f64;
        let pitch = 0.2_f64;
        let yaw = -0.3_f64;
        let m = mat3_from_euler_zyx(roll, pitch, yaw);
        let (r2, p2, y2) = mat3_to_euler_zyx(&m);
        assert!((r2 - roll).abs() < 1e-8, "roll: {} vs {}", r2, roll);
        assert!((p2 - pitch).abs() < 1e-8, "pitch: {} vs {}", p2, pitch);
        assert!((y2 - yaw).abs() < 1e-8, "yaw: {} vs {}", y2, yaw);
    }
    #[test]
    fn test_mat3_cofactor_identity() {
        let m = Mat3::identity();
        let c = mat3_cofactor(&m);
        let diff = (c - Mat3::identity()).norm();
        assert!(diff < 1e-10, "cofactor(I) should be I: diff={}", diff);
    }
    #[test]
    fn test_mat3_adjugate_na_relation() {
        let m = Mat3::new(1.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 1.0);
        let adj = mat3_adjugate_na(&m);
        let det = m.determinant();
        let prod = m * adj;
        let expected = Mat3::identity() * det;
        let diff = (prod - expected).norm();
        assert!(diff < 1e-8, "M*adj(M) diff={}", diff);
    }
    #[test]
    fn test_quat_from_two_vectors_basic() {
        let from = Vec3::new(1.0, 0.0, 0.0);
        let to = Vec3::new(0.0, 1.0, 0.0);
        let q = quat_from_two_vectors(&from, &to);
        let rotated = q.transform_vector(&from);
        assert!((rotated.x).abs() < 1e-8, "x={}", rotated.x);
        assert!((rotated.y - 1.0).abs() < 1e-8, "y={}", rotated.y);
    }
    #[test]
    fn test_quat_from_two_vectors_same_direction() {
        let from = Vec3::new(0.0, 0.0, 1.0);
        let q = quat_from_two_vectors(&from, &from);
        let rotated = q.transform_vector(&from);
        assert!((rotated - from).norm() < 1e-8);
    }
    #[test]
    fn test_quat_log_exp_roundtrip_na() {
        let q = quat_from_axis_angle(&Vec3::new(1.0, 0.0, 0.0), 1.2);
        let log_q = quat_log(&q);
        let exp_q = quat_exp(&log_q);
        let diff = q.angle_to(&exp_q);
        assert!(diff < 1e-8, "log/exp roundtrip angle diff={}", diff);
    }
    #[test]
    fn test_quat_double_cover_fix_positive_w() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.5);
        let fixed = quat_double_cover_fix(q);
        assert!(fixed.w >= 0.0, "w should be non-negative after fix");
    }
    #[test]
    fn test_quat_double_cover_fix_negative_w() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.5);
        let neg = UnitQuaternion::new_normalize(-q.into_inner());
        let fixed = quat_double_cover_fix(neg);
        assert!(fixed.w >= 0.0, "w should be non-negative after fix");
    }
    #[test]
    fn test_quat_squad_na_endpoints() {
        let q1 = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let q2 = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let r0 = quat_squad_na(&q1, &q2, &q1, &q2, 0.0);
        let r1 = quat_squad_na(&q1, &q2, &q1, &q2, 1.0);
        assert!(r0.angle_to(&q1) < 1e-6, "t=0 should be q1");
        assert!(r1.angle_to(&q2) < 1e-6, "t=1 should be q2");
    }
    #[test]
    fn test_plane_intersect_ray_method() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let ray = Ray::new(Vec3::new(0.0, 2.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let t = plane.intersect_ray(&ray).expect("should intersect");
        assert!((t - 2.0).abs() < 1e-10, "t={}", t);
    }
    #[test]
    fn test_plane_intersect_ray_parallel() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let ray = Ray::new(Vec3::new(0.0, 1.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert!(plane.intersect_ray(&ray).is_none());
    }
    #[test]
    fn test_frustum_contains_sphere_small() {
        let vp = perspective(std::f64::consts::FRAC_PI_2, 1.0, 1.0, 100.0);
        let frustum = Frustum::from_view_projection(&vp);
        let _result = frustum.contains_sphere(&Vec3::new(0.0, 0.0, 0.0), 0.01);
    }
    #[test]
    fn test_frustum_intersects_aabb_basic() {
        let vp = perspective(std::f64::consts::FRAC_PI_2, 1.0, 1.0, 100.0);
        let frustum = Frustum::from_view_projection(&vp);
        let aabb_far = Aabb::new([1000.0, 1000.0, 1000.0], [1001.0, 1001.0, 1001.0]);
        let _ = frustum.intersects_aabb(&aabb_far);
    }
    #[test]
    fn test_frustum_extract_from_view_proj_same_as_constructor() {
        let vp = perspective(std::f64::consts::FRAC_PI_4, 1.0, 0.5, 50.0);
        let f1 = Frustum::from_view_projection(&vp);
        let f2 = Frustum::extract_from_view_proj(&vp);
        for (p1, p2) in f1.planes.iter().zip(f2.planes.iter()) {
            assert!((p1.normal - p2.normal).norm() < 1e-10);
            assert!((p1.distance - p2.distance).abs() < 1e-10);
        }
    }
    #[test]
    fn test_aabb_merge() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([-1.0, 2.0, -1.0], [0.5, 3.0, 2.0]);
        let m = a.merge(&b);
        assert!((m.min[0] - (-1.0)).abs() < 1e-12);
        assert!((m.min[1] - 0.0).abs() < 1e-12);
        assert!((m.max[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_expand_by() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let e = a.expand_by(1.0);
        for i in 0..3 {
            assert!((e.min[i] - (-1.0)).abs() < 1e-12);
            assert!((e.max[i] - 3.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_aabb_closest_point_inside() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
        let p = [2.0, 2.0, 2.0];
        let cp = aabb.closest_point(p);
        for &cp_i in cp.iter() {
            assert!((cp_i - 2.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_aabb_closest_point_outside() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let p = [-1.0, 0.5, 2.0];
        let cp = aabb.closest_point(p);
        assert!((cp[0] - 0.0).abs() < 1e-12, "x={}", cp[0]);
        assert!((cp[1] - 0.5).abs() < 1e-12, "y={}", cp[1]);
        assert!((cp[2] - 1.0).abs() < 1e-12, "z={}", cp[2]);
    }
    #[test]
    fn test_aabb_contains_point_strict() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(aabb.contains_point_strict([1.0, 1.0, 1.0]));
        assert!(!aabb.contains_point_strict([0.0, 1.0, 1.0]));
        assert!(!aabb.contains_point_strict([3.0, 1.0, 1.0]));
    }
    #[test]
    fn test_aabb_contains_point_on_boundary() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(aabb.contains([0.0, 0.5, 0.5]));
        assert!(aabb.contains([1.0, 0.5, 0.5]));
    }
    #[test]
    fn test_aabb_ray_intersect_extended() {
        let aabb = Aabb::new([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        let hit = aabb.intersect_ray([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_some(), "should hit unit cube");
        let (t0, t1) = hit.unwrap();
        assert!((t0 - 4.0).abs() < 1e-10, "t0={}", t0);
        assert!((t1 - 6.0).abs() < 1e-10, "t1={}", t1);
    }
}
#[cfg(test)]
mod new_math_tests {

    use crate::Mat3;
    use crate::Vec3;

    use crate::math::cartesian_to_cylindrical;
    use crate::math::cartesian_to_spherical;
    use crate::math::cross_matrix;
    use crate::math::cylindrical_to_cartesian;

    use crate::math::gram_schmidt;

    use crate::math::mat3_exp_skew;
    use crate::math::mat3_from_axis_angle;

    use crate::math::mat3_log_rotation;

    use crate::math::quat_from_axis_angle;
    use crate::math::quat_geodesic_distance;

    use crate::math::quat_nlerp;
    use crate::math::quat_squad_control;
    use crate::math::rodrigues_rotate;

    use crate::math::skew_symmetric;
    use crate::math::spherical_to_cartesian;

    use crate::math::vec3_outer_product;

    #[test]
    fn test_vec3_outer_product_basic() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let op = vec3_outer_product(&a, &b);
        assert!((op[(0, 0)] - 4.0).abs() < 1e-12);
        assert!((op[(0, 1)] - 5.0).abs() < 1e-12);
        assert!((op[(0, 2)] - 6.0).abs() < 1e-12);
        assert!((op[(1, 0)] - 8.0).abs() < 1e-12);
        assert!((op[(2, 2)] - 18.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_outer_product_unit_vectors() {
        let e0 = Vec3::new(1.0, 0.0, 0.0);
        let e1 = Vec3::new(0.0, 1.0, 0.0);
        let op = vec3_outer_product(&e0, &e1);
        assert!((op[(0, 1)] - 1.0).abs() < 1e-12);
        assert!(op[(0, 0)].abs() < 1e-12);
        assert!(op[(1, 0)].abs() < 1e-12);
    }
    #[test]
    fn test_vec3_outer_product_trace() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let op = vec3_outer_product(&a, &a);
        let trace = op.trace();
        let norm_sq = a.norm_squared();
        assert!(
            (trace - norm_sq).abs() < 1e-10,
            "trace={}, |a|²={}",
            trace,
            norm_sq
        );
    }
    #[test]
    fn test_cross_matrix_equals_skew_symmetric() {
        let v = Vec3::new(2.0, -3.0, 1.0);
        let cm = cross_matrix(&v);
        let ss = skew_symmetric(&v);
        let diff = (cm - ss).norm();
        assert!(
            diff < 1e-12,
            "cross_matrix should equal skew_symmetric: diff={}",
            diff
        );
    }
    #[test]
    fn test_cross_matrix_cross_product_equivalence() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-1.0, 0.5, 4.0);
        let expected = a.cross(&b);
        let actual = cross_matrix(&a) * b;
        let diff = (expected - actual).norm();
        assert!(
            diff < 1e-10,
            "cross matrix * b should equal a x b: diff={}",
            diff
        );
    }
    #[test]
    fn test_gram_schmidt_orthonormality() {
        let v0 = Vec3::new(1.0, 1.0, 0.0);
        let v1 = Vec3::new(0.0, 1.0, 1.0);
        let v2 = Vec3::new(1.0, 0.0, 1.0);
        let (e0, e1, e2) = gram_schmidt(&v0, &v1, &v2).expect("should succeed");
        assert!(
            (e0.norm() - 1.0).abs() < 1e-10,
            "e0 not unit: {}",
            e0.norm()
        );
        assert!(
            (e1.norm() - 1.0).abs() < 1e-10,
            "e1 not unit: {}",
            e1.norm()
        );
        assert!(
            (e2.norm() - 1.0).abs() < 1e-10,
            "e2 not unit: {}",
            e2.norm()
        );
        assert!(e0.dot(&e1).abs() < 1e-10, "e0·e1={}", e0.dot(&e1));
        assert!(e0.dot(&e2).abs() < 1e-10, "e0·e2={}", e0.dot(&e2));
        assert!(e1.dot(&e2).abs() < 1e-10, "e1·e2={}", e1.dot(&e2));
    }
    #[test]
    fn test_gram_schmidt_span_preserved() {
        let v0 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 1.0, 0.0);
        let v2 = Vec3::new(1.0, 1.0, 1.0);
        let (e0, e1, e2) = gram_schmidt(&v0, &v1, &v2).expect("should succeed");
        assert!((e0.x - 1.0).abs() < 1e-10);
        assert!(e0.y.abs() < 1e-10);
        assert!(e2.z.abs() > 0.99, "e2.z={}", e2.z);
        let _ = e1;
    }
    #[test]
    fn test_gram_schmidt_degenerate_returns_none() {
        let v0 = Vec3::new(1.0, 0.0, 0.0);
        let v1 = Vec3::new(2.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);
        let result = gram_schmidt(&v0, &v1, &v2);
        assert!(result.is_none(), "collinear vectors should return None");
    }
    #[test]
    fn test_cartesian_to_spherical_z_axis() {
        let v = Vec3::new(0.0, 0.0, 2.0);
        let (r, theta, _phi) = cartesian_to_spherical(&v);
        assert!((r - 2.0).abs() < 1e-10, "r={}", r);
        assert!(theta.abs() < 1e-10, "theta={}", theta);
    }
    #[test]
    fn test_cartesian_to_spherical_x_axis() {
        let v = Vec3::new(3.0, 0.0, 0.0);
        let (r, theta, phi) = cartesian_to_spherical(&v);
        assert!((r - 3.0).abs() < 1e-10, "r={}", r);
        assert!(
            (theta - std::f64::consts::FRAC_PI_2).abs() < 1e-10,
            "theta={}",
            theta
        );
        assert!(phi.abs() < 1e-10, "phi={}", phi);
    }
    #[test]
    fn test_spherical_cartesian_roundtrip() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let (r, theta, phi) = cartesian_to_spherical(&v);
        let v2 = spherical_to_cartesian(r, theta, phi);
        let diff = (v - v2).norm();
        assert!(diff < 1e-10, "roundtrip diff={}", diff);
    }
    #[test]
    fn test_spherical_to_cartesian_unit_radius() {
        let v = spherical_to_cartesian(1.0, std::f64::consts::FRAC_PI_2, 0.0);
        assert!((v.x - 1.0).abs() < 1e-10);
        assert!(v.y.abs() < 1e-10);
        assert!(v.z.abs() < 1e-10);
    }
    #[test]
    fn test_cylindrical_roundtrip() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let (rho, phi, z) = cartesian_to_cylindrical(&v);
        let v2 = cylindrical_to_cartesian(rho, phi, z);
        let diff = (v - v2).norm();
        assert!(diff < 1e-10, "cylindrical roundtrip diff={}", diff);
    }
    #[test]
    fn test_cylindrical_z_preserved() {
        let v = Vec3::new(1.0, 1.0, 5.0);
        let (_, _, z) = cartesian_to_cylindrical(&v);
        assert!((z - 5.0).abs() < 1e-10, "z={}", z);
    }
    #[test]
    fn test_rodrigues_rotate_90z() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let r = rodrigues_rotate(&v, &axis, std::f64::consts::FRAC_PI_2);
        assert!(r.x.abs() < 1e-10, "x={}", r.x);
        assert!((r.y - 1.0).abs() < 1e-10, "y={}", r.y);
        assert!(r.z.abs() < 1e-10, "z={}", r.z);
    }
    #[test]
    fn test_rodrigues_rotate_180() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let r = rodrigues_rotate(&v, &axis, std::f64::consts::PI);
        assert!((r.x + 1.0).abs() < 1e-10, "x={}", r.x);
        assert!(r.y.abs() < 1e-10, "y={}", r.y);
        assert!(r.z.abs() < 1e-10, "z={}", r.z);
    }
    #[test]
    fn test_rodrigues_rotate_zero_angle() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let r = rodrigues_rotate(&v, &axis, 0.0);
        let diff = (r - v).norm();
        assert!(
            diff < 1e-10,
            "zero rotation should leave vector unchanged: diff={}",
            diff
        );
    }
    #[test]
    fn test_rodrigues_rotate_preserves_magnitude() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let r = rodrigues_rotate(&v, &axis, 1.234);
        let diff = (r.norm() - v.norm()).abs();
        assert!(
            diff < 1e-10,
            "rotation should preserve magnitude: diff={}",
            diff
        );
    }
    #[test]
    fn test_mat3_exp_skew_zero() {
        let omega = Vec3::zeros();
        let r = mat3_exp_skew(&omega);
        let diff = (r - Mat3::identity()).norm();
        assert!(diff < 1e-10, "exp(0) should be identity: diff={}", diff);
    }
    #[test]
    fn test_mat3_exp_skew_90z() {
        let omega = Vec3::new(0.0, 0.0, std::f64::consts::FRAC_PI_2);
        let r = mat3_exp_skew(&omega);
        let v = r * Vec3::new(1.0, 0.0, 0.0);
        assert!(v.x.abs() < 1e-10, "x={}", v.x);
        assert!((v.y - 1.0).abs() < 1e-10, "y={}", v.y);
    }
    #[test]
    fn test_mat3_exp_skew_is_rotation() {
        let omega = Vec3::new(0.3, -0.5, 0.7);
        let r = mat3_exp_skew(&omega);
        let det = r.determinant();
        let rrt = r * r.transpose();
        let diff = (rrt - Mat3::identity()).norm();
        assert!((det - 1.0).abs() < 1e-8, "det={}", det);
        assert!(diff < 1e-8, "R R^T should be I: diff={}", diff);
    }
    #[test]
    fn test_mat3_log_rotation_identity() {
        let r = Mat3::identity();
        let log = mat3_log_rotation(&r);
        let diff = log.norm();
        assert!(diff < 1e-10, "log(I) should be zero: norm={}", diff);
    }
    #[test]
    fn test_mat3_exp_log_roundtrip() {
        let omega = Vec3::new(0.2, 0.4, -0.3);
        let r = mat3_exp_skew(&omega);
        let log_r = mat3_log_rotation(&r);
        let omega2 = Vec3::new(log_r[(2, 1)], log_r[(0, 2)], log_r[(1, 0)]);
        let r2 = mat3_exp_skew(&omega2);
        let diff = (r - r2).norm();
        assert!(diff < 1e-7, "exp/log roundtrip diff={}", diff);
    }
    #[test]
    fn test_quat_nlerp_at_endpoints() {
        let a = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let q0 = quat_nlerp(&a, &b, 0.0);
        let q1 = quat_nlerp(&a, &b, 1.0);
        assert!(q0.angle_to(&a) < 1e-8, "nlerp(t=0) should be a");
        assert!(q1.angle_to(&b) < 1e-8, "nlerp(t=1) should be b");
    }
    #[test]
    fn test_quat_nlerp_is_unit() {
        let a = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.3);
        let b = quat_from_axis_angle(&Vec3::new(1.0, 0.0, 0.0), 0.9);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let q = quat_nlerp(&a, &b, t);
            let norm = q.into_inner().norm();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "nlerp result not unit at t={}: norm={}",
                t,
                norm
            );
        }
    }
    #[test]
    fn test_quat_squad_control_produces_unit_quat() {
        let q0 = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let q1 = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.5);
        let q2 = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let s = quat_squad_control(&q0, &q1, &q2);
        let norm = s.into_inner().norm();
        assert!(
            (norm - 1.0).abs() < 1e-8,
            "squad control should be unit quat: norm={}",
            norm
        );
    }
    #[test]
    fn test_quat_geodesic_distance_same() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.7);
        let d = quat_geodesic_distance(&q, &q);
        assert!(
            d.abs() < 1e-10,
            "geodesic distance to self should be 0: d={}",
            d
        );
    }
    #[test]
    fn test_quat_geodesic_distance_known() {
        let a = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2);
        let d = quat_geodesic_distance(&a, &b);
        assert!((d - std::f64::consts::FRAC_PI_2).abs() < 1e-8, "d={}", d);
    }
    #[test]
    fn test_outer_product_rank_one() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(1.0, -1.0, 2.0);
        let op = vec3_outer_product(&a, &b);
        let det = op.determinant();
        assert!(
            det.abs() < 1e-10,
            "outer product should have det=0: det={}",
            det
        );
    }
    #[test]
    fn test_cross_matrix_antisymmetric() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        let cm = cross_matrix(&v);
        let diff = (cm + cm.transpose()).norm();
        assert!(
            diff < 1e-10,
            "cross matrix should be antisymmetric: diff={}",
            diff
        );
    }
    #[test]
    fn test_cross_matrix_zero_diagonal() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let cm = cross_matrix(&v);
        assert!(cm[(0, 0)].abs() < 1e-12);
        assert!(cm[(1, 1)].abs() < 1e-12);
        assert!(cm[(2, 2)].abs() < 1e-12);
    }
    #[test]
    fn test_cartesian_to_spherical_origin() {
        let v = Vec3::zeros();
        let (r, theta, phi) = cartesian_to_spherical(&v);
        assert!(r.abs() < 1e-10, "r should be 0 for origin");
        assert!(theta.abs() < 1e-10);
        assert!(phi.abs() < 1e-10);
    }
    #[test]
    fn test_spherical_coordinate_r_is_norm() {
        let v = Vec3::new(2.0, 3.0, 6.0);
        let (r, _, _) = cartesian_to_spherical(&v);
        assert!((r - v.norm()).abs() < 1e-10, "r should equal |v|");
    }
    #[test]
    fn test_gram_schmidt_axis_aligned() {
        let e0 = Vec3::new(1.0, 0.0, 0.0);
        let e1 = Vec3::new(0.0, 1.0, 0.0);
        let e2 = Vec3::new(0.0, 0.0, 1.0);
        let (g0, g1, g2) = gram_schmidt(&e0, &e1, &e2).expect("axis-aligned should succeed");
        assert!((g0.dot(&e0).abs() - 1.0).abs() < 1e-10);
        assert!((g1.dot(&e1).abs() - 1.0).abs() < 1e-10);
        assert!((g2.dot(&e2).abs() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_log_rotation_90z() {
        let omega = Vec3::new(0.0, 0.0, std::f64::consts::FRAC_PI_2);
        let r = mat3_exp_skew(&omega);
        let log_r = mat3_log_rotation(&r);
        let recovered_angle = log_r[(1, 0)];
        assert!(
            (recovered_angle - std::f64::consts::FRAC_PI_2).abs() < 1e-8,
            "recovered angle={}",
            recovered_angle
        );
    }
    #[test]
    fn test_rodrigues_matches_mat3_from_axis_angle() {
        let axis = Vec3::new(0.0, 1.0, 0.0);
        let angle = 0.8_f64;
        let v = Vec3::new(1.0, 0.5, -0.3);
        let r1 = rodrigues_rotate(&v, &axis, angle);
        let r2 = mat3_from_axis_angle(&axis, angle) * v;
        let diff = (r1 - r2).norm();
        assert!(
            diff < 1e-10,
            "rodrigues_rotate should match mat3_from_axis_angle: diff={}",
            diff
        );
    }
    #[test]
    fn test_quat_nlerp_symmetric() {
        let a = quat_from_axis_angle(&Vec3::new(1.0, 0.0, 0.0), 0.3);
        let b = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.6);
        let q1 = quat_nlerp(&a, &b, 0.3);
        let q2 = quat_nlerp(&b, &a, 0.7);
        let d = quat_geodesic_distance(&q1, &q2);
        assert!(d < 1e-8, "nlerp should be symmetric: angle diff={}", d);
    }
}
/// Project vector `v` onto unit direction `onto_unit`.
///
/// `onto_unit` is assumed to be normalised; the caller is responsible for
/// ensuring that is the case.
pub fn vec3_project_onto_unit(v: &Vec3, onto_unit: &Vec3) -> Vec3 {
    onto_unit * v.dot(onto_unit)
}
/// Reject (orthogonal complement) of `v` with respect to unit direction `u`:
/// `reject(v, u) = v - project(v, u)`.
pub fn vec3_reject(v: &Vec3, onto_unit: &Vec3) -> Vec3 {
    v - vec3_project_onto_unit(v, onto_unit)
}
/// Reflect vector `v` about unit normal `n`: `v - 2*(v·n)*n`.
pub fn vec3_reflect_about_normal(v: &Vec3, n: &Vec3) -> Vec3 {
    v - n * (2.0 * v.dot(n))
}
/// Refract vector `v` (unit incident) through a surface with unit normal `n`.
///
/// `eta_ratio = eta_i / eta_t`.  Returns `None` for total internal reflection.
/// This variant uses the `eta_ratio` parameter name convention.
pub fn vec3_refract_ratio(v: &Vec3, n: &Vec3, eta_ratio: Real) -> Option<Vec3> {
    let cos_i = -v.dot(n);
    let sin2_t = eta_ratio * eta_ratio * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        return None;
    }
    let cos_t = (1.0 - sin2_t).sqrt();
    Some(v * eta_ratio + n * (eta_ratio * cos_i - cos_t))
}
/// Component-wise minimum of two `Vec3` values.
pub fn vec3_min(a: &Vec3, b: &Vec3) -> Vec3 {
    Vec3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z))
}
/// Component-wise maximum of two `Vec3` values.
pub fn vec3_max(a: &Vec3, b: &Vec3) -> Vec3 {
    Vec3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z))
}
/// Component-wise clamp of `v` to `[lo, hi]`.
pub fn vec3_clamp(v: &Vec3, lo: Real, hi: Real) -> Vec3 {
    Vec3::new(v.x.clamp(lo, hi), v.y.clamp(lo, hi), v.z.clamp(lo, hi))
}
/// Absolute value of each component of `v`.
pub fn vec3_abs(v: &Vec3) -> Vec3 {
    Vec3::new(v.x.abs(), v.y.abs(), v.z.abs())
}
/// Signed angle (in radians) from `a` to `b` measured about `axis`.
///
/// All three vectors should be unit vectors; `axis` must be perpendicular to
/// both `a` and `b` for the result to be meaningful.
pub fn vec3_signed_angle(a: &Vec3, b: &Vec3, axis: &Vec3) -> Real {
    let cross = a.cross(b);
    let sin_angle = cross.dot(axis);
    let cos_angle = a.dot(b);
    sin_angle.atan2(cos_angle)
}
/// Build an orthonormal basis `(tangent, bitangent, normal)` given a surface
/// normal `n` (must be unit length).
///
/// Uses the Frisvad / Duff–Burgess–Christensen (2017) numerically stable
/// construction.
pub fn build_orthonormal_basis(n: &Vec3) -> (Vec3, Vec3, Vec3) {
    let sign = if n.z >= 0.0 { 1.0 } else { -1.0 };
    let a = -1.0 / (sign + n.z);
    let b = n.x * n.y * a;
    let t = Vec3::new(1.0 + sign * n.x * n.x * a, sign * b, -sign * n.x);
    let bt = Vec3::new(b, sign + n.y * n.y * a, -n.y);
    (t, bt, *n)
}
/// Decompose a `Vec3` into components parallel and perpendicular to `axis`
/// (which must be a unit vector).
///
/// Returns `(parallel, perpendicular)` such that `parallel + perpendicular == v`.
pub fn vec3_decompose(v: &Vec3, axis: &Vec3) -> (Vec3, Vec3) {
    let parallel = vec3_project_onto_unit(v, axis);
    let perp = v - parallel;
    (parallel, perp)
}
/// Scalar triple product `a · (b × c)`.
pub fn vec3_scalar_triple(a: &Vec3, b: &Vec3, c: &Vec3) -> Real {
    a.dot(&b.cross(c))
}
/// Vector triple product `a × (b × c) = b*(a·c) - c*(a·b)`.
pub fn vec3_vector_triple(a: &Vec3, b: &Vec3, c: &Vec3) -> Vec3 {
    b * a.dot(c) - c * a.dot(b)
}
/// Distance from point `p` to the line defined by `origin` and unit `dir`.
pub fn point_to_line_distance(p: &Vec3, origin: &Vec3, dir: &Vec3) -> Real {
    let dp = p - origin;
    vec3_reject(&dp, dir).norm()
}
/// Closest point on line (`origin` + t * `dir`) to point `p`.
/// `dir` should be a unit vector.
pub fn closest_point_on_line(p: &Vec3, origin: &Vec3, dir: &Vec3) -> Vec3 {
    let t = (p - origin).dot(dir);
    origin + dir * t
}
/// Closest point on segment `a`–`b` to point `p`.
pub fn closest_point_on_segment(p: &Vec3, a: &Vec3, b: &Vec3) -> Vec3 {
    let ab = b - a;
    let ab_len2 = ab.norm_squared();
    if ab_len2 < 1e-30 {
        return *a;
    }
    let t = ((p - a).dot(&ab) / ab_len2).clamp(0.0, 1.0);
    a + ab * t
}
/// Distance from point `p` to segment `a`–`b`.
pub fn point_to_segment_distance(p: &Vec3, a: &Vec3, b: &Vec3) -> Real {
    (p - closest_point_on_segment(p, a, b)).norm()
}
/// Build a rotation matrix that takes unit vector `from` to unit vector `to`.
///
/// Uses Rodrigues' formula via the cross-product axis.
/// For anti-parallel vectors it returns a 180° rotation about an arbitrary
/// perpendicular axis.
pub fn mat3_rotation_from_to(from: &Vec3, to: &Vec3) -> Mat3 {
    let dot = from.dot(to).clamp(-1.0, 1.0);
    if (dot - 1.0).abs() < 1e-12 {
        return Mat3::identity();
    }
    if (dot + 1.0).abs() < 1e-6 {
        let perp = if from.x.abs() < 0.9 {
            Vec3::new(0.0, -from.z, from.y).normalize()
        } else {
            Vec3::new(-from.z, 0.0, from.x).normalize()
        };
        return 2.0 * perp * perp.transpose() - Mat3::identity();
    }
    let axis = from.cross(to);
    let axis_len = axis.norm();
    if axis_len < 1e-10 {
        let perp = if from.x.abs() < 0.9 {
            Vec3::new(0.0, -from.z, from.y).normalize()
        } else {
            Vec3::new(-from.z, 0.0, from.x).normalize()
        };
        return 2.0 * perp * perp.transpose() - Mat3::identity();
    }
    let angle = dot.acos();
    mat3_from_axis_angle(&(axis / axis_len), angle)
}
/// Extract the rotation angle from a 3×3 rotation matrix.
///
/// Returns the angle in `[0, π]`.
pub fn mat3_rotation_angle(r: &Mat3) -> Real {
    let tr = r.trace().clamp(-1.0, 3.0);
    ((tr - 1.0) / 2.0).clamp(-1.0, 1.0).acos()
}
/// Extract the rotation axis (unit vector) from a 3×3 rotation matrix.
///
/// Returns `None` for identity or 180° rotations where the axis is
/// ill-defined; otherwise returns the normalised axis.
pub fn mat3_rotation_axis(r: &Mat3) -> Option<Vec3> {
    let s = Vec3::new(
        r[(2, 1)] - r[(1, 2)],
        r[(0, 2)] - r[(2, 0)],
        r[(1, 0)] - r[(0, 1)],
    );
    let len = s.norm();
    if len < 1e-10 {
        return None;
    }
    Some(s / len)
}
/// Compose two rotation matrices: `r_total = r2 * r1` (first apply r1, then r2).
pub fn mat3_compose_rotations(r1: &Mat3, r2: &Mat3) -> Mat3 {
    r2 * r1
}
/// Interpolate between two rotation matrices using matrix slerp.
///
/// Computes `r1 * exp(t * log(r1^T * r2))`.
pub fn mat3_slerp(r1: &Mat3, r2: &Mat3, t: Real) -> Mat3 {
    let r_rel = r1.transpose() * r2;
    let log_r_rel = mat3_log_rotation(&r_rel);
    let scaled_log = log_r_rel * t;
    let omega = Vec3::new(scaled_log[(2, 1)], scaled_log[(0, 2)], scaled_log[(1, 0)]);
    r1 * mat3_exp_skew(&omega)
}
/// Compute the average of multiple rotation matrices using the iterative
/// geodesic mean (Moakher 2002).  `weights` must be non-negative and sum to 1.
pub fn mat3_weighted_mean(rotations: &[Mat3], weights: &[f64], max_iter: usize) -> Mat3 {
    assert_eq!(rotations.len(), weights.len());
    if rotations.is_empty() {
        return Mat3::identity();
    }
    let mut r_mean = rotations[0];
    for _ in 0..max_iter {
        let mut delta_sum = Mat3::zeros();
        for (r, &w) in rotations.iter().zip(weights.iter()) {
            let log = mat3_log_rotation(&(r_mean.transpose() * r));
            delta_sum += log * w;
        }
        if delta_sum.norm() < 1e-12 {
            break;
        }
        let omega = Vec3::new(delta_sum[(2, 1)], delta_sum[(0, 2)], delta_sum[(1, 0)]);
        r_mean *= mat3_exp_skew(&omega);
    }
    r_mean
}
/// Convert a unit quaternion `[x,y,z,w]` to its conjugate.
pub fn quat_arr_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}
/// Normalize a quaternion `[x,y,z,w]` to unit length.
pub fn quat_arr_normalize(q: [f64; 4]) -> [f64; 4] {
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if len < 1e-30 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / len, q[1] / len, q[2] / len, q[3] / len]
}
/// Dot product of two quaternions (as 4-vectors).
pub fn quat_arr_dot(p: [f64; 4], q: [f64; 4]) -> f64 {
    p[0] * q[0] + p[1] * q[1] + p[2] * q[2] + p[3] * q[3]
}
/// Power of a unit quaternion: `q^t = exp(t * log(q))`.
pub fn quat_arr_pow(q: [f64; 4], t: f64) -> [f64; 4] {
    let log_q = quat_arr_log(q);
    let scaled = [log_q[0] * t, log_q[1] * t, log_q[2] * t, 0.0];
    quat_arr_exp(scaled)
}
/// Rotate a 3-vector by a unit quaternion `q = [x,y,z,w]` using
/// the sandwich product `q * v * q^*`.
pub fn quat_arr_rotate_vec3(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    let [vx, vy, vz] = v;
    let tx = 2.0 * (qy * vz - qz * vy);
    let ty = 2.0 * (qz * vx - qx * vz);
    let tz = 2.0 * (qx * vy - qy * vx);
    [
        vx + qw * tx + qy * tz - qz * ty,
        vy + qw * ty + qz * tx - qx * tz,
        vz + qw * tz + qx * ty - qy * tx,
    ]
}
/// Double-geodesic (SQUAD) interpolation at midpoint between two key frames.
///
/// Given quaternions `q_prev`, `q_curr`, `q_next`, compute the SQUAD inner
/// control point `s` such that the curve is C1 at `q_curr`.
///
/// `s = q_curr * exp(-0.25 * (log(conj(q_curr)*q_next) + log(conj(q_curr)*q_prev)))`
pub fn quat_arr_squad_control(q_prev: [f64; 4], q_curr: [f64; 4], q_next: [f64; 4]) -> [f64; 4] {
    let conj_curr = quat_arr_conjugate(q_curr);
    let log_a = quat_arr_log(quat_multiply(conj_curr, q_next));
    let log_b = quat_arr_log(quat_multiply(conj_curr, q_prev));
    let sum = [
        (log_a[0] + log_b[0]) * (-0.25),
        (log_a[1] + log_b[1]) * (-0.25),
        (log_a[2] + log_b[2]) * (-0.25),
        0.0,
    ];
    quat_multiply(q_curr, quat_arr_exp(sum))
}
/// Gram-Schmidt orthonormalization of up to `n` vectors stored in `vecs`.
///
/// Returns a new vector of orthonormal vectors.  Linearly-dependent vectors
/// produce a near-zero output and are replaced by a vector from the standard
/// basis so the output always spans the same subspace dimension.
pub fn gram_schmidt_vectors(vecs: &[Vec3]) -> Vec<Vec3> {
    let mut result: Vec<Vec3> = Vec::with_capacity(vecs.len());
    for v in vecs {
        let mut u = *v;
        for q in &result {
            u -= q * u.dot(q);
        }
        let len = u.norm();
        if len < 1e-10 {
            let fallbacks = [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ];
            for fb in &fallbacks {
                let mut candidate = *fb;
                for q in &result {
                    candidate -= q * candidate.dot(q);
                }
                if candidate.norm() > 1e-10 {
                    u = candidate.normalize();
                    break;
                }
            }
        } else {
            u /= len;
        }
        result.push(u);
    }
    result
}
/// Modified Gram-Schmidt applied to the columns of a `Mat3`.
///
/// Returns an orthonormal `Mat3` whose columns span the same space as the input
/// (assuming the input is full-rank).  Returns `None` if the matrix is rank-deficient.
pub fn mat3_gram_schmidt(m: &Mat3) -> Option<Mat3> {
    let c0 = Vec3::new(m[(0, 0)], m[(1, 0)], m[(2, 0)]);
    let c1 = Vec3::new(m[(0, 1)], m[(1, 1)], m[(2, 1)]);
    let c2 = Vec3::new(m[(0, 2)], m[(1, 2)], m[(2, 2)]);
    let (q0, q1, q2) = gram_schmidt(&c0, &c1, &c2)?;
    Some(Mat3::from_columns(&[q0, q1, q2]))
}
/// Project a `Vec3` onto a plane defined by its unit normal `n`.
///
/// Returns the vector with its component along `n` removed.
pub fn vec3_project_onto_plane(v: &Vec3, plane_normal: &Vec3) -> Vec3 {
    v - plane_normal * v.dot(plane_normal)
}
/// Compute the reflection of direction `d` about unit normal `n`.
///
/// Equivalent to `vec3_reflect_about_normal` but uses the physics
/// convention where `d` points *toward* the surface.
pub fn reflect_direction(d: &Vec3, n: &Vec3) -> Vec3 {
    d - n * (2.0 * d.dot(n))
}
/// Barycentric coordinates of point `p` with respect to triangle `(a, b, c)`.
///
/// Returns `(u, v, w)` such that `p = u*a + v*b + w*c` and `u+v+w = 1`.
/// Returns `None` if the triangle is degenerate.
pub fn barycentric_coords(p: &Vec3, a: &Vec3, b: &Vec3, c: &Vec3) -> Option<(Real, Real, Real)> {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let normal = ab.cross(&ac);
    let denom = normal.norm_squared();
    if denom < 1e-20 {
        return None;
    }
    let v = normal.dot(&ap.cross(&ac)) / denom;
    let w = normal.dot(&ab.cross(&ap)) / denom;
    let u = 1.0 - v - w;
    Some((u, v, w))
}
/// Test whether point `p` lies inside triangle `(a, b, c)`.
pub fn point_in_triangle(p: &Vec3, a: &Vec3, b: &Vec3, c: &Vec3) -> bool {
    match barycentric_coords(p, a, b, c) {
        None => false,
        Some((u, v, w)) => u >= -1e-10 && v >= -1e-10 && w >= -1e-10,
    }
}
/// Cross product matrix `[a]×` such that `[a]× * b = a × b`.
///
/// Equivalent to `skew_symmetric` but named for clarity in the context of the
/// cross-product-as-linear-map interpretation.
pub fn cross_product_matrix(a: &Vec3) -> Mat3 {
    skew_symmetric(a)
}
/// Double cross product `a × (a × b)` expressed as a matrix times `b`:
/// the matrix is `-[a]×² = a*aᵀ - (a·a)*I`.
pub fn double_cross_matrix(a: &Vec3) -> Mat3 {
    let aa = a.norm_squared();
    a * a.transpose() - Mat3::identity() * aa
}
/// Spin matrix for angular velocity `omega`: the time-derivative of a rotation
/// matrix R satisfies `Ṙ = omega_hat * R` where `omega_hat = skew(omega)`.
pub fn spin_matrix(omega: &Vec3) -> Mat3 {
    skew_symmetric(omega)
}
/// Runge–Kutta 4th-order integration of a quaternion ODE.
///
/// Given a unit quaternion `q` and angular velocity `omega` (body frame),
/// integrates `q̇ = 0.5 * q ⊗ [omega, 0]` forward by `dt`.
pub fn quat_integrate_rk4(q: &Quat, omega: &Vec3, dt: Real) -> Quat {
    let half_dt = dt * 0.5;
    let omega_quat = |q_in: &Quat| -> Quat {
        let ow = UnitQuaternion::from_quaternion(nalgebra::Quaternion::new(
            0.0, omega.x, omega.y, omega.z,
        ));
        UnitQuaternion::new_normalize(q_in.as_ref() * (ow.as_ref() * 0.5))
    };
    let k1 = omega_quat(q);
    let q2 = UnitQuaternion::new_normalize(q.as_ref() + k1.as_ref() * half_dt);
    let k2 = omega_quat(&q2);
    let q3 = UnitQuaternion::new_normalize(q.as_ref() + k2.as_ref() * half_dt);
    let k3 = omega_quat(&q3);
    let q4 = UnitQuaternion::new_normalize(q.as_ref() + k3.as_ref() * dt);
    let k4 = omega_quat(&q4);
    let sum = k1.as_ref() + k2.as_ref() * 2.0 + k3.as_ref() * 2.0 + k4.as_ref();
    UnitQuaternion::new_normalize(q.as_ref() + sum * (dt / 6.0))
}
/// Linear blending of two quaternions followed by renormalization (NLerp).
pub fn quat_blend(q1: &Quat, q2: &Quat, t: Real) -> Quat {
    q1.slerp(q2, t)
}
/// Convert a Cartesian 3-vector to cylindrical coordinates `(rho, phi, z)`.
///
/// `rho` ≥ 0, `phi ∈ (-π, π]`, `z` is unchanged.
pub fn vec3_to_cylindrical(v: &Vec3) -> (Real, Real, Real) {
    let rho = (v.x * v.x + v.y * v.y).sqrt();
    let phi = v.y.atan2(v.x);
    (rho, phi, v.z)
}
/// Convert cylindrical `(rho, phi, z)` to a Cartesian `Vec3`.
pub fn cylindrical_to_vec3(rho: Real, phi: Real, z: Real) -> Vec3 {
    Vec3::new(rho * phi.cos(), rho * phi.sin(), z)
}
/// Convert a Cartesian `Vec3` to spherical coordinates `(r, theta, phi)`.
///
/// `r` ≥ 0, `theta ∈ [0, π]` (polar/colatitude), `phi ∈ (-π, π]` (azimuth).
pub fn vec3_to_spherical(v: &Vec3) -> (Real, Real, Real) {
    let r = v.norm();
    if r < 1e-30 {
        return (0.0, 0.0, 0.0);
    }
    let theta = (v.z / r).clamp(-1.0, 1.0).acos();
    let phi = v.y.atan2(v.x);
    (r, theta, phi)
}
/// Convert spherical `(r, theta, phi)` to a Cartesian `Vec3`.
pub fn spherical_to_vec3(r: Real, theta: Real, phi: Real) -> Vec3 {
    Vec3::new(
        r * theta.sin() * phi.cos(),
        r * theta.sin() * phi.sin(),
        r * theta.cos(),
    )
}
/// Uniformly sample a direction on the unit hemisphere (z ≥ 0) given two
/// uniform random numbers `u1, u2 ∈ [0, 1)`.
///
/// Returns a unit vector `(x, y, z)` with `z ≥ 0`.
pub fn sample_hemisphere_uniform(u1: Real, u2: Real) -> Vec3 {
    let cos_theta = u1;
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let phi = 2.0 * std::f64::consts::PI * u2;
    Vec3::new(sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta)
}
/// Cosine-weighted sample on the unit hemisphere (z ≥ 0).
///
/// Returns a unit vector; the PDF is `cos(theta) / π`.
pub fn sample_hemisphere_cosine(u1: Real, u2: Real) -> Vec3 {
    let r = u1.sqrt();
    let phi = 2.0 * std::f64::consts::PI * u2;
    let x = r * phi.cos();
    let y = r * phi.sin();
    let z = (1.0 - u1).max(0.0).sqrt();
    Vec3::new(x, y, z)
}
/// Uniformly sample a point on the surface of a unit disk.
///
/// Uses Shirley's concentric-disk mapping for low distortion.
/// Returns `(x, y)` with `x² + y² ≤ 1`.
pub fn sample_unit_disk_concentric(u1: Real, u2: Real) -> (Real, Real) {
    let a = 2.0 * u1 - 1.0;
    let b = 2.0 * u2 - 1.0;
    if a == 0.0 && b == 0.0 {
        return (0.0, 0.0);
    }
    let (r, phi) = if a.abs() > b.abs() {
        (a, std::f64::consts::FRAC_PI_4 * (b / a))
    } else {
        (
            b,
            std::f64::consts::FRAC_PI_2 - std::f64::consts::FRAC_PI_4 * (a / b),
        )
    };
    (r * phi.cos(), r * phi.sin())
}
#[cfg(test)]
mod extended_math_tests_v2 {

    use crate::Mat3;
    use crate::Vec3;
    use crate::math::barycentric_coords;
    use crate::math::build_orthonormal_basis;

    use crate::math::cylindrical_to_vec3;
    use crate::math::double_cross_matrix;

    use crate::math::gram_schmidt_vectors;

    use crate::math::mat3_from_axis_angle;
    use crate::math::mat3_gram_schmidt;

    use crate::math::mat3_rotation_angle;
    use crate::math::mat3_rotation_axis;
    use crate::math::mat3_rotation_from_to;
    use crate::math::mat3_slerp;
    use crate::math::mat3_weighted_mean;
    use crate::math::point_in_triangle;
    use crate::math::point_to_segment_distance;
    use crate::math::quat_arr_conjugate;
    use crate::math::quat_arr_dot;
    use crate::math::quat_arr_from_axis_angle;
    use crate::math::quat_arr_normalize;
    use crate::math::quat_arr_pow;
    use crate::math::quat_arr_rotate_vec3;

    use crate::math::quat_multiply;

    use crate::math::sample_hemisphere_cosine;
    use crate::math::sample_hemisphere_uniform;
    use crate::math::sample_unit_disk_concentric;

    use crate::math::spherical_to_vec3;
    use crate::math::vec3_abs;
    use crate::math::vec3_clamp;
    use crate::math::vec3_decompose;
    use crate::math::vec3_max;
    use crate::math::vec3_min;

    use crate::math::vec3_project_onto_unit;
    use crate::math::vec3_reflect_about_normal;
    use crate::math::vec3_refract;
    use crate::math::vec3_reject;
    use crate::math::vec3_scalar_triple;
    use crate::math::vec3_to_cylindrical;
    use crate::math::vec3_to_spherical;
    #[test]
    fn test_project_onto_unit_along_axis() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let u = Vec3::new(1.0, 0.0, 0.0);
        let proj = vec3_project_onto_unit(&v, &u);
        assert!((proj.x - 3.0).abs() < 1e-12);
        assert!(proj.y.abs() < 1e-12);
        assert!(proj.z.abs() < 1e-12);
    }
    #[test]
    fn test_reject_perpendicular_to_axis() {
        let v = Vec3::new(3.0, 4.0, 5.0);
        let u = Vec3::new(0.0, 0.0, 1.0);
        let rej = vec3_reject(&v, &u);
        assert!(rej.dot(&u).abs() < 1e-12, "reject should be perp to axis");
    }
    #[test]
    fn test_project_plus_reject_equals_v() {
        let v = Vec3::new(1.0, -2.0, 3.0);
        let u = Vec3::new(0.0, 1.0, 0.0);
        let proj = vec3_project_onto_unit(&v, &u);
        let rej = vec3_reject(&v, &u);
        let diff = (proj + rej - v).norm();
        assert!(diff < 1e-12, "proj + rej should equal v, diff={}", diff);
    }
    #[test]
    fn test_reflect_normal_preserves_magnitude() {
        let v = Vec3::new(1.0, -1.0, 0.0).normalize();
        let n = Vec3::new(0.0, 1.0, 0.0);
        let r = vec3_reflect_about_normal(&v, &n);
        assert!((r.norm() - v.norm()).abs() < 1e-12);
    }
    #[test]
    fn test_reflect_vertical_normal() {
        let v = Vec3::new(1.0, -1.0, 0.0);
        let n = Vec3::new(0.0, 1.0, 0.0);
        let r = vec3_reflect_about_normal(&v, &n);
        assert!((r.x - 1.0).abs() < 1e-12);
        assert!((r.y - 1.0).abs() < 1e-12);
        assert!(r.z.abs() < 1e-12);
    }
    #[test]
    fn test_refract_normal_incidence() {
        let v = Vec3::new(0.0, -1.0, 0.0);
        let n = Vec3::new(0.0, 1.0, 0.0);
        let t = vec3_refract(&v, &n, 1.0).unwrap();
        let diff = (t - v).norm();
        assert!(diff < 1e-10, "normal incidence same-medium: diff={}", diff);
    }
    #[test]
    fn test_refract_total_internal_reflection() {
        let v = Vec3::new(0.99, -0.14, 0.0).normalize();
        let n = Vec3::new(0.0, 1.0, 0.0);
        let result = vec3_refract(&v, &n, 1.5);
        let _ = result;
    }
    #[test]
    fn test_vec3_min_max() {
        let a = Vec3::new(1.0, 5.0, 3.0);
        let b = Vec3::new(2.0, 3.0, 4.0);
        let mn = vec3_min(&a, &b);
        let mx = vec3_max(&a, &b);
        assert!((mn.x - 1.0).abs() < 1e-12);
        assert!((mn.y - 3.0).abs() < 1e-12);
        assert!((mn.z - 3.0).abs() < 1e-12);
        assert!((mx.x - 2.0).abs() < 1e-12);
        assert!((mx.y - 5.0).abs() < 1e-12);
        assert!((mx.z - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_abs_all_negative() {
        let v = Vec3::new(-1.0, -2.0, -3.0);
        let a = vec3_abs(&v);
        assert!((a.x - 1.0).abs() < 1e-12);
        assert!((a.y - 2.0).abs() < 1e-12);
        assert!((a.z - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_clamp() {
        let v = Vec3::new(-5.0, 0.5, 10.0);
        let c = vec3_clamp(&v, 0.0, 1.0);
        assert!((c.x - 0.0).abs() < 1e-12);
        assert!((c.y - 0.5).abs() < 1e-12);
        assert!((c.z - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_build_orthonormal_basis_orthogonal() {
        let n = Vec3::new(0.0, 1.0, 0.0);
        let (t, bt, nn) = build_orthonormal_basis(&n);
        assert!(t.dot(&bt).abs() < 1e-10, "t·bt should be 0");
        assert!(t.dot(&nn).abs() < 1e-10, "t·n should be 0");
        assert!(bt.dot(&nn).abs() < 1e-10, "bt·n should be 0");
    }
    #[test]
    fn test_build_orthonormal_basis_unit_vectors() {
        let n = Vec3::new(1.0, 2.0, 3.0).normalize();
        let (t, bt, _) = build_orthonormal_basis(&n);
        assert!((t.norm() - 1.0).abs() < 1e-10, "t should be unit");
        assert!((bt.norm() - 1.0).abs() < 1e-10, "bt should be unit");
    }
    #[test]
    fn test_vec3_decompose_sums_to_original() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let (par, perp) = vec3_decompose(&v, &axis);
        let diff = (par + perp - v).norm();
        assert!(diff < 1e-12, "parallel + perp should be v, diff={}", diff);
    }
    #[test]
    fn test_vec3_decompose_parallel_along_axis() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let (par, _) = vec3_decompose(&v, &axis);
        assert!(par.x.abs() < 1e-12);
        assert!(par.y.abs() < 1e-12);
        assert!((par.z - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_triple_basis_vectors() {
        let e0 = Vec3::new(1.0, 0.0, 0.0);
        let e1 = Vec3::new(0.0, 1.0, 0.0);
        let e2 = Vec3::new(0.0, 0.0, 1.0);
        let triple = vec3_scalar_triple(&e0, &e1, &e2);
        assert!((triple - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_triple_antisymmetry() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let c = Vec3::new(7.0, 8.0, 9.0);
        let abc = vec3_scalar_triple(&a, &b, &c);
        let bca = vec3_scalar_triple(&b, &c, &a);
        assert!(
            (abc - bca).abs() < 1e-10,
            "cyclic permutation should be equal"
        );
    }
    #[test]
    fn test_point_to_segment_endpoint() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let p = Vec3::new(2.0, 0.0, 0.0);
        let d = point_to_segment_distance(&p, &a, &b);
        assert!((d - 1.0).abs() < 1e-12, "dist should be 1.0, got {}", d);
    }
    #[test]
    fn test_point_to_segment_perpendicular() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(2.0, 0.0, 0.0);
        let p = Vec3::new(1.0, 3.0, 0.0);
        let d = point_to_segment_distance(&p, &a, &b);
        assert!((d - 3.0).abs() < 1e-12, "dist should be 3.0, got {}", d);
    }
    #[test]
    fn test_mat3_rotation_from_to_identity() {
        let v = Vec3::new(1.0, 0.0, 0.0);
        let r = mat3_rotation_from_to(&v, &v);
        let diff = (r - Mat3::identity()).norm();
        assert!(diff < 1e-10, "from_to same direction should be identity");
    }
    #[test]
    fn test_mat3_rotation_from_to_correct() {
        let from = Vec3::new(1.0, 0.0, 0.0);
        let to = Vec3::new(0.0, 1.0, 0.0);
        let r = mat3_rotation_from_to(&from, &to);
        let result = r * from;
        let diff = (result - to).norm();
        assert!(diff < 1e-10, "rotation_from_to mismatch: diff={}", diff);
    }
    #[test]
    fn test_mat3_rotation_from_to_antiparallel() {
        let from = Vec3::new(1.0, 0.0, 0.0);
        let to = Vec3::new(-1.0, 0.0, 0.0);
        let r = mat3_rotation_from_to(&from, &to);
        let result = r * from;
        let len = result.norm();
        assert!(
            (len - 1.0).abs() < 1e-8,
            "rotation preserves length: len={}",
            len
        );
    }
    #[test]
    fn test_mat3_slerp_at_endpoints() {
        let r1 = mat3_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), 0.3);
        let r2 = mat3_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), 0.9);
        let s0 = mat3_slerp(&r1, &r2, 0.0);
        let s1 = mat3_slerp(&r1, &r2, 1.0);
        assert!((s0 - r1).norm() < 1e-8, "slerp(t=0) should be r1");
        assert!((s1 - r2).norm() < 1e-8, "slerp(t=1) should be r2");
    }
    #[test]
    fn test_mat3_slerp_midpoint_is_rotation() {
        let r1 = mat3_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let r2 = mat3_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let mid = mat3_slerp(&r1, &r2, 0.5);
        let det = mid.determinant();
        let rrt = mid * mid.transpose();
        let diff = (rrt - Mat3::identity()).norm();
        assert!((det - 1.0).abs() < 1e-7, "slerp midpoint det={}", det);
        assert!(diff < 1e-7, "slerp midpoint R R^T diff={}", diff);
    }
    #[test]
    fn test_quat_arr_conjugate_product_is_identity() {
        let q = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 0.5);
        let conj = quat_arr_conjugate(q);
        let prod = quat_multiply(q, conj);
        assert!(prod[0].abs() < 1e-12);
        assert!(prod[1].abs() < 1e-12);
        assert!(prod[2].abs() < 1e-12);
        assert!((prod[3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_quat_arr_normalize_unit_input() {
        let q = [0.0, 0.0, 0.0, 1.0];
        let n = quat_arr_normalize(q);
        assert!((n[3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_quat_arr_rotate_vec3_z90() {
        let q = quat_arr_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let v = [1.0_f64, 0.0, 0.0];
        let r = quat_arr_rotate_vec3(q, v);
        assert!(r[0].abs() < 1e-10, "x={}", r[0]);
        assert!((r[1] - 1.0).abs() < 1e-10, "y={}", r[1]);
        assert!(r[2].abs() < 1e-10, "z={}", r[2]);
    }
    #[test]
    fn test_quat_arr_pow_half_is_sqrt() {
        let q = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 0.8);
        let q_half = quat_arr_pow(q, 0.5);
        let q_reconstructed = quat_multiply(q_half, q_half);
        let q_norm = quat_arr_normalize(q_reconstructed);
        let dot = quat_arr_dot(q, q_norm).abs();
        assert!(dot > 0.9999, "q^0.5 * q^0.5 should equal q: dot={}", dot);
    }
    #[test]
    fn test_gram_schmidt_vectors_orthonormal() {
        let vecs = vec![
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 1.0),
        ];
        let qs = gram_schmidt_vectors(&vecs);
        assert_eq!(qs.len(), 3);
        for i in 0..3 {
            assert!((qs[i].norm() - 1.0).abs() < 1e-10, "q[{}] not unit", i);
            for j in (i + 1)..3 {
                let dot = qs[i].dot(&qs[j]);
                assert!(dot.abs() < 1e-10, "q[{}]·q[{}]={} != 0", i, j, dot);
            }
        }
    }
    #[test]
    fn test_mat3_gram_schmidt_gives_orthogonal_result() {
        let m = Mat3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0);
        let q = mat3_gram_schmidt(&m).expect("should succeed for non-singular");
        let diff = (q * q.transpose() - Mat3::identity()).norm();
        assert!(
            diff < 1e-8,
            "gram-schmidt result not orthogonal: diff={}",
            diff
        );
    }
    #[test]
    fn test_barycentric_vertex_a() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 1.0, 0.0);
        let (u, v, w) = barycentric_coords(&a, &a, &b, &c).unwrap();
        assert!((u - 1.0).abs() < 1e-10);
        assert!(v.abs() < 1e-10);
        assert!(w.abs() < 1e-10);
    }
    #[test]
    fn test_barycentric_centroid() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(3.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 3.0, 0.0);
        let centroid = (a + b + c) / 3.0;
        let (u, v, w) = barycentric_coords(&centroid, &a, &b, &c).unwrap();
        let third = 1.0 / 3.0;
        assert!((u - third).abs() < 1e-10, "u={}", u);
        assert!((v - third).abs() < 1e-10, "v={}", v);
        assert!((w - third).abs() < 1e-10, "w={}", w);
    }
    #[test]
    fn test_point_in_triangle_inside() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(0.2, 0.2, 0.0);
        assert!(point_in_triangle(&p, &a, &b, &c));
    }
    #[test]
    fn test_point_in_triangle_outside() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(0.8, 0.8, 0.0);
        assert!(!point_in_triangle(&p, &a, &b, &c));
    }
    #[test]
    fn test_sample_hemisphere_uniform_unit_length() {
        let v = sample_hemisphere_uniform(0.5, 0.5);
        assert!((v.norm() - 1.0).abs() < 1e-12, "norm={}", v.norm());
    }
    #[test]
    fn test_sample_hemisphere_uniform_nonneg_z() {
        for i in 0..10 {
            let u1 = i as f64 / 10.0;
            let u2 = (i + 1) as f64 / 11.0;
            let v = sample_hemisphere_uniform(u1, u2);
            assert!(v.z >= -1e-12, "z={}", v.z);
        }
    }
    #[test]
    fn test_sample_hemisphere_cosine_unit_length() {
        let v = sample_hemisphere_cosine(0.3, 0.7);
        assert!((v.norm() - 1.0).abs() < 1e-12, "norm={}", v.norm());
    }
    #[test]
    fn test_sample_unit_disk_in_unit_circle() {
        for i in 0..10 {
            let u1 = i as f64 / 10.0;
            let u2 = (i + 3) as f64 / 13.0;
            let (x, y) = sample_unit_disk_concentric(u1, u2);
            assert!(
                x * x + y * y <= 1.0 + 1e-12,
                "outside unit circle: {},{}",
                x,
                y
            );
        }
    }
    #[test]
    fn test_mat3_rotation_angle_known() {
        let angle = std::f64::consts::FRAC_PI_3;
        let r = mat3_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), angle);
        let extracted = mat3_rotation_angle(&r);
        assert!((extracted - angle).abs() < 1e-8, "angle={}", extracted);
    }
    #[test]
    fn test_mat3_rotation_axis_known() {
        let r = mat3_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.5);
        let axis = mat3_rotation_axis(&r).expect("should have axis");
        assert!(axis.y.abs() > 0.99, "axis y={}", axis.y);
    }
    #[test]
    fn test_cylindrical_roundtrip_new() {
        let v = Vec3::new(2.0, 3.0, 4.0);
        let (rho, phi, z) = vec3_to_cylindrical(&v);
        let v2 = cylindrical_to_vec3(rho, phi, z);
        let diff = (v - v2).norm();
        assert!(diff < 1e-12, "cylindrical roundtrip diff={}", diff);
    }
    #[test]
    fn test_spherical_roundtrip_new() {
        let v = Vec3::new(1.0, 2.0, 2.0);
        let (r, theta, phi) = vec3_to_spherical(&v);
        let v2 = spherical_to_vec3(r, theta, phi);
        let diff = (v - v2).norm();
        assert!(diff < 1e-12, "spherical roundtrip diff={}", diff);
    }
    #[test]
    fn test_double_cross_matrix_bac_cab() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let dcm = double_cross_matrix(&a);
        let via_matrix = dcm * b;
        let direct = a * a.dot(&b) - b * a.norm_squared();
        let diff = (via_matrix - direct).norm();
        assert!(diff < 1e-10, "double cross matrix diff={}", diff);
    }
    #[test]
    fn test_mat3_weighted_mean_single_rotation() {
        let r = mat3_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.7);
        let mean = mat3_weighted_mean(&[r], &[1.0], 20);
        let diff = (mean - r).norm();
        assert!(
            diff < 1e-6,
            "mean of single rotation should be itself: diff={}",
            diff
        );
    }
    #[test]
    fn test_mat3_weighted_mean_two_equal_rotations() {
        let r = mat3_from_axis_angle(&Vec3::new(1.0, 0.0, 0.0), 0.4);
        let mean = mat3_weighted_mean(&[r, r], &[0.5, 0.5], 20);
        let diff = (mean - r).norm();
        assert!(diff < 1e-6, "mean of same rotation: diff={}", diff);
    }
}
