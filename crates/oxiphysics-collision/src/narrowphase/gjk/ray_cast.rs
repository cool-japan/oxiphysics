// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ray casting against a convex shape with GJK.
//!
//! This is the algorithm of G. van den Bergen, "Ray Casting against General
//! Convex Objects with Application to Continuous Collision Detection" (2004).
//! Only the support mapping of the shape is needed. The ray point `x` starts at
//! the ray origin; at every step the point `q` of the current support-point
//! simplex that is closest to `x` is found (the GJK sub-algorithm), the shape's
//! support point `p` in direction `x - q` is queried, and if the plane through
//! `p` with normal `x - q` separates `x` from the shape, `x` is advanced along
//! the ray until it reaches that plane. Because the plane supports the shape,
//! the ray cannot enter the shape before it gets there, so the ray parameter
//! only ever moves to a certified lower bound of the hit parameter; the cast
//! ends once `x` lies on the shape.

use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::Shape;

use super::types::RayCastResult;
use crate::narrowphase::convex_contact::closest_to_origin;

/// Maximum number of GJK/advance steps before the analytic fallback is used.
const MAX_RAY_ITERATIONS: usize = 128;
/// Squared relative distance (to the size of the simplex) below which the ray
/// point counts as lying on the shape.
const CONVERGENCE_SQ: f64 = 1e-22;

/// Outcome of the support-mapping ray cast.
enum Cast {
    /// The ray hits the shape.
    Hit(RayCastResult),
    /// The ray provably misses the shape within `max_t`.
    Miss,
    /// The iteration budget ran out before the cast converged.
    Unconverged,
}

/// Cast a ray against a convex shape using the GJK algorithm.
///
/// `ray_origin` and `ray_dir` define the ray in world space (`ray_dir` need not
/// be normalised) and `max_t` limits the ray length, measured along the
/// normalised direction.
///
/// Returns `Some(RayCastResult)` for the first point of the shape hit within
/// `[0, max_t]`: `t` is the distance along the normalised direction, `point` the
/// hit point and `normal` the outward unit surface normal there. A ray that
/// starts inside (or on) the shape hits at `t = 0` with the normal `-ray_dir`.
/// Rays that pass the shape at a distance, point away from it or would need more
/// than `max_t` to reach it return `None`.
///
/// Only the support mapping of the shape is used; should the iteration fail to
/// converge within its budget the shape's own [`Shape::ray_cast`] decides.
pub fn gjk_ray_cast(
    shape: &dyn Shape,
    transform: &Transform,
    ray_origin: Vec3,
    ray_dir: Vec3,
    max_t: f64,
) -> Option<RayCastResult> {
    let dir_len = ray_dir.norm();
    if !dir_len.is_finite() || dir_len < 1e-14 || max_t.is_nan() || max_t < 0.0 {
        return None;
    }
    let unit_dir = ray_dir / dir_len;
    match cast_by_support(shape, transform, ray_origin, unit_dir, max_t) {
        Cast::Hit(hit) => Some(hit),
        Cast::Miss => None,
        Cast::Unconverged => cast_analytic(shape, transform, ray_origin, unit_dir, max_t),
    }
}

/// van den Bergen's GJK ray cast in world space.
fn cast_by_support(
    shape: &dyn Shape,
    transform: &Transform,
    origin: Vec3,
    dir: Vec3,
    max_t: f64,
) -> Cast {
    let inverse_rotation = transform.rotation.inverse();
    let support =
        |d: Vec3| transform.transform_point(&shape.support_point(&(inverse_rotation * d)));
    let towards_shape = transform.position - origin;
    let seed = if towards_shape.norm() > 1e-12 {
        towards_shape.normalize()
    } else {
        dir
    };
    // Support points of the shape (not of `x - shape`): the simplex in the
    // Minkowski space is `x - p` for the current ray point `x`.
    let mut points: Vec<Vec3> = vec![support(seed)];
    let mut t = 0.0_f64;
    let mut x = origin;
    let mut normal: Option<Vec3> = None;
    for _ in 0..MAX_RAY_ITERATIONS {
        let offsets: Vec<[f64; 3]> = points
            .iter()
            .map(|p| {
                let w = x - p;
                [w.x, w.y, w.z]
            })
            .collect();
        let closest = closest_to_origin(&offsets);
        points = closest.keep[..closest.len]
            .iter()
            .map(|&i| points[i])
            .collect();
        let v = Vec3::new(closest.point[0], closest.point[1], closest.point[2]);
        let v_sq = v.norm_squared();
        let size_sq = offsets
            .iter()
            .map(|w| w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
            .fold(0.0_f64, f64::max);
        if v_sq <= CONVERGENCE_SQ * size_sq {
            return Cast::Hit(RayCastResult {
                t,
                point: x,
                normal: normal.unwrap_or(-dir),
            });
        }
        let v_len = v_sq.sqrt();
        let p = support(v / v_len);
        let v_dot_w = v.dot(&(x - p));
        if v_dot_w > 0.0 {
            // The plane through `p` with normal `v` separates `x` from the shape.
            let v_dot_r = v.dot(&dir);
            if v_dot_r >= 0.0 {
                return Cast::Miss;
            }
            t -= v_dot_w / v_dot_r;
            if t > max_t {
                return Cast::Miss;
            }
            x = origin + dir * t;
            normal = Some(v / v_len);
        }
        points.push(p);
    }
    Cast::Unconverged
}

/// Ray cast through the shape's analytic [`Shape::ray_cast`] in its local frame.
fn cast_analytic(
    shape: &dyn Shape,
    transform: &Transform,
    origin: Vec3,
    dir: Vec3,
    max_t: f64,
) -> Option<RayCastResult> {
    let local_origin = transform.inverse().transform_point(&origin);
    let local_dir = transform.rotation.inverse() * dir;
    let hit = shape.ray_cast(&local_origin, &local_dir, max_t)?;
    let world_normal = transform.transform_vector(&hit.normal);
    let length = world_normal.norm();
    let normal = if length > 1e-10 {
        world_normal / length
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    Some(RayCastResult {
        t: hit.toi,
        point: transform.transform_point(&hit.point),
        normal,
    })
}

#[cfg(test)]
mod tests {
    use super::super::Gjk;
    use super::*;
    use oxiphysics_geometry::{BoxShape, Capsule, Cone, ConvexHull, Cylinder, Sphere};

    /// Small deterministic generator for the randomised cross-checks.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }

        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.next_f64()
        }

        fn vector(&mut self, extent: f64) -> Vec3 {
            Vec3::new(
                self.range(-extent, extent),
                self.range(-extent, extent),
                self.range(-extent, extent),
            )
        }
    }

    fn at(x: f64, y: f64, z: f64) -> Transform {
        Transform::from_position(Vec3::new(x, y, z))
    }

    fn cube_hull() -> ConvexHull {
        let mut vertices = Vec::new();
        for sx in [-0.5, 0.5] {
            for sy in [-0.5, 0.5] {
                for sz in [-0.5, 0.5] {
                    vertices.push(Vec3::new(sx, sy, sz));
                }
            }
        }
        ConvexHull::new(vertices)
    }

    #[test]
    fn ray_passing_far_beside_the_shape_is_a_miss() {
        // The shape's projection on the ray is at z = 5 but it lies 100 units off.
        let sphere = Sphere::new(1.0);
        let origin = Vec3::zeros();
        let hit = gjk_ray_cast(&sphere, &at(0.0, 100.0, 5.0), origin, Vec3::z(), 20.0);
        assert!(hit.is_none(), "lateral miss reported as {hit:?}");
        let hit = gjk_ray_cast(&sphere, &at(5.0, 10.0, 0.0), origin, Vec3::x(), 20.0);
        assert!(hit.is_none());
        let miss_by_a_hair = gjk_ray_cast(&sphere, &at(5.0, 1.001, 0.0), origin, Vec3::x(), 20.0);
        assert!(miss_by_a_hair.is_none());
    }

    #[test]
    fn head_on_sphere_hit_reports_distance_point_and_normal() {
        let hit = gjk_ray_cast(
            &Sphere::new(1.0),
            &at(5.0, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::x(),
            20.0,
        )
        .expect("ray hits the sphere");
        assert!((hit.t - 4.0).abs() < 1e-8, "t = {}", hit.t);
        assert!((hit.point - Vec3::new(4.0, 0.0, 0.0)).norm() < 1e-8);
        assert!((hit.normal - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-4);
    }

    #[test]
    fn hit_distance_is_independent_of_the_direction_length() {
        let hit = gjk_ray_cast(
            &Sphere::new(1.0),
            &at(5.0, 0.0, 0.0),
            Vec3::zeros(),
            Vec3::new(7.0, 0.0, 0.0),
            20.0,
        )
        .expect("ray hits the sphere");
        assert!((hit.t - 4.0).abs() < 1e-8, "t = {}", hit.t);
    }

    #[test]
    fn ray_limit_is_respected() {
        let sphere = Sphere::new(1.0);
        let t = at(100.0, 0.0, 0.0);
        assert!(gjk_ray_cast(&sphere, &t, Vec3::zeros(), Vec3::x(), 5.0).is_none());
        assert!(gjk_ray_cast(&sphere, &t, Vec3::zeros(), Vec3::x(), 98.9).is_none());
        let hit = gjk_ray_cast(&sphere, &t, Vec3::zeros(), Vec3::x(), 99.1).expect("within reach");
        assert!((hit.t - 99.0).abs() < 1e-7, "t = {}", hit.t);
        assert!(gjk_ray_cast(&sphere, &t, Vec3::zeros(), Vec3::x(), -1.0).is_none());
    }

    #[test]
    fn ray_pointing_away_from_the_shape_misses() {
        let sphere = Sphere::new(1.0);
        let t = at(5.0, 0.0, 0.0);
        assert!(gjk_ray_cast(&sphere, &t, Vec3::zeros(), -Vec3::x(), 50.0).is_none());
        assert!(gjk_ray_cast(&sphere, &t, Vec3::zeros(), Vec3::zeros(), 50.0).is_none());
    }

    #[test]
    fn ray_starting_inside_hits_immediately() {
        let hit = gjk_ray_cast(
            &BoxShape::new(Vec3::new(1.0, 1.0, 1.0)),
            &at(0.0, 0.0, 0.0),
            Vec3::new(0.2, 0.1, -0.3),
            Vec3::y(),
            10.0,
        )
        .expect("origin inside the box");
        assert_eq!(hit.t, 0.0);
        assert!((hit.normal + Vec3::y()).norm() < 1e-12);
        assert!((hit.point - Vec3::new(0.2, 0.1, -0.3)).norm() < 1e-12);
    }

    #[test]
    fn rotated_box_face_hit_has_the_face_normal() {
        let diamond = Transform::from_axis_angle(
            Vec3::new(5.0, 0.0, 0.0),
            Vec3::z(),
            std::f64::consts::FRAC_PI_4,
        );
        let hit = gjk_ray_cast(
            &BoxShape::new(Vec3::new(1.0, 1.0, 1.0)),
            &diamond,
            Vec3::new(0.0, 0.3, 0.0),
            Vec3::x(),
            20.0,
        )
        .expect("ray hits the diamond");
        let expected_t = 5.0 - (std::f64::consts::SQRT_2 - 0.3);
        assert!((hit.t - expected_t).abs() < 1e-8, "t = {}", hit.t);
        let expected_normal = Vec3::new(-1.0, 1.0, 0.0).normalize();
        assert!(
            (hit.normal - expected_normal).norm() < 1e-8,
            "{:?}",
            hit.normal
        );
    }

    /// Compare the support-mapping cast with the shape's own ray cast.
    fn compare_with_analytic(shape: &dyn Shape, seed: u64) {
        let mut rng = Lcg(seed);
        let mut compared = 0;
        for _ in 0..400 {
            let rotation_axis = rng.vector(1.0);
            let transform = if rotation_axis.norm() > 0.1 {
                Transform::from_axis_angle(rng.vector(0.5), rotation_axis, rng.range(0.0, 3.0))
            } else {
                at(0.0, 0.0, 0.0)
            };
            let origin = rng.vector(4.0);
            let dir = rng.vector(1.0);
            if dir.norm() < 0.1 {
                continue;
            }
            let unit = dir.normalize();
            let local_origin = transform.inverse().transform_point(&origin);
            let local_dir = transform.rotation.inverse() * unit;
            // Skip rays that start inside: the analytic cast reports the exit.
            let expected = shape.ray_cast(&local_origin, &local_dir, 50.0);
            let started_inside = Gjk::intersect(
                shape,
                &transform,
                &Sphere::new(1e-9),
                &Transform::from_position(origin),
            );
            if started_inside {
                continue;
            }
            let got = gjk_ray_cast(shape, &transform, origin, dir, 50.0);
            match (expected, got) {
                (None, None) => {}
                (Some(e), Some(g)) => {
                    compared += 1;
                    assert!((e.toi - g.t).abs() < 1e-6, "t {} vs {}", e.toi, g.t);
                    let grazing = e.normal.dot(&local_dir).abs() < 1e-2;
                    if !grazing {
                        let expected_normal = transform.transform_vector(&e.normal).normalize();
                        assert!(
                            (expected_normal - g.normal).norm() < 2e-3,
                            "normal {expected_normal:?} vs {:?}",
                            g.normal
                        );
                    }
                }
                (e, g) => {
                    // A grazing ray may sit on either side of the numerical edge.
                    let near_edge = e
                        .as_ref()
                        .is_some_and(|hit| hit.normal.dot(&local_dir).abs() < 1e-3);
                    assert!(near_edge, "analytic {e:?} vs gjk {g:?}");
                }
            }
        }
        assert!(compared > 40, "only {compared} rays hit");
    }

    #[test]
    fn matches_the_analytic_cast_for_a_sphere() {
        compare_with_analytic(&Sphere::new(1.3), 11);
    }

    #[test]
    fn matches_the_analytic_cast_for_a_box() {
        compare_with_analytic(&BoxShape::new(Vec3::new(1.0, 0.6, 0.4)), 12);
    }

    #[test]
    fn matches_the_analytic_cast_for_a_capsule() {
        compare_with_analytic(&Capsule::new(0.5, 0.8), 13);
    }

    /// Verify hits and misses geometrically, without a reference ray cast.
    fn verify_geometrically(shape: &dyn Shape, seed: u64) {
        let mut rng = Lcg(seed);
        let probe = Sphere::new(1e-7);
        let (mut hits, mut misses) = (0, 0);
        for _ in 0..300 {
            let transform = Transform::from_axis_angle(
                rng.vector(0.5),
                Vec3::new(1.0, 0.7, 0.2),
                rng.range(0.0, 3.0),
            );
            let origin = rng.vector(4.0);
            let dir = rng.vector(1.0);
            if dir.norm() < 0.1 {
                continue;
            }
            let unit = dir.normalize();
            if Gjk::intersect(shape, &transform, &probe, &Transform::from_position(origin)) {
                continue;
            }
            let inverse_rotation = transform.rotation.inverse();
            match gjk_ray_cast(shape, &transform, origin, dir, 30.0) {
                Some(hit) => {
                    hits += 1;
                    assert!((hit.point - (origin + unit * hit.t)).norm() < 1e-9);
                    assert!((hit.normal.norm() - 1.0).abs() < 1e-9);
                    // The hit point lies on the shape ...
                    assert!(
                        Gjk::intersect(
                            shape,
                            &transform,
                            &probe,
                            &Transform::from_position(hit.point)
                        ),
                        "hit point {:?} is off the shape",
                        hit.point
                    );
                    // ... on a supporting plane that has the ray origin in front of it.
                    let top = transform
                        .transform_point(&shape.support_point(&(inverse_rotation * hit.normal)));
                    assert!(
                        (top.dot(&hit.normal) - hit.point.dot(&hit.normal)).abs() < 1e-4,
                        "plane {} vs {}",
                        top.dot(&hit.normal),
                        hit.point.dot(&hit.normal)
                    );
                    assert!(origin.dot(&hit.normal) > hit.point.dot(&hit.normal));
                }
                None => {
                    misses += 1;
                    for step in 0..=150 {
                        let sample = origin + unit * (30.0 * step as f64 / 150.0);
                        assert!(
                            !Gjk::intersect(
                                shape,
                                &transform,
                                &probe,
                                &Transform::from_position(sample)
                            ),
                            "ray reported as a miss passes through the shape at {sample:?}"
                        );
                    }
                }
            }
        }
        assert!(hits > 30 && misses > 30, "hits {hits}, misses {misses}");
    }

    #[test]
    fn hits_and_misses_are_geometrically_consistent_for_a_cylinder() {
        verify_geometrically(&Cylinder::new(0.7, 0.9), 21);
    }

    #[test]
    fn hits_and_misses_are_geometrically_consistent_for_a_cone() {
        verify_geometrically(&Cone::new(0.8, 1.0), 22);
    }

    #[test]
    fn hits_and_misses_are_geometrically_consistent_for_a_convex_hull() {
        verify_geometrically(&cube_hull(), 23);
        let mut rng = Lcg(24);
        let cloud = ConvexHull::new((0..14).map(|_| rng.vector(1.0)).collect());
        verify_geometrically(&cloud, 25);
    }
}
