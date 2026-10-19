// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Converged distance and witness points for separated GJK queries.
//!
//! The boolean GJK loop behind [`Gjk::query`](super::Gjk::query) stops as soon
//! as a separating axis is certified. At that point its simplex only bounds the
//! distance from above: the closest point of the simplex to the origin is in
//! general *not* the closest point of the Minkowski difference `A - B`, so a
//! distance or a pair of witness points read off it overestimates the gap.
//! Conservative advancement built on such a distance overshoots (tunnelling),
//! and speculative contacts get a wrong normal and anchor.
//!
//! [`converge_distance`] continues from that simplex with the distance form of
//! GJK. Each step takes the point `v` of the current simplex that is closest to
//! the origin, queries the support point `w` of `A - B` in direction `-v` and
//! stops once `|v|^2 - v.w` drops below a small fraction of `|v|^2` (van den
//! Bergen's termination criterion): `v.w / |v|` is then a proven lower bound on
//! the distance that agrees with the upper bound `|v|`. The witness points are
//! the barycentric mixes of the stored per-vertex support points of `A` and `B`.
//!
//! The test cannot be met below the rounding noise of the vertex coordinates:
//! `v` is a difference of vertices of magnitude `size`, so its direction (and
//! with it `v.w / |v|`) is only known to about `epsilon * size / |v|`. The
//! tolerance therefore has a floor of `epsilon * size^2 / |v|` ([`certifiable_gap`]),
//! and the iteration also ends when the support point repeats a vertex of the
//! simplex or the distance stops decreasing for several steps in a row, so that
//! nearly touching shapes do not burn the whole iteration budget.

use oxiphysics_core::Transform;
use oxiphysics_core::math::{Mat3, Real, Vec3};
use oxiphysics_geometry::Shape;

use super::functions::support;
use super::types::{GjkResult, Simplex, SupportPoint};
use crate::narrowphase::convex_contact::closest_to_origin;

/// Upper bound on the number of refinement steps. Polyhedral shapes terminate
/// finitely; smooth shapes (spheres, capsules) converge linearly and reach the
/// tolerance well within this budget unless they nearly touch.
const MAX_REFINEMENT_STEPS: usize = 128;
/// Relative gap `(|v|^2 - v.w) / |v|^2` below which the distance is converged.
const RELATIVE_GAP: Real = 1e-10;
/// Absolute floor of the gap, relative to the magnitude of the Minkowski
/// vertices.
const ABSOLUTE_GAP: Real = 1e-13;
/// Safety factor on the rounding-noise floor of the gap.
const NOISE_FACTOR: Real = 32.0;
/// Relative decrease of the distance per step below which a step counts as
/// making no progress (the distance is non-increasing, so a smaller decrease is
/// rounding). On smooth shapes the decrease alternates between a large and a
/// quadratically small step, so a single such step proves nothing.
pub(crate) const STALL_TOLERANCE: Real = 8.0 * Real::EPSILON;
/// Number of consecutive steps without progress that end the iteration.
pub(crate) const STALL_STEPS: usize = 3;
/// Two Minkowski vertices closer than this fraction of `size` count as one.
pub(crate) const DUPLICATE_TOLERANCE: Real = 1e-12;

/// The gap between the upper distance bound `distance = |v|` and the lower bound
/// `v.w / |v|` that counts as converged, for Minkowski vertices of magnitude
/// `size`.
///
/// It is the largest of a relative tolerance, an absolute tolerance scaled by
/// `size` and the rounding noise `epsilon * size^2 / distance` of the lower bound
/// (a direction error of `epsilon * size / distance` times a vertex of length
/// `size`). `distance` must be positive.
pub(crate) fn certifiable_gap(distance: Real, size: Real) -> Real {
    let rounding = NOISE_FACTOR * Real::EPSILON * size * size / distance;
    (RELATIVE_GAP * distance)
        .max(ABSOLUTE_GAP * size)
        .max(rounding)
}

/// Closest pair of points between two separated convex shapes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ClosestPair {
    /// Distance between the shapes: the length of `point_a - point_b`.
    pub(crate) distance: Real,
    /// Closest point on shape A (world space).
    pub(crate) point_a: Vec3,
    /// Closest point on shape B (world space).
    pub(crate) point_b: Vec3,
    /// Number of support queries spent on the refinement.
    pub(crate) iterations: usize,
}

impl ClosestPair {
    /// The pair as the `Separated` variant of a GJK result.
    pub(crate) fn into_result(self) -> GjkResult {
        GjkResult::Separated {
            point_a: self.point_a,
            point_b: self.point_b,
            distance: self.distance,
        }
    }
}

/// Barycentric weights of `target` with respect to the 1 to 4 vertices of the
/// closest feature, clamped to the feature and renormalised to sum to one. Four
/// vertices mean that the origin is enclosed; the weights then reproduce it, so
/// the two witnesses coincide.
fn feature_weights(vertices: &[SupportPoint], target: &Vec3) -> [Real; 4] {
    match vertices {
        [] => [0.0; 4],
        [_] => [1.0, 0.0, 0.0, 0.0],
        [a, b] => {
            let edge = b.point - a.point;
            let length_sq = edge.norm_squared();
            if length_sq <= 0.0 {
                return [1.0, 0.0, 0.0, 0.0];
            }
            let t = ((target - a.point).dot(&edge) / length_sq).clamp(0.0, 1.0);
            [1.0 - t, t, 0.0, 0.0]
        }
        [a, b, c] => {
            let v0 = b.point - a.point;
            let v1 = c.point - a.point;
            let v2 = target - a.point;
            let d00 = v0.dot(&v0);
            let d01 = v0.dot(&v1);
            let d11 = v1.dot(&v1);
            let d20 = v2.dot(&v0);
            let d21 = v2.dot(&v1);
            let denominator = d00 * d11 - d01 * d01;
            if denominator.abs() <= Real::MIN_POSITIVE {
                return [1.0, 0.0, 0.0, 0.0];
            }
            let v = (d11 * d20 - d01 * d21) / denominator;
            let w = (d00 * d21 - d01 * d20) / denominator;
            let u = 1.0 - v - w;
            let (u, v, w) = (u.max(0.0), v.max(0.0), w.max(0.0));
            let total = u + v + w;
            if total <= 0.0 {
                return [1.0, 0.0, 0.0, 0.0];
            }
            [u / total, v / total, w / total, 0.0]
        }
        [a, b, c, d, ..] => {
            let edges =
                Mat3::from_columns(&[b.point - a.point, c.point - a.point, d.point - a.point]);
            let Some(inverse) = edges.try_inverse() else {
                return [0.25; 4];
            };
            let local = inverse * (target - a.point);
            let raw = [1.0 - local.x - local.y - local.z, local.x, local.y, local.z];
            let clamped = raw.map(|weight| weight.max(0.0));
            let total: Real = clamped.iter().sum();
            if total <= 0.0 {
                return [0.25; 4];
            }
            clamped.map(|weight| weight / total)
        }
    }
}

/// Witness points of the closest feature `vertices` to the origin.
fn witnesses(vertices: &[SupportPoint], closest: &Vec3) -> (Vec3, Vec3) {
    let weights = feature_weights(vertices, closest);
    let mut point_a = Vec3::zeros();
    let mut point_b = Vec3::zeros();
    for (vertex, weight) in vertices.iter().zip(weights) {
        point_a += vertex.support_a * weight;
        point_b += vertex.support_b * weight;
    }
    (point_a, point_b)
}

/// Refine the distance of two separated convex shapes starting from the
/// terminating simplex of the boolean GJK search.
///
/// `start` must hold between one and four vertices of the Minkowski difference
/// `A - B` (as produced by [`support`]). The returned distance never undercuts
/// the true distance by more than rounding and exceeds it by at most
/// [`certifiable_gap`] (a relative `1e-10`, or the rounding noise of the
/// vertices for nearly touching shapes). That bound is exact for polyhedra; on
/// curved shapes the simplex can start to cycle once its closest point is only
/// determined to about `1e-8` laterally, which limits the accuracy to roughly
/// `1e-9` relative. `point_a - point_b` is the closest-point vector of length
/// `distance`.
pub(crate) fn converge_distance(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    start: &Simplex,
) -> ClosestPair {
    let seed = match start.points.first() {
        Some(point) => *point,
        None => support(shape_a, transform_a, shape_b, transform_b, &Vec3::x()),
    };
    let mut vertices = [seed; 4];
    let mut count = start.points.len().clamp(1, 4);
    for (slot, point) in vertices.iter_mut().zip(&start.points).take(count) {
        *slot = *point;
    }
    let mut size = vertices[..count]
        .iter()
        .map(|vertex| vertex.point.norm())
        .fold(0.0, Real::max);
    let mut iterations = 0;
    let mut previous_length = Real::INFINITY;
    let mut stalled_steps = 0;
    let closest = loop {
        // Closest point of the current simplex, and the vertices that carry it.
        let mut points = [[0.0; 3]; 4];
        for (slot, vertex) in points.iter_mut().zip(&vertices[..count]) {
            *slot = [vertex.point.x, vertex.point.y, vertex.point.z];
        }
        let nearest = closest_to_origin(&points[..count]);
        let mut carried = [seed; 4];
        for (slot, &index) in carried.iter_mut().zip(&nearest.keep[..nearest.len]) {
            *slot = vertices[index];
        }
        vertices = carried;
        count = nearest.len;
        let closest = Vec3::new(nearest.point[0], nearest.point[1], nearest.point[2]);
        let closest_length = closest.norm();
        // The distance never grows as vertices are added; decreases within the
        // rounding error for several steps in a row mean that the iteration
        // cannot improve any more.
        if iterations > 0 && previous_length - closest_length <= STALL_TOLERANCE * previous_length {
            stalled_steps += 1;
        } else {
            stalled_steps = 0;
        }
        // An enclosed origin, a vanishing distance, a stall or a non-finite state
        // leave nothing to refine.
        if count == 4
            || !closest_length.is_finite()
            || closest_length <= ABSOLUTE_GAP * size
            || stalled_steps >= STALL_STEPS
            || iterations >= MAX_REFINEMENT_STEPS
        {
            break closest;
        }
        previous_length = closest_length;
        iterations += 1;
        let direction = -closest / closest_length;
        let extreme = support(shape_a, transform_a, shape_b, transform_b, &direction);
        size = size.max(extreme.point.norm());
        // `extreme` is the point of A - B furthest along `direction`, so the
        // plane through it with normal `direction` bounds A - B: the distance is
        // at least `lower`, while `closest_length` is an upper bound.
        let lower = -extreme.point.dot(&direction);
        let repeated = vertices[..count]
            .iter()
            .any(|vertex| (vertex.point - extreme.point).norm() <= DUPLICATE_TOLERANCE * size);
        if closest_length - lower <= certifiable_gap(closest_length, size) || repeated {
            break closest;
        }
        vertices[count] = extreme;
        count += 1;
    };
    let (point_a, point_b) = witnesses(&vertices[..count], &closest);
    ClosestPair {
        distance: (point_a - point_b).norm(),
        point_a,
        point_b,
        iterations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::narrowphase::gjk::{Gjk, WarmStartGjk, gjk_distance_query};
    use oxiphysics_geometry::{BoxShape, Capsule, ConvexHull, Sphere};

    const TOLERANCE: Real = 1e-9;
    /// How far a witness may sit behind the supporting plane of its shape. On
    /// curved shapes the witnesses are barycentric mixes of surface points, which
    /// lie inside the shape by the sagitta of the (tiny) simplex.
    const SUPPORT_TOLERANCE: Real = 1e-7;

    fn at(x: Real, y: Real, z: Real) -> Transform {
        Transform::from_position(Vec3::new(x, y, z))
    }

    fn separated(result: GjkResult) -> (Real, Vec3, Vec3) {
        match result {
            GjkResult::Separated {
                point_a,
                point_b,
                distance,
            } => (distance, point_a, point_b),
            GjkResult::Intersecting(_) => panic!("the shapes must be separated"),
        }
    }

    fn assert_close(actual: Vec3, expected: Vec3) {
        assert!(
            (actual - expected).norm() <= TOLERANCE,
            "{actual:?} differs from {expected:?}"
        );
    }

    /// World-space point of `shape` that is furthest along `direction`.
    fn world_support(shape: &dyn Shape, transform: &Transform, direction: &Vec3) -> Vec3 {
        let local = shape.support_point(&(transform.rotation.inverse() * direction));
        transform.transform_point(&local)
    }

    /// Prove that `(point_a, point_b)` is the closest pair of the two shapes.
    ///
    /// Both witnesses lie on their shapes, and the planes through them with the
    /// normal `n = (point_b - point_a) / distance` support the shapes from
    /// opposite sides, so no pair of points can be closer than `distance`.
    fn assert_closest_pair(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
        pair: (Real, Vec3, Vec3),
    ) {
        let (distance, point_a, point_b) = pair;
        assert!(distance > 0.0, "distance {distance}");
        assert!(
            ((point_a - point_b).norm() - distance).abs() <= 1e-9 * (1.0 + distance),
            "witness gap {} vs distance {distance}",
            (point_a - point_b).norm()
        );
        let probe = Sphere::new(1e-7);
        assert!(
            Gjk::intersect(
                shape_a,
                transform_a,
                &probe,
                &at(point_a.x, point_a.y, point_a.z)
            ),
            "witness {point_a:?} is off shape A"
        );
        assert!(
            Gjk::intersect(
                shape_b,
                transform_b,
                &probe,
                &at(point_b.x, point_b.y, point_b.z)
            ),
            "witness {point_b:?} is off shape B"
        );
        let n = (point_b - point_a) / distance;
        let top_a = n.dot(&world_support(shape_a, transform_a, &n));
        let bottom_b = n.dot(&world_support(shape_b, transform_b, &(-n)));
        assert!(
            (top_a - n.dot(&point_a)).abs() <= SUPPORT_TOLERANCE,
            "A extends {} beyond its witness along the closest-point normal",
            top_a - n.dot(&point_a)
        );
        assert!(
            (bottom_b - n.dot(&point_b)).abs() <= SUPPORT_TOLERANCE,
            "B extends {} beyond its witness along the closest-point normal",
            n.dot(&point_b) - bottom_b
        );
    }

    #[test]
    fn unit_spheres_five_apart_have_the_true_gap_and_witnesses() {
        let sphere = Sphere::new(1.0);
        let (distance, point_a, point_b) = separated(Gjk::query(
            &sphere,
            &at(0.0, 0.0, 0.0),
            &sphere,
            &at(5.0, 0.0, 0.0),
        ));
        assert!((distance - 3.0).abs() <= TOLERANCE, "distance {distance}");
        assert_close(point_a, Vec3::new(1.0, 0.0, 0.0));
        assert_close(point_b, Vec3::new(4.0, 0.0, 0.0));
    }

    #[test]
    fn spheres_on_a_diagonal_have_witnesses_on_the_line_of_centres() {
        let big = Sphere::new(1.0);
        let small = Sphere::new(0.5);
        let centre = Vec3::new(3.0, 3.0, 3.0);
        let (distance, point_a, point_b) = separated(Gjk::query(
            &big,
            &at(0.0, 0.0, 0.0),
            &small,
            &Transform::from_position(centre),
        ));
        let axis = centre.normalize();
        assert!(
            (distance - (centre.norm() - 1.5)).abs() <= TOLERANCE,
            "distance {distance}"
        );
        assert_close(point_a, axis);
        assert_close(point_b, centre - axis * 0.5);
    }

    #[test]
    fn staggered_boxes_are_two_apart_with_witnesses_on_the_facing_sides() {
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let (distance, point_a, point_b) = separated(Gjk::query(
            &cube,
            &at(0.0, 0.0, 0.0),
            &cube,
            &at(4.0, 0.5, 0.0),
        ));
        assert!((distance - 2.0).abs() <= TOLERANCE, "distance {distance}");
        assert!((point_a.x - 1.0).abs() <= TOLERANCE, "{point_a:?}");
        assert!((point_b.x - 3.0).abs() <= TOLERANCE, "{point_b:?}");
        // The facing sides overlap in y in [-0.5, 1] and in z in [-1, 1].
        assert!((point_a.y - point_b.y).abs() <= TOLERANCE);
        assert!((point_a.z - point_b.z).abs() <= TOLERANCE);
        assert!(point_a.y >= -0.5 - TOLERANCE && point_a.y <= 1.0 + TOLERANCE);
        assert!(point_a.z.abs() <= 1.0 + TOLERANCE);
    }

    #[test]
    fn diamond_against_a_sphere_touches_with_its_corner() {
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let diamond = Transform::from_axis_angle(
            Vec3::zeros(),
            Vec3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_4,
        );
        let sphere = Sphere::new(0.5);
        let (distance, point_a, point_b) =
            separated(Gjk::query(&cube, &diamond, &sphere, &at(4.0, 0.0, 0.0)));
        let corner = std::f64::consts::SQRT_2;
        assert!(
            (distance - (4.0 - corner - 0.5)).abs() <= TOLERANCE,
            "distance {distance}"
        );
        assert_close(point_a, Vec3::new(corner, 0.0, 0.0));
        assert_close(point_b, Vec3::new(3.5, 0.0, 0.0));
    }

    #[test]
    fn boxes_meeting_edge_to_edge_are_separated_by_the_diagonal_gap() {
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let (distance, point_a, point_b) = separated(Gjk::query(
            &cube,
            &at(0.0, 0.0, 0.0),
            &cube,
            &at(3.5, 3.5, 0.0),
        ));
        assert!(
            (distance - 1.5 * std::f64::consts::SQRT_2).abs() <= TOLERANCE,
            "distance {distance}"
        );
        assert!((point_a.x - 1.0).abs() <= TOLERANCE && (point_a.y - 1.0).abs() <= TOLERANCE);
        assert!((point_b.x - 2.5).abs() <= TOLERANCE && (point_b.y - 2.5).abs() <= TOLERANCE);
        assert!((point_a.z - point_b.z).abs() <= TOLERANCE);
    }

    /// Small deterministic generator for the randomised sweep.
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

        fn shape(&mut self, kind: usize) -> Box<dyn Shape> {
            match kind % 4 {
                0 => Box::new(BoxShape::new(Vec3::new(
                    self.range(0.3, 1.0),
                    self.range(0.3, 1.0),
                    self.range(0.3, 1.0),
                ))),
                1 => Box::new(Capsule::new(self.range(0.2, 0.6), self.range(0.3, 0.9))),
                2 => Box::new(Sphere::new(self.range(0.3, 0.9))),
                _ => Box::new(ConvexHull::new(
                    (0..9).map(|_| self.vector(0.9)).collect::<Vec<_>>(),
                )),
            }
        }
    }

    #[test]
    fn random_pairs_satisfy_the_closest_pair_certificate() {
        let mut rng = Lcg(0x5EED_0178);
        let mut certified = 0;
        for case in 0..240 {
            let shape_a = rng.shape(case);
            let shape_b = rng.shape(case / 4 + 1);
            let transform_a = Transform::from_axis_angle(
                Vec3::zeros(),
                rng.vector(1.0) + Vec3::new(0.0, 0.0, 0.3),
                rng.range(0.0, 3.0),
            );
            let transform_b = Transform::from_axis_angle(
                rng.vector(1.0) + Vec3::new(3.0, 0.0, 0.0),
                rng.vector(1.0) + Vec3::new(0.0, 0.3, 0.0),
                rng.range(0.0, 3.0),
            );
            let result = Gjk::query(
                shape_a.as_ref(),
                &transform_a,
                shape_b.as_ref(),
                &transform_b,
            );
            let GjkResult::Separated {
                point_a,
                point_b,
                distance,
            } = result
            else {
                continue;
            };
            certified += 1;
            assert_closest_pair(
                shape_a.as_ref(),
                &transform_a,
                shape_b.as_ref(),
                &transform_b,
                (distance, point_a, point_b),
            );
        }
        assert!(certified >= 150, "only {certified} separated pairs");
    }

    #[test]
    fn distance_query_returns_the_same_converged_pair() {
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let capsule = Capsule::new(0.4, 0.7);
        let rotated =
            Transform::from_axis_angle(Vec3::new(2.5, 1.5, -0.5), Vec3::new(0.3, 1.0, 0.2), 0.8);
        let query = separated(Gjk::query(&cube, &at(0.0, 0.0, 0.0), &capsule, &rotated));
        let detailed = gjk_distance_query(&cube, &at(0.0, 0.0, 0.0), &capsule, &rotated);
        assert!((detailed.distance - query.0).abs() <= 1e-12);
        assert_close(detailed.closest_a, query.1);
        assert_close(detailed.closest_b, query.2);
        assert!(detailed.iterations > 0);
        assert_closest_pair(
            &cube,
            &at(0.0, 0.0, 0.0),
            &capsule,
            &rotated,
            (detailed.distance, detailed.closest_a, detailed.closest_b),
        );
    }

    #[test]
    fn warm_started_queries_stay_converged() {
        let mut warm = WarmStartGjk::new();
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        let sphere = Sphere::new(0.75);
        let moved =
            Transform::from_axis_angle(Vec3::new(3.0, 2.0, 1.0), Vec3::new(1.0, 1.0, 0.0), 0.5);
        // The sphere centre (3, 2, 1) sits above the cube edge x = y = 1, so the
        // gap is |(2, 1, 0)| minus the radius.
        let exact = 5.0_f64.sqrt() - 0.75;
        for _ in 0..3 {
            let pair = separated(warm.query(&cube, &at(0.0, 0.0, 0.0), &sphere, &moved));
            assert!((pair.0 - exact).abs() < 1e-8, "{} vs {exact}", pair.0);
            assert!(pair.0 >= exact - 1e-12, "{} undercuts {exact}", pair.0);
            assert_closest_pair(&cube, &at(0.0, 0.0, 0.0), &sphere, &moved, pair);
        }
    }

    #[test]
    fn overlap_simplex_is_some_only_for_overlapping_pairs() {
        let cube = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
        assert!(
            Gjk::overlap_simplex(&cube, &at(0.0, 0.0, 0.0), &cube, &at(1.5, 0.0, 0.0)).is_some()
        );
        assert!(
            Gjk::overlap_simplex(&cube, &at(0.0, 0.0, 0.0), &cube, &at(2.5, 0.0, 0.0)).is_none()
        );
    }

    #[test]
    fn empty_and_enclosing_start_simplices_are_handled() {
        let sphere = Sphere::new(1.0);
        let far = at(5.0, 0.0, 0.0);
        let pair = converge_distance(&sphere, &at(0.0, 0.0, 0.0), &sphere, &far, &Simplex::new());
        assert!(
            (pair.distance - 3.0).abs() <= TOLERANCE,
            "distance {}",
            pair.distance
        );

        // A simplex that encloses the origin leaves nothing to refine: the
        // witnesses coincide.
        let near = at(1.0, 0.0, 0.0);
        let GjkResult::Intersecting(enclosing) =
            Gjk::query(&sphere, &at(0.0, 0.0, 0.0), &sphere, &near)
        else {
            panic!("spheres one apart overlap");
        };
        let pair = converge_distance(&sphere, &at(0.0, 0.0, 0.0), &sphere, &near, &enclosing);
        assert!(pair.distance <= 1e-9, "distance {}", pair.distance);
    }
}
