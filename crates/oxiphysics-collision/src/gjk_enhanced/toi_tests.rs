//! Regression tests for the GJK distance accuracy and the conservative
//! advancement time of impact.

use super::functions::{
    ConvexShape, closest_point_on_box, closest_point_on_segment, gjk_enhanced, gjk_enhanced_core,
    gjk_time_of_impact, vadd, vdot, vlen, vscale, vsub,
};
use super::types::{ConvexHull, EnhBox, EnhCapsule, EnhSphere, ToiResult, TranslatedShape};

/// Small deterministic linear congruential generator for reproducible scenes.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Uniform sample in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as f64 / (1u64 << 31) as f64
    }

    /// Uniform sample in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn point(&mut self, extent: f64) -> [f64; 3] {
        [
            self.range(-extent, extent),
            self.range(-extent, extent),
            self.range(-extent, extent),
        ]
    }
}

/// Exact signed distance between a sphere and an axis-aligned box (negative if
/// they overlap).
fn sphere_box_distance(sphere: &EnhSphere, bx: &EnhBox) -> f64 {
    let closest = closest_point_on_box(bx, sphere.centre);
    vlen(vsub(sphere.centre, closest)) - sphere.radius
}

/// Exact distance between two axis-aligned boxes (zero if they overlap).
fn box_box_distance(a: &EnhBox, b: &EnhBox) -> f64 {
    let mut sum = 0.0;
    for k in 0..3 {
        let gap = ((a.centre[k] - b.centre[k]).abs() - a.half[k] - b.half[k]).max(0.0);
        sum += gap * gap;
    }
    sum.sqrt()
}

/// Exact signed distance between a capsule and a sphere.
fn capsule_sphere_distance(capsule: &EnhCapsule, sphere: &EnhSphere) -> f64 {
    let (closest, _) = closest_point_on_segment(capsule.a, capsule.b, sphere.centre);
    vlen(vsub(sphere.centre, closest)) - capsule.radius - sphere.radius
}

#[test]
fn off_axis_boxes_report_the_true_distance() {
    // The y ranges overlap, so the gap is purely along x: exactly 3.0. A GJK
    // that stops at the first support point reports 3.45 here.
    let a = EnhBox::new([0.0; 3], [1.0; 3]);
    let b = EnhBox::new([5.0, 0.3, 0.0], [1.0; 3]);
    let result = gjk_enhanced(&a, &b, None, 64);
    assert!(!result.intersecting);
    assert!(
        (result.distance() - 3.0).abs() < 1e-9,
        "distance {}",
        result.distance()
    );
    // Diagonal separation: gaps 3 (x) and 1 (y).
    let c = EnhBox::new([5.0, 3.0, 0.0], [1.0; 3]);
    let result = gjk_enhanced(&a, &c, None, 64);
    assert!((result.distance() - 10.0_f64.sqrt()).abs() < 1e-9);
}

#[test]
fn sphere_box_distance_matches_closed_form() {
    let mut rng = Lcg::new(99);
    let mut separated = 0;
    for _ in 0..300 {
        let sphere = EnhSphere::new(rng.point(5.0), rng.range(0.2, 1.0));
        let bx = EnhBox::new(
            rng.point(2.0),
            [
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
            ],
        );
        let exact = sphere_box_distance(&sphere, &bx);
        if exact.abs() < 1e-4 {
            continue;
        }
        let result = gjk_enhanced(&sphere, &bx, None, 64);
        assert_eq!(result.intersecting, exact < 0.0, "exact distance {exact}");
        if exact > 0.0 {
            separated += 1;
            assert!(
                (result.distance() - exact).abs() < 1e-6,
                "gjk {} vs exact {exact}",
                result.distance()
            );
        }
    }
    assert!(separated > 100, "only {separated} separated samples");
}

#[test]
fn box_box_and_capsule_sphere_distances_match_closed_form() {
    let mut rng = Lcg::new(7);
    for _ in 0..300 {
        let a = EnhBox::new(
            rng.point(3.0),
            [
                rng.range(0.3, 1.2),
                rng.range(0.3, 1.2),
                rng.range(0.3, 1.2),
            ],
        );
        let b = EnhBox::new(
            rng.point(3.0),
            [
                rng.range(0.3, 1.2),
                rng.range(0.3, 1.2),
                rng.range(0.3, 1.2),
            ],
        );
        let exact = box_box_distance(&a, &b);
        let result = gjk_enhanced(&a, &b, None, 64);
        if exact > 1e-4 {
            assert!(!result.intersecting);
            assert!((result.distance() - exact).abs() < 1e-6);
        } else if exact == 0.0 {
            assert!(result.intersecting || result.distance() < 1e-6);
        }

        let capsule = EnhCapsule::new(rng.point(3.0), rng.point(3.0), rng.range(0.1, 0.6));
        let sphere = EnhSphere::new(rng.point(4.0), rng.range(0.2, 0.8));
        let exact = capsule_sphere_distance(&capsule, &sphere);
        if exact.abs() < 1e-4 {
            continue;
        }
        let result = gjk_enhanced(&capsule, &sphere, None, 64);
        assert_eq!(result.intersecting, exact < 0.0, "exact {exact}");
        if exact > 0.0 {
            assert!(
                (result.distance() - exact).abs() < 1e-6,
                "gjk {} vs exact {exact}",
                result.distance()
            );
        }
    }
}

#[test]
fn separation_certificate_is_a_sound_lower_bound() {
    let mut rng = Lcg::new(31);
    let hull = ConvexHull::new(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    let mut checked = 0;
    for _ in 0..300 {
        let sphere = EnhSphere::new(rng.point(4.0), rng.range(0.1, 0.6));
        let core = gjk_enhanced_core(&sphere, &hull, None, 64);
        if core.result.intersecting {
            assert_eq!(core.gap, 0.0);
            continue;
        }
        checked += 1;
        let upper = core.result.distance();
        assert!(core.gap > 0.0 || upper < 1e-6);
        // Lower bound: never above the distance between the closest points.
        assert!(core.gap <= upper + 1e-9, "gap {} > dist {upper}", core.gap);
        // Converged: the bound is tight.
        assert!(core.gap >= upper - 1e-6, "gap {} << dist {upper}", core.gap);
        assert!((vlen(core.axis) - 1.0).abs() < 1e-9);
        // The axis really separates the shapes: all hull vertices are at least
        // `gap` ahead of every sphere point along the axis.
        let sphere_reach = vdot(core.axis, sphere.centre) + sphere.radius;
        for vertex in &hull.vertices {
            assert!(vdot(core.axis, *vertex) >= sphere_reach + core.gap - 1e-9);
        }
    }
    assert!(checked > 100);
}

#[test]
fn fast_sphere_does_not_tunnel_through_box() {
    // At t = 0 the sphere is 4.4 left of the box, at t = 1 it is 4.4 right of
    // it: both end points are separated, the path crosses the box.
    let sphere = EnhSphere::new([-5.0, 0.0, 0.0], 0.1);
    let wall = EnhBox::new([0.0; 3], [0.5; 3]);
    let result = gjk_time_of_impact(&sphere, &wall, [10.0, 0.0, 0.0], 64, 1e-3);
    assert!(result.hit, "tunnelled through the box: {result:?}");
    // Contact when the sphere centre reaches x = -0.6, i.e. t = 0.44.
    assert!(result.toi <= 0.44 + 1e-9, "toi {} is late", result.toi);
    assert!((result.toi - 0.44).abs() < 2e-4, "toi {}", result.toi);
    assert!(result.normal[0] < -0.99, "normal {:?}", result.normal);
}

#[test]
fn thin_wall_is_not_tunnelled() {
    let mover = EnhBox::new([-3.0, 0.0, 0.0], [0.1; 3]);
    let wall = EnhBox::new([0.0; 3], [0.01, 5.0, 5.0]);
    let result = gjk_time_of_impact(&mover, &wall, [6.0, 0.0, 0.0], 64, 1e-4);
    assert!(result.hit);
    let exact = (3.0 - 0.1 - 0.01) / 6.0;
    assert!(result.toi <= exact + 1e-9);
    assert!((result.toi - exact).abs() < 1e-4, "toi {}", result.toi);
}

#[test]
fn passing_beside_the_obstacle_is_a_miss() {
    // Sphere passing 3 units beside a unit-radius sphere.
    let mover = EnhSphere::new([-5.0, 3.0, 0.0], 0.5);
    let target = EnhSphere::new([0.0; 3], 0.5);
    let result = gjk_time_of_impact(&mover, &target, [10.0, 0.0, 0.0], 64, 1e-3);
    assert!(!result.hit);
    assert_eq!(result.toi, 1.0);

    // Diagonal path that clears the corner of a box; both end points are far
    // from the box as well.
    let mover = EnhSphere::new([-5.0, -8.0, 0.0], 0.2);
    let target = EnhBox::new([0.0; 3], [1.0; 3]);
    let result = gjk_time_of_impact(&mover, &target, [10.0, 10.0, 0.0], 64, 1e-3);
    assert!(!result.hit, "{result:?}");
}

#[test]
fn diverging_stationary_and_invalid_motion_miss() {
    let mover = EnhSphere::new([0.0; 3], 0.5);
    let target = EnhSphere::new([5.0, 0.0, 0.0], 0.5);
    for velocity in [[-1.0, 0.0, 0.0], [0.0; 3], [0.0, 3.0, 0.0]] {
        let result = gjk_time_of_impact(&mover, &target, velocity, 32, 1e-3);
        assert!(!result.hit, "velocity {velocity:?}");
        assert_eq!(result.toi, 1.0);
    }
    let result = gjk_time_of_impact(&mover, &target, [f64::NAN, 0.0, 0.0], 32, 1e-3);
    assert!(!result.hit);
}

#[test]
fn too_short_motion_is_a_miss() {
    // The gap is 4; moving 3 units never reaches the target.
    let mover = EnhSphere::new([0.0; 3], 0.5);
    let target = EnhSphere::new([5.0, 0.0, 0.0], 0.5);
    let result = gjk_time_of_impact(&mover, &target, [3.0, 0.0, 0.0], 32, 1e-3);
    assert!(!result.hit);
}

#[test]
fn initially_overlapping_or_touching_shapes_hit_at_time_zero() {
    let a = EnhSphere::new([0.0; 3], 1.0);
    let overlapping = EnhSphere::new([1.0, 0.0, 0.0], 1.0);
    let result = gjk_time_of_impact(&a, &overlapping, [1.0, 0.0, 0.0], 32, 1e-3);
    assert!(result.hit);
    assert_eq!(result.toi, 0.0);
    // Within the tolerance of touching.
    let near = EnhSphere::new([2.0005, 0.0, 0.0], 1.0);
    let result = gjk_time_of_impact(&a, &near, [1.0, 0.0, 0.0], 32, 1e-3);
    assert!(result.hit);
    assert_eq!(result.toi, 0.0);
}

#[test]
fn exhausted_iteration_budget_stays_conservative() {
    let sphere = EnhSphere::new([-5.0, 0.0, 0.0], 0.1);
    let wall = EnhBox::new([0.0; 3], [0.5; 3]);
    // One iteration only advances the time once; the reported time must still
    // not be later than the true time of impact.
    let result = gjk_time_of_impact(&sphere, &wall, [10.0, 0.0, 0.0], 1, 1e-3);
    assert!(result.hit);
    assert!(result.toi <= 0.44 + 1e-9);
}

/// First time on a fine grid at which `a` swept by `velocity` intersects `b`.
fn first_grid_contact(
    a: &dyn ConvexShape,
    b: &dyn ConvexShape,
    velocity: [f64; 3],
    steps: u32,
) -> Option<f64> {
    (0..=steps).find_map(|k| {
        let t = k as f64 / steps as f64;
        let moved = TranslatedShape {
            shape: a,
            offset: vscale(velocity, t),
        };
        gjk_enhanced(&moved, b, None, 64).intersecting.then_some(t)
    })
}

fn distance_at(a: &dyn ConvexShape, b: &dyn ConvexShape, velocity: [f64; 3], t: f64) -> f64 {
    let moved = TranslatedShape {
        shape: a,
        offset: vscale(velocity, t),
    };
    gjk_enhanced(&moved, b, None, 64).distance()
}

#[test]
fn time_of_impact_is_conservative_and_tight_on_random_sweeps() {
    let mut rng = Lcg::new(2025);
    let tolerance = 1e-3;
    let steps = 2000;
    let target_box = EnhBox::new([0.0; 3], [1.0, 0.7, 1.3]);
    let target_hull = ConvexHull::new(vec![
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.2],
        [0.7, 0.7, 0.7],
    ]);
    let mut hits = 0;
    let mut misses = 0;
    for scenario in 0..120 {
        let mover = EnhSphere::new(rng.point(5.0), rng.range(0.1, 0.5));
        // Aim at a random point around the target and fall short of it or
        // overshoot it, so the sweeps are a mix of clean hits, grazing passes,
        // near misses and sweeps that are too short to arrive.
        let aim = rng.point(2.2);
        let velocity = vscale(vsub(aim, mover.centre), rng.range(0.6, 2.2));
        let target: &dyn ConvexShape = if scenario % 2 == 0 {
            &target_box
        } else {
            &target_hull
        };
        let result = gjk_time_of_impact(&mover, target, velocity, 64, tolerance);
        let reference = first_grid_contact(&mover, target, velocity, steps);
        if result.hit {
            hits += 1;
            // Within the tolerance of the target at the reported time ...
            assert!(
                distance_at(&mover, target, velocity, result.toi) <= tolerance + 1e-6,
                "scenario {scenario}: not within tolerance at {}",
                result.toi
            );
            // ... and never later than the first actual contact.
            if let Some(first) = reference {
                assert!(
                    result.toi <= first + 1e-9,
                    "scenario {scenario}: toi {} after contact {first}",
                    result.toi
                );
            }
            // Nothing in front of the reported time touches already.
            for k in 0..=steps {
                let t = k as f64 / steps as f64;
                if t >= result.toi {
                    break;
                }
                assert!(
                    distance_at(&mover, target, velocity, t) > 0.0,
                    "scenario {scenario}: contact at {t} precedes toi {}",
                    result.toi
                );
            }
        } else {
            misses += 1;
            assert!(
                reference.is_none(),
                "scenario {scenario}: missed a contact at {reference:?}"
            );
            assert_eq!(result.toi, 1.0);
        }
    }
    assert!(hits > 15 && misses > 15, "hits {hits}, misses {misses}");
}

#[test]
fn returns_the_tolerance_neighbourhood_of_the_impact_for_moving_pairs() {
    // Both shapes move; only the relative displacement matters.
    let a = EnhSphere::new([0.0; 3], 0.5);
    let b = EnhSphere::new([6.0, 0.0, 0.0], 0.5);
    let vel_a = [4.0, 0.0, 0.0];
    let vel_b = [-2.0, 0.0, 0.0];
    let relative = vsub(vel_a, vel_b);
    let result: ToiResult = gjk_time_of_impact(&a, &b, relative, 64, 1e-6);
    // Gap 5 closed at speed 6.
    assert!(result.hit);
    assert!((result.toi - 5.0 / 6.0).abs() < 1e-5, "toi {}", result.toi);
    let at_impact = vadd(a.centre, vscale(vel_a, result.toi));
    let b_at_impact = vadd(b.centre, vscale(vel_b, result.toi));
    assert!((vlen(vsub(b_at_impact, at_impact)) - 1.0).abs() < 1e-5);
}
