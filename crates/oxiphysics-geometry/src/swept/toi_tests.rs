//! Regression tests for the swept time-of-impact queries: initial overlap,
//! contact windows that a sampled/bisected search tunnels through, degenerate
//! (zero-length) capsules, and a brute-force cross-check of the closed-form
//! sweeps.

use super::functions::*;
use super::types::{SweptCapsule, SweptSphere};

/// Minimal deterministic generator (64-bit LCG) so the randomised cross-checks
/// are reproducible without pulling in an RNG dependency.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn point(&mut self, half_extent: f64) -> [f64; 3] {
        [
            self.range(-half_extent, half_extent),
            self.range(-half_extent, half_extent),
            self.range(-half_extent, half_extent),
        ]
    }
}

// ---------------------------------------------------------------------------
// F1687: initially overlapping spheres report t = 0, never the exit time
// ---------------------------------------------------------------------------

#[test]
fn overlapping_sphere_moving_apart_reports_zero() {
    // Starts 0.1 from the centre of a radius-1 combined ball and leaves it at
    // t ~ 0.0909; the old code returned that exit time.
    let swept = SweptSphere::new([0.1, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5);
    assert_eq!(swept_sphere_vs_sphere(&swept, [0.0; 3], 0.5), Some(0.0));
}

#[test]
fn overlapping_sphere_exiting_after_the_sweep_reports_zero() {
    // The exit time (t = 5) is far beyond the sweep: the old code returned
    // `None` for a pair that is penetrating at t = 0.
    let swept = SweptSphere::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0);
    assert_eq!(swept_sphere_vs_sphere(&swept, [0.0; 3], 2.0), Some(0.0));
}

#[test]
fn separated_sphere_sweeps_are_unchanged() {
    let hit = SweptSphere::new([-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5);
    let t = swept_sphere_vs_sphere(&hit, [0.0; 3], 0.5).expect("head-on sweep must hit");
    assert!((t - 0.4).abs() < 1e-12, "t = {t}");

    // Moving away from the sphere and a pure miss.
    let away = SweptSphere::new([-5.0, 0.0, 0.0], [-9.0, 0.0, 0.0], 0.5);
    assert_eq!(swept_sphere_vs_sphere(&away, [0.0; 3], 0.5), None);
    let miss = SweptSphere::new([-5.0, 3.0, 0.0], [5.0, 3.0, 0.0], 0.5);
    assert_eq!(swept_sphere_vs_sphere(&miss, [0.0; 3], 0.5), None);

    // Contact beyond the end of the sweep.
    let short = SweptSphere::new([-5.0, 0.0, 0.0], [-3.0, 0.0, 0.0], 0.5);
    assert_eq!(swept_sphere_vs_sphere(&short, [0.0; 3], 0.5), None);

    // Grazing contact: passes at exactly the combined radius.
    let graze = SweptSphere::new([-5.0, 1.0, 0.0], [5.0, 1.0, 0.0], 0.5);
    let t = swept_sphere_vs_sphere(&graze, [0.0; 3], 0.5).expect("tangent sweep touches");
    assert!((t - 0.5).abs() < 1e-12, "t = {t}");
}

#[test]
fn stationary_sphere_overlap_and_separation() {
    let inside = SweptSphere::new([0.2; 3], [0.2; 3], 1.0);
    assert_eq!(swept_sphere_vs_sphere(&inside, [0.0; 3], 1.0), Some(0.0));
    let apart = SweptSphere::new([9.0; 3], [9.0; 3], 1.0);
    assert_eq!(swept_sphere_vs_sphere(&apart, [0.0; 3], 1.0), None);
}

#[test]
fn swept_sphere_vs_swept_sphere_initial_overlap_is_zero() {
    let t = swept_sphere_vs_swept_sphere(
        [0.0; 3],
        [10.0, 0.0, 0.0],
        1.0,
        [0.5, 0.0, 0.0],
        [-10.0, 0.0, 0.0],
        1.0,
    );
    assert_eq!(t, Some(0.0));
}

#[test]
fn swept_sphere_vs_swept_sphere_head_on_contact_time() {
    // Gap closes at 20 units per unit time from a 20-unit separation; contact
    // when the centres are 1 apart: t = 19 / 20.
    let t = swept_sphere_vs_swept_sphere(
        [-10.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        0.5,
        [10.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        0.5,
    )
    .expect("head-on spheres collide");
    assert!((t - 0.95).abs() < 1e-12, "t = {t}");
}

#[test]
fn swept_capsule_overlap_at_start_reports_zero() {
    let capsule = SweptCapsule::new([0.0; 3], [10.0, 0.0, 0.0], 0.5, 0.5);
    assert_eq!(capsule.toi_vs_sphere([0.3, 0.2, 0.0], 0.5), Some(0.0));
}

// ---------------------------------------------------------------------------
// F1689: capsule_capsule_toi no longer tunnels through contact windows
// ---------------------------------------------------------------------------

#[test]
fn capsule_capsule_toi_finds_contact_that_bisection_missed() {
    // d(0.5) = 10 and d(1) = 30 (both outside), yet the spheres touch at
    // x = -2, i.e. t = 0.2.
    let a = SweptSphere::new([-10.0, 0.0, 0.0], [30.0, 0.0, 0.0], 1.0);
    let b = SweptSphere::new([0.0; 3], [0.0; 3], 1.0);
    let t = capsule_capsule_toi(&a, &b).expect("sweep passes through the target");
    assert!((t - 0.2).abs() < 1e-9, "t = {t}");
}

#[test]
fn capsule_capsule_toi_two_moving_spheres_crossing() {
    // The centres start 40 apart and close at 40 units per unit time; with a
    // combined radius of 4 the spheres first touch at t = 36 / 40.
    let a = SweptSphere::new([-20.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0);
    let b = SweptSphere::new([20.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0);
    let t = capsule_capsule_toi(&a, &b).expect("spheres collide");
    assert!((t - 0.9).abs() < 1e-12, "t = {t}");
}

#[test]
fn capsule_capsule_toi_crossing_paths_at_different_times_do_not_collide() {
    // The two sweep paths cross at the origin, but A is at the origin at
    // t = 0.5 while B passed it at t ~ 0.09, so the spheres never touch.
    let a = SweptSphere::new([-10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5);
    let b = SweptSphere::new([0.0, -3.0, 0.0], [0.0, 30.0, 0.0], 0.5);
    assert_eq!(capsule_capsule_toi(&a, &b), None);
}

#[test]
fn capsule_capsule_toi_matches_swept_sphere_vs_swept_sphere() {
    let mut rng = Lcg(0x5eed_1234_abcd_0001);
    for _ in 0..500 {
        let a = SweptSphere::new(rng.point(6.0), rng.point(6.0), rng.range(0.2, 1.5));
        let b = SweptSphere::new(rng.point(6.0), rng.point(6.0), rng.range(0.2, 1.5));
        let direct = swept_sphere_vs_swept_sphere(
            a.center_start,
            a.center_end,
            a.radius,
            b.center_start,
            b.center_end,
            b.radius,
        );
        assert_eq!(capsule_capsule_toi(&a, &b), direct);
    }
}

#[test]
fn moving_sphere_pair_matches_dense_sampling() {
    let mut rng = Lcg(0x0bad_cafe_0000_0042);
    let samples = 4000usize;
    for case in 0..1500 {
        let a = SweptSphere::new(rng.point(8.0), rng.point(8.0), rng.range(0.2, 1.5));
        let b = SweptSphere::new(rng.point(8.0), rng.point(8.0), rng.range(0.2, 1.5));
        let combined = a.radius + b.radius;
        let dist_at = |t: f64| len3(sub3(a.center_at(t), b.center_at(t)));
        let first_inside = (0..=samples).find(|&k| dist_at(k as f64 / samples as f64) <= combined);
        let toi = capsule_capsule_toi(&a, &b);
        match (toi, first_inside) {
            (None, Some(k)) => panic!("case {case}: sample {k} is inside but no contact reported"),
            (Some(t), None) => {
                assert!(
                    (dist_at(t) - combined).abs() < 1e-7,
                    "case {case}: grazing contact at t = {t} has d = {}",
                    dist_at(t)
                );
            }
            (Some(t), Some(k)) => {
                let t_k = k as f64 / samples as f64;
                assert!(t <= t_k + 1e-9, "case {case}: t = {t} after sample {t_k}");
                if k > 0 {
                    assert!(
                        t >= (k - 1) as f64 / samples as f64 - 1e-9,
                        "case {case}: t = {t} before the last outside sample"
                    );
                    assert!((dist_at(t) - combined).abs() < 1e-7);
                } else {
                    assert_eq!(t, 0.0);
                }
            }
            (None, None) => {}
        }
    }
}

// ---------------------------------------------------------------------------
// F1690: ccd_sphere_vs_capsule is exact (no t = 0.5 / 1 sampling)
// ---------------------------------------------------------------------------

#[test]
fn ccd_sphere_vs_capsule_finds_short_contact_window() {
    // Contact window t in [0.225, 0.275]; d(0.5) = 10 and d(1) = 30.
    let hit = ccd_sphere_vs_capsule(
        [-10.0, 0.5, 0.0],
        [30.0, 0.5, 0.0],
        0.5,
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        0.5,
    )
    .expect("the sweep crosses the capsule");
    assert!((hit.toi - 0.225).abs() < 1e-12, "toi = {}", hit.toi);
    // The sphere centre is at (-1, 0.5, 0): the closest axis point is (0, 0.5, 0).
    assert!(hit.contact_point[0].abs() < 1e-12);
    assert!((hit.contact_point[1] - 0.5).abs() < 1e-12);
    assert!(hit.contact_point[2].abs() < 1e-12);
}

#[test]
fn ccd_sphere_vs_capsule_hits_the_end_cap() {
    // Sweeping along the axis the sphere meets the lower hemispherical cap at
    // y = -2 (distance 1 = sphere + capsule radius from the axis end -1).
    let hit = ccd_sphere_vs_capsule(
        [0.0, -10.0, 0.0],
        [0.0, 10.0, 0.0],
        0.5,
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
        0.5,
    )
    .expect("axial sweep hits the cap");
    assert!((hit.toi - 0.4).abs() < 1e-12, "toi = {}", hit.toi);
    assert!((hit.contact_point[1] + 1.0).abs() < 1e-12);
}

#[test]
fn ccd_sphere_vs_capsule_degenerate_axis_is_a_sphere() {
    // A zero-length capsule used to be missed (division by the segment
    // length); it is a sphere of radius 0.5 at the origin.
    let hit = ccd_sphere_vs_capsule(
        [-5.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
        0.5,
        [0.0; 3],
        [0.0; 3],
        0.5,
    )
    .expect("sweep through a point capsule");
    assert!((hit.toi - 0.4).abs() < 1e-12, "toi = {}", hit.toi);
    assert_eq!(hit.contact_point, [0.0; 3]);
}

#[test]
fn ccd_sphere_vs_capsule_misses_and_stationary_cases() {
    let cap_a = [0.0, -1.0, 0.0];
    let cap_b = [0.0, 1.0, 0.0];
    // Moving away.
    assert!(
        ccd_sphere_vs_capsule([-3.0, 0.0, 0.0], [-9.0, 0.0, 0.0], 0.5, cap_a, cap_b, 0.5).is_none()
    );
    // Contact only beyond the end of the sweep.
    assert!(
        ccd_sphere_vs_capsule([-10.0, 0.0, 0.0], [-5.0, 0.0, 0.0], 0.5, cap_a, cap_b, 0.5)
            .is_none()
    );
    // Passes beside the capsule.
    assert!(
        ccd_sphere_vs_capsule([-10.0, 0.0, 3.0], [10.0, 0.0, 3.0], 0.5, cap_a, cap_b, 0.5)
            .is_none()
    );
    // Stationary and separated / overlapping.
    assert!(ccd_sphere_vs_capsule([4.0; 3], [4.0; 3], 0.5, cap_a, cap_b, 0.5).is_none());
    let inside = ccd_sphere_vs_capsule([0.5, 0.0, 0.0], [0.5, 0.0, 0.0], 0.5, cap_a, cap_b, 0.5)
        .expect("overlapping stationary sphere");
    assert_eq!(inside.toi, 0.0);
}

#[test]
fn ccd_sphere_vs_capsule_matches_dense_sampling() {
    let mut rng = Lcg(0xfeed_beef_dead_0007);
    let samples = 3000usize;
    for case in 0..1500 {
        let cap_a = rng.point(3.0);
        let cap_b = rng.point(3.0);
        let cap_radius = rng.range(0.1, 1.5);
        let start = rng.point(9.0);
        let end = rng.point(9.0);
        let sphere_radius = rng.range(0.1, 1.0);
        let combined = sphere_radius + cap_radius;
        let dist_at = |t: f64| point_segment_distance(lerp3(start, end, t), cap_a, cap_b);
        let first_inside = (0..=samples).find(|&k| dist_at(k as f64 / samples as f64) <= combined);
        let result = ccd_sphere_vs_capsule(start, end, sphere_radius, cap_a, cap_b, cap_radius);
        match (&result, first_inside) {
            (None, Some(k)) => panic!("case {case}: sample {k} is inside but no contact reported"),
            (Some(hit), None) => {
                assert!(
                    (dist_at(hit.toi) - combined).abs() < 1e-6,
                    "case {case}: grazing contact at t = {} has d = {}",
                    hit.toi,
                    dist_at(hit.toi)
                );
            }
            (Some(hit), Some(k)) => {
                let t_k = k as f64 / samples as f64;
                assert!(
                    hit.toi <= t_k + 1e-9,
                    "case {case}: t = {} after {t_k}",
                    hit.toi
                );
                if k > 0 {
                    assert!(
                        hit.toi >= (k - 1) as f64 / samples as f64 - 1e-9,
                        "case {case}: t = {} before last outside sample",
                        hit.toi
                    );
                    assert!((dist_at(hit.toi) - combined).abs() < 1e-6);
                } else {
                    assert_eq!(hit.toi, 0.0);
                }
            }
            (None, None) => {}
        }
    }
}

#[test]
fn ccd_sphere_vs_capsule_parallel_sweeps() {
    let cap_a = [0.0, -1.0, 0.0];
    let cap_b = [0.0, 1.0, 0.0];
    // Offset 0.3 from the axis, travelling along it: the sphere centre comes
    // within 1.0 of the lower axis end at y = -1 - sqrt(1 - 0.3^2).
    let hit = ccd_sphere_vs_capsule([0.3, -10.0, 0.0], [0.3, 10.0, 0.0], 0.5, cap_a, cap_b, 0.5)
        .expect("axial sweep inside the capsule radius");
    let expected = (-1.0 - (1.0f64 - 0.09).sqrt() + 10.0) / 20.0;
    assert!((hit.toi - expected).abs() < 1e-12, "toi = {}", hit.toi);
    // Offset beyond the combined radius: parallel and never touching.
    assert!(
        ccd_sphere_vs_capsule([1.2, -10.0, 0.0], [1.2, 10.0, 0.0], 0.5, cap_a, cap_b, 0.5)
            .is_none()
    );
}

#[test]
fn linear_cast_sphere_vs_sphere_overlap_reports_zero() {
    let cast = linear_cast_sphere_vs_sphere([0.0; 3], 1.0, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0)
        .expect("overlapping spheres");
    assert_eq!(cast.toi, 0.0);
    assert_eq!(cast.normal, [-1.0, 0.0, 0.0]);
}

// ---------------------------------------------------------------------------
// SweptCapsule::toi_vs_sphere is exact, not a bounding-sphere test (F1692)
// ---------------------------------------------------------------------------

#[test]
fn swept_capsule_toi_vs_sphere_matches_dense_sampling() {
    let mut rng = Lcg(0xc0ff_ee00_1234_5678);
    let samples = 3000usize;
    for case in 0..1500 {
        let start = rng.point(6.0);
        let end = rng.point(6.0);
        let radius = rng.range(0.1, 1.2);
        let half_height = rng.range(0.0, 1.5);
        let capsule = SweptCapsule::new(start, end, radius, half_height);
        let centre = rng.point(6.0);
        let sphere_radius = rng.range(0.1, 1.0);
        let combined = radius + sphere_radius;
        let axis = [0.0, half_height, 0.0];
        let dist_at = |t: f64| {
            let position = lerp3(start, end, t);
            point_segment_distance(centre, sub3(position, axis), add3(position, axis))
        };
        let first_inside = (0..=samples).find(|&k| dist_at(k as f64 / samples as f64) <= combined);
        match (capsule.toi_vs_sphere(centre, sphere_radius), first_inside) {
            (None, Some(k)) => panic!("case {case}: sample {k} touches but no TOI reported"),
            (Some(t), None) => assert!(
                (dist_at(t) - combined).abs() < 1e-6,
                "case {case}: grazing TOI {t} has d = {}",
                dist_at(t)
            ),
            (Some(t), Some(k)) => {
                let t_k = k as f64 / samples as f64;
                assert!(t <= t_k + 1e-9, "case {case}: t = {t} after sample {t_k}");
                if k > 0 {
                    assert!(
                        t >= (k - 1) as f64 / samples as f64 - 1e-9,
                        "case {case}: t = {t}"
                    );
                    assert!((dist_at(t) - combined).abs() < 1e-6, "case {case}");
                } else {
                    assert_eq!(t, 0.0);
                }
            }
            (None, None) => {}
        }
    }
}

#[test]
fn swept_capsule_toi_vs_sphere_ignores_bounding_sphere_slack() {
    // Capsule radius 0.5 and half height 1.0 (bounding-sphere radius 1.5).
    // The sphere sits 1.2 away in z: farther than radius + sphere radius (1.0)
    // from the axis but inside the old bounding-sphere distance (2.0).
    let capsule = SweptCapsule::new([-10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5, 1.0);
    assert_eq!(capsule.toi_vs_sphere([0.0, 0.0, 1.2], 0.5), None);
}

#[test]
fn swept_capsule_toi_vs_sphere_hits_the_cylinder_side() {
    let capsule = SweptCapsule::new([-10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5, 1.0);
    // Sphere on the line y = 0.3: axis distance is |dx| until the capsule
    // arrives; contact at dx = -1 -> capsule centre x = -1, t = 9/20.
    let t = capsule
        .toi_vs_sphere([0.0, 0.3, 0.0], 0.5)
        .expect("capsule sweeps through the sphere");
    assert!((t - 0.45).abs() < 1e-12, "t = {t}");
}

#[test]
fn swept_capsule_toi_vs_sphere_hits_the_top_cap() {
    // Sphere centre 0.4 above the axis end (y = 1.4): horizontal offset at
    // contact is sqrt(1 - 0.4^2).
    let capsule = SweptCapsule::new([-10.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.5, 1.0);
    let t = capsule
        .toi_vs_sphere([0.0, 1.4, 0.0], 0.5)
        .expect("capsule cap touches the sphere");
    let expected = (10.0 - (1.0f64 - 0.16).sqrt()) / 20.0;
    assert!((t - expected).abs() < 1e-12, "t = {t}, expected {expected}");
}

#[test]
fn swept_capsule_with_zero_half_height_is_a_swept_sphere() {
    let capsule = SweptCapsule::new([-5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.5, 0.0);
    let t = capsule
        .toi_vs_sphere([0.0; 3], 0.5)
        .expect("a sphere-shaped capsule hits head-on");
    assert!((t - 0.4).abs() < 1e-12, "t = {t}");
}

// ---------------------------------------------------------------------------
// SweptSphere::ray_intersect with a zero-length sweep (a stationary sphere)
// ---------------------------------------------------------------------------

#[test]
fn ray_intersect_handles_a_stationary_sphere() {
    let ball = SweptSphere::new([0.0; 3], [0.0; 3], 1.0);
    let t = ball
        .ray_intersect([0.0, 5.0, 0.0], [0.0, -1.0, 0.0])
        .expect("ray towards the sphere hits it");
    assert!((t - 4.0).abs() < 1e-12, "t = {t}");
    assert!(
        ball.ray_intersect([3.0, 5.0, 0.0], [0.0, -1.0, 0.0])
            .is_none()
    );
    assert!(
        ball.ray_intersect([0.0, 5.0, 0.0], [0.0, 1.0, 0.0])
            .is_none()
    );
    // Origin inside the sphere: the first surface crossing is the exit point.
    let exit = ball
        .ray_intersect([0.0; 3], [1.0, 0.0, 0.0])
        .expect("ray leaves the sphere");
    assert!((exit - 1.0).abs() < 1e-12, "exit = {exit}");
    assert!(ball.ray_intersect([0.0, 5.0, 0.0], [0.0; 3]).is_none());
}
