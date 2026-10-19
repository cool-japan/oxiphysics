// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the exact swept sphere queries (triangle, capsule) and
//! the conservative advancement of rotating bodies in [`super`].

use super::point_cast::{closest_point_on_segment, cylinder_entry, sphere_entry};
use super::sphere_triangle::closest_point_on_triangle;
use super::*;

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn vec_close(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    len(sub(a, b)) <= tol
}

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    add(scale(a, 1.0 - t), scale(b, t))
}

/// Deterministic xorshift generator for the randomized tests.
struct TestRng(u64);

impl TestRng {
    fn new(seed: u64) -> Self {
        TestRng(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

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

    /// Uniformly distributed unit vector (rejection sampling in the unit ball).
    fn direction(&mut self) -> [f64; 3] {
        loop {
            let p = self.point(1.0);
            let l = len_sq(p);
            if l > 0.01 && l <= 1.0 {
                return normalize(p);
            }
        }
    }

    fn pick(&mut self, values: &[f64]) -> f64 {
        values[(self.next_u64() % values.len() as u64) as usize]
    }

    fn lattice_point(&mut self, values: &[f64]) -> [f64; 3] {
        [self.pick(values), self.pick(values), self.pick(values)]
    }
}

// ---------------------------------------------------------------------------
// Reference distance functions and time-of-impact verification
// ---------------------------------------------------------------------------

/// Distance from `p` to the segment `a`-`b` by clamped projection.
fn reference_segment_distance(a: [f64; 3], b: [f64; 3], p: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let denom = dot(ab, ab);
    let u = if denom > 0.0 {
        (dot(sub(p, a), ab) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    len(sub(p, add(a, scale(ab, u))))
}

/// Distance from `p` to a triangle: the perpendicular distance when the
/// projection falls inside, the distance to the nearest edge otherwise.
fn reference_triangle_distance(tri: &[[f64; 3]; 3], p: [f64; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..3 {
        best = best.min(reference_segment_distance(tri[i], tri[(i + 1) % 3], p));
    }
    let n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
    let n_len = len(n);
    if n_len > 0.0 {
        let n = scale(n, 1.0 / n_len);
        let s = dot(n, sub(p, tri[0]));
        let q = sub(p, scale(n, s));
        let inside =
            (0..3).all(|i| dot(cross(sub(tri[(i + 1) % 3], tri[i]), sub(q, tri[i])), n) >= 0.0);
        if inside {
            best = best.min(s.abs());
        }
    }
    best
}

/// Check a computed time of impact against a distance function `dist(t)` (the
/// distance between the sphere centre at time `t` and the target body).
///
/// * a reported contact must be real: `dist(toi) == radius`;
/// * the sphere must be outside the body at every sampled time before it;
/// * no contact that penetrates by more than a small margin may be missed, and
///   the reported time may not be later than the first such sampled time.
fn verify_time_of_impact(toi: Option<f64>, radius: f64, dist: &dyn Fn(f64) -> f64, label: &str) {
    const SAMPLES: usize = 400;
    const MARGIN: f64 = 1e-7;
    // A path through a feature with zero radius is a double root of the
    // quadratic (the squared distance is a perfect square), which resolves the
    // time only to about the square root of the rounding error.
    const VALIDITY_TOLERANCE: f64 = 1e-6;
    let sample = |k: usize| k as f64 / SAMPLES as f64;
    let penetrating = (0..=SAMPLES)
        .map(sample)
        .find(|&t| dist(t) <= radius - MARGIN);
    let Some(t) = toi else {
        assert!(
            penetrating.is_none(),
            "{label}: missed a contact that penetrates at t = {penetrating:?}"
        );
        return;
    };
    assert!((0.0..=1.0).contains(&t), "{label}: toi {t} outside [0, 1]");
    if t > 0.0 {
        assert!(
            (dist(t) - radius).abs() <= VALIDITY_TOLERANCE,
            "{label}: no contact at the reported toi {t}: distance {} radius {radius}",
            dist(t)
        );
    } else {
        assert!(
            dist(0.0) <= radius + VALIDITY_TOLERANCE,
            "{label}: toi 0 but the sphere starts apart: distance {} radius {radius}",
            dist(0.0)
        );
    }
    if let Some(t_ref) = penetrating {
        assert!(
            t <= t_ref + 1e-9,
            "{label}: toi {t} later than the penetration at {t_ref}"
        );
    }
    for k in 0..=SAMPLES {
        let s = sample(k);
        if s < t - 1e-9 {
            assert!(
                dist(s) >= radius - 1e-8,
                "{label}: the sphere penetrates at t = {s}, before the reported toi {t}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Point-cast primitives
// ---------------------------------------------------------------------------

#[test]
fn sphere_entry_returns_the_first_crossing() {
    // Unit sphere at the origin, point from x = -3 to x = 1 (|dir| = 4).
    let t = sphere_entry([-3.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0; 3], 1.0).expect("hit");
    assert!(approx(t, 0.5, 1e-12), "t={t}");
    // Starting inside reports 0, moving away or short of the sphere reports none.
    assert_eq!(
        sphere_entry([0.5, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3], 1.0),
        Some(0.0)
    );
    assert_eq!(
        sphere_entry([-3.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0; 3], 1.0),
        None
    );
    assert_eq!(
        sphere_entry([-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3], 1.0),
        None
    );
    // No motion at all must not divide by zero.
    assert_eq!(
        sphere_entry([-3.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0),
        None
    );
}

#[test]
fn cylinder_entry_requires_the_axis_position_inside_the_segment() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    // Straight down onto the middle of the segment: touches at y = 1.
    let t = cylinder_entry([0.5, 3.0, 0.0], [0.0, -6.0, 0.0], a, b, 1.0).expect("hit");
    assert!(approx(t, 1.0 / 3.0, 1e-12), "t={t}");
    // The same descent beyond the end of the segment is not a lateral hit.
    assert_eq!(
        cylinder_entry([1.5, 3.0, 0.0], [0.0, -6.0, 0.0], a, b, 1.0),
        None
    );
    // Motion parallel to the axis never reaches the lateral surface.
    assert_eq!(
        cylinder_entry([-2.0, 3.0, 0.0], [4.0, 0.0, 0.0], a, b, 1.0),
        None
    );
    // A zero-length segment has no lateral surface.
    assert_eq!(
        cylinder_entry([0.5, 3.0, 0.0], [0.0, -6.0, 0.0], a, a, 1.0),
        None
    );
}

#[test]
fn closest_point_on_segment_clamps_to_the_end_points() {
    let a = [0.0, 0.0, 0.0];
    let b = [2.0, 0.0, 0.0];
    let (p, u) = closest_point_on_segment(a, b, [1.0, 5.0, 0.0]);
    assert!(vec_close(p, [1.0, 0.0, 0.0], 1e-15) && approx(u, 0.5, 1e-15));
    let (p, u) = closest_point_on_segment(a, b, [-4.0, 1.0, 0.0]);
    assert!(vec_close(p, a, 1e-15) && u == 0.0);
    let (p, u) = closest_point_on_segment(a, b, [9.0, 1.0, 0.0]);
    assert!(vec_close(p, b, 1e-15) && u == 1.0);
    let (p, u) = closest_point_on_segment(a, a, [9.0, 1.0, 0.0]);
    assert!(vec_close(p, a, 1e-15) && u == 0.0);
}

#[test]
fn closest_point_on_triangle_handles_every_region_and_degenerate_triangles() {
    let (v0, v1, v2) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
    let cases: [([f64; 3], [f64; 3]); 7] = [
        ([0.25, 2.0, 0.25], [0.25, 0.0, 0.25]), // face
        ([-1.0, 1.0, -1.0], v0),                // vertex 0
        ([3.0, 1.0, -1.0], v1),                 // vertex 1
        ([-1.0, 1.0, 3.0], v2),                 // vertex 2
        ([0.5, 1.0, -1.0], [0.5, 0.0, 0.0]),    // edge 0-1
        ([-1.0, 1.0, 0.5], [0.0, 0.0, 0.5]),    // edge 0-2
        ([1.0, 1.0, 1.0], [0.5, 0.0, 0.5]),     // edge 1-2
    ];
    for (p, expected) in cases {
        let q = closest_point_on_triangle(v0, v1, v2, p);
        assert!(
            vec_close(q, expected, 1e-12),
            "p={p:?}: {q:?} != {expected:?}"
        );
    }
    // Collinear and coincident vertices give finite points on the segment.
    let q = closest_point_on_triangle([0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 0.5, 0.0]);
    assert!(vec_close(q, [1.5, 0.0, 0.0], 1e-12), "{q:?}");
    let q = closest_point_on_triangle([1.0; 3], [1.0; 3], [1.0; 3], [4.0, 5.0, 6.0]);
    assert!(vec_close(q, [1.0; 3], 1e-12), "{q:?}");
}

// ---------------------------------------------------------------------------
// conservative_advancement_rotating
// ---------------------------------------------------------------------------

/// Sphere A at the origin moving along +X at `speed`, static sphere B so that
/// the bounding-sphere gap is `gap` (both radii 0.5).
fn approaching_pair(
    speed: f64,
    gap: f64,
    angular_speed_a: f64,
) -> (RotatingBodyState, RotatingBodyState) {
    (
        RotatingBodyState::new([0.0; 3], [speed, 0.0, 0.0], angular_speed_a, 0.5),
        RotatingBodyState::new([1.0 + gap, 0.0, 0.0], [0.0; 3], 0.0, 0.5),
    )
}

#[test]
fn advancement_step_is_gap_over_rate_not_gap_over_distance() {
    // Gap 0.1, closing speed 10, interval [0, 0.016]: the contact is at 0.01.
    // The step used to be a dimensionless ratio that was clamped to the whole
    // interval, so the result was 0.016 with the spheres 6 cm inside each other.
    let (a, b) = approaching_pair(10.0, 0.1, 0.0);
    let toi = conservative_advancement_rotating(&a, &b, 0.1, 0.0, 0.016, 1e-6, 100)
        .expect("the spheres meet within the interval");
    assert!(approx(toi, 0.01, 1e-4), "toi={toi}");
    assert!(toi <= 0.01 + 1e-12, "advanced past the contact: toi={toi}");
}

#[test]
fn advancement_returns_absolute_time_not_a_fraction() {
    // The same approach with the interval starting at t = 1: the positions are
    // those at t_start and the contact is 0.01 after it.
    let (a, b) = approaching_pair(10.0, 0.1, 0.0);
    let toi = conservative_advancement_rotating(&a, &b, 0.1, 1.0, 1.016, 1e-6, 100)
        .expect("the spheres meet within the interval");
    assert!(approx(toi, 1.01, 1e-4), "toi={toi}");
    assert!(toi <= 1.01 + 1e-12);
    // And with a long interval the contact time is not a fraction of it.
    let toi = conservative_advancement_rotating(&a, &b, 0.1, 0.0, 100.0, 1e-6, 100).expect("hit");
    assert!(approx(toi, 0.01, 1e-4), "toi={toi}");
}

#[test]
fn advancement_with_rotation_never_passes_the_contact() {
    // Rotation enlarges the speed bound, which only makes the steps smaller.
    let (a, b) = approaching_pair(10.0, 0.1, 20.0);
    let toi = conservative_advancement_rotating(&a, &b, 0.1, 0.0, 0.016, 1e-6, 200).expect("hit");
    assert!(toi <= 0.01 + 1e-12, "toi={toi}");
    assert!(toi > 0.0099, "toi={toi}");
}

#[test]
fn advancement_reports_none_when_the_interval_ends_before_contact() {
    // Gap 1.0, closing speed 10: 0.05 s covers only 0.5.
    let (a, b) = approaching_pair(10.0, 1.0, 0.0);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 1.0, 0.0, 0.05, 1e-6, 100),
        None
    );
}

#[test]
fn advancement_reports_contact_exactly_at_the_end_of_the_interval() {
    let (a, b) = approaching_pair(10.0, 0.16, 0.0);
    let toi = conservative_advancement_rotating(&a, &b, 0.16, 0.0, 0.016, 1e-6, 100)
        .expect("touching at the end of the interval");
    assert!(approx(toi, 0.016, 1e-9), "toi={toi}");
}

#[test]
fn advancement_does_not_trust_a_separation_larger_than_the_bounding_gap() {
    // The bounding spheres are 0.1 apart. A caller-supplied separation of 5 can
    // only describe shapes that do not fill their spheres; stepping by it would
    // carry the bounding spheres 0.4 into each other, so the first step uses the
    // smaller of the two gaps and the contact time is still found exactly.
    let (a, b) = approaching_pair(10.0, 0.1, 0.0);
    let toi = conservative_advancement_rotating(&a, &b, 5.0, 0.0, 1.0, 1e-9, 100).expect("hit");
    assert!(approx(toi, 0.01, 1e-8), "toi={toi}");
    assert!(toi <= 0.01 + 1e-12, "toi={toi}");
    // A smaller (more cautious) separation only slows the first step down.
    let toi = conservative_advancement_rotating(&a, &b, 0.05, 0.0, 1.0, 1e-9, 100).expect("hit");
    assert!(approx(toi, 0.01, 1e-8), "toi={toi}");
}

#[test]
fn advancement_with_no_iterations_reports_the_start_time_as_a_lower_bound() {
    let (a, b) = approaching_pair(10.0, 0.1, 0.0);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 0.1, 0.5, 1.0, 1e-9, 0),
        Some(0.5)
    );
}

#[test]
fn advancement_handles_bodies_without_relative_motion() {
    let a = RotatingBodyState::new([0.0; 3], [3.0, 1.0, 0.0], 0.0, 0.5);
    let b = RotatingBodyState::new([5.0, 0.0, 0.0], [3.0, 1.0, 0.0], 0.0, 0.5);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 4.0, 0.0, 10.0, 1e-6, 100),
        None
    );
}

#[test]
fn advancement_uses_the_relative_speed_when_one_body_chases_the_other() {
    // A chases B at 10 versus 2: the closing speed is 8, not 12.
    let a = RotatingBodyState::new([0.0; 3], [10.0, 0.0, 0.0], 0.0, 0.5);
    let b = RotatingBodyState::new([1.8, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0, 0.5);
    let toi = conservative_advancement_rotating(&a, &b, 0.8, 0.0, 1.0, 1e-9, 100).expect("hit");
    assert!(approx(toi, 0.8 / 8.0, 1e-8), "toi={toi}");
    assert!(toi <= 0.1 + 1e-12);
}

#[test]
fn advancement_reports_overlap_at_the_start_time() {
    let a = RotatingBodyState::new([0.0; 3], [0.0; 3], 0.0, 1.0);
    let b = RotatingBodyState::new([0.5, 0.0, 0.0], [0.0; 3], 0.0, 1.0);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 0.0, 2.0, 3.0, 1e-4, 100),
        Some(2.0)
    );
    // A caller that knows the shapes already penetrate gets the start time too.
    let (a, b) = approaching_pair(1.0, 5.0, 0.0);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, -0.1, 0.5, 3.0, 1e-4, 100),
        Some(0.5)
    );
}

#[test]
fn advancement_rejects_empty_or_invalid_intervals() {
    let (a, b) = approaching_pair(10.0, 0.1, 0.0);
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 0.1, 1.0, 1.0, 1e-6, 100),
        None
    );
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 0.1, 2.0, 1.0, 1e-6, 100),
        None
    );
    assert_eq!(
        conservative_advancement_rotating(&a, &b, 0.1, 0.0, f64::NAN, 1e-6, 100),
        None
    );
}

#[test]
fn advancement_exhausted_budget_returns_a_lower_bound_of_the_contact_time() {
    // An oblique approach needs more than two steps to reach the tolerance.
    let a = RotatingBodyState::new([0.0; 3], [10.0, 0.0, 0.0], 0.0, 0.5);
    let b = RotatingBodyState::new([1.2, 0.5, 0.0], [0.0; 3], 0.0, 0.5);
    let exact = swept_sphere_sphere(
        [0.0; 3],
        [10.0, 0.0, 0.0],
        0.5,
        [1.2, 0.5, 0.0],
        [1.2, 0.5, 0.0],
        0.5,
    )
    .expect("the spheres meet")
    .toi;
    let gap = len(sub(b.position, a.position)) - 1.0;
    let toi = conservative_advancement_rotating(&a, &b, gap, 0.0, 1.0, 1e-9, 2)
        .expect("a lower bound is reported, not a miss");
    assert!(
        toi < exact - 1e-4,
        "not converged yet: toi={toi}, exact={exact}"
    );
    assert!(toi > 0.0, "toi={toi}");
    // With a sufficient budget the contact time is reached from below.
    let toi = conservative_advancement_rotating(&a, &b, gap, 0.0, 1.0, 1e-9, 200).expect("hit");
    assert!(
        toi <= exact + 1e-12 && exact - toi < 1e-6,
        "toi={toi}, exact={exact}"
    );
}

#[test]
fn advancement_matches_the_analytic_sphere_sphere_time_of_impact() {
    let mut rng = TestRng::new(0x5EED_0A11);
    let (mut hits, mut misses) = (0, 0);
    for case in 0..1500 {
        let radius_a = rng.range(0.2, 1.0);
        let radius_b = rng.range(0.2, 1.0);
        let vel_a = rng.point(4.0);
        let vel_b = rng.point(2.0);
        let start_a = rng.point(2.0);
        // B starts near the path A follows, so that about half of the pairs meet.
        let near_path = add(start_a, scale(vel_a, rng.range(0.0, 1.0)));
        let start_b = add(near_path, scale(rng.direction(), rng.range(0.0, 2.0)));
        let a = RotatingBodyState::new(start_a, vel_a, 0.0, radius_a);
        let b = RotatingBodyState::new(start_b, vel_b, 0.0, radius_b);
        let exact = swept_sphere_sphere(
            a.position,
            add(a.position, vel_a),
            radius_a,
            b.position,
            add(b.position, vel_b),
            radius_b,
        );
        let gap0 = len(sub(b.position, a.position)) - radius_a - radius_b;
        let advanced =
            conservative_advancement_rotating(&a, &b, gap0.max(0.0), 0.0, 1.0, 1e-9, 100_000);
        match exact {
            Some(hit) => {
                hits += 1;
                let t = advanced.unwrap_or(f64::NAN);
                assert!(
                    t <= hit.toi + 1e-9,
                    "case {case}: advanced {t} past exact {}",
                    hit.toi
                );
                // Transversal contacts converge to the contact time itself.
                let rel = sub(vel_a, vel_b);
                let dp = sub(sub(b.position, a.position), scale(rel, hit.toi));
                let closing = dot(rel, dp) / len(dp).max(1e-300);
                if hit.toi > 0.0 && closing > 0.3 * len(rel) {
                    assert!(hit.toi - t < 1e-6, "case {case}: t={t}, exact {}", hit.toi);
                }
            }
            None => {
                misses += 1;
                let rel = sub(vel_a, vel_b);
                let dp = sub(b.position, a.position);
                let rel_sq = dot(rel, rel);
                let closest_t = if rel_sq > 0.0 {
                    (dot(dp, rel) / rel_sq).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let min_gap = len(sub(dp, scale(rel, closest_t))) - radius_a - radius_b;
                if min_gap > 1e-6 {
                    assert_eq!(advanced, None, "case {case}: reported a contact for a miss");
                }
            }
        }
    }
    assert!(hits > 200 && misses > 200, "hits {hits}, misses {misses}");
}

#[test]
fn advancement_with_rotation_never_passes_the_analytic_contact() {
    // Rotation only enlarges the speed bound, so the steps shrink but still
    // never step past the contact of the (rotation-invariant) bounding spheres.
    let mut rng = TestRng::new(0x0707_A7E5);
    let mut hits = 0;
    for case in 0..800 {
        let radius_a = rng.range(0.2, 1.0);
        let radius_b = rng.range(0.2, 1.0);
        let vel_a = rng.point(4.0);
        let vel_b = rng.point(2.0);
        let start_a = rng.point(2.0);
        let near_path = add(start_a, scale(vel_a, rng.range(0.0, 1.0)));
        let start_b = add(near_path, scale(rng.direction(), rng.range(0.0, 2.0)));
        let a = RotatingBodyState::new(start_a, vel_a, rng.range(0.0, 6.0), radius_a);
        let b = RotatingBodyState::new(start_b, vel_b, rng.range(0.0, 6.0), radius_b);
        let exact = swept_sphere_sphere(
            a.position,
            add(a.position, vel_a),
            radius_a,
            b.position,
            add(b.position, vel_b),
            radius_b,
        );
        let gap0 = len(sub(b.position, a.position)) - radius_a - radius_b;
        let advanced =
            conservative_advancement_rotating(&a, &b, gap0.max(0.0), 0.0, 1.0, 1e-6, 10_000);
        if let Some(hit) = exact {
            hits += 1;
            let t = advanced.unwrap_or(f64::NAN);
            assert!(
                t <= hit.toi + 1e-9,
                "case {case}: advanced {t} past exact {}",
                hit.toi
            );
        }
    }
    assert!(hits > 100, "hits {hits}");
}

#[test]
fn advancement_with_distance_uses_the_shape_separation() {
    // Both bodies have bounding spheres of radius 2 but the shapes inside are
    // points of radius 0.05. A approaches B at 10 from 10 apart: the bounding
    // spheres touch at 0.6, the shapes at 0.99.
    let a = RotatingBodyState::new([0.0; 3], [10.0, 0.0, 0.0], 0.0, 2.0);
    let b = RotatingBodyState::new([10.0, 0.0, 0.0], [0.0; 3], 0.0, 2.0);
    let shape_gap = |t: f64| 10.0 - 10.0 * t - 0.1;
    let shapes =
        conservative_advancement_rotating_with_distance(&a, &b, 0.0, 1.0, 1e-9, 100, shape_gap)
            .expect("the shapes meet");
    assert!(approx(shapes, 0.99, 1e-6), "toi={shapes}");
    let spheres = conservative_advancement_rotating(&a, &b, 9.9, 0.0, 1.0, 1e-9, 100)
        .expect("the bounding spheres meet");
    assert!(approx(spheres, 0.6, 1e-6), "toi={spheres}");
    assert!(spheres < shapes);
}

#[test]
fn advancement_with_distance_reports_start_overlap_and_misses() {
    let a = RotatingBodyState::new([0.0; 3], [1.0, 0.0, 0.0], 0.0, 1.0);
    let b = RotatingBodyState::new([10.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1.0);
    // Penetrating at the start.
    assert_eq!(
        conservative_advancement_rotating_with_distance(&a, &b, 3.0, 4.0, 1e-6, 50, |_| -0.2),
        Some(3.0)
    );
    // Equal velocities never close a gap of 8.
    assert_eq!(
        conservative_advancement_rotating_with_distance(&a, &b, 0.0, 4.0, 1e-6, 50, |_| 8.0),
        None
    );
    // Invalid interval.
    assert_eq!(
        conservative_advancement_rotating_with_distance(&a, &b, 1.0, 1.0, 1e-6, 50, |_| 8.0),
        None
    );
}

#[test]
fn rotating_body_speed_bound_uses_the_angular_speed_magnitude() {
    let body = RotatingBodyState::new([0.0; 3], [3.0, 4.0, 0.0], -2.0, 0.5);
    assert!(approx(body.max_point_speed(), 6.0, 1e-12));
    assert!(approx(body.motion_bound(0.5), 3.0, 1e-12));
}

// ---------------------------------------------------------------------------
// swept_sphere_triangle / swept_sphere_triangle_mesh
// ---------------------------------------------------------------------------

/// Triangle in the XZ plane with its normal along -Y.
const FLOOR_TRI: [[f64; 3]; 3] = [[-1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 0.0, 1.0]];

/// Right triangle `(0,0,0) (1,0,0) (0,0,1)` in the XZ plane.
const CORNER_TRI: [[f64; 3]; 3] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

#[test]
fn triangle_face_hit_from_above_and_below() {
    let above = swept_sphere_triangle([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &FLOOR_TRI)
        .expect("a face hit");
    assert!(approx(above.toi, 0.625, 1e-12), "toi={}", above.toi);
    assert!(
        vec_close(above.normal, [0.0, 1.0, 0.0], 1e-12),
        "{:?}",
        above.normal
    );
    assert!(vec_close(above.point_a, [0.0, 0.0, 0.0], 1e-12));
    assert!(vec_close(above.point_b, [0.0, 0.0, 0.0], 1e-12));

    let below = swept_sphere_triangle([0.0, -3.0, 0.0], [0.0, 1.0, 0.0], 0.5, &FLOOR_TRI)
        .expect("a back-face hit");
    assert!(approx(below.toi, 0.625, 1e-12), "toi={}", below.toi);
    assert!(
        vec_close(below.normal, [0.0, -1.0, 0.0], 1e-12),
        "{:?}",
        below.normal
    );
}

#[test]
fn triangle_parallel_slide_into_an_edge_is_detected() {
    // The sphere moves parallel to the triangle plane, 0.3 above it, and runs
    // into the edge along Z at x = 0: the old plane-crossing test skipped it.
    let hit = swept_sphere_triangle([-2.0, 0.3, 0.25], [2.0, 0.3, 0.25], 0.5, &CORNER_TRI)
        .expect("an edge hit");
    assert!(approx(hit.toi, 0.4, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [-0.8, 0.6, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    assert!(
        vec_close(hit.point_b, [0.0, 0.0, 0.25], 1e-12),
        "{:?}",
        hit.point_b
    );
    assert!(
        vec_close(hit.point_a, [0.0, 0.0, 0.25], 1e-12),
        "{:?}",
        hit.point_a
    );
}

#[test]
fn triangle_vertex_hit_uses_the_vertex_sphere() {
    // Approach the corner at the origin along the diagonal, outside the triangle.
    let hit = swept_sphere_triangle([-3.0, 0.0, -3.0], [0.0, 0.0, 0.0], 0.5, &CORNER_TRI)
        .expect("a vertex hit");
    let expected = 1.0 - 0.5 / (3.0 * 2.0_f64.sqrt());
    assert!(approx(hit.toi, expected, 1e-12), "toi={}", hit.toi);
    let s = std::f64::consts::FRAC_1_SQRT_2;
    assert!(
        vec_close(hit.normal, [-s, 0.0, -s], 1e-12),
        "{:?}",
        hit.normal
    );
    assert!(vec_close(hit.point_b, [0.0; 3], 1e-12));
}

#[test]
fn triangle_edge_hit_ignores_the_infinite_cylinder_beyond_the_edge() {
    // Along the diagonal the sphere crosses the cylinders of both edges through
    // the corner before it touches the corner itself, but only outside the edge
    // segments (negative edge parameter): the contact is with the vertex.
    let hit =
        swept_sphere_triangle([-3.0, 0.0, -3.0], [0.0, 0.0, 0.0], 0.5, &CORNER_TRI).expect("hit");
    assert!(hit.toi > 0.88, "toi={}", hit.toi);
}

#[test]
fn triangle_initial_contact_reports_time_zero() {
    // Resting on the face (distance exactly the radius), moving away or into it.
    for end in [[0.25, 5.0, 0.25], [0.25, -1.0, 0.25]] {
        let hit =
            swept_sphere_triangle([0.25, 0.5, 0.25], end, 0.5, &CORNER_TRI).expect("touching");
        assert_eq!(hit.toi, 0.0);
        assert!(
            vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12),
            "{:?}",
            hit.normal
        );
    }
    // Deeply overlapping: the centre is inside the slab but over the triangle.
    let hit = swept_sphere_triangle([0.25, 0.1, 0.25], [0.25, 0.1, 0.25], 0.5, &CORNER_TRI)
        .expect("overlapping");
    assert_eq!(hit.toi, 0.0);
    assert!(vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12));
    assert!(vec_close(hit.point_b, [0.25, 0.0, 0.25], 1e-12));
    // Centre exactly on the triangle: the normal faces against the motion.
    let hit = swept_sphere_triangle([0.25, 0.0, 0.25], [0.25, 2.0, 0.25], 0.5, &CORNER_TRI)
        .expect("overlapping");
    assert_eq!(hit.toi, 0.0);
    assert!(
        vec_close(hit.normal, [0.0, -1.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
}

#[test]
fn triangle_misses_report_none() {
    // Sliding parallel above the face, farther than the radius.
    assert!(swept_sphere_triangle([-2.0, 1.0, 0.25], [2.0, 1.0, 0.25], 0.5, &CORNER_TRI).is_none());
    // Moving away.
    assert!(
        swept_sphere_triangle([0.25, 1.0, 0.25], [0.25, 4.0, 0.25], 0.5, &CORNER_TRI).is_none()
    );
    // Stopping short of the plane offset.
    assert!(
        swept_sphere_triangle([0.25, 3.0, 0.25], [0.25, 0.9, 0.25], 0.5, &CORNER_TRI).is_none()
    );
    // Passing beside the triangle.
    assert!(swept_sphere_triangle([5.0, 3.0, 0.0], [5.0, -3.0, 0.0], 0.5, &CORNER_TRI).is_none());
    // A static sphere that is apart.
    assert!(swept_sphere_triangle([0.0, 3.0, 0.0], [0.0, 3.0, 0.0], 0.5, &CORNER_TRI).is_none());
}

#[test]
fn triangle_point_sweep_crosses_the_face() {
    let hit = swept_sphere_triangle([0.25, 2.0, 0.25], [0.25, -2.0, 0.25], 0.0, &CORNER_TRI)
        .expect("a ray through the face");
    assert!(approx(hit.toi, 0.5, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    assert!(vec_close(hit.point_b, [0.25, 0.0, 0.25], 1e-12));
}

#[test]
fn triangle_without_area_still_hits_through_its_edges() {
    let collinear = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let hit = swept_sphere_triangle([1.0, 3.0, 0.0], [1.0, -3.0, 0.0], 0.5, &collinear)
        .expect("the sliver is hit");
    assert!(approx(hit.toi, 2.5 / 6.0, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    // Coincident vertices form a point.
    let point = [[1.0, 1.0, 1.0]; 3];
    let hit = swept_sphere_triangle([-2.0, 1.0, 1.0], [4.0, 1.0, 1.0], 0.5, &point).expect("hit");
    assert!(approx(hit.toi, 2.5 / 6.0, 1e-12), "toi={}", hit.toi);
    assert!(vec_close(hit.normal, [-1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn triangle_rejects_invalid_input() {
    assert!(swept_sphere_triangle([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], -0.5, &FLOOR_TRI).is_none());
    assert!(
        swept_sphere_triangle([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], f64::NAN, &FLOOR_TRI).is_none()
    );
    assert!(
        swept_sphere_triangle([f64::NAN, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &FLOOR_TRI).is_none()
    );
    let mut broken = FLOOR_TRI;
    broken[1][0] = f64::INFINITY;
    assert!(swept_sphere_triangle([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &broken).is_none());
}

#[test]
fn triangle_mesh_reports_the_earliest_triangle() {
    // Sphere from y = 5 to y = -5, radius 0.3: the triangle at y = 1 is reached
    // when the centre is at 1.3 (t = 0.37), the one at y = -5 much later.
    let triangles = [
        [[-1.0, -5.0, -1.0], [1.0, -5.0, -1.0], [0.0, -5.0, 1.0]],
        [[-1.0, 1.0, -1.0], [1.0, 1.0, -1.0], [0.0, 1.0, 1.0]],
    ];
    let hit = swept_sphere_triangle_mesh([0.0, 5.0, 0.0], [0.0, -5.0, 0.0], 0.3, &triangles)
        .expect("hit");
    assert_eq!(hit.tri_index, 1);
    assert!(approx(hit.toi, 0.37, 1e-12), "toi={}", hit.toi);
    assert!(vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12));
    assert!(vec_close(hit.point_b, [0.0, 1.0, 0.0], 1e-12));
    assert!(vec_close(hit.point_a, [0.0, 1.0, 0.0], 1e-12));
}

#[test]
fn triangle_mesh_prefers_an_earlier_edge_contact_over_a_later_face_hit() {
    // The sphere slides along +X, 0.3 above the corner triangle, and first runs
    // into its edge at x = 0 (t = 0.4); the wall listed before it is reached
    // only when the centre arrives at x = 1.0 (t = 0.75).
    let wall = [[1.5, -1.0, -1.0], [1.5, 1.0, -1.0], [1.5, 0.0, 1.0]];
    let triangles = [wall, CORNER_TRI];
    let hit = swept_sphere_triangle_mesh([-2.0, 0.3, 0.25], [2.0, 0.3, 0.25], 0.5, &triangles)
        .expect("hit");
    assert_eq!(hit.tri_index, 1);
    assert!(approx(hit.toi, 0.4, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [-0.8, 0.6, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    // Without the corner triangle the wall is hit, face on.
    let hit =
        swept_sphere_triangle_mesh([-2.0, 0.3, 0.25], [2.0, 0.3, 0.25], 0.5, &[wall]).expect("hit");
    assert_eq!(hit.tri_index, 0);
    assert!(approx(hit.toi, 0.75, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [-1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
}

#[test]
fn triangle_mesh_ties_go_to_the_lower_index_and_empty_meshes_miss() {
    let triangles = [FLOOR_TRI, FLOOR_TRI];
    let hit = swept_sphere_triangle_mesh([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &triangles)
        .expect("hit");
    assert_eq!(hit.tri_index, 0);
    assert!(swept_sphere_triangle_mesh([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &[]).is_none());
    // A triangle the swept sphere cannot reach is skipped without a hit.
    let far = [[[50.0, 0.0, 0.0], [51.0, 0.0, 0.0], [50.0, 0.0, 1.0]]];
    assert!(swept_sphere_triangle_mesh([0.0, 3.0, 0.0], [0.0, -1.0, 0.0], 0.5, &far).is_none());
}

/// Which feature of the triangle a contact point lies on.
fn contact_feature(tri: &[[f64; 3]; 3], point: [f64; 3]) -> &'static str {
    if tri.iter().any(|&v| len(sub(point, v)) < 1e-7) {
        "vertex"
    } else if (0..3).any(|i| reference_segment_distance(tri[i], tri[(i + 1) % 3], point) < 1e-7) {
        "edge"
    } else {
        "face"
    }
}

#[test]
fn triangle_sweep_matches_a_sampled_reference_on_random_input() {
    let mut rng = TestRng::new(0xC0DE_F00D);
    let (mut vertex_hits, mut edge_hits, mut face_hits, mut misses) = (0, 0, 0, 0);
    for case in 0..600 {
        let tri = [rng.point(1.0), rng.point(1.0), rng.point(1.0)];
        let (c0, c1) = if case % 3 == 0 {
            (rng.point(2.5), rng.point(2.5))
        } else {
            // A path aimed at a point near the triangle.
            let target = rng.point(0.8);
            let heading = rng.direction();
            (
                sub(target, scale(heading, rng.range(0.3, 3.0))),
                add(target, scale(heading, rng.range(0.0, 3.0))),
            )
        };
        let radius = if case % 11 == 0 {
            0.0
        } else {
            rng.range(0.05, 0.8)
        };
        let hit = swept_sphere_triangle(c0, c1, radius, &tri);
        let dist = |t: f64| reference_triangle_distance(&tri, lerp3(c0, c1, t));
        verify_time_of_impact(
            hit.as_ref().map(|h| h.toi),
            radius,
            &dist,
            &format!("case {case}"),
        );
        match hit {
            Some(h) => {
                assert!(
                    approx(len(h.normal), 1.0, 1e-9),
                    "case {case}: normal {:?}",
                    h.normal
                );
                let centre = lerp3(c0, c1, h.toi);
                // The contact point lies on the triangle and the sphere point
                // on the sphere surface, on the line through the centre.
                assert!(
                    reference_triangle_distance(&tri, h.point_b) < 1e-9,
                    "case {case}"
                );
                if h.toi > 0.0 {
                    assert!(
                        vec_close(h.point_a, sub(centre, scale(h.normal, radius)), 1e-9),
                        "case {case}"
                    );
                    assert!(
                        vec_close(centre, add(h.point_b, scale(h.normal, radius)), 1e-8),
                        "case {case}: the normal does not join the contact points"
                    );
                    match contact_feature(&tri, h.point_b) {
                        "vertex" => vertex_hits += 1,
                        "edge" => edge_hits += 1,
                        _ => face_hits += 1,
                    }
                }
            }
            None => misses += 1,
        }
    }
    assert!(
        vertex_hits >= 5 && edge_hits >= 20 && face_hits >= 40 && misses >= 20,
        "feature coverage: vertex {vertex_hits}, edge {edge_hits}, face {face_hits}, miss {misses}"
    );
}

#[test]
fn triangle_sweep_matches_a_sampled_reference_on_lattice_input() {
    // Coordinates on a coarse lattice make parallel motion, exact tangency,
    // zero-area triangles and exactly touching starts common.
    let mut rng = TestRng::new(0x1A77_1CE5);
    let tri_values = [-1.0, 0.0, 1.0];
    let centre_values = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0];
    let radii = [0.0, 0.25, 0.5, 1.0];
    for case in 0..800 {
        let tri = [
            rng.lattice_point(&tri_values),
            rng.lattice_point(&tri_values),
            rng.lattice_point(&tri_values),
        ];
        let c0 = rng.lattice_point(&centre_values);
        let c1 = rng.lattice_point(&centre_values);
        let radius = rng.pick(&radii);
        let hit = swept_sphere_triangle(c0, c1, radius, &tri);
        let dist = |t: f64| reference_triangle_distance(&tri, lerp3(c0, c1, t));
        verify_time_of_impact(
            hit.as_ref().map(|h| h.toi),
            radius,
            &dist,
            &format!("lattice case {case}: {tri:?} {c0:?} -> {c1:?} r={radius}"),
        );
        if let Some(h) = hit {
            assert!(
                approx(len(h.normal), 1.0, 1e-9),
                "lattice case {case}: {:?}",
                h.normal
            );
        }
    }
}

// ---------------------------------------------------------------------------
// swept_sphere_capsule
// ---------------------------------------------------------------------------

#[test]
fn capsule_sweep_detects_the_crossing_far_from_the_axis_midpoint() {
    // Sphere (-10,5,0) -> (10,0.5,0) against the capsule on y = 0, x in [-10, 10]
    // with r_s + r_c = 1: the centre reaches y = 1 at t = 8/9, x = 7.78. The old
    // code replaced the capsule by one static sphere at the axis point nearest the
    // sweep midpoint (the origin) and missed the contact entirely.
    let hit = swept_sphere_capsule(
        [-10.0, 5.0, 0.0],
        [10.0, 0.5, 0.0],
        0.5,
        [-10.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        0.5,
    )
    .expect("the sphere hits the capsule body");
    assert!(approx(hit.toi, 8.0 / 9.0, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [0.0, 1.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    let x = -10.0 + 20.0 * 8.0 / 9.0;
    assert!(
        vec_close(hit.point_a, [x, 0.5, 0.0], 1e-12),
        "{:?}",
        hit.point_a
    );
    assert!(
        vec_close(hit.point_b, [x, 0.5, 0.0], 1e-12),
        "{:?}",
        hit.point_b
    );
}

#[test]
fn capsule_sweep_hits_the_end_cap_head_on() {
    // Along the axis towards the cap at (1, 0, 0): contact at x = 2.
    let hit = swept_sphere_capsule(
        [5.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        0.5,
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        0.5,
    )
    .expect("hit");
    assert!(approx(hit.toi, 0.6, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    // From the other side the first cap is the one at the origin.
    let hit = swept_sphere_capsule(
        [-5.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        0.5,
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        0.5,
    )
    .expect("hit");
    assert!(approx(hit.toi, 0.8, 1e-12), "toi={}", hit.toi);
    assert!(vec_close(hit.normal, [-1.0, 0.0, 0.0], 1e-12));
}

#[test]
fn capsule_sweep_ignores_the_cylinder_beyond_the_end_point() {
    // Coming down at x = 1.5, past the end point at x = 1: the contact is with
    // the end cap sphere, at y = 0.866 where the centre is 1 away from it.
    let hit = swept_sphere_capsule(
        [1.5, 3.0, 0.0],
        [1.5, -3.0, 0.0],
        0.5,
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        0.5,
    )
    .expect("hit");
    let y = 0.75_f64.sqrt();
    assert!(approx(hit.toi, (3.0 - y) / 6.0, 1e-12), "toi={}", hit.toi);
    assert!(
        vec_close(hit.normal, [0.5, y, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
}

#[test]
fn capsule_sweep_parallel_to_the_axis_enters_through_the_cap() {
    let hit = swept_sphere_capsule(
        [0.5, 5.0, 0.0],
        [0.5, -5.0, 0.0],
        0.3,
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        0.5,
    )
    .expect("hit");
    let dy = (0.8_f64 * 0.8 - 0.25).sqrt();
    assert!(
        approx(hit.toi, (5.0 - 1.0 - dy) / 10.0, 1e-12),
        "toi={}",
        hit.toi
    );
    assert!(
        vec_close(hit.normal, [0.5 / 0.8, dy / 0.8, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
}

#[test]
fn capsule_sweep_misses_and_initial_overlap() {
    let (cap0, cap1) = ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]);
    // Moving away, and passing beside the capsule.
    assert!(
        swept_sphere_capsule([5.0, 0.0, 0.0], [10.0, 0.0, 0.0], 0.3, cap0, cap1, 0.5).is_none()
    );
    assert!(
        swept_sphere_capsule([5.0, 5.0, 0.0], [-5.0, 5.0, 0.0], 0.3, cap0, cap1, 0.5).is_none()
    );
    // Stopping short.
    assert!(swept_sphere_capsule([5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3, cap0, cap1, 0.5).is_none());
    // A static sphere that is apart.
    assert!(swept_sphere_capsule([5.0, 0.0, 0.0], [5.0, 0.0, 0.0], 0.3, cap0, cap1, 0.5).is_none());
    // Overlapping at the start: toi 0, the normal points from the axis to the sphere.
    let hit = swept_sphere_capsule([0.3, 0.2, 0.0], [0.5, 0.2, 0.0], 0.3, cap0, cap1, 0.5)
        .expect("overlap");
    assert_eq!(hit.toi, 0.0);
    assert!(
        vec_close(hit.normal, [1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    // Centre exactly on the axis: the normal is perpendicular to it and opposes the motion.
    let hit = swept_sphere_capsule([0.0, 0.2, 0.0], [1.0, 0.2, 0.0], 0.3, cap0, cap1, 0.5)
        .expect("overlap");
    assert_eq!(hit.toi, 0.0);
    assert!(
        vec_close(hit.normal, [-1.0, 0.0, 0.0], 1e-12),
        "{:?}",
        hit.normal
    );
    // Motion along the axis from a centre on it: still a perpendicular normal.
    let hit = swept_sphere_capsule([0.0, 0.2, 0.0], [0.0, 0.9, 0.0], 0.3, cap0, cap1, 0.5)
        .expect("overlap");
    assert!(
        approx(dot(hit.normal, [0.0, 1.0, 0.0]), 0.0, 1e-12),
        "{:?}",
        hit.normal
    );
    assert!(approx(len(hit.normal), 1.0, 1e-12));
}

#[test]
fn capsule_sweep_rejects_invalid_input() {
    let (cap0, cap1) = ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(
        swept_sphere_capsule([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], -0.3, cap0, cap1, 0.5).is_none()
    );
    assert!(
        swept_sphere_capsule([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.3, cap0, cap1, f64::NAN).is_none()
    );
    assert!(
        swept_sphere_capsule([f64::NAN, 0.0, 0.0], [0.0, 0.0, 0.0], 0.3, cap0, cap1, 0.5).is_none()
    );
}

#[test]
fn capsule_with_a_zero_length_axis_is_a_sphere() {
    let mut rng = TestRng::new(0xB10C_5EED);
    let mut hits = 0;
    for case in 0..400 {
        let centre = rng.point(1.0);
        let (r_s, r_c) = (rng.range(0.1, 0.8), rng.range(0.1, 0.8));
        let (c0, c1) = (rng.point(2.5), rng.point(2.5));
        let capsule = swept_sphere_capsule(c0, c1, r_s, centre, centre, r_c);
        let sphere = swept_sphere_sphere(c0, c1, r_s, centre, centre, r_c);
        match (capsule, sphere) {
            (Some(a), Some(b)) => {
                hits += 1;
                assert!(
                    approx(a.toi, b.toi, 1e-9),
                    "case {case}: {} vs {}",
                    a.toi,
                    b.toi
                );
            }
            (None, None) => {}
            (a, b) => panic!("case {case}: capsule {a:?} vs sphere {b:?}"),
        }
    }
    assert!(hits > 50, "hits {hits}");
}

#[test]
fn capsule_sweep_matches_a_sampled_reference_on_random_input() {
    let mut rng = TestRng::new(0xCA95_01E5);
    let (mut lateral, mut caps, mut misses) = (0, 0, 0);
    for case in 0..600 {
        let (cap0, cap1) = (rng.point(1.5), rng.point(1.5));
        let (r_s, r_c) = (rng.range(0.05, 0.6), rng.range(0.05, 0.6));
        let (c0, c1) = (rng.point(4.0), rng.point(4.0));
        let hit = swept_sphere_capsule(c0, c1, r_s, cap0, cap1, r_c);
        let dist = |t: f64| reference_segment_distance(cap0, cap1, lerp3(c0, c1, t));
        verify_time_of_impact(
            hit.as_ref().map(|h| h.toi),
            r_s + r_c,
            &dist,
            &format!("case {case}"),
        );
        match hit {
            Some(h) => {
                assert!(
                    approx(len(h.normal), 1.0, 1e-9),
                    "case {case}: normal {:?}",
                    h.normal
                );
                if h.toi > 0.0 {
                    let centre = lerp3(c0, c1, h.toi);
                    assert!(
                        vec_close(h.point_a, sub(centre, scale(h.normal, r_s)), 1e-9),
                        "case {case}"
                    );
                    assert!(
                        reference_segment_distance(
                            cap0,
                            cap1,
                            sub(h.point_b, scale(h.normal, r_c))
                        ) < 1e-8,
                        "case {case}: point_b is not on the capsule surface"
                    );
                    let (_, u) = closest_point_on_segment(cap0, cap1, centre);
                    if u > 0.0 && u < 1.0 {
                        lateral += 1;
                    } else {
                        caps += 1;
                    }
                }
            }
            None => misses += 1,
        }
    }
    assert!(
        lateral >= 20 && caps >= 20 && misses >= 20,
        "lateral {lateral}, caps {caps}, misses {misses}"
    );
}
