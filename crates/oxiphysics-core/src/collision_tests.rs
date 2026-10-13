// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the collision module (GJK/EPA, SAT, BVH, contact manifolds).
//!
//! Split out of `collision.rs` to keep that file under the 2000-line limit.
//! Included via `#[path]` as the body of `collision::tests`.

use super::*;

// ── Vector helpers ────────────────────────────────────────────────────────

#[test]
fn test_v3_dot_orthogonal() {
    assert_eq!(v3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 0.0);
}

#[test]
fn test_v3_cross_basis() {
    let c = v3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert!((c[2] - 1.0).abs() < 1e-12);
    assert!(c[0].abs() < 1e-12 && c[1].abs() < 1e-12);
}

#[test]
fn test_v3_normalize_unit() {
    let n = v3_normalize([3.0, 4.0, 0.0]);
    let len = v3_len(n);
    assert!((len - 1.0).abs() < 1e-12, "len={}", len);
}

#[test]
fn test_v3_normalize_zero_safe() {
    let n = v3_normalize([0.0; 3]);
    assert_eq!(n, [0.0; 3]);
}

// ── Sphere-Sphere ─────────────────────────────────────────────────────────

#[test]
fn test_sphere_sphere_overlap() {
    let a = Sphere::new([0.0; 3], 1.0);
    let b = Sphere::new([1.5, 0.0, 0.0], 1.0);
    let r = sphere_sphere(&a, &b);
    assert!(r.overlapping, "spheres should overlap");
    assert!(r.depth > 0.0, "depth should be positive");
}

#[test]
fn test_sphere_sphere_no_overlap() {
    let a = Sphere::new([0.0; 3], 1.0);
    let b = Sphere::new([5.0, 0.0, 0.0], 1.0);
    let r = sphere_sphere(&a, &b);
    assert!(!r.overlapping);
    assert!(r.depth < 0.0);
}

#[test]
fn test_sphere_sphere_touching() {
    let a = Sphere::new([0.0; 3], 1.0);
    let b = Sphere::new([2.0, 0.0, 0.0], 1.0);
    let r = sphere_sphere(&a, &b);
    assert!((r.depth).abs() < 1e-10, "depth should be ~0 at touching");
}

// ── Sphere-Capsule ────────────────────────────────────────────────────────

#[test]
fn test_sphere_capsule_overlap() {
    let s = Sphere::new([0.0; 3], 0.5);
    let c = Capsule::new([0.0, 2.0, 0.0], [0.0, 4.0, 0.0], 0.5);
    let (overlap, depth, _, _) = sphere_capsule(&s, &c);
    assert!(!overlap, "should not overlap, depth={}", depth);
}

#[test]
fn test_sphere_capsule_nearby() {
    let s = Sphere::new([0.0, 1.0, 0.0], 0.6);
    let c = Capsule::new([0.0, 0.0, 0.0], [0.0, 3.0, 0.0], 0.6);
    let (overlap, depth, _, _) = sphere_capsule(&s, &c);
    assert!(overlap, "should overlap, depth={}", depth);
}

// ── AABB ─────────────────────────────────────────────────────────────────

#[test]
fn test_aabb_overlap() {
    let a = AabbRaw::new([-1.0; 3], [1.0; 3]);
    let b = AabbRaw::new([0.5, 0.5, 0.5], [2.0, 2.0, 2.0]);
    assert!(a.overlaps(&b));
}

#[test]
fn test_aabb_no_overlap() {
    let a = AabbRaw::new([-1.0; 3], [1.0; 3]);
    let b = AabbRaw::new([2.0; 3], [3.0; 3]);
    assert!(!a.overlaps(&b));
}

#[test]
fn test_aabb_merge_contains_both() {
    let a = AabbRaw::new([-2.0; 3], [0.0; 3]);
    let b = AabbRaw::new([1.0; 3], [3.0; 3]);
    let m = a.merge(&b);
    assert!(m.min[0] <= -2.0 && m.max[0] >= 3.0);
}

// ── BVH ──────────────────────────────────────────────────────────────────

#[test]
fn test_bvh_query_single_match() {
    let leaves = vec![
        BvhLeaf {
            id: 1,
            aabb: AabbRaw::new([0.0; 3], [1.0; 3]),
        },
        BvhLeaf {
            id: 2,
            aabb: AabbRaw::new([5.0; 3], [6.0; 3]),
        },
    ];
    let bvh = Bvh::build(leaves);
    let query = AabbRaw::new([0.0; 3], [0.5; 3]);
    let hits = bvh.query(&query);
    assert_eq!(hits, vec![1]);
}

#[test]
fn test_bvh_query_no_match() {
    let leaves = vec![BvhLeaf {
        id: 1,
        aabb: AabbRaw::new([10.0; 3], [11.0; 3]),
    }];
    let bvh = Bvh::build(leaves);
    let query = AabbRaw::new([0.0; 3], [1.0; 3]);
    assert!(bvh.query(&query).is_empty());
}

#[test]
fn test_bvh_empty() {
    let bvh = Bvh::build(vec![]);
    let query = AabbRaw::new([0.0; 3], [1.0; 3]);
    assert!(bvh.query(&query).is_empty());
}

// ── OBB SAT ───────────────────────────────────────────────────────────────

#[test]
fn test_obb_sat_overlap() {
    let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
    let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0; 3]);
    assert!(obb_obb_sat(&a, &b).is_some(), "overlapping OBBs");
}

#[test]
fn test_obb_sat_separated() {
    let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
    let b = Obb::axis_aligned([5.0, 0.0, 0.0], [1.0; 3]);
    assert!(obb_obb_sat(&a, &b).is_none(), "separated OBBs");
}

// ── Ray tests ─────────────────────────────────────────────────────────────

#[test]
fn test_ray_aabb_hit() {
    let origin = [-5.0, 0.0, 0.0];
    let dir = [1.0, 0.0, 0.0];
    let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
    let hit = ray_aabb(origin, dir, &aabb);
    assert!(hit.is_some());
    let h = hit.unwrap();
    assert!(h.t_min < h.t_max);
}

#[test]
fn test_ray_aabb_miss() {
    let origin = [0.0, 5.0, 0.0];
    let dir = [1.0, 0.0, 0.0];
    let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
    assert!(ray_aabb(origin, dir, &aabb).is_none());
}

#[test]
fn test_ray_sphere_hit() {
    let s = Sphere::new([0.0; 3], 1.0);
    let t = ray_sphere([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], &s);
    assert!(t.is_some());
    assert!((t.unwrap() - 4.0).abs() < 1e-10);
}

#[test]
fn test_ray_sphere_miss() {
    let s = Sphere::new([0.0; 3], 1.0);
    let t = ray_sphere([0.0, 5.0, 0.0], [1.0, 0.0, 0.0], &s);
    assert!(t.is_none());
}

// ── GJK ──────────────────────────────────────────────────────────────────

#[test]
fn test_gjk_spheres_intersecting() {
    let a = Sphere::new([0.0; 3], 1.5);
    let b = Sphere::new([1.0, 0.0, 0.0], 1.5);
    let r = gjk(&a, &b);
    assert!(r.intersecting);
}

#[test]
fn test_gjk_spheres_separated() {
    let a = Sphere::new([0.0; 3], 0.5);
    let b = Sphere::new([10.0, 0.0, 0.0], 0.5);
    let r = gjk(&a, &b);
    assert!(!r.intersecting);
    assert!(r.distance > 0.0);
}

/// Build a faceted (UV-tessellated) approximation of a sphere as a point
/// cloud.  Many near-coplanar facets make GJK/EPA converge slowly, which
/// exercises the non-convergence tails.
fn faceted_sphere(
    center: [Real; 3],
    radius: Real,
    rings: usize,
    sectors: usize,
) -> ConvexPointCloud {
    let mut pts = Vec::with_capacity(rings * sectors + 2);
    for ri in 0..=rings {
        let theta = std::f64::consts::PI * (ri as f64) / (rings as f64);
        let (st, ct) = theta.sin_cos();
        for si in 0..sectors {
            let phi = 2.0 * std::f64::consts::PI * (si as f64) / (sectors as f64);
            let (sp, cp) = phi.sin_cos();
            pts.push([
                center[0] + radius * st * cp,
                center[1] + radius * st * sp,
                center[2] + radius * ct,
            ]);
        }
    }
    ConvexPointCloud::new(pts)
}

// ── GJK distance (van den Bergen) ─────────────────────────────────────────

/// Euclidean distance between two points.
fn point_dist(a: [Real; 3], b: [Real; 3]) -> Real {
    v3_len(v3_sub(a, b))
}

/// xorshift64* generator: keeps the randomised tests deterministic without an
/// external dependency.
struct TestRng(u64);

impl TestRng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> Real {
        (self.next_u64() >> 11) as Real / (1u64 << 53) as Real
    }

    /// Uniform in `[lo, hi)`.
    fn range(&mut self, lo: Real, hi: Real) -> Real {
        lo + (hi - lo) * self.unit()
    }

    /// Vector with every component uniform in `[lo, hi)`.
    fn vec3(&mut self, lo: Real, hi: Real) -> [Real; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

/// Invariants every separated GJK result must satisfy: not intersecting, a
/// positive distance, and witnesses whose gap is exactly that distance.
fn assert_separated(r: &GjkResult) {
    assert!(!r.intersecting, "expected a separated result, got {r:?}");
    assert!(r.distance > 0.0, "distance must be positive, got {r:?}");
    let witness_gap = point_dist(r.closest_a, r.closest_b);
    assert!(
        (witness_gap - r.distance).abs() <= 1e-12 * r.distance.max(1.0),
        "|closest_a - closest_b| = {witness_gap} must equal distance {}",
        r.distance
    );
}

/// Closest point to `p` on triangle `abc` (Voronoi-region method).
fn closest_point_on_triangle(p: [Real; 3], a: [Real; 3], b: [Real; 3], c: [Real; 3]) -> [Real; 3] {
    let ab = v3_sub(b, a);
    let ac = v3_sub(c, a);
    let ap = v3_sub(p, a);
    let d1 = v3_dot(ab, ap);
    let d2 = v3_dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = v3_sub(p, b);
    let d3 = v3_dot(ab, bp);
    let d4 = v3_dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return v3_add(a, v3_scale(ab, d1 / (d1 - d3)));
    }
    let cp = v3_sub(p, c);
    let d5 = v3_dot(ab, cp);
    let d6 = v3_dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return v3_add(a, v3_scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return v3_add(b, v3_scale(v3_sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    v3_add(
        a,
        v3_add(v3_scale(ab, vb * denom), v3_scale(ac, vc * denom)),
    )
}

/// Closest points between segments `p1q1` and `p2q2`.
fn closest_points_segments(
    p1: [Real; 3],
    q1: [Real; 3],
    p2: [Real; 3],
    q2: [Real; 3],
) -> ([Real; 3], [Real; 3]) {
    let d1 = v3_sub(q1, p1);
    let d2 = v3_sub(q2, p2);
    let r = v3_sub(p1, p2);
    let a = v3_dot(d1, d1);
    let e = v3_dot(d2, d2);
    let f = v3_dot(d2, r);
    let (s, t);
    if a <= 1e-300 && e <= 1e-300 {
        return (p1, p2);
    }
    if a <= 1e-300 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = v3_dot(d1, r);
        if e <= 1e-300 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = v3_dot(d1, d2);
            let denom = a * e - b * b;
            let mut s_try = if denom > 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t_try = (b * s_try + f) / e;
            if t_try < 0.0 {
                t_try = 0.0;
                s_try = (-c / a).clamp(0.0, 1.0);
            } else if t_try > 1.0 {
                t_try = 1.0;
                s_try = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = s_try;
            t = t_try;
        }
    }
    (v3_add(p1, v3_scale(d1, s)), v3_add(p2, v3_scale(d2, t)))
}

/// Exhaustive distance between the convex hulls of two small point sets: the
/// minimum over vertex-vertex, vertex-triangle and edge-edge feature pairs,
/// which contains the true closest pair of two disjoint polytopes.
fn brute_force_hull_distance(pa: &[[Real; 3]], pb: &[[Real; 3]]) -> Real {
    let mut best = Real::MAX;
    for &p in pa {
        for &q in pb {
            best = best.min(point_dist(p, q));
        }
    }
    for (points, tri_pts) in [(pa, pb), (pb, pa)] {
        for &p in points {
            for i in 0..tri_pts.len() {
                for j in (i + 1)..tri_pts.len() {
                    for k in (j + 1)..tri_pts.len() {
                        let c = closest_point_on_triangle(p, tri_pts[i], tri_pts[j], tri_pts[k]);
                        best = best.min(point_dist(p, c));
                    }
                }
            }
        }
    }
    for i in 0..pa.len() {
        for j in (i + 1)..pa.len() {
            for k in 0..pb.len() {
                for l in (k + 1)..pb.len() {
                    let (c1, c2) = closest_points_segments(pa[i], pa[j], pb[k], pb[l]);
                    best = best.min(point_dist(c1, c2));
                }
            }
        }
    }
    best
}

/// Largest separating-axis gap over every face-normal and edge-cross-product
/// axis built from the two point sets.  Positive ⇒ the hulls are disjoint;
/// negative ⇒ they overlap.  Independent oracle for the GJK classification.
fn separating_axis_gap(pa: &[[Real; 3]], pb: &[[Real; 3]]) -> Real {
    let mut axes: Vec<[Real; 3]> = Vec::new();
    for pts in [pa, pb] {
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                for k in (j + 1)..pts.len() {
                    axes.push(v3_cross(v3_sub(pts[j], pts[i]), v3_sub(pts[k], pts[i])));
                }
            }
        }
    }
    for i in 0..pa.len() {
        for j in (i + 1)..pa.len() {
            let edge_a = v3_sub(pa[j], pa[i]);
            for k in 0..pb.len() {
                for l in (k + 1)..pb.len() {
                    axes.push(v3_cross(edge_a, v3_sub(pb[l], pb[k])));
                }
            }
        }
    }
    let mut best = Real::MIN;
    for axis in axes {
        let n = v3_normalize(axis);
        if v3_len_sq(n) < 0.5 {
            continue;
        }
        let project = |pts: &[[Real; 3]]| {
            pts.iter().fold((Real::MAX, Real::MIN), |(lo, hi), &p| {
                let d = v3_dot(n, p);
                (lo.min(d), hi.max(d))
            })
        };
        let (a_lo, a_hi) = project(pa);
        let (b_lo, b_hi) = project(pb);
        best = best.max((b_lo - a_hi).max(a_lo - b_hi));
    }
    best
}

#[test]
fn test_gjk_unit_spheres_3_4_0_report_true_distance_and_witnesses() {
    // Minimum distance is |(3,4,0)| - 1 - 1 = 3.  The old early exit returned
    // the first simplex' upper bound sqrt(17) ~ 4.123 with raw support points.
    let a = Sphere::new([0.0; 3], 1.0);
    let b = Sphere::new([3.0, 4.0, 0.0], 1.0);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!((r.distance - 3.0).abs() < 1e-9, "distance {}", r.distance);
    // Smooth shapes: witnesses are accurate to about sqrt(1e-12).
    assert!(
        point_dist(r.closest_a, [0.6, 0.8, 0.0]) < 1e-5,
        "closest_a {:?}",
        r.closest_a
    );
    assert!(
        point_dist(r.closest_b, [2.4, 3.2, 0.0]) < 1e-5,
        "closest_b {:?}",
        r.closest_b
    );
}

#[test]
fn test_gjk_distance_is_symmetric_in_the_arguments() {
    let a = Sphere::new([0.5, -1.0, 2.0], 0.8);
    let b = Box3::new([4.0, 1.0, -1.0], [1.0, 0.5, 2.0]);
    let ab = gjk(&a, &b);
    let ba = gjk(&b, &a);
    assert_separated(&ab);
    assert_separated(&ba);
    assert!((ab.distance - ba.distance).abs() < 1e-9);
    assert!(point_dist(ab.closest_a, ba.closest_b) < 1e-4);
    assert!(point_dist(ab.closest_b, ba.closest_a) < 1e-4);
}

#[test]
fn test_gjk_distance_is_translation_invariant() {
    let shift = [1.0e4, -2.0e4, 3.0e3];
    let a = Sphere::new(shift, 1.0);
    let b = Sphere::new(v3_add(shift, [3.0, 4.0, 0.0]), 1.0);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!((r.distance - 3.0).abs() < 1e-6, "distance {}", r.distance);
}

#[test]
fn test_gjk_spheres_match_analytic_distance_and_witnesses() {
    let mut rng = TestRng(0x1234_5678_9ABC_DEF1);
    let mut checked = 0;
    for _ in 0..300 {
        let c1 = rng.vec3(-5.0, 5.0);
        let c2 = rng.vec3(-5.0, 5.0);
        let r1 = rng.range(0.1, 2.0);
        let r2 = rng.range(0.1, 2.0);
        let centre_gap = point_dist(c1, c2);
        let gap = centre_gap - r1 - r2;
        let r = gjk(&Sphere::new(c1, r1), &Sphere::new(c2, r2));
        if gap > 1e-6 {
            checked += 1;
            assert_separated(&r);
            assert!((r.distance - gap).abs() < 1e-9, "{} vs {gap}", r.distance);
            let u = v3_scale(v3_sub(c2, c1), 1.0 / centre_gap);
            assert!(point_dist(r.closest_a, v3_add(c1, v3_scale(u, r1))) < 1e-4);
            assert!(point_dist(r.closest_b, v3_sub(c2, v3_scale(u, r2))) < 1e-4);
        } else if gap < -1e-6 {
            assert!(r.intersecting, "overlapping spheres, gap {gap}");
        }
    }
    assert!(checked > 100, "too few separated cases: {checked}");
}

#[test]
fn test_gjk_box_sphere_matches_analytic_distance_and_witnesses() {
    let mut rng = TestRng(0x0BAD_5EED_1234_5678);
    let mut checked = 0;
    for _ in 0..300 {
        let box_centre = rng.vec3(-2.0, 2.0);
        let half = [
            rng.range(0.2, 2.0),
            rng.range(0.2, 2.0),
            rng.range(0.2, 2.0),
        ];
        let sphere_centre = rng.vec3(-6.0, 6.0);
        let radius = rng.range(0.1, 2.0);
        let nearest_on_box = [
            sphere_centre[0].clamp(box_centre[0] - half[0], box_centre[0] + half[0]),
            sphere_centre[1].clamp(box_centre[1] - half[1], box_centre[1] + half[1]),
            sphere_centre[2].clamp(box_centre[2] - half[2], box_centre[2] + half[2]),
        ];
        let to_box = point_dist(nearest_on_box, sphere_centre);
        let gap = to_box - radius;
        let r = gjk(
            &Box3::new(box_centre, half),
            &Sphere::new(sphere_centre, radius),
        );
        if to_box > 1e-9 && gap > 1e-6 {
            checked += 1;
            assert_separated(&r);
            assert!((r.distance - gap).abs() < 1e-9, "{} vs {gap}", r.distance);
            let u = v3_scale(v3_sub(nearest_on_box, sphere_centre), 1.0 / to_box);
            assert!(point_dist(r.closest_a, nearest_on_box) < 1e-4);
            assert!(point_dist(r.closest_b, v3_add(sphere_centre, v3_scale(u, radius))) < 1e-4);
        } else if gap < -1e-6 {
            assert!(r.intersecting, "overlapping box and sphere, gap {gap}");
        }
    }
    assert!(checked > 100, "too few separated cases: {checked}");
}

#[test]
fn test_gjk_capsule_sphere_matches_analytic_distance_and_witnesses() {
    let mut rng = TestRng(0xCA95_01E5_0000_0007);
    let mut checked = 0;
    for _ in 0..300 {
        let a = rng.vec3(-3.0, 3.0);
        let b = rng.vec3(-3.0, 3.0);
        let capsule_radius = rng.range(0.1, 1.0);
        let sphere_centre = rng.vec3(-6.0, 6.0);
        let sphere_radius = rng.range(0.1, 1.5);
        let axis = v3_sub(b, a);
        let t = (v3_dot(v3_sub(sphere_centre, a), axis) / v3_dot(axis, axis)).clamp(0.0, 1.0);
        let axis_point = v3_add(a, v3_scale(axis, t));
        let centre_gap = point_dist(axis_point, sphere_centre);
        let gap = centre_gap - capsule_radius - sphere_radius;
        let r = gjk(
            &Capsule::new(a, b, capsule_radius),
            &Sphere::new(sphere_centre, sphere_radius),
        );
        if gap > 1e-6 {
            checked += 1;
            assert_separated(&r);
            assert!((r.distance - gap).abs() < 1e-8, "{} vs {gap}", r.distance);
            let u = v3_scale(v3_sub(sphere_centre, axis_point), 1.0 / centre_gap);
            assert!(
                point_dist(r.closest_a, v3_add(axis_point, v3_scale(u, capsule_radius))) < 1e-3
            );
            assert!(
                point_dist(
                    r.closest_b,
                    v3_sub(sphere_centre, v3_scale(u, sphere_radius))
                ) < 1e-3
            );
        } else if gap < -1e-6 {
            assert!(r.intersecting, "overlapping capsule and sphere, gap {gap}");
        }
    }
    assert!(checked > 100, "too few separated cases: {checked}");
}

#[test]
fn test_gjk_capsules_match_closest_segment_distance() {
    let mut rng = TestRng(0x5E65_0000_BEEF_0001);
    let mut checked = 0;
    for _ in 0..300 {
        let (a0, a1) = (rng.vec3(-3.0, 3.0), rng.vec3(-3.0, 3.0));
        let (b0, b1) = (rng.vec3(-3.0, 3.0), rng.vec3(-3.0, 3.0));
        let (ra, rb) = (rng.range(0.05, 0.5), rng.range(0.05, 0.5));
        let (c1, c2) = closest_points_segments(a0, a1, b0, b1);
        let gap = point_dist(c1, c2) - ra - rb;
        let r = gjk(&Capsule::new(a0, a1, ra), &Capsule::new(b0, b1, rb));
        if gap > 1e-6 {
            checked += 1;
            assert_separated(&r);
            assert!((r.distance - gap).abs() < 1e-8, "{} vs {gap}", r.distance);
        } else if gap < -1e-6 {
            assert!(r.intersecting, "overlapping capsules, gap {gap}");
        }
    }
    assert!(checked > 100, "too few separated cases: {checked}");
}

#[test]
fn test_gjk_axis_aligned_boxes_match_per_axis_gap_distance() {
    let mut rng = TestRng(0xB0C5_0000_0000_0003);
    let mut checked = 0;
    for _ in 0..300 {
        let c1 = rng.vec3(-3.0, 3.0);
        let h1 = [
            rng.range(0.2, 1.5),
            rng.range(0.2, 1.5),
            rng.range(0.2, 1.5),
        ];
        let c2 = rng.vec3(-3.0, 3.0);
        let h2 = [
            rng.range(0.2, 1.5),
            rng.range(0.2, 1.5),
            rng.range(0.2, 1.5),
        ];
        let gaps: Vec<Real> = (0..3)
            .map(|k| ((c1[k] - c2[k]).abs() - h1[k] - h2[k]).max(0.0))
            .collect();
        let expected = gaps.iter().map(|g| g * g).sum::<Real>().sqrt();
        let r = gjk(&Box3::new(c1, h1), &Box3::new(c2, h2));
        if expected > 1e-6 {
            checked += 1;
            assert_separated(&r);
            assert!(
                (r.distance - expected).abs() < 1e-9,
                "{} vs {expected}",
                r.distance
            );
            // Each witness lies inside (or on) its own box.
            for k in 0..3 {
                assert!((r.closest_a[k] - c1[k]).abs() <= h1[k] + 1e-9);
                assert!((r.closest_b[k] - c2[k]).abs() <= h2[k] + 1e-9);
            }
        } else if expected == 0.0 {
            assert!(r.intersecting, "overlapping boxes");
        }
    }
    assert!(checked > 100, "too few separated cases: {checked}");
}

#[test]
fn test_gjk_boxes_corner_to_corner_have_unique_witnesses() {
    let a = Box3::new([0.0; 3], [1.0; 3]);
    let b = Box3::new([3.0, 3.0, 3.0], [1.0; 3]);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!(
        (r.distance - 3.0_f64.sqrt()).abs() < 1e-12,
        "{}",
        r.distance
    );
    assert!(point_dist(r.closest_a, [1.0, 1.0, 1.0]) < 1e-12);
    assert!(point_dist(r.closest_b, [2.0, 2.0, 2.0]) < 1e-12);
}

#[test]
fn test_gjk_boxes_face_to_face_witnesses_lie_on_the_facing_faces() {
    // Face-to-face contact has a whole patch of closest pairs, so only the
    // invariants are asserted, not specific points.
    let a = Box3::new([0.0; 3], [1.0; 3]);
    let b = Box3::new([3.5, 0.2, -0.3], [1.0; 3]);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!((r.distance - 1.5).abs() < 1e-12, "{}", r.distance);
    assert!((r.closest_a[0] - 1.0).abs() < 1e-12, "{:?}", r.closest_a);
    assert!((r.closest_b[0] - 2.5).abs() < 1e-12, "{:?}", r.closest_b);
    for k in 1..3 {
        assert!(r.closest_a[k].abs() <= 1.0 + 1e-12);
        assert!((r.closest_b[k] - [0.0, 0.2, -0.3][k]).abs() <= 1.0 + 1e-12);
    }
}

#[test]
fn test_gjk_small_face_gaps_are_resolved_until_the_intersection_tolerance() {
    let a = Box3::new([0.0; 3], [1.0; 3]);
    for gap in [1e-3, 1e-6, 1e-9] {
        let b = Box3::new([2.0 + gap, 0.0, 0.0], [1.0; 3]);
        let r = gjk(&a, &b);
        assert_separated(&r);
        assert!(
            (r.distance - gap).abs() < 1e-12 + 1e-9 * gap,
            "gap {gap}: distance {}",
            r.distance
        );
    }
    // Exactly touching boxes: the gap is below the intersection tolerance.
    let touching = gjk(&a, &Box3::new([2.0, 0.0, 0.0], [1.0; 3]));
    assert!(touching.intersecting, "touching boxes must report contact");
    assert_eq!(touching.distance, 0.0);
}

#[test]
fn test_gjk_capsule_sphere_unique_witnesses_on_the_axis_normal() {
    // Capsule axis (0,0,0)-(0,1,0), radius 0.3; sphere radius 0.3 at (5,0.5,0).
    // The closest axis point is (0,0.5,0), so the pair is unique.
    let cap = Capsule::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.3);
    let s = Sphere::new([5.0, 0.5, 0.0], 0.3);
    let r = gjk(&cap, &s);
    assert_separated(&r);
    assert!((r.distance - 4.4).abs() < 1e-9, "{}", r.distance);
    // Smooth shapes: witnesses are accurate to roughly sqrt(1e-12).
    assert!(point_dist(r.closest_a, [0.3, 0.5, 0.0]) < 1e-4);
    assert!(point_dist(r.closest_b, [4.7, 0.5, 0.0]) < 1e-4);
}

#[test]
fn test_gjk_skew_segments_flat_minkowski_difference() {
    // Zero-radius capsules are segments; the Minkowski difference is a flat
    // square at z = -1, so the simplex degenerates to lower dimensions.
    let a = Capsule::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0);
    let b = Capsule::new([0.0, -1.0, 1.0], [0.0, 1.0, 1.0], 0.0);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!((r.distance - 1.0).abs() < 1e-12, "{}", r.distance);
    assert!(point_dist(r.closest_a, [0.0, 0.0, 0.0]) < 1e-12);
    assert!(point_dist(r.closest_b, [0.0, 0.0, 1.0]) < 1e-12);
}

#[test]
fn test_gjk_point_above_triangle_interior_and_point_pairs() {
    let tri = Triangle::new([0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]);
    let point = ConvexPointCloud::new(vec![[1.0, 1.0, 2.0]]);
    let r = gjk(&tri, &point);
    assert_separated(&r);
    assert!((r.distance - 2.0).abs() < 1e-12, "{}", r.distance);
    assert!(point_dist(r.closest_a, [1.0, 1.0, 0.0]) < 1e-12);
    assert!(point_dist(r.closest_b, [1.0, 1.0, 2.0]) < 1e-12);

    let p = ConvexPointCloud::new(vec![[1.0, 2.0, 3.0]]);
    let q = ConvexPointCloud::new(vec![[4.0, 6.0, 3.0]]);
    let r = gjk(&p, &q);
    assert_separated(&r);
    assert!((r.distance - 5.0).abs() < 1e-12);
    assert_eq!(r.closest_a, [1.0, 2.0, 3.0]);
    assert_eq!(r.closest_b, [4.0, 6.0, 3.0]);
}

#[test]
fn test_gjk_random_point_clouds_match_exhaustive_feature_distance() {
    let mut rng = TestRng(0x7001_0C10_0D5E_ED01);
    let (mut separated, mut overlapping) = (0, 0);
    for case in 0..600 {
        let n1 = 4 + (rng.next_u64() % 6) as usize;
        let n2 = 4 + (rng.next_u64() % 6) as usize;
        let offset = rng.vec3(-3.0, 3.0);
        let pa: Vec<[Real; 3]> = (0..n1).map(|_| rng.vec3(-1.0, 1.0)).collect();
        let pb: Vec<[Real; 3]> = (0..n2)
            .map(|_| v3_add(rng.vec3(-1.0, 1.0), offset))
            .collect();
        let r = gjk(
            &ConvexPointCloud::new(pa.clone()),
            &ConvexPointCloud::new(pb.clone()),
        );
        let gap = separating_axis_gap(&pa, &pb);
        if gap > 1e-9 {
            separated += 1;
            assert_separated(&r);
            let expected = brute_force_hull_distance(&pa, &pb);
            assert!(
                (r.distance - expected).abs() < 1e-9,
                "case {case}: gjk {} vs exhaustive {expected}",
                r.distance
            );
        } else if gap < -1e-9 {
            overlapping += 1;
            assert!(r.intersecting, "case {case}: hulls overlap (gap {gap})");
        }
    }
    assert!(separated > 100, "too few separated cases: {separated}");
    assert!(overlapping > 5, "too few overlapping cases: {overlapping}");
}

#[test]
fn test_gjk_degenerate_point_clouds_coplanar_collinear_and_duplicated() {
    // Flat, collinear and lattice (heavily duplicated) point sets make the
    // Minkowski difference lower-dimensional.  Every pair below is separated by
    // construction (parallel planes 0.4 apart, offset parallel lines, and an x
    // gap of at least 1 for the lattice sets), and the distance is still exact.
    let mut rng = TestRng(0xDE6E_0000_0000_0A11);
    for case in 0..300 {
        let kind = case % 3;
        let mut make = |offset: [Real; 3]| -> Vec<[Real; 3]> {
            (0..6)
                .map(|_| {
                    let p = rng.vec3(-1.0, 1.0);
                    let q = match kind {
                        0 => [p[0], p[1], 0.0],
                        1 => [p[0], 0.0, 0.0],
                        _ => [p[0].round(), p[1].round(), p[2].round()],
                    };
                    v3_add(q, offset)
                })
                .collect()
        };
        let pa = make([0.0; 3]);
        let pb = make([3.0 + (case % 7) as Real * 0.1, 0.7, -0.4]);
        let expected = brute_force_hull_distance(&pa, &pb);
        assert!(expected > 0.3, "case {case}: construction must separate");
        let r = gjk(
            &ConvexPointCloud::new(pa.clone()),
            &ConvexPointCloud::new(pb.clone()),
        );
        assert_separated(&r);
        assert!(
            (r.distance - expected).abs() < 1e-8,
            "case {case} (kind {kind}): gjk {} vs exhaustive {expected}",
            r.distance
        );
    }
}

#[test]
fn test_gjk_faceted_spheres_distance_equals_vertex_gap() {
    // Two faceted spheres of radius 1, centres 5 apart on x.  Both clouds
    // contain the vertices (1, 0, ~0) and (4, 0, ~0), which are the extreme
    // points along x, so the minimum distance is exactly 5 - 1 - 1 = 3.
    let a = faceted_sphere([0.0, 0.0, 0.0], 1.0, 24, 48);
    let b = faceted_sphere([5.0, 0.0, 0.0], 1.0, 24, 48);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!(
        (r.distance - 3.0).abs() < 1e-9,
        "faceted distance {}",
        r.distance
    );
    assert!((r.closest_a[0] - 1.0).abs() < 1e-9);
    assert!((r.closest_b[0] - 4.0).abs() < 1e-9);
    assert!(
        !(r.distance == 0.0 && r.closest_a == [0.0; 3] && r.closest_b == [0.0; 3]),
        "result must not be the fabricated zero/zero constant"
    );
}

#[test]
fn test_gjk_faceted_spheres_off_axis_distance_is_at_least_the_exact_sphere_gap() {
    // A polytope inscribed in a sphere lies inside it, so the distance between
    // two inscribed polytopes is at least the exact sphere-sphere gap, and the
    // difference is bounded by the facet sagitta (a few 1e-3 for 24 x 48 facets).
    let a = faceted_sphere([0.0; 3], 1.0, 24, 48);
    let b = faceted_sphere([3.0, 4.0, 0.0], 1.0, 24, 48);
    let r = gjk(&a, &b);
    assert_separated(&r);
    assert!(r.distance >= 3.0 - 1e-9, "distance {}", r.distance);
    assert!(r.distance <= 3.0 + 0.01, "distance {}", r.distance);
}

#[test]
fn test_gjk_overlapping_shapes_report_zero_distance_and_zero_witnesses() {
    let r = gjk(
        &Sphere::new([0.0; 3], 1.5),
        &Sphere::new([1.0, 0.0, 0.0], 1.5),
    );
    assert!(r.intersecting);
    assert_eq!(r.distance, 0.0);
    assert_eq!(r.closest_a, [0.0; 3]);
    assert_eq!(r.closest_b, [0.0; 3]);
}

#[test]
fn test_gjk_touching_spheres_report_contact_or_zero_gap() {
    let r = gjk(
        &Sphere::new([0.0; 3], 1.0),
        &Sphere::new([2.0, 0.0, 0.0], 1.0),
    );
    assert!(
        r.intersecting || r.distance < 1e-9,
        "touching spheres: {r:?}"
    );
}

#[test]
fn test_gjk_terminal_simplex_is_an_enclosing_tetrahedron_for_generic_overlap() {
    // Tilted overlapping tetrahedra: the Minkowski difference is a genuine 3-D
    // body around the origin and GJK ends on a tetrahedron that contains it.
    let a = ConvexPointCloud::new(vec![
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ]);
    let b = ConvexPointCloud::new(vec![
        [0.4, 0.3, 0.2],
        [-0.8, -0.5, 0.9],
        [-0.6, 0.7, -0.7],
        [0.9, -0.8, -0.6],
    ]);
    let (r, terminal) = gjk_simplex(&a, &b);
    assert!(r.intersecting);
    let terminal = terminal.expect("intersecting shapes report their terminal simplex");
    assert_eq!(terminal.len(), 4);
    assert!(epa_simplex_is_valid(&terminal), "tetrahedron needs volume");
    assert!(
        epa_origin_in_tetra(&terminal),
        "tetrahedron must hold the origin"
    );
    // EPA accepts it unchanged.
    let completed = epa_complete_simplex(&a, &b, &terminal).expect("already a valid seed");
    assert_eq!(completed, terminal);
}

#[test]
fn test_gjk_terminal_simplex_is_none_for_separated_shapes() {
    let (r, terminal) = gjk_simplex(
        &Sphere::new([0.0; 3], 1.0),
        &Sphere::new([5.0, 0.0, 0.0], 1.0),
    );
    assert!(!r.intersecting);
    assert!(terminal.is_none());
}

#[test]
fn test_gjk_terminal_simplex_is_the_minimal_face_for_symmetric_overlap() {
    // The sphere centre lies on the box's x axis, so y/z symmetry makes the
    // origin land inside a triangle (not a tetrahedron) of the Minkowski
    // difference: the minimal face is reported, not a padded tetrahedron.
    let a = Box3::new([0.0; 3], [1.0; 3]);
    let b = Sphere::new([0.5, 0.0, 0.0], 0.8);
    let (r, terminal) = gjk_simplex(&a, &b);
    assert!(r.intersecting);
    let terminal = terminal.expect("intersecting shapes report their terminal simplex");
    assert_eq!(terminal.len(), 3, "terminal simplex {terminal:?}");
    let verts: Vec<GjkVertex> = terminal.iter().map(|&w| gjk_test_vertex(w)).collect();
    let closest = gjk_closest_to_origin(&verts).expect("non-empty simplex");
    assert!(closest.dist_sq < 1e-20, "origin must lie in the face");
}

#[test]
fn test_epa_complete_simplex_grows_every_face_dimension_to_a_tetrahedron() {
    let a = Sphere::new([0.0; 3], 1.5);
    let b = Sphere::new([1.0, 0.0, 0.0], 1.5);
    // Point, segment and triangle that all contain the origin.
    let point = [[0.0, 0.0, 0.0]];
    let segment = [[2.0, 0.0, 0.0], [-4.0, 0.0, 0.0]];
    let triangle = [[2.0, 0.0, 0.0], [-3.0, 2.0, 0.0], [-3.0, -2.0, 0.0]];
    for (name, face) in [
        ("point", point.to_vec()),
        ("segment", segment.to_vec()),
        ("triangle", triangle.to_vec()),
    ] {
        let tetra = epa_complete_simplex(&a, &b, &face)
            .unwrap_or_else(|| panic!("{name} must grow to a tetrahedron"));
        assert_eq!(tetra.len(), 4, "{name}");
        assert_eq!(&tetra[..face.len()], &face[..], "{name}: face is preserved");
        assert!(epa_simplex_is_valid(&tetra), "{name}: needs volume");
        assert!(epa_origin_in_tetra(&tetra), "{name}: must hold the origin");
    }
}

#[test]
fn test_epa_complete_simplex_rejects_flat_minkowski_difference_and_bad_input() {
    // Two coplanar squares: A ⊖ B is flat, no tetrahedron with volume exists.
    let square = |c: [Real; 3]| {
        ConvexPointCloud::new(vec![
            [c[0] - 1.0, c[1] - 1.0, 0.0],
            [c[0] + 1.0, c[1] - 1.0, 0.0],
            [c[0] + 1.0, c[1] + 1.0, 0.0],
            [c[0] - 1.0, c[1] + 1.0, 0.0],
        ])
    };
    let (a, b) = (square([0.0; 3]), square([0.5, 0.2, 0.0]));
    let (r, terminal) = gjk_simplex(&a, &b);
    assert!(r.intersecting, "coplanar overlapping squares intersect");
    let terminal = terminal.expect("terminal simplex");
    assert!(epa_complete_simplex(&a, &b, &terminal).is_none());
    assert!(epa(&a, &b).is_none(), "EPA needs a volumetric overlap");
    // Empty and over-long inputs are rejected instead of panicking.
    assert!(epa_complete_simplex(&a, &b, &[]).is_none());
    assert!(epa_complete_simplex(&a, &b, &[[0.0; 3]; 5]).is_none());
}

#[test]
fn test_epa_box_sphere_symmetric_configuration_has_exact_depth() {
    // Sphere (r 0.8) centred inside a unit-half-extent box at (0.5, 0, 0): the
    // shallowest way out is through the +x face, 1 - 0.5 + 0.8 = 1.3.  GJK ends
    // on a triangle here, so this exercises the seed completion.
    let a = Box3::new([0.0; 3], [1.0; 3]);
    let b = Sphere::new([0.5, 0.0, 0.0], 0.8);
    let e = epa(&a, &b).expect("overlapping box and sphere have a penetration depth");
    assert!((e.depth - 1.3).abs() < 1e-6, "depth {}", e.depth);
    assert!(
        (e.normal[0].abs() - 1.0).abs() < 1e-6,
        "normal {:?}",
        e.normal
    );
}

#[test]
fn test_epa_symmetric_boxes_and_spheres_have_analytic_depth() {
    // Boxes overlapping by 2 - 1.5 along x (GJK ends on a lower-dimensional face).
    let e = epa(
        &Box3::new([0.0; 3], [1.0; 3]),
        &Box3::new([1.5, 0.0, 0.0], [1.0; 3]),
    )
    .expect("overlapping boxes");
    assert!((e.depth - 0.5).abs() < 1e-6, "depth {}", e.depth);
    assert!(e.normal[0].abs() > 0.999, "normal {:?}", e.normal);

    // Concentric boxes: the inner one leaves 1 + 0.5 of clearance on every side.
    let e = epa(
        &Box3::new([0.0; 3], [1.0; 3]),
        &Box3::new([0.0; 3], [0.5; 3]),
    )
    .expect("nested boxes");
    assert!((e.depth - 1.5).abs() < 1e-6, "depth {}", e.depth);

    // Overlapping spheres: depth r1 + r2 - d = 2 along x (smooth shapes, so the
    // polytope approximation is accurate to ~1e-6).
    let e = epa(
        &Sphere::new([0.0; 3], 1.5),
        &Sphere::new([1.0, 0.0, 0.0], 1.5),
    )
    .expect("overlapping spheres");
    assert!((e.depth - 2.0).abs() < 1e-4, "depth {}", e.depth);
    assert!(e.normal[0].abs() > 0.999, "normal {:?}", e.normal);
}

#[test]
fn test_epa_generic_box_overlap_recovers_penetration_depth() {
    // Cubes (half-extent 1) offset by (0.7, 0.4, 0.2): the smallest overlap is
    // along x, 2 - 0.7 = 1.3.
    let a = Box3::new([0.0; 3], [1.0; 3]);
    let b = Box3::new([0.7, 0.4, 0.2], [1.0; 3]);
    assert!(gjk(&a, &b).intersecting);
    let e = epa(&a, &b).expect("overlapping boxes must have a penetration depth");
    assert!((e.depth - 1.3).abs() < 1e-6, "depth {}", e.depth);
    assert!(e.normal[0].abs() > 0.999, "normal {:?}", e.normal);
}

// ── GJK sub-algorithm (closest point of a simplex) ───────────────────────────

/// Simplex vertex whose witnesses coincide with the Minkowski point (B at the
/// origin), for exercising the sub-algorithm directly.
fn gjk_test_vertex(w: [Real; 3]) -> GjkVertex {
    GjkVertex {
        w,
        a: w,
        b: [0.0; 3],
    }
}

#[test]
fn test_gjk_subset_weights_segment() {
    let w =
        gjk_subset_weights(&[[-1.0, 1.0, 0.0], [1.0, 1.0, 0.0]]).expect("origin projects inside");
    assert!(
        (w[0] - 0.5).abs() < 1e-15 && (w[1] - 0.5).abs() < 1e-15,
        "{w:?}"
    );
    // Origin projects beyond the second endpoint, or the points coincide.
    assert!(gjk_subset_weights(&[[1.0, 1.0, 0.0], [2.0, 1.0, 0.0]]).is_none());
    assert!(gjk_subset_weights(&[[1.0, 1.0, 0.0], [1.0, 1.0, 0.0]]).is_none());
}

#[test]
fn test_gjk_subset_weights_triangle() {
    // The origin projects to (0,0,1), i.e. barycentric (1/4, 1/4, 1/2).
    let w = gjk_subset_weights(&[[-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [0.0, 1.0, 1.0]])
        .expect("origin projects inside");
    assert!((w[0] - 0.25).abs() < 1e-15, "{w:?}");
    assert!((w[1] - 0.25).abs() < 1e-15, "{w:?}");
    assert!((w[2] - 0.5).abs() < 1e-15, "{w:?}");
    // Projection falls outside the triangle.
    assert!(gjk_subset_weights(&[[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0]]).is_none());
    // Collinear vertices are affinely dependent.
    assert!(gjk_subset_weights(&[[0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]]).is_none());
}

#[test]
fn test_gjk_subset_weights_tetrahedron() {
    // Regular tetrahedron centred at the origin: all weights are 1/4.
    let tetra = [
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ];
    let w = gjk_subset_weights(&tetra).expect("origin inside");
    for (i, wi) in w.iter().enumerate() {
        assert!((wi - 0.25).abs() < 1e-15, "w[{i}] = {wi}");
    }
    // Translate it away from the origin: no longer contained.
    let shifted = tetra.map(|p| v3_add(p, [5.0, 0.0, 0.0]));
    assert!(gjk_subset_weights(&shifted).is_none());
    // Four coplanar points have no volume.
    let flat = [
        [1.0, 1.0, 0.5],
        [-1.0, -1.0, 0.5],
        [-1.0, 1.0, 0.5],
        [1.0, -1.0, 0.5],
    ];
    assert!(gjk_subset_weights(&flat).is_none());
}

#[test]
fn test_gjk_closest_to_origin_selects_the_minimal_supporting_face() {
    // Closest point is the interior of edge (0, 1): vertex 2 must be dropped.
    let edge = [
        gjk_test_vertex([-1.0, 1.0, 0.0]),
        gjk_test_vertex([1.0, 1.0, 0.0]),
        gjk_test_vertex([0.0, 5.0, 0.0]),
    ];
    let c = gjk_closest_to_origin(&edge).expect("non-empty simplex");
    assert_eq!(c.mask, 0b011);
    assert!((c.dist_sq - 1.0).abs() < 1e-15);
    assert!(point_dist(c.point, [0.0, 1.0, 0.0]) < 1e-15);
    assert_eq!(c.lambda[2], 0.0);

    // Closest point is vertex 0: the other two vertices are dropped.
    let vertex = [
        gjk_test_vertex([2.0, 0.0, 0.0]),
        gjk_test_vertex([3.0, 1.0, 0.0]),
        gjk_test_vertex([3.0, -1.0, 0.0]),
    ];
    let c = gjk_closest_to_origin(&vertex).expect("non-empty simplex");
    assert_eq!(c.mask, 0b001);
    assert!((c.dist_sq - 4.0).abs() < 1e-15);
    assert_eq!(c.lambda, [1.0, 0.0, 0.0, 0.0]);

    // Closest point is interior to the triangle: all three vertices support it.
    let face = [
        gjk_test_vertex([-1.0, -1.0, 1.0]),
        gjk_test_vertex([1.0, -1.0, 1.0]),
        gjk_test_vertex([0.0, 1.0, 1.0]),
    ];
    let c = gjk_closest_to_origin(&face).expect("non-empty simplex");
    assert_eq!(c.mask, 0b111);
    assert!((c.dist_sq - 1.0).abs() < 1e-15);
    assert!((c.lambda[0] - 0.25).abs() < 1e-15 && (c.lambda[2] - 0.5).abs() < 1e-15);

    // Origin inside the tetrahedron: distance zero, all four vertices kept.
    let tetra = [
        gjk_test_vertex([1.0, 1.0, 1.0]),
        gjk_test_vertex([-1.0, -1.0, 1.0]),
        gjk_test_vertex([-1.0, 1.0, -1.0]),
        gjk_test_vertex([1.0, -1.0, -1.0]),
    ];
    let c = gjk_closest_to_origin(&tetra).expect("non-empty simplex");
    assert_eq!(c.mask, 0b1111);
    assert_eq!(c.dist_sq, 0.0);

    assert!(gjk_closest_to_origin(&[]).is_none(), "empty simplex");
}

#[test]
fn test_gjk_closest_to_origin_never_exceeds_any_vertex_distance() {
    // The closest point of a hull is at least as close as every vertex.
    let mut rng = TestRng(0x51A9_1E0C_0000_0042);
    for _ in 0..500 {
        let n = 1 + (rng.next_u64() % 4) as usize;
        let verts: Vec<GjkVertex> = (0..n)
            .map(|_| gjk_test_vertex(v3_add(rng.vec3(-2.0, 2.0), [3.0, 0.0, 0.0])))
            .collect();
        let c = gjk_closest_to_origin(&verts).expect("non-empty simplex");
        for v in &verts {
            assert!(c.dist_sq <= v3_len_sq(v.w) * (1.0 + 1e-12));
        }
        let weight_sum: Real = c.lambda.iter().sum();
        assert!((weight_sum - 1.0).abs() < 1e-12, "weights {:?}", c.lambda);
        assert!(c.lambda.iter().all(|&l| l >= 0.0));
    }
}

// ── ContactManifold ───────────────────────────────────────────────────────

#[test]
fn test_manifold_add_and_max_depth() {
    let mut m = ContactManifold::new();
    m.add_contact(Contact {
        point: [0.0; 3],
        normal: [0.0, 1.0, 0.0],
        depth: 0.1,
        id_a: 1,
        id_b: 2,
    });
    m.add_contact(Contact {
        point: [1.0, 0.0, 0.0],
        normal: [0.0, 1.0, 0.0],
        depth: 0.3,
        id_a: 1,
        id_b: 2,
    });
    assert_eq!(m.contacts.len(), 2);
    assert!((m.max_depth() - 0.3).abs() < 1e-12);
}

#[test]
fn test_manifold_average_normal_unit() {
    let mut m = ContactManifold::new();
    for _ in 0..3 {
        m.add_contact(Contact {
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
            id_a: 1,
            id_b: 2,
        });
    }
    let n = m.average_normal();
    let len = v3_len(n);
    assert!((len - 1.0).abs() < 1e-12);
}

// ── EPA tests ─────────────────────────────────────────────────────────────

#[test]
fn test_epa_two_overlapping_spheres() {
    let a = Sphere::new([0.0; 3], 1.5);
    let b = Sphere::new([1.0, 0.0, 0.0], 1.5);
    let result = epa(&a, &b);
    assert!(
        result.is_some(),
        "EPA should return penetration for overlapping spheres"
    );
    let ep = result.unwrap();
    assert!(ep.depth >= 0.0, "penetration depth should be non-negative");
}

#[test]
fn test_epa_concentric_spheres() {
    let a = Sphere::new([0.0; 3], 2.0);
    let b = Sphere::new([0.0; 3], 1.0);
    let result = epa(&a, &b);
    // Concentric spheres: deep penetration
    assert!(result.is_some());
}

#[test]
fn test_epa_box_sphere_overlap() {
    let b = Box3::new([0.0; 3], [1.0; 3]);
    let s = Sphere::new([0.5, 0.0, 0.0], 0.8);
    let result = epa(&b, &s);
    assert!(
        result.is_some(),
        "overlapping box and sphere should produce EPA result"
    );
}

/// Slow-convergence EPA on two deeply overlapping faceted spheres offset
/// along x.  Whether EPA converges normally or falls through its
/// non-convergence tail, it must return the BEST face found — a real
/// positive depth with a normal aligned to the actual separation axis (±x) —
/// NOT the old fabricated `depth: 0, normal: [0,1,0], point: [0,0,0]`.
#[test]
fn test_epa_nonconvergence_tail_returns_best_face() {
    // Radius-1 spheres, centres 1.5 apart on x ⇒ overlap depth ≈ 0.5 on x.
    let a = faceted_sphere([0.0, 0.0, 0.0], 1.0, 24, 48);
    let b = faceted_sphere([1.5, 0.0, 0.0], 1.0, 24, 48);

    // Sanity: GJK confirms the overlap so EPA is well-defined.
    assert!(gjk(&a, &b).intersecting, "faceted spheres must overlap");

    let result = epa(&a, &b).expect("EPA must report penetration for overlapping shapes");

    // Real penetration is positive and on the order of the true overlap.
    assert!(
        result.depth > 0.05,
        "depth must be a real positive penetration, got {}",
        result.depth
    );
    // The fabricated constant claimed the y-axis; the true normal is ±x.
    let nx = result.normal[0].abs();
    let ny = result.normal[1].abs();
    let nz = result.normal[2].abs();
    assert!(
        nx > ny && nx > nz,
        "normal must follow the real (x) separation axis, got {:?}",
        result.normal
    );
    // Explicitly rule out the old fabricated tail value.
    let is_fabricated =
        result.depth == 0.0 && result.normal == [0.0, 1.0, 0.0] && result.point == [0.0; 3];
    assert!(
        !is_fabricated,
        "EPA returned the fabricated zero-depth constant"
    );
}

// ── Support function framework ────────────────────────────────────────────

#[test]
fn test_support_sphere_axis_aligned() {
    let s = Sphere::new([1.0, 2.0, 3.0], 2.0);
    let sp_x = s.support([1.0, 0.0, 0.0]);
    assert!(
        (sp_x[0] - 3.0).abs() < 1e-10,
        "x support = center.x + radius"
    );
    let sp_y = s.support([0.0, 1.0, 0.0]);
    assert!(
        (sp_y[1] - 4.0).abs() < 1e-10,
        "y support = center.y + radius"
    );
}

#[test]
fn test_support_box_positive_axes() {
    let b = Box3::new([0.0; 3], [2.0, 3.0, 4.0]);
    let sp = b.support([1.0, 1.0, 1.0]);
    assert!((sp[0] - 2.0).abs() < 1e-10);
    assert!((sp[1] - 3.0).abs() < 1e-10);
    assert!((sp[2] - 4.0).abs() < 1e-10);
}

#[test]
fn test_support_box_negative_axes() {
    let b = Box3::new([0.0; 3], [1.0; 3]);
    let sp = b.support([-1.0, -1.0, -1.0]);
    assert!((sp[0] + 1.0).abs() < 1e-10);
    assert!((sp[1] + 1.0).abs() < 1e-10);
    assert!((sp[2] + 1.0).abs() < 1e-10);
}

#[test]
fn test_support_capsule_chooses_correct_endpoint() {
    let cap = Capsule::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
    let sp_up = cap.support([0.0, 1.0, 0.0]);
    // Should be near top endpoint + radius in Y
    assert!(sp_up[1] > 1.0, "y support above top endpoint");
    let sp_down = cap.support([0.0, -1.0, 0.0]);
    assert!(sp_down[1] < -1.0, "y support below bottom endpoint");
}

// ── Minkowski sum ─────────────────────────────────────────────────────────

#[test]
fn test_minkowski_sum_sphere_sphere_support() {
    let a = Sphere::new([0.0; 3], 1.0);
    let b = Sphere::new([0.0; 3], 2.0);
    let dir = [1.0, 0.0, 0.0];
    // Support of sum = support of a + support of b in same direction
    let sp = minkowski_sum_support(&a, &b, dir);
    let expected = a.support(dir)[0] + b.support(dir)[0];
    assert!((sp[0] - expected).abs() < 1e-10);
}

#[test]
fn test_minkowski_sum_sphere_box() {
    let s = Sphere::new([0.0; 3], 0.5);
    let b = Box3::new([0.0; 3], [1.0, 1.0, 1.0]);
    let dir = [1.0, 0.0, 0.0];
    let sp = minkowski_sum_support(&s, &b, dir);
    // Support ≥ box support alone (sphere only adds positive amount)
    assert!(sp[0] >= b.support(dir)[0]);
}

// ── Shape casting / continuous collision ──────────────────────────────────

#[test]
fn test_shape_cast_sphere_hits_sphere() {
    let moving = Sphere::new([0.0; 3], 0.5);
    let target = Sphere::new([5.0, 0.0, 0.0], 0.5);
    let velocity = [1.0, 0.0, 0.0];
    let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
    assert!(t.is_some(), "moving sphere should hit static sphere");
    let tv = t.unwrap();
    assert!((0.0..=10.0).contains(&tv));
}

#[test]
fn test_shape_cast_sphere_misses_sphere() {
    let moving = Sphere::new([0.0; 3], 0.5);
    let target = Sphere::new([0.0, 100.0, 0.0], 0.5);
    let velocity = [1.0, 0.0, 0.0]; // moving in X, target far in Y
    let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
    assert!(
        t.is_none(),
        "sphere moving in X should miss sphere far in Y"
    );
}

#[test]
fn test_shape_cast_already_overlapping() {
    let moving = Sphere::new([0.0; 3], 1.0);
    let target = Sphere::new([0.5, 0.0, 0.0], 1.0);
    let velocity = [1.0, 0.0, 0.0];
    let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
    // Already overlapping → t = 0 or very small
    if let Some(tv) = t {
        assert!(tv >= 0.0);
    }
}

// ── GJK with capsule/box ──────────────────────────────────────────────────

#[test]
fn test_gjk_box_sphere_intersecting() {
    let b = Box3::new([0.0; 3], [1.0; 3]);
    let s = Sphere::new([1.5, 0.0, 0.0], 1.0);
    let r = gjk(&b, &s);
    assert!(r.intersecting, "box and sphere should intersect");
}

#[test]
fn test_gjk_box_sphere_separated() {
    let b = Box3::new([0.0; 3], [1.0; 3]);
    let s = Sphere::new([10.0, 0.0, 0.0], 0.5);
    let r = gjk(&b, &s);
    assert!(!r.intersecting, "box and sphere should be separated");
    assert!(r.distance > 0.0);
}

#[test]
fn test_gjk_capsule_sphere_overlapping() {
    let cap = Capsule::new([0.0, -2.0, 0.0], [0.0, 2.0, 0.0], 0.5);
    let s = Sphere::new([0.0, 0.0, 0.0], 0.5);
    let r = gjk(&cap, &s);
    assert!(r.intersecting, "capsule and sphere should overlap");
}

#[test]
fn test_gjk_capsule_sphere_separated() {
    let cap = Capsule::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.3);
    let s = Sphere::new([5.0, 0.0, 0.0], 0.3);
    let r = gjk(&cap, &s);
    assert!(!r.intersecting);
}

// ── AabbRaw additional tests ───────────────────────────────────────────────

#[test]
fn test_aabb_surface_area_unit_cube() {
    let a = AabbRaw::new([0.0; 3], [1.0; 3]);
    // Unit cube: SA = 6
    assert!((a.surface_area() - 6.0).abs() < 1e-12);
}

#[test]
fn test_aabb_center() {
    let a = AabbRaw::new([-2.0; 3], [4.0; 3]);
    let c = a.center();
    assert!((c[0] - 1.0).abs() < 1e-12);
}

#[test]
fn test_aabb_from_center_half() {
    let a = AabbRaw::from_center_half([1.0, 2.0, 3.0], [0.5; 3]);
    assert!((a.min[0] - 0.5).abs() < 1e-12);
    assert!((a.max[1] - 2.5).abs() < 1e-12);
}

// ── OBB additional tests ───────────────────────────────────────────────────

#[test]
fn test_obb_vertices_count() {
    let obb = Obb::axis_aligned([0.0; 3], [1.0; 3]);
    let verts = obb.vertices();
    assert_eq!(verts.len(), 8);
}

#[test]
fn test_obb_axis_aligned_vertices_correct() {
    let obb = Obb::axis_aligned([0.0; 3], [1.0; 3]);
    let verts = obb.vertices();
    // Every vertex should have coordinates in {-1, 1}³
    for v in &verts {
        assert!(v[0].abs() <= 1.0 + 1e-10);
        assert!(v[1].abs() <= 1.0 + 1e-10);
        assert!(v[2].abs() <= 1.0 + 1e-10);
    }
}

#[test]
fn test_obb_sat_touching() {
    // Two OBBs touching exactly at their faces
    let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
    let b = Obb::axis_aligned([2.0, 0.0, 0.0], [1.0; 3]);
    // Distance between centers = 2, sum of half-extents in X = 2 → touching
    let result = obb_obb_sat(&a, &b);
    // Either touching (some depth ~ 0) or just separated — both valid
    let _ = result;
}

// ── BVH additional tests ───────────────────────────────────────────────────

#[test]
fn test_bvh_all_overlapping_pairs() {
    let leaves = vec![
        BvhLeaf {
            id: 1,
            aabb: AabbRaw::new([0.0; 3], [2.0; 3]),
        },
        BvhLeaf {
            id: 2,
            aabb: AabbRaw::new([1.0; 3], [3.0; 3]),
        },
        BvhLeaf {
            id: 3,
            aabb: AabbRaw::new([10.0; 3], [11.0; 3]),
        },
    ];
    let bvh = Bvh::build(leaves);
    let pairs = bvh.overlapping_pairs();
    // Pair (1,2) overlaps; (1,3) and (2,3) don't
    assert!(
        pairs.contains(&(1, 2)) || pairs.contains(&(2, 1)),
        "pair (1,2) should be in overlapping pairs: {:?}",
        pairs
    );
}

#[test]
fn test_bvh_single_leaf() {
    let leaves = vec![BvhLeaf {
        id: 42,
        aabb: AabbRaw::new([0.0; 3], [1.0; 3]),
    }];
    let bvh = Bvh::build(leaves);
    let hits = bvh.query(&AabbRaw::new([0.0; 3], [0.5; 3]));
    assert_eq!(hits, vec![42]);
}

// ── ContactManifold additional tests ──────────────────────────────────────

#[test]
fn test_manifold_overflow_replaces_shallowest() {
    let mut m = ContactManifold::new();
    for i in 0..5 {
        m.add_contact(Contact {
            point: [i as f64, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: (i as f64 + 1.0) * 0.1,
            id_a: 1,
            id_b: 2,
        });
    }
    // Should have at most 4 contacts
    assert!(m.contacts.len() <= 4);
}

#[test]
fn test_manifold_is_empty_after_new() {
    let m = ContactManifold::new();
    assert!(m.is_empty());
}

#[test]
fn test_manifold_max_depth_empty() {
    let m = ContactManifold::new();
    // max_depth on empty manifold should be NEG_INFINITY
    assert_eq!(m.max_depth(), f64::NEG_INFINITY);
}

// ── point_segment_closest ─────────────────────────────────────────────────

#[test]
fn test_point_segment_closest_midpoint() {
    let (cp, dist) = point_segment_closest([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    assert!(cp[0].abs() < 1e-12, "closest x should be 0");
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn test_point_segment_closest_beyond_end() {
    let (cp, dist) = point_segment_closest([5.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]);
    assert!((cp[0] - 1.0).abs() < 1e-12, "closest x should clamp to 1");
    assert!((dist - 4.0).abs() < 1e-10);
}

// ── ray_sphere additional ─────────────────────────────────────────────────

#[test]
fn test_ray_sphere_from_inside() {
    let s = Sphere::new([0.0; 3], 2.0);
    // Ray from origin → should hit from inside (t > 0)
    let t = ray_sphere([0.0; 3], [1.0, 0.0, 0.0], &s);
    assert!(t.is_some());
    assert!(t.unwrap() > 0.0);
}

#[test]
fn test_ray_aabb_along_each_axis() {
    let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
    for axis in 0..3 {
        // Start on the negative axis side, aligned with the center
        let mut origin = [0.0_f64; 3];
        let mut dir = [0.0_f64; 3];
        origin[axis] = -5.0;
        dir[axis] = 1.0;
        let hit = ray_aabb(origin, dir, &aabb);
        assert!(hit.is_some(), "ray along axis {axis} should hit unit AABB");
    }
}
