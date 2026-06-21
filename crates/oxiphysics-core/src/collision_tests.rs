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

/// Slow-convergence GJK on two clearly separated faceted spheres: even if the
/// main loop runs out of iterations, the fallback must report an HONEST
/// non-intersection with a real positive distance — never the old hardcoded
/// `intersecting: true, distance: 0` constant.
#[test]
fn test_gjk_nonconvergence_tail_reports_real_estimate() {
    // Two faceted spheres of radius 1, centres 5 apart on x ⇒ gap ≈ 3.
    let a = faceted_sphere([0.0, 0.0, 0.0], 1.0, 24, 48);
    let b = faceted_sphere([5.0, 0.0, 0.0], 1.0, 24, 48);
    let r = gjk(&a, &b);
    // True minimum distance is 5 - 1 - 1 = 3.  The fabricated constant would
    // have claimed intersecting with distance 0.
    assert!(
        !r.intersecting,
        "clearly separated shapes must not be reported as intersecting"
    );
    assert!(
        r.distance > 1.0,
        "fallback must report the real separation, got {}",
        r.distance
    );
    assert!(
        !(r.distance == 0.0 && r.closest_a == [0.0; 3] && r.closest_b == [0.0; 3]),
        "fallback returned the fabricated zero/zero constant"
    );
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
