// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the persistent contact manifolds: manifolds of pairs
//! that separated must expire, and a contact that is detected again at the
//! place of a stored point must refresh that point instead of piling up
//! duplicates.

use super::{
    AabbGpu, CollisionGpuPipeline, ContactResult, DEFAULT_MANIFOLD_MERGE_RADIUS,
    GpuCollisionPipeline, ManifoldPoint, PersistentManifoldGpu,
};

/// Contact of the body pair (0, 1) at `(x, 0, 0)` with the given depth.
fn manifold_contact(x: f32, depth: f32) -> ContactResult {
    ContactResult {
        prim_a: 0,
        prim_b: 1,
        contact_point: [x, 0.0, 0.0],
        normal: [0.0, 1.0, 0.0],
        depth,
    }
}
/// Two boxes that overlap by 0.5 along x.
fn touching_pair() -> Vec<AabbGpu> {
    vec![
        AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
    ]
}
/// Two boxes 4 units apart.
fn separated_pair() -> Vec<AabbGpu> {
    vec![
        AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        AabbGpu::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
    ]
}
/// A chain: box 1 touches box 0 and box 2, boxes 0 and 2 do not touch.
fn touching_chain() -> Vec<AabbGpu> {
    vec![
        AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        AabbGpu::new([1.2, 0.0, 0.0], [2.2, 1.0, 1.0]),
    ]
}
fn body_key(m: &PersistentManifoldGpu) -> (u32, u32) {
    (m.body_a.min(m.body_b), m.body_a.max(m.body_b))
}
#[test]
fn test_manifold_merges_a_redetected_contact_and_keeps_its_warm_impulses() {
    let mut m = PersistentManifoldGpu::new(0, 1);
    let mut first = ManifoldPoint::from_contact(&manifold_contact(0.0, 0.10));
    first.warm_impulse_normal = 3.0;
    first.warm_impulse_t1 = 1.0;
    first.warm_impulse_t2 = 2.0;
    m.add_contact(first);
    // The same contact one frame later: 0.01 away, a little deeper.
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.01, 0.12)));
    assert_eq!(m.num_points, 1, "within the merge radius: no duplicate");
    let kept = m.active_points().next().expect("one point");
    assert!(
        (kept.position[0] - 0.01).abs() < 1.0e-7,
        "geometry refreshed"
    );
    assert!((kept.depth - 0.12).abs() < 1.0e-7);
    assert_eq!(
        (
            kept.warm_impulse_normal,
            kept.warm_impulse_t1,
            kept.warm_impulse_t2
        ),
        (3.0, 1.0, 2.0),
        "the replaced point's warm impulses must survive"
    );
}
#[test]
fn test_manifold_incoming_warm_start_wins_over_the_stored_impulses() {
    let mut m = PersistentManifoldGpu::new(0, 1);
    let mut stored = ManifoldPoint::from_contact(&manifold_contact(0.0, 0.1));
    stored.warm_impulse_normal = 3.0;
    stored.warm_impulse_t1 = 1.0;
    m.add_contact(stored);
    let mut incoming = ManifoldPoint::from_contact(&manifold_contact(0.005, 0.1));
    incoming.warm_impulse_normal = 7.0;
    m.add_contact(incoming);
    assert_eq!(m.num_points, 1);
    let kept = m.active_points().next().expect("one point");
    assert_eq!(
        (
            kept.warm_impulse_normal,
            kept.warm_impulse_t1,
            kept.warm_impulse_t2
        ),
        (7.0, 0.0, 0.0)
    );
}
#[test]
fn test_manifold_merge_radius_decides_between_merge_and_append() {
    let mut default_radius = PersistentManifoldGpu::new(0, 1);
    assert_eq!(default_radius.merge_radius, DEFAULT_MANIFOLD_MERGE_RADIUS);
    assert_eq!(DEFAULT_MANIFOLD_MERGE_RADIUS, 0.02);
    default_radius.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.1)));
    default_radius.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.05, 0.1)));
    assert_eq!(default_radius.num_points, 2, "0.05 is beyond the radius");

    let mut wide = PersistentManifoldGpu::new(0, 1).with_merge_radius(0.1);
    wide.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.1)));
    wide.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.05, 0.1)));
    assert_eq!(wide.num_points, 1, "0.05 is within a radius of 0.1");

    for radius in [0.0_f32, -1.0] {
        let mut exact = PersistentManifoldGpu::new(0, 1).with_merge_radius(radius);
        exact.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.1)));
        exact.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.2)));
        assert_eq!(
            exact.num_points, 1,
            "identical positions merge (radius {radius})"
        );
        exact.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.001, 0.2)));
        assert_eq!(exact.num_points, 2, "radius {radius} merges nothing else");
    }
}
#[test]
fn test_manifold_merges_with_the_nearest_of_two_candidates() {
    let mut m = PersistentManifoldGpu::new(0, 1).with_merge_radius(0.5);
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.1)));
    // 0.6 is farther than the radius from the first point: a second point.
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.6, 0.1)));
    assert_eq!(m.num_points, 2);
    // 0.4 is within the radius of both; it replaces the nearer one (0.6).
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.4, 0.3)));
    assert_eq!(m.num_points, 2);
    let mut xs: Vec<f32> = m.active_points().map(|p| p.position[0]).collect();
    xs.sort_by(f32::total_cmp);
    assert_eq!(xs, vec![0.0, 0.4]);
}
#[test]
fn test_manifold_add_after_remove_stale_never_overwrites_a_live_point() {
    let mut m = PersistentManifoldGpu::new(0, 1);
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(0.0, 0.5)));
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(1.0, 0.01)));
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(2.0, 0.5)));
    // Dropping the shallow middle point leaves a hole in the slot array.
    m.remove_stale(0.05);
    assert_eq!(m.num_points, 2);
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(5.0, 0.4)));
    assert_eq!(m.num_points, 3);
    let mut xs: Vec<f32> = m.active_points().map(|p| p.position[0]).collect();
    xs.sort_by(f32::total_cmp);
    assert_eq!(xs, vec![0.0, 2.0, 5.0], "no live point may be overwritten");
}
#[test]
fn test_manifold_full_replaces_the_shallowest_point_only_if_the_new_one_is_deeper() {
    let mut m = PersistentManifoldGpu::new(0, 1);
    for (i, depth) in [0.1_f32, 0.2, 0.3, 0.4].into_iter().enumerate() {
        m.add_contact(ManifoldPoint::from_contact(&manifold_contact(
            i as f32, depth,
        )));
    }
    assert_eq!(m.num_points, 4);
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(10.0, 0.05)));
    assert!(
        m.active_points().all(|p| p.position[0] != 10.0),
        "a shallower point must not displace a deeper one"
    );
    m.add_contact(ManifoldPoint::from_contact(&manifold_contact(11.0, 0.35)));
    assert_eq!(m.num_points, 4);
    let depths: Vec<f32> = m.active_points().map(|p| p.depth).collect();
    assert!(!depths.contains(&0.1), "the shallowest point was replaced");
    assert!(m.active_points().any(|p| p.position[0] == 11.0));
}
#[test]
fn test_pipeline_redetected_contact_does_not_pile_up_duplicate_points() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    let aabbs = touching_pair();
    for frame in 0..5 {
        pipeline.dispatch(&aabbs);
        assert_eq!(pipeline.manifolds.len(), 1, "frame {frame}");
        assert_eq!(pipeline.manifolds[0].num_points, 1, "frame {frame}");
        assert_eq!(pipeline.total_contact_points(), 1, "frame {frame}");
    }
}
#[test]
fn test_pipeline_manifold_expires_when_the_pair_separates() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    pipeline.dispatch(&touching_pair());
    assert_eq!(pipeline.total_contact_points(), 1);
    pipeline.dispatch(&separated_pair());
    assert!(pipeline.contacts.is_empty());
    assert!(
        pipeline.manifolds.is_empty(),
        "a separated pair must not keep its manifold"
    );
    assert_eq!(pipeline.total_contact_points(), 0, "no phantom contacts");
    // The BVH broadphase path behaves the same way.
    pipeline.dispatch_bvh(&touching_pair());
    assert_eq!(pipeline.total_contact_points(), 1);
    pipeline.dispatch_bvh(&separated_pair());
    assert!(pipeline.manifolds.is_empty());
    // A pair that touches again starts a fresh manifold with a single point.
    pipeline.dispatch(&touching_pair());
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(pipeline.total_contact_points(), 1);
}
#[test]
fn test_pipeline_only_the_separated_pair_loses_its_manifold() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    let mut aabbs = touching_chain();
    pipeline.dispatch(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 2);
    aabbs[0] = AabbGpu::new([-5.0, 0.0, 0.0], [-4.0, 1.0, 1.0]);
    pipeline.dispatch(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(body_key(&pipeline.manifolds[0]), (1, 2));
    assert_eq!(pipeline.manifolds[0].num_points, 1);
}
#[test]
fn test_pipeline_matches_manifolds_regardless_of_body_order() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    pipeline.dispatch(&touching_pair());
    for m in pipeline.manifolds.iter_mut() {
        std::mem::swap(&mut m.body_a, &mut m.body_b);
    }
    pipeline.dispatch(&touching_pair());
    assert_eq!(
        pipeline.manifolds.len(),
        1,
        "the swapped manifold is reused"
    );
    assert_eq!(pipeline.manifolds[0].num_points, 1);
    assert_eq!(body_key(&pipeline.manifolds[0]), (0, 1));
}
#[test]
fn test_pipeline_contact_that_moved_beyond_the_merge_radius_adds_a_point() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    let mut aabbs = touching_pair();
    pipeline.dispatch(&aabbs);
    // The contact point is the midpoint of the box centres: it moves by 0.05.
    aabbs[1] = AabbGpu::new([0.6, 0.0, 0.0], [1.6, 1.0, 1.0]);
    pipeline.dispatch(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(pipeline.manifolds[0].num_points, 2);
    // A sub-radius drift (0.005) refreshes the nearest point instead.
    aabbs[1] = AabbGpu::new([0.61, 0.0, 0.0], [1.61, 1.0, 1.0]);
    pipeline.dispatch(&aabbs);
    assert_eq!(pipeline.manifolds[0].num_points, 2);
}
#[test]
fn test_pipeline_redetected_contact_keeps_the_warm_impulses_of_its_point() {
    let mut pipeline = CollisionGpuPipeline::new(0.0);
    pipeline.dispatch(&touching_pair());
    // A solver stores the impulse it applied in the manifold point.
    let point = pipeline.manifolds[0].points[0]
        .as_mut()
        .expect("the first contact point");
    point.warm_impulse_normal = 3.0;
    pipeline.dispatch(&touching_pair());
    assert_eq!(pipeline.manifolds[0].num_points, 1);
    let kept = pipeline.manifolds[0].active_points().next().expect("point");
    assert_eq!(kept.warm_impulse_normal, 3.0);
}
#[test]
fn test_gpu_pipeline_redetected_contact_does_not_pile_up_duplicate_points() {
    let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
    let aabbs = touching_pair();
    for frame in 0..5 {
        pipeline.step(&aabbs);
        assert_eq!(pipeline.manifolds.len(), 1, "frame {frame}");
        assert_eq!(pipeline.manifolds[0].num_points, 1, "frame {frame}");
        assert_eq!(pipeline.total_contact_points(), 1, "frame {frame}");
    }
}
#[test]
fn test_gpu_pipeline_manifold_expires_when_the_pair_separates() {
    let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
    pipeline.step(&touching_pair());
    assert_eq!(pipeline.total_contact_points(), 1);
    pipeline.step(&separated_pair());
    assert!(pipeline.contacts.is_empty());
    assert!(
        pipeline.manifolds.is_empty(),
        "a separated pair must not keep its manifold"
    );
    assert_eq!(pipeline.total_contact_points(), 0, "no phantom contacts");
    pipeline.step(&touching_pair());
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(pipeline.total_contact_points(), 1);
}
#[test]
fn test_gpu_pipeline_only_the_separated_pair_loses_its_manifold() {
    let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
    let mut aabbs = touching_chain();
    pipeline.step(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 2);
    aabbs[0] = AabbGpu::new([-5.0, 0.0, 0.0], [-4.0, 1.0, 1.0]);
    pipeline.step(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(body_key(&pipeline.manifolds[0]), (1, 2));
    assert_eq!(pipeline.manifolds[0].num_points, 1);
}
#[test]
fn test_gpu_pipeline_solver_update_of_the_cache_reaches_the_persistent_point() {
    let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
    let aabbs = touching_pair();
    pipeline.step(&aabbs);
    assert_eq!(pipeline.manifolds[0].num_points, 1);
    pipeline
        .contact_cache
        .find_mut(0, 1)
        .expect("pair (0, 1) is in contact")
        .update(7.0_f32, 0.5_f32, 0.25_f32);
    // The manifold is kept (not cleared): the re-detected point replaces the
    // stored one and must carry the impulses the solver just cached.
    pipeline.step(&aabbs);
    assert_eq!(pipeline.manifolds.len(), 1);
    assert_eq!(pipeline.manifolds[0].num_points, 1);
    let kept = pipeline.manifolds[0].active_points().next().expect("point");
    assert_eq!(kept.warm_impulse_normal, 7.0);
    assert_eq!(kept.warm_impulse_t1, 0.5);
    assert_eq!(kept.warm_impulse_t2, 0.25);
}
