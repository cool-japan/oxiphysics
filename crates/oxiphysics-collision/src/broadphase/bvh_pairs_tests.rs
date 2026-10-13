//! Regression tests for [`BvhBroadphase`] pair queries.
//!
//! The hierarchy must (a) answer `find_pairs` even when it was never built,
//! (b) agree exactly with the brute-force reference, and (c) stay correct when
//! the hierarchy is stale with respect to the AABBs that are queried.

use super::types::{BruteForceBroadPhase, BvhBroadphase};
use super::{BroadPhase, BroadphaseProfiler, BroadphaseSceneGraph};
use crate::CollisionPair;
use oxiphysics_core::{Aabb, Vec3};

/// Small deterministic linear congruential generator for reproducible scenes.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_bits(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    /// Uniform sample in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.next_bits() as f64 / (1u64 << 31) as f64
    }
}

fn boxed(cx: f64, cy: f64, cz: f64, half: f64) -> Aabb {
    Aabb::new(
        Vec3::new(cx - half, cy - half, cz - half),
        Vec3::new(cx + half, cy + half, cz + half),
    )
}

/// `n` boxes with random centres in `[0, extent)^3` and random half sizes in
/// `[0.05, 0.05 + max_half)`.
fn random_boxes(n: usize, extent: f64, max_half: f64, seed: u64) -> Vec<Aabb> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|_| {
            let cx = rng.unit() * extent;
            let cy = rng.unit() * extent;
            let cz = rng.unit() * extent;
            boxed(cx, cy, cz, 0.05 + rng.unit() * max_half)
        })
        .collect()
}

/// Boxes on an integer lattice so that many of them touch exactly at faces,
/// edges and corners (the AABB test is inclusive).
fn lattice_boxes(n: usize, cells: u64, seed: u64) -> Vec<Aabb> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|_| {
            let cx = (rng.next_bits() % cells) as f64;
            let cy = (rng.next_bits() % cells) as f64;
            let cz = (rng.next_bits() % cells) as f64;
            boxed(cx, cy, cz, 0.5)
        })
        .collect()
}

#[test]
fn default_bvh_find_pairs_reports_overlap() {
    let aabbs = vec![boxed(0.0, 0.0, 0.0, 1.0), boxed(1.5, 0.0, 0.0, 1.0)];
    let pairs = BvhBroadphase::default().find_pairs(&aabbs);
    assert_eq!(pairs, vec![CollisionPair::new(0, 1)]);
}

#[test]
fn unbuilt_bvh_find_pairs_matches_brute_force() {
    let aabbs = random_boxes(120, 12.0, 1.0, 7);
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    assert!(!expected.is_empty(), "scene must contain overlaps");
    let pairs = BvhBroadphase::new().find_pairs(&aabbs);
    assert_eq!(pairs, expected);
}

#[test]
fn built_bvh_matches_brute_force_on_500_random_boxes() {
    let aabbs = random_boxes(500, 25.0, 1.2, 11);
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    assert!(expected.len() > 20, "got only {} pairs", expected.len());
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    assert!(bvh.validate());
    assert_eq!(bvh.find_pairs(&aabbs), expected);
}

#[test]
fn built_bvh_matches_brute_force_on_2000_random_boxes() {
    let aabbs = random_boxes(2000, 40.0, 1.0, 23);
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    assert!(expected.len() > 100, "got only {} pairs", expected.len());
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    assert_eq!(bvh.find_all_pairs(&aabbs), expected);
}

#[test]
fn touching_lattice_boxes_match_brute_force() {
    let aabbs = lattice_boxes(300, 6, 5);
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    assert_eq!(bvh.find_pairs(&aabbs), expected);
}

#[test]
fn identical_boxes_report_every_pair_once() {
    let aabbs: Vec<Aabb> = (0..40).map(|_| boxed(1.0, 2.0, 3.0, 0.5)).collect();
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    let pairs = bvh.find_pairs(&aabbs);
    assert_eq!(pairs.len(), 40 * 39 / 2);
    assert_eq!(pairs, BruteForceBroadPhase.find_pairs(&aabbs));
}

#[test]
fn pairs_are_canonical_and_sorted() {
    let aabbs = random_boxes(300, 15.0, 1.5, 31);
    let mut bvh = BvhBroadphase::new();
    bvh.build(&aabbs);
    let pairs = bvh.find_pairs(&aabbs);
    assert!(!pairs.is_empty());
    assert!(pairs.iter().all(|p| p.a < p.b));
    assert!(
        pairs
            .windows(2)
            .all(|w| (w[0].a, w[0].b) < (w[1].a, w[1].b))
    );
}

#[test]
fn degenerate_inputs_have_no_pairs() {
    let mut bvh = BvhBroadphase::new();
    assert!(bvh.find_pairs(&[]).is_empty());
    let single = vec![boxed(0.0, 0.0, 0.0, 1.0)];
    assert!(bvh.find_pairs(&single).is_empty());
    bvh.build(&single);
    assert!(bvh.find_pairs(&single).is_empty());
}

#[test]
fn stale_hierarchy_is_not_trusted() {
    // Tree built for the old positions; the objects have since moved, so some
    // AABBs are no longer enclosed by their leaves.
    let before = random_boxes(150, 20.0, 1.0, 3);
    let after = random_boxes(150, 20.0, 1.0, 4);
    let mut bvh = BvhBroadphase::new();
    bvh.build(&before);
    assert_eq!(
        bvh.find_pairs(&after),
        BruteForceBroadPhase.find_pairs(&after)
    );
}

#[test]
fn hierarchy_built_for_other_object_count_does_not_panic() {
    let big = random_boxes(60, 10.0, 1.0, 9);
    let small: Vec<Aabb> = big[..7].to_vec();
    let mut bvh = BvhBroadphase::new();
    bvh.build(&big);
    assert_eq!(
        bvh.find_pairs(&small),
        BruteForceBroadPhase.find_pairs(&small)
    );
    let bigger = random_boxes(90, 10.0, 1.0, 10);
    assert_eq!(
        bvh.find_pairs(&bigger),
        BruteForceBroadPhase.find_pairs(&bigger)
    );
}

#[test]
fn fat_hierarchy_answers_with_tight_boxes() {
    // A hierarchy built from fattened AABBs may be queried with the tight ones:
    // only the tight overlaps are reported.
    let tight = random_boxes(200, 18.0, 0.8, 17);
    let fat: Vec<Aabb> = tight.iter().map(|a| a.expand(0.4)).collect();
    let mut bvh = BvhBroadphase::new();
    bvh.build(&fat);
    let tight_pairs = BruteForceBroadPhase.find_pairs(&tight);
    let fat_pairs = BruteForceBroadPhase.find_pairs(&fat);
    assert!(fat_pairs.len() > tight_pairs.len());
    assert_eq!(bvh.find_pairs(&tight), tight_pairs);
    assert_eq!(bvh.find_pairs(&fat), fat_pairs);
}

#[test]
fn profiler_and_scene_graph_use_the_pruned_query() {
    let aabbs = random_boxes(80, 10.0, 1.0, 41);
    let expected = BruteForceBroadPhase.find_pairs(&aabbs);
    let mut profiler = BroadphaseProfiler::new();
    assert_eq!(profiler.profile_bvh(&aabbs), expected);

    let mut graph = BroadphaseSceneGraph::new();
    for aabb in &aabbs {
        graph.add_node(aabb.clone(), None);
    }
    assert_eq!(graph.find_pairs_bvh(), expected);
}
