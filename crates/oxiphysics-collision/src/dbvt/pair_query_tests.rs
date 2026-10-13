//! Regression tests for [`DynamicBvh::query_overlapping_pairs`].
//!
//! The pruned self-collision traversal must report exactly the leaf pairs a
//! brute-force scan over the (fat) leaf AABBs reports, for trees that have been
//! shaped by arbitrary insert / remove / update sequences.

use super::functions::NodeIdx;
use super::types::{BvhAabb, DynamicBvh};
use oxiphysics_core::Vec3;

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

    /// Uniform index in `0..n` (`n > 0`).
    fn below(&mut self, n: usize) -> usize {
        (self.next_bits() % n as u64) as usize
    }
}

fn random_aabb(rng: &mut Lcg, extent: f64) -> BvhAabb {
    let c = Vec3::new(
        rng.unit() * extent,
        rng.unit() * extent,
        rng.unit() * extent,
    );
    let h = 0.2 + rng.unit();
    BvhAabb::new(c - Vec3::new(h, h, h), c + Vec3::new(h, h, h))
}

/// Reference answer: test every pair of live leaves against each other.
fn brute_force_pairs(bvh: &DynamicBvh) -> Vec<(u32, u32)> {
    let leaves: Vec<(u32, BvhAabb)> = bvh
        .nodes
        .iter()
        .filter_map(|node| node.data.map(|data| (data, node.aabb)))
        .collect();
    let mut pairs = Vec::new();
    for (i, (data_a, aabb_a)) in leaves.iter().enumerate() {
        for (data_b, aabb_b) in &leaves[(i + 1)..] {
            if data_a != data_b && aabb_a.intersects(aabb_b) {
                pairs.push((*data_a.min(data_b), *data_a.max(data_b)));
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

#[test]
fn empty_and_single_leaf_trees_have_no_pairs() {
    let mut bvh = DynamicBvh::new();
    assert!(bvh.query_overlapping_pairs().is_empty());
    bvh.insert(random_aabb(&mut Lcg::new(1), 5.0), 7);
    assert!(bvh.query_overlapping_pairs().is_empty());
}

#[test]
fn matches_brute_force_on_1000_leaves_after_mixed_edits() {
    let mut rng = Lcg::new(2024);
    let extent = 30.0;
    let mut bvh = DynamicBvh::new();
    let mut live: Vec<(NodeIdx, u32)> = Vec::new();
    let mut next_data = 0u32;

    for _ in 0..1000 {
        let leaf = bvh.insert(random_aabb(&mut rng, extent), next_data);
        live.push((leaf, next_data));
        next_data += 1;
    }
    assert!(bvh.validate());

    // Remove roughly a third of the leaves.
    for _ in 0..330 {
        let victim = rng.below(live.len());
        let (leaf, _) = live.swap_remove(victim);
        bvh.remove(leaf);
    }
    assert!(bvh.validate());

    // Move some leaves far enough to force a re-insertion.
    let mut moved = 0;
    for _ in 0..250 {
        let which = rng.below(live.len());
        let (leaf, data) = live[which];
        let target = random_aabb(&mut rng, extent);
        if bvh.update(leaf, target) {
            moved += 1;
            let new_leaf = bvh
                .nodes
                .iter()
                .position(|node| node.data == Some(data))
                .expect("re-inserted leaf must exist");
            live[which] = (new_leaf, data);
        }
    }
    assert!(moved > 100, "only {moved} leaves were re-inserted");

    // A few more inserts on top of the edited tree.
    for _ in 0..100 {
        let leaf = bvh.insert(random_aabb(&mut rng, extent), next_data);
        live.push((leaf, next_data));
        next_data += 1;
    }
    assert!(bvh.validate());
    assert_eq!(bvh.n_leaves(), live.len());

    let expected = brute_force_pairs(&bvh);
    assert!(expected.len() > 50, "got only {} pairs", expected.len());
    assert_eq!(bvh.query_overlapping_pairs(), expected);
}

#[test]
fn result_is_sorted_canonical_and_duplicate_free() {
    let mut rng = Lcg::new(77);
    let mut bvh = DynamicBvh::new();
    for data in 0..200u32 {
        bvh.insert(random_aabb(&mut rng, 12.0), data);
    }
    let pairs = bvh.query_overlapping_pairs();
    assert!(pairs.len() > 50);
    assert!(pairs.iter().all(|&(a, b)| a < b));
    assert!(pairs.windows(2).all(|w| w[0] < w[1]));
}

#[test]
fn leaves_sharing_a_payload_are_never_paired() {
    let mut bvh = DynamicBvh::new();
    let aabb = BvhAabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
    bvh.insert(aabb, 5);
    bvh.insert(aabb, 5);
    bvh.insert(aabb, 6);
    assert_eq!(bvh.query_overlapping_pairs(), vec![(5, 6)]);
}

#[test]
fn identical_leaves_report_every_pair_once() {
    let mut bvh = DynamicBvh::new();
    let aabb = BvhAabb::new(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
    for data in 0..30u32 {
        bvh.insert(aabb, data);
    }
    assert_eq!(bvh.query_overlapping_pairs().len(), 30 * 29 / 2);
}

#[test]
fn far_apart_clusters_only_pair_within_themselves() {
    let mut bvh = DynamicBvh::new();
    let mut rng = Lcg::new(5);
    // Two dense clusters 1000 units apart: no pair may ever cross them.
    for data in 0..60u32 {
        let local = random_aabb(&mut rng, 4.0);
        let offset = if data < 30 { 0.0 } else { 1000.0 };
        let shift = Vec3::new(offset, 0.0, 0.0);
        bvh.insert(BvhAabb::new(local.min + shift, local.max + shift), data);
    }
    let pairs = bvh.query_overlapping_pairs();
    assert!(!pairs.is_empty());
    assert!(pairs.iter().all(|&(a, b)| (a < 30) == (b < 30)));
    assert_eq!(pairs, brute_force_pairs(&bvh));
}

#[test]
fn pairs_with_stats_reports_the_pair_count() {
    let mut rng = Lcg::new(13);
    let mut bvh = DynamicBvh::new();
    for data in 0..150u32 {
        bvh.insert(random_aabb(&mut rng, 10.0), data);
    }
    let (pairs, stats) = bvh.query_pairs_with_stats();
    assert_eq!(stats.overlap_pairs, pairs.len());
    assert_eq!(stats.leaf_count, 150);
    assert_eq!(pairs, brute_force_pairs(&bvh));
}
