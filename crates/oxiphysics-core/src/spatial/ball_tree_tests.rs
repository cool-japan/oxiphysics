//! Regression and property tests for [`BallTree`].
//!
//! The defect behind these tests: nodes were pruned with a ball centred on the
//! pivot point, while the stored radius had been measured from the centroid of
//! the subtree.  Such a ball does not bound the subtree, so nearest-neighbour
//! and range queries on skewed point clouds silently missed points, and leaves
//! reported their pivot twice because it was also listed among their items.

use super::BallTree;
use crate::math::Vec3;
use proptest::prelude::*;
use std::f64::consts::TAU;

/// SplitMix64: a tiny deterministic generator for reproducible test clouds.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform sample in the open interval (0, 1).
    fn next_open_unit(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    /// Standard normal sample (Box-Muller).
    fn next_gaussian(&mut self) -> f64 {
        let radius = (-2.0 * self.next_open_unit().ln()).sqrt();
        radius * (TAU * self.next_open_unit()).cos()
    }

    /// A point whose coordinates are cubed standard normals: heavily skewed, so
    /// the centroid of a subtree sits far from its median (pivot) point.
    fn next_cubed_gaussian_point(&mut self) -> Vec3 {
        let mut cube = || {
            let g = self.next_gaussian();
            g * g * g
        };
        Vec3::new(cube(), cube(), cube())
    }
}

fn cubed_gaussian_cloud(rng: &mut SplitMix64, n: usize) -> Vec<Vec3> {
    (0..n).map(|_| rng.next_cubed_gaussian_point()).collect()
}

fn brute_force_nearest_distance(points: &[Vec3], query: Vec3) -> f64 {
    points
        .iter()
        .map(|p| (p - query).norm_squared())
        .fold(f64::INFINITY, f64::min)
        .sqrt()
}

fn brute_force_range(points: &[Vec3], query: Vec3, radius: f64) -> Vec<usize> {
    let r2 = radius * radius;
    points
        .iter()
        .enumerate()
        .filter(|(_, p)| (*p - query).norm_squared() <= r2)
        .map(|(i, _)| i)
        .collect()
}

fn sorted(mut ids: Vec<usize>) -> Vec<usize> {
    ids.sort_unstable();
    ids
}

fn has_duplicates(sorted_ids: &[usize]) -> bool {
    sorted_ids.windows(2).any(|pair| pair[0] == pair[1])
}

/// Collects the indices of every point stored in the subtree below `node_idx`.
fn collect_subtree(tree: &BallTree, node_idx: usize, out: &mut Vec<usize>) {
    let node = &tree.nodes[node_idx];
    out.extend(node.pivot);
    out.extend_from_slice(&node.items);
    if let Some(left) = node.left {
        collect_subtree(tree, left, out);
    }
    if let Some(right) = node.right {
        collect_subtree(tree, right, out);
    }
}

#[test]
fn every_bounding_ball_contains_its_whole_subtree() {
    let mut rng = SplitMix64::new(0x5EED_0001);
    let points = cubed_gaussian_cloud(&mut rng, 400);
    let tree = BallTree::build(points.clone());
    for node_idx in 0..tree.nodes.len() {
        let mut members = Vec::new();
        collect_subtree(&tree, node_idx, &mut members);
        let node = &tree.nodes[node_idx];
        let centroid = members
            .iter()
            .fold(Vec3::zeros(), |acc, &m| acc + points[m])
            / members.len() as f64;
        assert!(
            (centroid - node.center).norm() <= 1e-12 * (1.0 + centroid.norm()),
            "node {node_idx}: centre is not the subtree centroid"
        );
        let farthest = members
            .iter()
            .map(|&m| (points[m] - node.center).norm())
            .fold(0.0_f64, f64::max);
        assert!(
            farthest <= node.radius * (1.0 + 1e-12),
            "node {node_idx}: a subtree point lies at {farthest}, outside radius {}",
            node.radius
        );
        assert!(
            (farthest - node.radius).abs() <= 1e-12 * (1.0 + node.radius),
            "node {node_idx}: radius {} is not tight (farthest member at {farthest})",
            node.radius
        );
    }
}

#[test]
fn every_point_is_stored_exactly_once() {
    for n in [1usize, 2, 7, 8, 9, 10, 16, 17, 33, 100, 257] {
        let mut rng = SplitMix64::new(0x5EED_0002 + n as u64);
        let tree = BallTree::build(cubed_gaussian_cloud(&mut rng, n));
        let mut stored = Vec::new();
        collect_subtree(&tree, 0, &mut stored);
        assert_eq!(sorted(stored), (0..n).collect::<Vec<_>>(), "n = {n}");
    }
}

#[test]
fn leaf_pivot_is_not_reported_twice_by_range_queries() {
    for n in [1usize, 2, 5, 8, 9, 17, 100] {
        let mut rng = SplitMix64::new(0x5EED_0003 + n as u64);
        let tree = BallTree::build(cubed_gaussian_cloud(&mut rng, n));
        let hits = sorted(tree.range_query(Vec3::zeros(), 1.0e6));
        assert_eq!(hits, (0..n).collect::<Vec<_>>(), "n = {n}");
    }
}

#[test]
fn nearest_matches_brute_force_on_skewed_cloud() {
    let mut rng = SplitMix64::new(0x5EED_0004);
    let points = cubed_gaussian_cloud(&mut rng, 300);
    let tree = BallTree::build(points.clone());
    for query_no in 0..200 {
        let query = rng.next_cubed_gaussian_point();
        let (idx, dist) = tree
            .nearest(query)
            .expect("a non-empty tree always has a nearest point");
        let expected = brute_force_nearest_distance(&points, query);
        assert_eq!(dist, expected, "query {query_no}: wrong nearest distance");
        let actual = (points[idx] - query).norm();
        assert!(
            (actual - dist).abs() <= 1e-12 * (1.0 + dist),
            "query {query_no}: returned index {idx} is at {actual}, not {dist}"
        );
    }
}

#[test]
fn range_query_matches_brute_force_on_skewed_cloud() {
    let mut rng = SplitMix64::new(0x5EED_0005);
    let points = cubed_gaussian_cloud(&mut rng, 300);
    let tree = BallTree::build(points.clone());
    for query_no in 0..200 {
        let query = rng.next_cubed_gaussian_point();
        for radius in [0.25, 1.0, 3.0, 10.0] {
            let hits = sorted(tree.range_query(query, radius));
            assert!(
                !has_duplicates(&hits),
                "query {query_no}, radius {radius}: duplicate indices in {hits:?}"
            );
            assert_eq!(
                hits,
                brute_force_range(&points, query, radius),
                "query {query_no}, radius {radius}"
            );
        }
    }
}

#[test]
fn queries_from_stored_points_find_themselves() {
    let mut rng = SplitMix64::new(0x5EED_0006);
    let points = cubed_gaussian_cloud(&mut rng, 120);
    let tree = BallTree::build(points.clone());
    for (i, &p) in points.iter().enumerate() {
        let (_, dist) = tree.nearest(p).expect("tree is not empty");
        assert_eq!(dist, 0.0, "point {i} must be its own nearest neighbour");
        assert!(tree.range_query(p, 0.0).contains(&i), "point {i}");
    }
}

#[test]
fn range_query_includes_points_exactly_on_the_boundary() {
    let mut points = Vec::new();
    for x in 0..6 {
        for y in 0..6 {
            for z in 0..6 {
                points.push(Vec3::new(f64::from(x), f64::from(y), f64::from(z)));
            }
        }
    }
    let tree = BallTree::build(points.clone());
    let queries = [
        Vec3::new(2.0, 2.0, 2.0),
        Vec3::new(0.5, 0.5, 0.5),
        Vec3::new(5.0, 0.0, 3.0),
        Vec3::new(-1.0, 7.0, 2.5),
    ];
    for query in queries {
        for radius in [0.0, 1.0, 2.0, std::f64::consts::SQRT_2, 2.5, 3.0] {
            assert_eq!(
                sorted(tree.range_query(query, radius)),
                brute_force_range(&points, query, radius),
                "query {query:?}, radius {radius}"
            );
        }
    }
}

#[test]
fn coincident_points_are_each_reported_once() {
    let points = vec![Vec3::new(1.0, -2.0, 3.0); 40];
    let tree = BallTree::build(points);
    let (_, dist) = tree
        .nearest(Vec3::new(1.0, -2.0, 3.0))
        .expect("tree is not empty");
    assert_eq!(dist, 0.0);
    assert_eq!(
        sorted(tree.range_query(Vec3::new(1.0, -2.0, 3.0), 0.0)),
        (0..40).collect::<Vec<_>>()
    );
}

#[test]
fn degenerate_radii_match_nothing() {
    let points: Vec<Vec3> = (0..20).map(|i| Vec3::new(f64::from(i), 0.0, 0.0)).collect();
    let tree = BallTree::build(points);
    assert!(tree.range_query(Vec3::zeros(), -1.0).is_empty());
    assert!(tree.range_query(Vec3::zeros(), f64::NAN).is_empty());
    assert!(
        BallTree::build(Vec::new())
            .range_query(Vec3::zeros(), 1.0)
            .is_empty()
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Nearest distances and range results agree with brute force over random
    /// cubed-Gaussian clouds, and no range query lists an index twice.
    #[test]
    fn prop_ball_tree_matches_brute_force(
        seed in any::<u64>(),
        n in 1usize..500,
        radius in 0.0f64..40.0,
    ) {
        let mut rng = SplitMix64::new(seed);
        let points = cubed_gaussian_cloud(&mut rng, n);
        let tree = BallTree::build(points.clone());
        for _ in 0..12 {
            let query = rng.next_cubed_gaussian_point();
            let Some((idx, dist)) = tree.nearest(query) else {
                return Err(TestCaseError::fail("non-empty tree returned no neighbour"));
            };
            prop_assert_eq!(dist, brute_force_nearest_distance(&points, query));
            prop_assert!(idx < n);
            let hits = sorted(tree.range_query(query, radius));
            prop_assert!(!has_duplicates(&hits), "duplicate indices in {hits:?}");
            prop_assert_eq!(hits, brute_force_range(&points, query, radius));
        }
    }
}
