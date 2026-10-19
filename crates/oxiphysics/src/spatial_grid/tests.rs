// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for [`SpatialGrid`]: behaviour, differential tests against a
//! brute-force reference, serialisation, and checks that queries only touch the
//! entries near them (the grid used to be a flat list scanned linearly).

use super::*;
use std::collections::{BTreeMap, HashSet};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Minimal xorshift64* generator so that the randomised tests are deterministic.
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// Uniform in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform in `0..n`.
    fn below(&mut self, n: usize) -> usize {
        (self.unit() * n as f64) as usize % n
    }

    fn point(&mut self, extent: f64) -> [f64; 3] {
        [
            self.range(-extent, extent),
            self.range(-extent, extent),
            self.range(-extent, extent),
        ]
    }
}

impl SpatialGrid {
    /// Verify every internal invariant of the index.
    fn assert_invariants(&self) {
        assert_eq!(self.slot.len(), self.entries.len(), "slot map size");
        for (index, entry) in self.entries.iter().enumerate() {
            assert_eq!(
                self.slot.get(&entry.id),
                Some(&index),
                "slot of id {}",
                entry.id
            );
            assert_eq!(
                entry.cell,
                CellCoord::from_pos(entry.position, self.inv_cell),
                "stale cell of id {}",
                entry.id
            );
            let bucket = self
                .cells
                .get(&entry.cell)
                .unwrap_or_else(|| panic!("no bucket for the cell of id {}", entry.id));
            let item = bucket
                .get(entry.bucket_pos)
                .unwrap_or_else(|| panic!("bucket position of id {} is out of range", entry.id));
            assert_eq!(item.id, entry.id, "bucket item id");
            assert_eq!(
                item.index as usize, index,
                "bucket item index of id {}",
                entry.id
            );
            // Compare bit patterns: NaN positions are legal and never equal.
            assert_eq!(
                item.position.map(f64::to_bits),
                entry.position.map(f64::to_bits),
                "bucket item position of id {}",
                entry.id
            );
        }
        let bucketed: usize = self.cells.values().map(Vec::len).sum();
        assert_eq!(bucketed, self.entries.len(), "every entry is in one bucket");
        assert!(
            self.cells.values().all(|bucket| !bucket.is_empty()),
            "empty buckets must be dropped"
        );
        assert!((self.cell_size * self.inv_cell - 1.0).abs() < 1e-12);
    }
}

/// Brute-force reference model: a plain ID → position map with linear scans.
#[derive(Default, Clone)]
struct Reference {
    points: BTreeMap<u32, [f64; 3]>,
}

impl Reference {
    fn radius(&self, center: [f64; 3], radius: f64) -> Vec<u32> {
        if radius.is_nan() || radius < 0.0 {
            return Vec::new();
        }
        let r2 = radius * radius;
        self.points
            .iter()
            .filter(|&(_, &p)| dist_sq(p, center) <= r2)
            .map(|(&id, _)| id)
            .collect()
    }

    fn aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<u32> {
        self.points
            .iter()
            .filter(|&(_, p)| (0..3).all(|a| p[a] >= min[a] && p[a] <= max[a]))
            .map(|(&id, _)| id)
            .collect()
    }

    fn sorted_by_distance(&self, pos: [f64; 3]) -> Vec<(u32, f64)> {
        let mut all: Vec<(u32, f64)> = self
            .points
            .iter()
            .map(|(&id, &p)| (id, dist_sq(p, pos)))
            .collect();
        // Order by squared distance, then by ID (the documented tie-break).
        all.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        all.into_iter().map(|(id, d2)| (id, d2.sqrt())).collect()
    }

    fn nearest(&self, pos: [f64; 3]) -> Option<(u32, f64)> {
        self.sorted_by_distance(pos).into_iter().next()
    }

    fn k_nearest(&self, pos: [f64; 3], k: usize) -> Vec<(u32, f64)> {
        let mut all = self.sorted_by_distance(pos);
        all.truncate(k);
        all
    }

    fn pairs(&self, radius: f64) -> Vec<(u32, u32)> {
        if radius.is_nan() || radius < 0.0 {
            return Vec::new();
        }
        let r2 = radius * radius;
        let items: Vec<(u32, [f64; 3])> = self.points.iter().map(|(&id, &p)| (id, p)).collect();
        let mut pairs = Vec::new();
        for (i, &(a, pa)) in items.iter().enumerate() {
            for &(b, pb) in items.iter().skip(i + 1) {
                if dist_sq(pa, pb) <= r2 {
                    pairs.push((a.min(b), a.max(b)));
                }
            }
        }
        pairs.sort_unstable();
        pairs
    }
}

/// Assert that `grid` answers every kind of query exactly like `reference`.
fn assert_matches_reference(grid: &SpatialGrid, reference: &Reference, rng: &mut XorShift) {
    grid.assert_invariants();
    assert_eq!(grid.len(), reference.points.len());
    assert_eq!(grid.is_empty(), reference.points.is_empty());
    for (&id, &position) in &reference.points {
        assert_eq!(grid.position(id), Some(position));
        assert!(grid.contains(id));
    }

    for _ in 0..6 {
        let center = rng.point(12.0);
        let radius = rng.range(0.0, 9.0);
        assert_eq!(
            grid.query_radius(center, radius),
            reference.radius(center, radius),
            "query_radius({center:?}, {radius})"
        );

        let a = rng.point(12.0);
        let b = rng.point(12.0);
        let min = [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])];
        let max = [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])];
        assert_eq!(
            grid.query_aabb(min, max),
            reference.aabb(min, max),
            "query_aabb({min:?}, {max:?})"
        );

        let probe = rng.point(15.0);
        assert_eq!(grid.nearest(probe), reference.nearest(probe), "nearest");
        for k in [1, 3, 8] {
            assert_eq!(
                grid.k_nearest(probe, k),
                reference.k_nearest(probe, k),
                "k_nearest({probe:?}, {k})"
            );
        }
    }

    for radius in [0.0, 0.7, 2.5, 6.0] {
        assert_eq!(
            grid.pairs_within_radius(radius),
            reference.pairs(radius),
            "pairs_within_radius({radius})"
        );
    }
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

#[test]
fn try_new_validates_the_cell_size() {
    for bad in [
        0.0,
        -0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        5e-324, // the reciprocal overflows
    ] {
        match SpatialGrid::try_new(bad) {
            Err(GridError::InvalidCellSize(reported)) => {
                assert!(reported.to_bits() == bad.to_bits() || (reported.is_nan() && bad.is_nan()));
            }
            Ok(_) => panic!("cell size {bad} must be rejected"),
        }
    }
    for good in [1e-300, 0.001, 1.0, 5.0, 1e300] {
        let grid = SpatialGrid::try_new(good).expect("valid cell size");
        assert_eq!(grid.cell_size(), good);
        assert!(grid.is_empty());
    }
}

#[test]
#[should_panic(expected = "cell_size must be positive and finite")]
fn new_panics_on_a_non_positive_cell_size() {
    let _ = SpatialGrid::new(0.0);
}

#[test]
fn grid_error_displays_the_offending_value() {
    let message = GridError::InvalidCellSize(-2.5).to_string();
    assert!(message.contains("-2.5"), "{message}");
    assert!(message.contains("positive"), "{message}");
}

// ---------------------------------------------------------------------------
// Basic behaviour
// ---------------------------------------------------------------------------

#[test]
fn insert_update_remove_and_lookup() {
    let mut grid = SpatialGrid::new(2.0);
    assert!(grid.is_empty());
    grid.insert(7, [1.0, 2.0, 3.0]);
    grid.insert(9, [-1.0, 0.5, 8.0]);
    assert_eq!(grid.len(), 2);
    assert!(grid.contains(7) && grid.contains(9) && !grid.contains(8));
    assert_eq!(grid.position(7), Some([1.0, 2.0, 3.0]));
    assert_eq!(grid.position(8), None);
    grid.assert_invariants();

    // Inserting an existing ID replaces the entry (no duplicate).
    grid.insert(7, [10.0, 10.0, 10.0]);
    assert_eq!(grid.len(), 2);
    assert_eq!(grid.position(7), Some([10.0, 10.0, 10.0]));
    assert_eq!(grid.query_radius([1.0, 2.0, 3.0], 0.5), Vec::<u32>::new());
    grid.assert_invariants();

    // `update` moves existing entries and ignores unknown IDs.
    grid.update(9, [0.0, 0.0, 0.0]);
    grid.update(1234, [5.0, 5.0, 5.0]);
    assert_eq!(grid.position(9), Some([0.0, 0.0, 0.0]));
    assert!(!grid.contains(1234));
    assert_eq!(grid.len(), 2);
    grid.assert_invariants();

    // `remove` drops entries and ignores unknown IDs.
    grid.remove(7);
    grid.remove(7);
    grid.remove(555);
    assert_eq!(grid.len(), 1);
    assert!(!grid.contains(7));
    assert_eq!(grid.position(7), None);
    grid.assert_invariants();

    grid.clear();
    assert!(grid.is_empty());
    assert_eq!(grid.query_radius([0.0; 3], 100.0), Vec::<u32>::new());
    grid.assert_invariants();
}

#[test]
fn removal_keeps_the_index_consistent_in_every_position() {
    // Many entries per cell, removed from the front, the middle and the back of
    // both the entry array and the buckets.
    for victim in [0_u32, 1, 4, 9, 18, 19] {
        let mut grid = SpatialGrid::new(10.0);
        for id in 0..20_u32 {
            // Two cells, interleaved.
            let x = if id.is_multiple_of(2) { 1.0 } else { 11.0 };
            grid.insert(id, [x, 0.0, f64::from(id) * 0.01]);
        }
        grid.assert_invariants();
        grid.remove(victim);
        grid.assert_invariants();
        assert!(!grid.contains(victim));
        assert_eq!(grid.len(), 19);
        for id in (0..20).filter(|&id| id != victim) {
            assert!(grid.contains(id), "id {id} lost after removing {victim}");
        }
    }

    // Remove everything in a scrambled order.
    let mut grid = SpatialGrid::new(1.0);
    for id in 0..64_u32 {
        grid.insert(id, [f64::from(id % 4), 0.0, f64::from(id % 3)]);
    }
    let mut order: Vec<u32> = (0..64).collect();
    let mut rng = XorShift(0xDEAD_BEEF_1234_5678);
    for i in (1..order.len()).rev() {
        order.swap(i, rng.below(i + 1));
    }
    for id in order {
        grid.remove(id);
        grid.assert_invariants();
    }
    assert!(grid.is_empty());
    assert!(grid.cells.is_empty(), "no empty buckets may linger");
}

#[test]
fn moving_an_entry_between_cells_updates_the_cell_index() {
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(1, [0.5, 0.5, 0.5]);
    grid.insert(2, [0.6, 0.5, 0.5]);
    assert_eq!(grid.entries_in_cell(0, 0, 0), vec![1, 2]);

    grid.update(1, [3.5, 0.5, 0.5]);
    assert_eq!(grid.entries_in_cell(0, 0, 0), vec![2]);
    assert_eq!(grid.entries_in_cell(3, 0, 0), vec![1]);
    assert_eq!(grid.cell_for([3.5, 0.5, 0.5]), (3, 0, 0));
    grid.assert_invariants();

    // Moving within the same cell only changes the position.
    grid.update(1, [3.9, 0.1, 0.9]);
    assert_eq!(grid.entries_in_cell(3, 0, 0), vec![1]);
    assert_eq!(grid.position(1), Some([3.9, 0.1, 0.9]));
    grid.assert_invariants();

    // The old cell disappears once its last entry has left.
    grid.update(2, [3.2, 0.5, 0.5]);
    assert!(grid.entries_in_cell(0, 0, 0).is_empty());
    assert_eq!(grid.entries_in_cell(3, 0, 0), vec![1, 2]);
    assert_eq!(grid.cells.len(), 1);
    grid.assert_invariants();
}

#[test]
fn cell_for_uses_floor_for_negative_coordinates() {
    let grid = SpatialGrid::new(2.0);
    assert_eq!(grid.cell_for([0.0, 1.999, 2.0]), (0, 0, 1));
    assert_eq!(grid.cell_for([-0.001, -2.0, -2.001]), (-1, -1, -2));
    assert_eq!(grid.cell_for([-4.0, 4.0, 5.9]), (-2, 2, 2));
}

#[test]
fn iter_lists_every_entry_in_insertion_order() {
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(30, [0.0, 0.0, 0.0]);
    grid.insert(10, [1.0, 0.0, 0.0]);
    grid.insert(20, [2.0, 0.0, 0.0]);
    let ids: Vec<u32> = grid.iter().map(|(id, _)| id).collect();
    assert_eq!(ids, vec![30, 10, 20]);

    // Removing the first entry hands its slot to the last one.
    grid.remove(30);
    let ids: Vec<u32> = grid.iter().map(|(id, _)| id).collect();
    assert_eq!(ids, vec![20, 10]);
    assert_eq!(grid.position(20), Some([2.0, 0.0, 0.0]));
    assert_eq!(grid.position(10), Some([1.0, 0.0, 0.0]));
    grid.assert_invariants();
}

// ---------------------------------------------------------------------------
// Radius and AABB queries
// ---------------------------------------------------------------------------

#[test]
fn query_radius_example_from_the_docs() {
    let mut grid = SpatialGrid::new(5.0);
    grid.insert(0, [1.0, 0.0, 2.0]);
    grid.insert(1, [3.0, 0.0, -1.0]);
    grid.insert(2, [100.0, 0.0, 100.0]);
    let hits = grid.query_radius([0.0, 0.0, 0.0], 5.0);
    assert_eq!(hits, vec![0, 1]);
}

#[test]
fn query_radius_includes_the_boundary_and_sorts_by_id() {
    // Lattice points (exactly representable): distance 2 from the origin.
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(9, [2.0, 0.0, 0.0]);
    grid.insert(3, [0.0, -2.0, 0.0]);
    grid.insert(5, [0.0, 0.0, 2.0]);
    grid.insert(1, [2.0, 1.0, 0.0]); // sqrt(5) > 2: outside
    grid.insert(2, [0.0, 0.0, 0.0]);
    assert_eq!(grid.query_radius([0.0; 3], 2.0), vec![2, 3, 5, 9]);
    assert_eq!(grid.query_radius([0.0; 3], 0.0), vec![2]);
    // Straddling a cell boundary at x = 1.
    assert_eq!(grid.query_radius([1.0, 0.0, 0.0], 1.0), vec![2, 9]);
}

#[test]
fn query_radius_rejects_negative_and_nan_radii() {
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(1, [0.0; 3]);
    assert!(grid.query_radius([0.0; 3], -1.0).is_empty());
    assert!(grid.query_radius([0.0; 3], f64::NAN).is_empty());
    assert!(grid.query_radius([f64::NAN, 0.0, 0.0], 5.0).is_empty());
    assert_eq!(grid.query_radius([0.0; 3], f64::INFINITY), vec![1]);
}

#[test]
fn query_aabb_includes_its_boundary() {
    let mut grid = SpatialGrid::new(2.0);
    grid.insert(1, [0.0, 0.0, 0.0]);
    grid.insert(2, [4.0, 4.0, 4.0]);
    grid.insert(3, [4.0001, 0.0, 0.0]);
    grid.insert(4, [2.0, 2.0, 2.0]);
    assert_eq!(grid.query_aabb([0.0; 3], [4.0; 3]), vec![1, 2, 4]);
    assert_eq!(grid.query_aabb([2.0; 3], [2.0; 3]), vec![4]);
    // Inverted boxes are empty.
    assert!(grid.query_aabb([4.0; 3], [0.0; 3]).is_empty());
    assert!(
        grid.query_aabb([f64::NAN; 3], [1.0; 3]).is_empty(),
        "NaN bounds match nothing"
    );
}

#[test]
fn queries_work_across_the_origin() {
    let mut grid = SpatialGrid::new(1.0);
    for ix in -3..=3 {
        for iz in -3..=3 {
            let id = u32::try_from((ix + 3) * 7 + (iz + 3)).expect("non-negative");
            grid.insert(id, [f64::from(ix), 0.0, f64::from(iz)]);
        }
    }
    assert_eq!(grid.len(), 49);
    // Radius 1 around the origin: itself and its four axis neighbours.
    assert_eq!(grid.query_radius([0.0; 3], 1.0).len(), 5);
    // A box spanning the origin: 3 x 3 points.
    assert_eq!(
        grid.query_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]).len(),
        9
    );
}

// ---------------------------------------------------------------------------
// Nearest-neighbour queries
// ---------------------------------------------------------------------------

#[test]
fn nearest_and_k_nearest_basics() {
    let mut grid = SpatialGrid::new(1.0);
    assert_eq!(grid.nearest([0.0; 3]), None);
    assert!(grid.k_nearest([0.0; 3], 3).is_empty());

    for i in 0..5_u32 {
        grid.insert(i, [f64::from(i), 0.0, 0.0]);
    }
    let (id, distance) = grid.nearest([0.2, 0.0, 0.0]).expect("non-empty");
    assert_eq!(id, 0);
    assert!((distance - 0.2).abs() < 1e-12);
    let (id, distance) = grid.nearest([3.9, 0.0, 0.0]).expect("non-empty");
    assert_eq!(id, 4);
    assert!((distance - 0.1).abs() < 1e-12);
    assert_eq!(
        grid.k_nearest([0.0; 3], 3),
        vec![(0, 0.0), (1, 1.0), (2, 2.0)]
    );
    assert!(grid.k_nearest([0.0; 3], 0).is_empty());
    // Asking for more than there is returns everything, nearest first.
    let all = grid.k_nearest([10.0, 0.0, 0.0], 100);
    assert_eq!(
        all.iter().map(|&(id, _)| id).collect::<Vec<_>>(),
        vec![4, 3, 2, 1, 0]
    );
}

#[test]
fn nearest_neighbours_break_ties_by_id() {
    let mut grid = SpatialGrid::new(1.0);
    // Four entries at distance exactly 1 from the origin, in different cells.
    grid.insert(40, [1.0, 0.0, 0.0]);
    grid.insert(20, [-1.0, 0.0, 0.0]);
    grid.insert(30, [0.0, 1.0, 0.0]);
    grid.insert(10, [0.0, 0.0, -1.0]);
    assert_eq!(grid.nearest([0.0; 3]), Some((10, 1.0)));
    assert_eq!(
        grid.k_nearest([0.0; 3], 3),
        vec![(10, 1.0), (20, 1.0), (30, 1.0)]
    );
    // A coincident entry always wins, whatever its ID.
    grid.insert(99, [0.0; 3]);
    assert_eq!(grid.nearest([0.0; 3]), Some((99, 0.0)));
}

#[test]
fn nearest_finds_far_away_entries_in_a_sparse_grid() {
    let mut grid = SpatialGrid::new(0.5);
    grid.insert(1, [1.0e6, 0.0, 0.0]);
    grid.insert(2, [-2.0e6, 0.0, 0.0]);
    assert_eq!(grid.nearest([0.0; 3]), Some((1, 1.0e6)));
    assert_eq!(grid.k_nearest([0.0; 3], 2), vec![(1, 1.0e6), (2, 2.0e6)]);
}

#[test]
fn non_finite_entries_never_match_range_queries_and_come_last() {
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(1, [f64::NAN, 0.0, 0.0]);
    grid.insert(2, [f64::INFINITY, 0.0, 0.0]);
    grid.insert(3, [1.0e150, -1.0e150, 0.0]);
    grid.insert(4, [2.0, 0.0, 0.0]);
    grid.assert_invariants();

    // Finite queries: only the finite entries inside them match.
    assert_eq!(grid.query_radius([2.0, 0.0, 0.0], 1.0), vec![4]);
    assert_eq!(grid.query_radius([0.0; 3], 1.0e100), vec![4]);
    assert_eq!(grid.query_aabb([-10.0; 3], [10.0; 3]), vec![4]);
    assert_eq!(grid.pairs_within_radius(1.0e6), Vec::<(u32, u32)>::new());
    // Only an unbounded query reaches the infinite entry; NaN is never matched.
    assert_eq!(grid.query_radius([0.0; 3], f64::INFINITY), vec![2, 3, 4]);
    assert_eq!(
        grid.query_aabb([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]),
        vec![2, 3, 4]
    );
    assert_eq!(grid.nearest([0.0; 3]), Some((4, 2.0)));

    let ranked: Vec<u32> = grid
        .k_nearest([0.0; 3], 4)
        .into_iter()
        .map(|(id, _)| id)
        .collect();
    assert_eq!(ranked[0], 4);
    assert_eq!(ranked[1], 3);
    let tail: HashSet<u32> = ranked[2..].iter().copied().collect();
    assert_eq!(tail, HashSet::from([1, 2]));

    // A grid of nothing but NaN still answers, at infinite distance.
    let mut only_nan = SpatialGrid::new(1.0);
    only_nan.insert(5, [f64::NAN; 3]);
    let (id, distance) = only_nan.nearest([0.0; 3]).expect("one entry");
    assert_eq!(id, 5);
    assert!(distance.is_infinite());
}

// ---------------------------------------------------------------------------
// Pairs
// ---------------------------------------------------------------------------

#[test]
fn pairs_within_radius_returns_sorted_unique_pairs() {
    let mut grid = SpatialGrid::new(1.0);
    grid.insert(5, [0.0, 0.0, 0.0]);
    grid.insert(2, [1.0, 0.0, 0.0]);
    grid.insert(8, [2.0, 0.0, 0.0]);
    grid.insert(1, [10.0, 0.0, 0.0]);
    // Distance exactly 1 counts.
    assert_eq!(grid.pairs_within_radius(1.0), vec![(2, 5), (2, 8)]);
    assert_eq!(grid.pairs_within_radius(2.0), vec![(2, 5), (2, 8), (5, 8)]);
    assert_eq!(grid.pairs_within_radius(0.5), Vec::<(u32, u32)>::new());
    assert!(grid.pairs_within_radius(-1.0).is_empty());
    assert!(grid.pairs_within_radius(f64::NAN).is_empty());
    // Coincident points are a pair at radius zero.
    grid.insert(7, [10.0, 0.0, 0.0]);
    assert_eq!(grid.pairs_within_radius(0.0), vec![(1, 7)]);
}

#[test]
fn pairs_are_found_across_cell_boundaries_and_for_radii_beyond_one_cell() {
    let mut grid = SpatialGrid::new(1.0);
    // 5 x 5 x 5 lattice with spacing 1.
    let mut id = 0_u32;
    let mut reference = Reference::default();
    for x in 0..5 {
        for y in 0..5 {
            for z in 0..5 {
                let p = [f64::from(x), f64::from(y), f64::from(z)];
                grid.insert(id, p);
                reference.points.insert(id, p);
                id += 1;
            }
        }
    }
    for radius in [1.0, 1.5, 2.0, 2.9, 3.0, 4.5, 100.0] {
        assert_eq!(
            grid.pairs_within_radius(radius),
            reference.pairs(radius),
            "radius {radius}"
        );
    }
    // Unit-distance pairs of an n^3 lattice: 3 n^2 (n - 1).
    assert_eq!(grid.pairs_within_radius(1.0).len(), 3 * 25 * 4);
}

#[test]
fn pairs_with_a_cell_size_far_below_the_radius_fall_back_to_all_pairs() {
    let mut grid = SpatialGrid::new(1.0e-3);
    let mut reference = Reference::default();
    let mut rng = XorShift(7);
    for id in 0..60_u32 {
        let p = rng.point(1.0);
        grid.insert(id, p);
        reference.points.insert(id, p);
    }
    for radius in [0.1, 0.5, 1.0, 3.0] {
        assert_eq!(grid.pairs_within_radius(radius), reference.pairs(radius));
    }
}

// ---------------------------------------------------------------------------
// Differential tests against the brute-force reference
// ---------------------------------------------------------------------------

#[test]
fn random_operations_match_the_brute_force_reference() {
    for (seed, cell_size) in [(1_u64, 1.0), (2, 0.37), (3, 4.0), (4, 25.0), (5, 0.05)] {
        let mut rng = XorShift(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let mut grid = SpatialGrid::new(cell_size);
        let mut reference = Reference::default();

        for step in 0..400_usize {
            let id = u32::try_from(rng.below(48)).expect("small");
            match rng.below(10) {
                0..=4 => {
                    let p = rng.point(12.0);
                    grid.insert(id, p);
                    reference.points.insert(id, p);
                }
                5..=6 => {
                    let p = rng.point(12.0);
                    grid.update(id, p);
                    if let Some(slot) = reference.points.get_mut(&id) {
                        *slot = p;
                    }
                }
                7..=8 => {
                    grid.remove(id);
                    reference.points.remove(&id);
                }
                _ => {
                    if step.is_multiple_of(3) {
                        let new_size = rng.range(0.1, 6.0);
                        grid.set_cell_size(new_size).expect("valid cell size");
                        assert_eq!(grid.cell_size(), new_size);
                    }
                }
            }
            if step.is_multiple_of(20) {
                assert_matches_reference(&grid, &reference, &mut rng);
            }
        }
        assert_matches_reference(&grid, &reference, &mut rng);
    }
}

#[test]
fn dense_clusters_match_the_brute_force_reference() {
    // Many entries per cell and many coincident points.
    let mut rng = XorShift(0xC0FFEE);
    let mut grid = SpatialGrid::new(3.0);
    let mut reference = Reference::default();
    for id in 0..150_u32 {
        let centre = [
            f64::from(u32::try_from(rng.below(3)).expect("small")) * 3.0,
            0.0,
            f64::from(u32::try_from(rng.below(2)).expect("small")) * 3.0,
        ];
        let p = [
            centre[0] + rng.range(0.0, 0.5),
            centre[1] + rng.range(0.0, 0.5),
            centre[2] + rng.range(0.0, 0.5),
        ];
        grid.insert(id, p);
        reference.points.insert(id, p);
    }
    for id in (0..150_u32).step_by(7) {
        let duplicate = reference.points[&id]; // exact duplicates of positions
        grid.insert(id + 1000, duplicate);
        reference.points.insert(id + 1000, duplicate);
    }
    assert_matches_reference(&grid, &reference, &mut rng);
}

#[test]
fn results_do_not_depend_on_insertion_order() {
    let mut rng = XorShift(99);
    let points: Vec<(u32, [f64; 3])> = (0..120_u32).map(|id| (id, rng.point(8.0))).collect();
    let mut shuffled = points.clone();
    for i in (1..shuffled.len()).rev() {
        shuffled.swap(i, rng.below(i + 1));
    }

    let build = |items: &[(u32, [f64; 3])]| {
        let mut grid = SpatialGrid::new(1.5);
        for &(id, p) in items {
            grid.insert(id, p);
        }
        grid
    };
    let a = build(&points);
    let b = build(&shuffled);

    for _ in 0..25 {
        let center = rng.point(8.0);
        assert_eq!(a.query_radius(center, 3.0), b.query_radius(center, 3.0));
        assert_eq!(a.nearest(center), b.nearest(center));
        assert_eq!(a.k_nearest(center, 6), b.k_nearest(center, 6));
    }
    assert_eq!(a.pairs_within_radius(1.2), b.pairs_within_radius(1.2));
    assert_eq!(
        a.query_aabb([-3.0; 3], [3.0; 3]),
        b.query_aabb([-3.0; 3], [3.0; 3])
    );
}

// ---------------------------------------------------------------------------
// Changing the cell size (regression: `cell_size` was a public field, so
// changing it left every cached cell stale)
// ---------------------------------------------------------------------------

#[test]
fn set_cell_size_rebuckets_every_entry() {
    let mut grid = SpatialGrid::new(10.0);
    grid.insert(1, [1.0, 1.0, 1.0]);
    grid.insert(2, [3.0, 1.0, 1.0]);
    grid.insert(3, [-5.0, 1.0, 1.0]);
    assert_eq!(grid.entries_in_cell(0, 0, 0), vec![1, 2]);
    assert_eq!(grid.entries_in_cell(-1, 0, 0), vec![3]);

    grid.set_cell_size(2.0).expect("valid");
    assert_eq!(grid.cell_size(), 2.0);
    assert_eq!(grid.cell_for([3.0, 1.0, 1.0]), (1, 0, 0));
    assert_eq!(grid.entries_in_cell(0, 0, 0), vec![1]);
    assert_eq!(grid.entries_in_cell(1, 0, 0), vec![2]);
    assert_eq!(grid.entries_in_cell(-3, 0, 0), vec![3]);
    grid.assert_invariants();

    // Queries see the same entries as before.
    assert_eq!(grid.query_radius([1.0, 1.0, 1.0], 2.5), vec![1, 2]);
    assert_eq!(grid.nearest([-4.0, 1.0, 1.0]), Some((3, 1.0)));

    // Entries inserted afterwards use the new cell size.
    grid.insert(4, [5.0, 1.0, 1.0]);
    assert_eq!(grid.entries_in_cell(2, 0, 0), vec![4]);
    grid.assert_invariants();
}

#[test]
fn set_cell_size_rejects_invalid_values_and_keeps_the_grid() {
    let mut grid = SpatialGrid::new(4.0);
    grid.insert(1, [5.0, 0.0, 0.0]);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            grid.set_cell_size(bad),
            Err(GridError::InvalidCellSize(_))
        ));
        assert_eq!(grid.cell_size(), 4.0);
        assert_eq!(grid.entries_in_cell(1, 0, 0), vec![1]);
        grid.assert_invariants();
    }
}

// ---------------------------------------------------------------------------
// Serialisation
// ---------------------------------------------------------------------------

#[test]
fn serde_round_trip_preserves_entries_and_answers() {
    let mut rng = XorShift(2024);
    let mut grid = SpatialGrid::new(2.5);
    for id in 0..80_u32 {
        grid.insert(id * 3, rng.point(10.0));
    }
    grid.remove(9);
    grid.remove(30);

    let json = serde_json::to_string(&grid).expect("serialise");
    let restored: SpatialGrid = serde_json::from_str(&json).expect("deserialise");
    restored.assert_invariants();
    assert_eq!(restored.cell_size(), grid.cell_size());
    assert_eq!(restored.len(), grid.len());
    // Entry order survives the round trip.  (JSON does not promise to
    // reproduce f64 values bit for bit, hence the tolerance on positions.)
    let original: Vec<(u32, [f64; 3])> = grid.iter().collect();
    let loaded: Vec<(u32, [f64; 3])> = restored.iter().collect();
    assert_eq!(
        original.iter().map(|&(id, _)| id).collect::<Vec<_>>(),
        loaded.iter().map(|&(id, _)| id).collect::<Vec<_>>(),
        "entry order survives the round trip"
    );
    for (&(id, a), &(_, b)) in original.iter().zip(&loaded) {
        assert!(
            (0..3).all(|axis| (a[axis] - b[axis]).abs() < 1e-12),
            "position of id {id} changed: {a:?} -> {b:?}"
        );
    }
    for _ in 0..10 {
        let c = rng.point(10.0);
        assert_eq!(restored.query_radius(c, 4.0), grid.query_radius(c, 4.0));
        let before = grid.k_nearest(c, 5);
        let after = restored.k_nearest(c, 5);
        assert_eq!(
            before.iter().map(|&(id, _)| id).collect::<Vec<_>>(),
            after.iter().map(|&(id, _)| id).collect::<Vec<_>>()
        );
        for (&(_, d_before), &(_, d_after)) in before.iter().zip(&after) {
            assert!((d_before - d_after).abs() < 1e-12);
        }
    }
    assert_eq!(
        restored.pairs_within_radius(2.0),
        grid.pairs_within_radius(2.0)
    );
}

#[test]
fn serialised_form_keeps_the_established_wire_format() {
    let mut grid = SpatialGrid::new(2.0);
    grid.insert(7, [3.0, -1.0, 0.5]);
    let value: serde_json::Value = serde_json::to_value(&grid).expect("serialise");
    assert_eq!(value["cell_size"], 2.0);
    let entries = value["entries"].as_array().expect("entries array");
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0]["id"], 7);
    assert_eq!(entries[0]["position"], serde_json::json!([3.0, -1.0, 0.5]));
    assert_eq!(
        entries[0]["cell"],
        serde_json::json!({"x": 1, "y": -1, "z": 0})
    );
    // Runtime-only bookkeeping must not leak into the format.
    assert!(entries[0].get("bucket_pos").is_none());
    assert!(value.get("slot").is_none() && value.get("cells").is_none());
}

#[test]
fn deserialising_recomputes_cells_instead_of_trusting_them() {
    // The stored cell of id 1 is wrong (and the grid was "edited" by hand).
    let json = r#"{
        "cell_size": 1.0,
        "entries": [
            {"id": 1, "position": [5.5, 0.5, 0.5], "cell": {"x": 0, "y": 0, "z": 0}},
            {"id": 2, "position": [0.5, 0.5, 0.5], "cell": {"x": 0, "y": 0, "z": 0}}
        ]
    }"#;
    let grid: SpatialGrid = serde_json::from_str(json).expect("deserialise");
    grid.assert_invariants();
    assert_eq!(grid.entries_in_cell(5, 0, 0), vec![1]);
    assert_eq!(grid.entries_in_cell(0, 0, 0), vec![2]);
    assert_eq!(grid.query_radius([5.5, 0.5, 0.5], 0.1), vec![1]);
}

#[test]
fn deserialising_accepts_entries_without_a_cell() {
    let json = r#"{"cell_size": 2.0, "entries": [{"id": 3, "position": [5.0, -1.0, 0.0]}]}"#;
    let grid: SpatialGrid = serde_json::from_str(json).expect("deserialise");
    grid.assert_invariants();
    assert_eq!(grid.entries_in_cell(2, -1, 0), vec![3]);
}

#[test]
fn deserialising_validates_the_cell_size() {
    for bad in ["0.0", "-3.0"] {
        let json = format!(r#"{{"cell_size": {bad}, "entries": []}}"#);
        let error = serde_json::from_str::<SpatialGrid>(&json).expect_err("must be rejected");
        assert!(error.to_string().contains("cell_size"), "{error}");
    }
}

#[test]
fn deserialising_duplicate_ids_keeps_the_last_entry() {
    let json = r#"{
        "cell_size": 1.0,
        "entries": [
            {"id": 4, "position": [0.0, 0.0, 0.0], "cell": {"x": 0, "y": 0, "z": 0}},
            {"id": 4, "position": [9.0, 0.0, 0.0], "cell": {"x": 9, "y": 0, "z": 0}}
        ]
    }"#;
    let grid: SpatialGrid = serde_json::from_str(json).expect("deserialise");
    assert_eq!(grid.len(), 1);
    assert_eq!(grid.position(4), Some([9.0, 0.0, 0.0]));
    grid.assert_invariants();
}

// ---------------------------------------------------------------------------
// Scalability: queries only touch the entries near them
// ---------------------------------------------------------------------------

/// `n x n x n` lattice with unit spacing, IDs in x-major order.
fn lattice(n: u32, cell_size: f64) -> SpatialGrid {
    let mut grid = SpatialGrid::new(cell_size);
    let mut id = 0;
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                grid.insert(id, [f64::from(x), f64::from(y), f64::from(z)]);
                id += 1;
            }
        }
    }
    grid
}

#[test]
fn a_range_query_only_visits_entries_in_nearby_cells() {
    // 64 000 entries; a radius-3 query must not look at (nearly) all of them.
    let grid = lattice(40, 2.0);
    assert_eq!(grid.len(), 64_000);

    let center = [20.0, 20.0, 20.0];
    let radius = 3.0;
    let inv = grid.inv_cell;
    let min = CellCoord::from_pos([17.0, 17.0, 17.0], inv);
    let max = CellCoord::from_pos([23.0, 23.0, 23.0], inv);
    let mut visited = 0_usize;
    grid.for_each_in_cell_box(min, max, |_| visited += 1);
    assert!(
        visited <= 64 * 8 + 8,
        "visited {visited} of {} entries",
        grid.len()
    );

    // The answer is the ball of lattice points: 123 points within radius 3.
    let hits = grid.query_radius(center, radius);
    let expected = (-3_i32..=3)
        .flat_map(|dx| (-3_i32..=3).flat_map(move |dy| (-3_i32..=3).map(move |dz| (dx, dy, dz))))
        .filter(|&(dx, dy, dz)| dx * dx + dy * dy + dz * dz <= 9)
        .count();
    assert_eq!(expected, 123);
    assert_eq!(hits.len(), expected);
}

/// Counts the entries a nearest-neighbour search looks at.
struct CountingSearch<'a> {
    inner: &'a mut KNearestSearch,
    visited: usize,
}

impl RingSearch for CountingSearch<'_> {
    fn visit(&mut self, item: &BucketItem) {
        self.visited += 1;
        self.inner.visit(item);
    }

    fn satisfied(&self, unvisited_min_dist_sq: f64) -> bool {
        self.inner.satisfied(unvisited_min_dist_sq)
    }
}

#[test]
fn nearest_neighbour_search_only_visits_entries_in_nearby_shells() {
    let grid = lattice(40, 2.0);
    let probe = [20.2, 20.2, 20.2];
    let mut knn = KNearestSearch {
        pos: probe,
        k: 5,
        heap: BinaryHeap::new(),
    };
    let mut counting = CountingSearch {
        inner: &mut knn,
        visited: 0,
    };
    grid.search_rings(probe, &mut counting);
    let visited = counting.visited;
    assert!(
        visited < 64_000 / 20,
        "a 5-NN query visited {visited} of 64000 entries"
    );
    assert_eq!(knn.heap.len(), 5);

    let nearest = grid.nearest(probe).expect("non-empty");
    assert_eq!(nearest.0, 20 * 1600 + 20 * 40 + 20);
    assert!((nearest.1 - (3.0_f64 * 0.04).sqrt()).abs() < 1e-12);
}

#[test]
fn pair_search_compares_only_entries_of_nearby_cells() {
    // 13 824 entries: all pairs would be ~10^8 comparisons.
    let grid = lattice(24, 1.0);
    let reach = grid.reach_cells(1.0);
    assert_eq!(reach, 1);
    let mut candidates = 0_usize;
    grid.for_each_candidate_pair(reach, |_, _| candidates += 1);
    assert!(
        candidates < 13_824 * 27,
        "visited {candidates} candidate pairs"
    );

    let pairs = grid.pairs_within_radius(1.0);
    assert_eq!(pairs.len(), 3 * 24 * 24 * 23);
    // Spot-check a few of the pairs (the list is sorted).
    assert!(pairs.binary_search(&(0, 1)).is_ok());
    assert!(pairs.binary_search(&(0, 24)).is_ok());
    assert!(pairs.binary_search(&(0, 576)).is_ok());
    assert!(pairs.binary_search(&(0, 25)).is_err());
}

#[test]
fn every_candidate_pair_is_visited_exactly_once() {
    let mut rng = XorShift(31337);
    let mut grid = SpatialGrid::new(1.0);
    for id in 0..200_u32 {
        grid.insert(id, rng.point(4.0));
    }
    for reach in [0_u32, 1, 2, 3] {
        let mut seen: HashSet<(u32, u32)> = HashSet::new();
        grid.for_each_candidate_pair(reach, |a, b| {
            let key = (a.id.min(b.id), a.id.max(b.id));
            assert!(
                seen.insert(key),
                "pair {key:?} visited twice (reach {reach})"
            );
        });
        // Every pair of entries within `reach` cells of each other (Chebyshev)
        // must be among the candidates.
        for a in &grid.entries {
            for b in &grid.entries {
                if a.id < b.id && a.cell.chebyshev_distance(b.cell) <= reach {
                    assert!(
                        seen.contains(&(a.id, b.id)),
                        "pair ({}, {}) missing at reach {reach}",
                        a.id,
                        b.id
                    );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Hasher
// ---------------------------------------------------------------------------

#[test]
fn cell_hasher_spreads_neighbouring_cells_over_the_low_bits() {
    use std::hash::BuildHasher;
    let build = BuildHasherDefault::<CellHasher>::default();
    let mut buckets: HashSet<u64> = HashSet::new();
    let mut total = 0_usize;
    for x in 0..16 {
        for y in 0..16 {
            for z in 0..16 {
                total += 1;
                buckets.insert(build.hash_one(CellCoord { x, y, z }) & 0xFFF);
            }
        }
    }
    // 4096 keys into 4096 buckets: a good hash fills ~63 % of them.
    assert_eq!(total, 4096);
    assert!(
        buckets.len() > 2300,
        "only {} distinct low-bit patterns for 4096 neighbouring cells",
        buckets.len()
    );
}

#[test]
fn cell_hasher_is_deterministic() {
    use std::hash::BuildHasher;
    let a = BuildHasherDefault::<CellHasher>::default();
    let b = BuildHasherDefault::<CellHasher>::default();
    let cell = CellCoord {
        x: -4,
        y: 9,
        z: 1 << 20,
    };
    assert_eq!(a.hash_one(cell), b.hash_one(cell));
    assert_ne!(
        a.hash_one(cell),
        a.hash_one(CellCoord {
            x: 9,
            y: -4,
            z: 1 << 20
        })
    );
}
