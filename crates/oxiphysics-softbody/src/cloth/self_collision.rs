// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uniform-grid broadphase for cloth vertex-vertex self-collision.
//!
//! [`SelfCollisionGrid`] hashes every vertex into a cell of edge length
//! `min_dist`, counting-sorts the vertex indices by cell (no `HashMap`) and
//! answers neighbourhood queries by visiting the 3 x 3 x 3 block of cells around
//! a point.  Two vertices closer than `min_dist` differ by less than one cell
//! in every coordinate, so they always sit in neighbouring cells: the grid
//! reports exactly the brute-force set of overlapping pairs while only touching
//! a number of candidate pairs that grows linearly with the vertex count for a
//! cloth of bounded vertex density.
//!
//! The grid keeps its scratch buffers between calls, so a simulation that
//! resolves self-collision every step allocates only on the first call.

use super::functions::{add, dot, scale, sub};
use super::types::ClothVertex;

/// Relative slack applied to the cell size.
///
/// It guarantees that two coordinates less than `min_dist` apart land in
/// neighbouring cells even when the division by the cell size rounds.
const CELL_SLACK: f64 = 1.0 + 1e-9;

/// Separation (m) at or below which two vertices count as coincident and no
/// separation direction can be derived.
const COINCIDENT_DISTANCE: f64 = 1e-14;

/// Combined inverse mass below which a pair is treated as immovable.
const IMMOVABLE_INVERSE_MASS: f64 = 1e-30;

/// Smallest number of hash buckets.
const MIN_BUCKETS: usize = 16;

/// Largest number of hash buckets (the largest power of two that fits).
const MAX_BUCKETS: usize = (usize::MAX >> 2) + 1;

/// Returns `true` when `min_dist` is a usable, strictly positive, finite length.
fn is_usable_distance(min_dist: f64) -> bool {
    min_dist.is_finite() && min_dist > 0.0
}

/// Reusable uniform-grid broadphase over cloth vertices.
///
/// Build it with [`SelfCollisionGrid::rebuild`] (or let
/// [`SelfCollisionGrid::overlapping_pairs`] / [`SelfCollisionGrid::resolve`]
/// do so), then query neighbourhoods with [`SelfCollisionGrid::neighbours_of`].
/// The grid is rebuilt from scratch on every call, which costs O(n); the
/// buffers it owns are reused so repeated calls do not allocate.
#[derive(Debug, Clone, Default)]
pub struct SelfCollisionGrid {
    /// Reciprocal of the cell edge length (`0.0` collapses everything into one cell).
    inv_cell: f64,
    /// `bucket_count - 1`; the bucket count is a power of two.
    bucket_mask: usize,
    /// `sorted[bucket_start[b]..bucket_start[b + 1]]` are the vertices of bucket `b`.
    bucket_start: Vec<usize>,
    /// Scratch write cursor of the counting sort.
    cursor: Vec<usize>,
    /// Vertex indices ordered by bucket, ascending within a bucket.
    sorted: Vec<usize>,
    /// Bucket of every vertex (scratch of the counting sort).
    vertex_bucket: Vec<usize>,
    /// Scratch list of the candidate vertices of one query.
    candidates: Vec<usize>,
}

impl SelfCollisionGrid {
    /// Create an empty grid.
    pub fn new() -> Self {
        Self::default()
    }

    /// Integer cell coordinate of a world coordinate.
    ///
    /// The float-to-integer cast saturates, and maps `NaN` to cell `0`, so
    /// degenerate coordinates never panic.
    fn cell_coord(&self, coordinate: f64) -> i64 {
        (coordinate * self.inv_cell).floor() as i64
    }

    /// Hash bucket of the cell `(cx, cy, cz)`.
    fn bucket(&self, cx: i64, cy: i64, cz: i64) -> usize {
        let mut hash = (cx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ (cy as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
            ^ (cz as u64).wrapping_mul(0x1656_67B1_9E37_79F9);
        hash ^= hash >> 32;
        hash = hash.wrapping_mul(0xD6E8_FEB8_6659_FD93);
        hash ^= hash >> 32;
        (hash as usize) & self.bucket_mask
    }

    /// Hash bucket of the cell that contains `position`.
    fn bucket_of(&self, position: [f64; 3]) -> usize {
        self.bucket(
            self.cell_coord(position[0]),
            self.cell_coord(position[1]),
            self.cell_coord(position[2]),
        )
    }

    /// Rebuild the grid from the current vertex positions.
    ///
    /// `cell_size` is the edge length of one cell and should equal the
    /// collision distance that the queries are used for.  A non-positive or
    /// non-finite `cell_size` collapses all vertices into a single cell, which
    /// stays correct (every vertex is a candidate of every other) but is slow.
    pub fn rebuild(&mut self, vertices: &[ClothVertex], cell_size: f64) {
        self.inv_cell = if is_usable_distance(cell_size) {
            1.0 / (cell_size * CELL_SLACK)
        } else {
            0.0
        };
        let vertex_count = vertices.len();
        let bucket_count = vertex_count
            .clamp(MIN_BUCKETS, MAX_BUCKETS)
            .next_power_of_two();
        self.bucket_mask = bucket_count - 1;

        self.bucket_start.clear();
        self.bucket_start.resize(bucket_count + 1, 0);
        self.vertex_bucket.clear();
        self.vertex_bucket.reserve(vertex_count);
        for vertex in vertices {
            let bucket = self.bucket_of(vertex.pos);
            self.vertex_bucket.push(bucket);
            self.bucket_start[bucket] += 1;
        }

        // Exclusive prefix sum: counts become the start offset of every bucket,
        // and the trailing entry becomes the total number of vertices.
        let mut running = 0;
        for entry in self.bucket_start.iter_mut() {
            let count = *entry;
            *entry = running;
            running += count;
        }

        self.cursor.clear();
        self.cursor
            .extend_from_slice(&self.bucket_start[..bucket_count]);
        self.sorted.clear();
        self.sorted.resize(vertex_count, 0);
        for (vertex, &bucket) in self.vertex_bucket.iter().enumerate() {
            let slot = self.cursor[bucket];
            self.sorted[slot] = vertex;
            self.cursor[bucket] = slot + 1;
        }
    }

    /// Collect into `out` the candidate neighbours of a point.
    ///
    /// `out` is cleared and then filled with the indices (ascending, without
    /// duplicates) of every vertex stored in the 3 x 3 x 3 block of cells
    /// around `position`.  This is a superset of all vertices within one cell
    /// size of `position`: hash collisions may add vertices from unrelated
    /// cells, so callers must still apply an exact distance test.
    pub fn neighbours_of(&self, position: [f64; 3], out: &mut Vec<usize>) {
        out.clear();
        if self.sorted.is_empty() {
            return;
        }
        let cx = self.cell_coord(position[0]);
        let cy = self.cell_coord(position[1]);
        let cz = self.cell_coord(position[2]);
        let mut buckets = [0_usize; 27];
        let mut count = 0;
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    buckets[count] = self.bucket(
                        cx.wrapping_add(dx),
                        cy.wrapping_add(dy),
                        cz.wrapping_add(dz),
                    );
                    count += 1;
                }
            }
        }
        buckets.sort_unstable();
        let mut previous = usize::MAX;
        for &bucket in &buckets {
            if bucket == previous {
                continue;
            }
            previous = bucket;
            out.extend_from_slice(
                &self.sorted[self.bucket_start[bucket]..self.bucket_start[bucket + 1]],
            );
        }
        out.sort_unstable();
    }

    /// All vertex pairs `(i, j)` with `i < j` that are closer than `min_dist`.
    ///
    /// The pairs are ordered by `i`, then `j`, and are exactly the pairs a
    /// brute-force double loop with the strict test `distance < min_dist`
    /// finds.  A non-positive or non-finite `min_dist` yields no pairs.
    pub fn overlapping_pairs(
        &mut self,
        vertices: &[ClothVertex],
        min_dist: f64,
    ) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if !is_usable_distance(min_dist) || vertices.len() < 2 {
            return pairs;
        }
        self.rebuild(vertices, min_dist);
        let min_dist_sq = min_dist * min_dist;
        let mut candidates = std::mem::take(&mut self.candidates);
        for (i, vertex) in vertices.iter().enumerate() {
            self.neighbours_of(vertex.pos, &mut candidates);
            let first = candidates.partition_point(|&j| j <= i);
            for &j in &candidates[first..] {
                let separation = sub(vertices[j].pos, vertex.pos);
                if dot(separation, separation) < min_dist_sq {
                    pairs.push((i, j));
                }
            }
        }
        self.candidates = candidates;
        pairs
    }

    /// Push apart every pair of vertices closer than `min_dist`.
    ///
    /// The grid is built once from the positions at the start of the call.
    /// Candidate pairs are then visited in ascending `(i, j)` order and
    /// corrected immediately (Gauss-Seidel), re-measuring the live distance of
    /// each pair, so earlier corrections influence later ones exactly as in a
    /// brute-force sweep.  A pair that only starts to overlap because of
    /// corrections made during this very sweep, and that was more than one
    /// cell apart when the grid was built, is picked up by the next call.
    ///
    /// Each overlapping pair is moved along its connecting line until it is
    /// `min_dist` apart, split by inverse mass; vertices flagged
    /// [`ClothVertex::fixed`] never move, and vertices that coincide exactly
    /// have no separation direction and are left alone.  A non-positive or
    /// non-finite `min_dist` is a no-op.
    pub fn resolve(&mut self, vertices: &mut [ClothVertex], min_dist: f64) {
        if !is_usable_distance(min_dist) || vertices.len() < 2 {
            return;
        }
        self.rebuild(vertices, min_dist);
        let min_dist_sq = min_dist * min_dist;
        let mut candidates = std::mem::take(&mut self.candidates);
        for i in 0..vertices.len() {
            self.neighbours_of(vertices[i].pos, &mut candidates);
            let first = candidates.partition_point(|&j| j <= i);
            for &j in &candidates[first..] {
                separate_pair(vertices, i, j, min_dist, min_dist_sq);
            }
        }
        self.candidates = candidates;
    }
}

/// Find all vertex pairs `(i, j)` with `i < j` closer than `min_dist`.
///
/// Convenience wrapper around [`SelfCollisionGrid::overlapping_pairs`] that
/// uses a temporary grid.
pub fn find_self_collision_pairs(vertices: &[ClothVertex], min_dist: f64) -> Vec<(usize, usize)> {
    SelfCollisionGrid::new().overlapping_pairs(vertices, min_dist)
}

/// Push apart vertices `i` and `j` if they are closer than `min_dist`.
fn separate_pair(
    vertices: &mut [ClothVertex],
    i: usize,
    j: usize,
    min_dist: f64,
    min_dist_sq: f64,
) {
    let d = sub(vertices[j].pos, vertices[i].pos);
    let dist_sq = dot(d, d);
    if dist_sq >= min_dist_sq {
        return;
    }
    let dist = dist_sq.sqrt();
    if dist <= COINCIDENT_DISTANCE {
        return;
    }
    let overlap = min_dist - dist;
    let dir = scale(d, 1.0 / dist);
    let wi = vertices[i].inv_mass;
    let wj = vertices[j].inv_mass;
    let w_sum = wi + wj;
    if w_sum < IMMOVABLE_INVERSE_MASS {
        return;
    }
    let corr = scale(dir, overlap / w_sum);
    if !vertices[i].fixed {
        vertices[i].pos = sub(vertices[i].pos, scale(corr, wi));
    }
    if !vertices[j].fixed {
        vertices[j].pos = add(vertices[j].pos, scale(corr, wj));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic 64-bit linear congruential generator for test clouds.
    struct Lcg(u64);

    impl Lcg {
        fn next_unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }

        fn point_in_box(&mut self, lo: f64, hi: f64) -> [f64; 3] {
            [
                lo + (hi - lo) * self.next_unit(),
                lo + (hi - lo) * self.next_unit(),
                lo + (hi - lo) * self.next_unit(),
            ]
        }
    }

    fn vertex(pos: [f64; 3], inv_mass: f64, fixed: bool) -> ClothVertex {
        ClothVertex {
            pos,
            prev_pos: pos,
            vel: [0.0; 3],
            mass: if inv_mass > 0.0 { 1.0 / inv_mass } else { 0.0 },
            inv_mass,
            fixed,
            normal: [0.0, 1.0, 0.0],
        }
    }

    fn random_cloud(count: usize, lo: f64, hi: f64, seed: u64) -> Vec<ClothVertex> {
        let mut rng = Lcg(seed);
        (0..count)
            .map(|_| vertex(rng.point_in_box(lo, hi), 1.0, false))
            .collect()
    }

    fn brute_force_pairs(vertices: &[ClothVertex], min_dist: f64) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                let d = sub(vertices[j].pos, vertices[i].pos);
                if dot(d, d).sqrt() < min_dist {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }

    /// The O(n^2) Gauss-Seidel sweep that the grid replaced.
    fn brute_force_resolve(vertices: &mut [ClothVertex], min_dist: f64) {
        let n = vertices.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = sub(vertices[j].pos, vertices[i].pos);
                let dist = dot(d, d).sqrt();
                if dist < min_dist && dist > 1e-14 {
                    let overlap = min_dist - dist;
                    let dir = scale(d, 1.0 / dist);
                    let wi = vertices[i].inv_mass;
                    let wj = vertices[j].inv_mass;
                    let w_sum = wi + wj;
                    if w_sum < 1e-30 {
                        continue;
                    }
                    let corr = scale(dir, overlap / w_sum);
                    if !vertices[i].fixed {
                        vertices[i].pos = sub(vertices[i].pos, scale(corr, wi));
                    }
                    if !vertices[j].fixed {
                        vertices[j].pos = add(vertices[j].pos, scale(corr, wj));
                    }
                }
            }
        }
    }

    #[test]
    fn grid_pairs_equal_brute_force_on_random_cloud() {
        for (min_dist, seed) in [(0.05, 11_u64), (0.1, 23), (0.3, 37)] {
            let cloud = random_cloud(500, 0.0, 2.0, seed);
            let expected = brute_force_pairs(&cloud, min_dist);
            let found = find_self_collision_pairs(&cloud, min_dist);
            assert!(
                !expected.is_empty(),
                "the test cloud must contain overlaps (min_dist = {min_dist})"
            );
            assert_eq!(
                found, expected,
                "grid and brute force disagree for min_dist = {min_dist}"
            );
        }
    }

    #[test]
    fn grid_pairs_equal_brute_force_on_regular_lattices_and_stacked_layers() {
        // A cloth is a regular lattice: the worst case for cell boundaries, since
        // neighbouring vertices sit at (almost) exactly one cell from each other.
        let count = 15;
        let spacing = 0.1;
        let lattice: Vec<ClothVertex> = (0..count)
            .flat_map(|j| (0..count).map(move |i| [i as f64 * spacing, 0.0, j as f64 * spacing]))
            .map(|pos| vertex(pos, 1.0, false))
            .collect();
        assert!(
            find_self_collision_pairs(&lattice, 0.0999).is_empty(),
            "vertices one spacing apart must not collide with a smaller distance"
        );
        for min_dist in [0.1001, 0.1415, 0.2001] {
            let expected = brute_force_pairs(&lattice, min_dist);
            assert!(!expected.is_empty(), "min_dist = {min_dist}");
            assert_eq!(
                find_self_collision_pairs(&lattice, min_dist),
                expected,
                "single layer, min_dist = {min_dist}"
            );
        }
        // A folded cloth: a second copy of the lattice lies 0.05 above and
        // 0.05 beside the first one.
        let mut stacked = lattice.clone();
        stacked.extend(
            lattice
                .iter()
                .map(|v| vertex(add(v.pos, [0.05, 0.05, 0.05]), 1.0, false)),
        );
        assert!(
            find_self_collision_pairs(&stacked, 0.0866).is_empty(),
            "the layers are sqrt(3) * 0.05 = 0.0866 apart"
        );
        for min_dist in [0.0867, 0.12, 0.2001] {
            let expected = brute_force_pairs(&stacked, min_dist);
            assert!(!expected.is_empty(), "min_dist = {min_dist}");
            assert_eq!(
                find_self_collision_pairs(&stacked, min_dist),
                expected,
                "stacked layers, min_dist = {min_dist}"
            );
        }
    }

    #[test]
    fn grid_pairs_equal_brute_force_with_negative_and_distant_coordinates() {
        // Cells straddle the origin and the cloud is centred far from it.
        let around_origin = random_cloud(400, -1.0, 1.0, 5);
        assert_eq!(
            find_self_collision_pairs(&around_origin, 0.12),
            brute_force_pairs(&around_origin, 0.12)
        );
        let far_away: Vec<ClothVertex> = random_cloud(400, 0.0, 1.5, 9)
            .into_iter()
            .map(|v| vertex(add(v.pos, [-2.5e4, 7.0e3, -1.0e5]), 1.0, false))
            .collect();
        assert_eq!(
            find_self_collision_pairs(&far_away, 0.1),
            brute_force_pairs(&far_away, 0.1)
        );
    }

    #[test]
    fn pair_exactly_one_cell_apart_is_not_reported_but_just_inside_is() {
        let min_dist = 0.25;
        let inside = [
            vertex([0.0, 0.0, 0.0], 1.0, false),
            vertex([min_dist * 0.999, 0.0, 0.0], 1.0, false),
        ];
        assert_eq!(find_self_collision_pairs(&inside, min_dist), vec![(0, 1)]);
        let outside = [
            vertex([0.0, 0.0, 0.0], 1.0, false),
            vertex([min_dist * 1.001, 0.0, 0.0], 1.0, false),
        ];
        assert!(find_self_collision_pairs(&outside, min_dist).is_empty());
    }

    #[test]
    fn candidate_count_is_far_below_all_pairs_for_a_spread_cloud() {
        let cloud = random_cloud(3000, 0.0, 30.0, 77);
        let mut grid = SelfCollisionGrid::new();
        grid.rebuild(&cloud, 0.5);
        let mut candidates = Vec::new();
        let mut total = 0_usize;
        for vertex in &cloud {
            grid.neighbours_of(vertex.pos, &mut candidates);
            assert!(
                candidates.windows(2).all(|w| w[0] < w[1]),
                "candidates must be strictly ascending (no duplicates)"
            );
            total += candidates.len();
        }
        let all_pairs = cloud.len() * cloud.len();
        assert!(
            total < 50 * cloud.len(),
            "grid visited {total} candidates for {} vertices ({all_pairs} all-pairs)",
            cloud.len()
        );
    }

    #[test]
    fn neighbours_include_the_query_vertex_and_all_vertices_within_one_cell() {
        let cloud = random_cloud(300, 0.0, 3.0, 3);
        let cell = 0.2;
        let mut grid = SelfCollisionGrid::new();
        grid.rebuild(&cloud, cell);
        let mut candidates = Vec::new();
        for (i, v) in cloud.iter().enumerate() {
            grid.neighbours_of(v.pos, &mut candidates);
            assert!(candidates.contains(&i), "vertex {i} must find itself");
            for (j, w) in cloud.iter().enumerate() {
                let d = sub(w.pos, v.pos);
                if dot(d, d).sqrt() < cell {
                    assert!(candidates.contains(&j), "vertex {j} within {cell} of {i}");
                }
            }
        }
    }

    #[test]
    fn resolve_matches_brute_force_sweep_for_isolated_overlaps() {
        // Place 120 well separated pairs (3 cm apart) on a coarse lattice so
        // that no vertex overlaps more than its partner.
        let min_dist = 0.1;
        let mut vertices = Vec::new();
        let mut rng = Lcg(2024);
        for ix in 0..6 {
            for iy in 0..5 {
                for iz in 0..4 {
                    let base = [ix as f64, iy as f64, iz as f64];
                    let angle = rng.next_unit() * std::f64::consts::TAU;
                    let offset = [0.03 * angle.cos(), 0.03 * angle.sin(), 0.0];
                    vertices.push(vertex(base, 1.0, false));
                    vertices.push(vertex(add(base, offset), 2.0, false));
                }
            }
        }
        let mut expected = vertices.clone();
        brute_force_resolve(&mut expected, min_dist);
        let mut resolved = vertices.clone();
        SelfCollisionGrid::new().resolve(&mut resolved, min_dist);
        for (k, (got, want)) in resolved.iter().zip(expected.iter()).enumerate() {
            assert_eq!(
                got.pos, want.pos,
                "vertex {k} differs from the brute-force sweep"
            );
        }
        assert!(find_self_collision_pairs(&resolved, min_dist * (1.0 - 1e-9)).is_empty());
    }

    #[test]
    fn resolve_separates_a_pair_to_min_dist_split_by_inverse_mass() {
        let min_dist = 0.2;
        // Vertex 1 is four times lighter, so it moves four times as far.
        let mut vertices = vec![
            vertex([0.0, 0.0, 0.0], 1.0, false),
            vertex([0.05, 0.0, 0.0], 4.0, false),
        ];
        resolve_with_fresh_grid(&mut vertices, min_dist);
        let d = sub(vertices[1].pos, vertices[0].pos);
        assert!((dot(d, d).sqrt() - min_dist).abs() < 1e-12);
        let moved_0 = vertices[0].pos[0].abs();
        let moved_1 = (vertices[1].pos[0] - 0.05).abs();
        assert!(
            (moved_1 / moved_0 - 4.0).abs() < 1e-9,
            "{moved_0} vs {moved_1}"
        );
        // Linear momentum (mass-weighted displacement) is conserved.
        assert!((moved_0 * 1.0 - moved_1 * 0.25).abs() < 1e-12);
    }

    fn resolve_with_fresh_grid(vertices: &mut [ClothVertex], min_dist: f64) {
        SelfCollisionGrid::new().resolve(vertices, min_dist);
    }

    #[test]
    fn resolve_never_moves_fixed_vertices() {
        let min_dist = 0.2;
        let mut vertices = vec![
            vertex([0.0, 0.0, 0.0], 0.0, true),
            vertex([0.05, 0.0, 0.0], 1.0, false),
            vertex([5.0, 5.0, 5.0], 0.0, true),
            vertex([5.0, 5.0, 5.05], 0.0, true),
        ];
        resolve_with_fresh_grid(&mut vertices, min_dist);
        assert_eq!(vertices[0].pos, [0.0, 0.0, 0.0]);
        assert!((vertices[1].pos[0] - min_dist).abs() < 1e-12);
        assert_eq!(vertices[2].pos, [5.0, 5.0, 5.0]);
        assert_eq!(
            vertices[3].pos,
            [5.0, 5.0, 5.05],
            "two pinned vertices cannot move"
        );
    }

    #[test]
    fn resolve_ignores_exactly_coincident_vertices() {
        let mut vertices = vec![
            vertex([1.0, 2.0, 3.0], 1.0, false),
            vertex([1.0, 2.0, 3.0], 1.0, false),
        ];
        resolve_with_fresh_grid(&mut vertices, 0.5);
        assert_eq!(vertices[0].pos, [1.0, 2.0, 3.0]);
        assert_eq!(vertices[1].pos, [1.0, 2.0, 3.0]);
        assert_eq!(find_self_collision_pairs(&vertices, 0.5), vec![(0, 1)]);
    }

    #[test]
    fn degenerate_inputs_are_safe() {
        let mut empty: Vec<ClothVertex> = Vec::new();
        resolve_with_fresh_grid(&mut empty, 0.1);
        assert!(find_self_collision_pairs(&empty, 0.1).is_empty());

        let mut single = vec![vertex([0.0; 3], 1.0, false)];
        resolve_with_fresh_grid(&mut single, 0.1);
        assert_eq!(single[0].pos, [0.0; 3]);

        let close = vec![
            vertex([0.0, 0.0, 0.0], 1.0, false),
            vertex([0.01, 0.0, 0.0], 1.0, false),
        ];
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                find_self_collision_pairs(&close, bad).is_empty(),
                "min_dist = {bad}"
            );
            let mut copy = close.clone();
            resolve_with_fresh_grid(&mut copy, bad);
            assert_eq!(copy[0].pos, close[0].pos);
            assert_eq!(copy[1].pos, close[1].pos);
        }

        let with_nan = vec![
            vertex([f64::NAN, 0.0, 0.0], 1.0, false),
            vertex([0.0, f64::INFINITY, 0.0], 1.0, false),
            vertex([0.0, 0.0, 0.0], 1.0, false),
            vertex([0.05, 0.0, 0.0], 1.0, false),
        ];
        assert_eq!(find_self_collision_pairs(&with_nan, 0.1), vec![(2, 3)]);
    }

    #[test]
    fn grid_is_reusable_across_calls_with_different_sizes() {
        let mut grid = SelfCollisionGrid::new();
        let big = random_cloud(600, 0.0, 2.0, 1);
        let small = random_cloud(40, 0.0, 1.0, 2);
        for _ in 0..3 {
            assert_eq!(
                grid.overlapping_pairs(&big, 0.1),
                brute_force_pairs(&big, 0.1)
            );
            assert_eq!(
                grid.overlapping_pairs(&small, 0.2),
                brute_force_pairs(&small, 0.2)
            );
        }
    }

    #[test]
    fn cloth_mesh_reuses_its_grid_and_matches_the_stateless_function() {
        use super::super::functions::resolve_self_collision;
        use super::super::types::ClothMesh;

        // Crumple a flat 12 x 12 sheet so that many vertices overlap.
        let mut mesh = ClothMesh::new_grid(12, 12, 1.0, 1.0, 1.0);
        let mut rng = Lcg(99);
        for vertex in mesh.vertices.iter_mut() {
            vertex.pos[0] *= 0.5;
            vertex.pos[1] = 0.03 * rng.next_unit();
            vertex.pos[2] *= 0.5;
        }
        mesh.vertices[0].fixed = true;
        mesh.vertices[0].inv_mass = 0.0;
        let min_dist = 0.1;
        let initial_overlaps = find_self_collision_pairs(&mesh.vertices, min_dist).len();
        assert!(initial_overlaps > 100, "{initial_overlaps} overlaps");

        let original = mesh.vertices.clone();
        let mut reference = mesh.vertices.clone();
        let pinned = mesh.vertices[0].pos;
        for _ in 0..5 {
            mesh.resolve_self_collision(min_dist);
            resolve_self_collision(&mut reference, min_dist);
            for (k, (got, want)) in mesh.vertices.iter().zip(&reference).enumerate() {
                assert_eq!(
                    got.pos, want.pos,
                    "vertex {k} differs from the stateless call"
                );
            }
        }
        assert_eq!(
            mesh.vertices[0].pos, pinned,
            "a pinned vertex must not move"
        );
        let moved = mesh
            .vertices
            .iter()
            .zip(&original)
            .filter(|(now, before)| now.pos != before.pos)
            .count();
        assert!(
            moved > mesh.vertices.len() / 2,
            "most of the crowded vertices must be pushed apart, only {moved} moved"
        );
    }
}
