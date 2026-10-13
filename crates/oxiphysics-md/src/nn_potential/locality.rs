// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cutoff neighbor lists for finite (non-periodic) atomic clusters.
//!
//! The neural-network potentials in this module are local: an atomic energy
//! depends only on the atoms within a cutoff radius of its centre.  Building
//! the neighbor lists once with a spatial hash grid turns every energy and
//! force evaluation from `O(N²)` pair scans into `O(N · k)` work for `k`
//! neighbors per atom.  (The periodic [`crate::neighbor_list::CellList`] is not
//! used here because these potentials act on isolated clusters, and a
//! periodic grid sized to a sparse or very extended cluster would be huge.)

use std::collections::HashMap;

/// One neighbor `j` of a centre atom `i`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Neighbor {
    /// Index of the neighbor atom.
    pub(super) j: usize,
    /// Displacement `r_j - r_i`.
    pub(super) dr: [f64; 3],
    /// Distance `|r_j - r_i|`.
    pub(super) r: f64,
}

/// Integer grid coordinates of a hash-grid cell.
type CellKey = [i64; 3];

/// Neighbors within `radius` (strictly closer than `radius`) of every atom.
///
/// `lists[i]` holds the neighbors of atom `i` sorted by ascending atom index,
/// which fixes the order in which descriptors accumulate their terms.  A
/// non-positive or NaN `radius` yields empty lists; an infinite `radius`
/// lists every other atom.
pub(super) fn neighbor_lists(positions: &[[f64; 3]], radius: f64) -> Vec<Vec<Neighbor>> {
    let n = positions.len();
    let mut lists: Vec<Vec<Neighbor>> = vec![Vec::new(); n];
    if radius.is_nan() || radius <= 0.0 || n < 2 {
        return lists;
    }
    if !radius.is_finite() {
        for (i, list) in lists.iter_mut().enumerate() {
            for (j, pos_j) in positions.iter().enumerate() {
                if j != i {
                    push_if_within(list, positions[i], j, *pos_j, radius);
                }
            }
        }
        return lists;
    }

    let origin = positions.iter().fold([f64::INFINITY; 3], |lo, p| {
        [lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])]
    });
    let key_of = |p: &[f64; 3]| -> CellKey {
        [
            ((p[0] - origin[0]) / radius).floor() as i64,
            ((p[1] - origin[1]) / radius).floor() as i64,
            ((p[2] - origin[2]) / radius).floor() as i64,
        ]
    };
    let mut grid: HashMap<CellKey, Vec<usize>> = HashMap::new();
    for (i, p) in positions.iter().enumerate() {
        grid.entry(key_of(p)).or_default().push(i);
    }

    for (i, p) in positions.iter().enumerate() {
        let centre = key_of(p);
        let mut keys: Vec<CellKey> = Vec::with_capacity(27);
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    keys.push([
                        centre[0].saturating_add(dx),
                        centre[1].saturating_add(dy),
                        centre[2].saturating_add(dz),
                    ]);
                }
            }
        }
        // Saturating arithmetic can only collapse keys for absurdly distant
        // atoms; deduplicate so no cell is scanned twice.
        keys.sort_unstable();
        keys.dedup();
        let list = &mut lists[i];
        for key in &keys {
            if let Some(members) = grid.get(key) {
                for &j in members {
                    if j != i {
                        push_if_within(list, *p, j, positions[j], radius);
                    }
                }
            }
        }
        list.sort_unstable_by_key(|nb| nb.j);
    }
    lists
}

/// Neighbors of the single atom `center` within `radius`, by a linear scan.
///
/// `O(N)`, in ascending atom-index order; use [`neighbor_lists`] when the
/// neighbors of every atom are needed.
///
/// # Panics
/// Panics if `center` is not a valid atom index.
pub(super) fn neighbors_of_atom(
    positions: &[[f64; 3]],
    center: usize,
    radius: f64,
) -> Vec<Neighbor> {
    let mut list = Vec::new();
    if radius.is_nan() || radius <= 0.0 {
        return list;
    }
    let centre = positions[center];
    for (j, pos_j) in positions.iter().enumerate() {
        if j != center {
            push_if_within(&mut list, centre, j, *pos_j, radius);
        }
    }
    list
}

/// Append `j` to `list` if it lies strictly closer than `radius` to `pos_i`.
fn push_if_within(
    list: &mut Vec<Neighbor>,
    pos_i: [f64; 3],
    j: usize,
    pos_j: [f64; 3],
    radius: f64,
) {
    let dr = [
        pos_j[0] - pos_i[0],
        pos_j[1] - pos_i[1],
        pos_j[2] - pos_i[2],
    ];
    let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
    if r < radius {
        list.push(Neighbor { j, dr, r });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{RngExt, SeedableRng};

    fn random_cluster(n: usize, extent: f64, seed: u64) -> Vec<[f64; 3]> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        (0..n)
            .map(|_| {
                [
                    rng.random_range(-extent..extent),
                    rng.random_range(-extent..extent),
                    rng.random_range(-extent..extent),
                ]
            })
            .collect()
    }

    fn brute_force(positions: &[[f64; 3]], radius: f64) -> Vec<Vec<usize>> {
        positions
            .iter()
            .enumerate()
            .map(|(i, pi)| {
                positions
                    .iter()
                    .enumerate()
                    .filter(|&(j, pj)| {
                        let dr = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
                        j != i && (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt() < radius
                    })
                    .map(|(j, _)| j)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn grid_neighbor_lists_match_brute_force() {
        for &(n, extent, radius) in &[(120, 4.0, 1.3), (200, 10.0, 2.5), (60, 1.0, 5.0)] {
            let positions = random_cluster(n, extent, 3);
            let lists = neighbor_lists(&positions, radius);
            let reference = brute_force(&positions, radius);
            for i in 0..n {
                let got: Vec<usize> = lists[i].iter().map(|nb| nb.j).collect();
                assert_eq!(got, reference[i], "atom {i}, radius {radius}");
                for nb in &lists[i] {
                    let dr = [
                        positions[nb.j][0] - positions[i][0],
                        positions[nb.j][1] - positions[i][1],
                        positions[nb.j][2] - positions[i][2],
                    ];
                    assert_eq!(nb.dr, dr);
                    assert!(nb.r < radius);
                }
            }
        }
    }

    #[test]
    fn lists_are_sorted_by_index_and_symmetric() {
        let positions = random_cluster(80, 3.0, 8);
        let lists = neighbor_lists(&positions, 1.7);
        for (i, list) in lists.iter().enumerate() {
            assert!(list.windows(2).all(|w| w[0].j < w[1].j));
            for nb in list {
                assert!(lists[nb.j].iter().any(|back| back.j == i));
            }
        }
    }

    #[test]
    fn degenerate_radii_and_inputs_are_handled() {
        let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [9.0, 0.0, 0.0]];
        assert!(neighbor_lists(&positions, 0.0).iter().all(|l| l.is_empty()));
        assert!(
            neighbor_lists(&positions, -1.0)
                .iter()
                .all(|l| l.is_empty())
        );
        assert!(
            neighbor_lists(&positions, f64::NAN)
                .iter()
                .all(|l| l.is_empty())
        );
        let all = neighbor_lists(&positions, f64::INFINITY);
        assert!(all.iter().all(|l| l.len() == 2));
        assert!(neighbor_lists(&[], 1.0).is_empty());
        assert!(neighbor_lists(&[[0.0; 3]], 1.0)[0].is_empty());
    }

    #[test]
    fn far_apart_clusters_do_not_allocate_a_dense_grid() {
        // A periodic cell grid spanning 1e12 cutoffs would be hopeless; the hash grid is fine.
        let positions = [
            [0.0; 3],
            [0.4, 0.0, 0.0],
            [1.0e12, 0.0, 0.0],
            [1.0e12, 0.3, 0.0],
        ];
        let lists = neighbor_lists(&positions, 1.0);
        let got: Vec<Vec<usize>> = lists
            .iter()
            .map(|l| l.iter().map(|nb| nb.j).collect())
            .collect();
        assert_eq!(got, vec![vec![1], vec![0], vec![3], vec![2]]);
    }

    #[test]
    fn nan_coordinates_do_not_panic_and_have_no_neighbors() {
        let positions = [[f64::NAN, 0.0, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let lists = neighbor_lists(&positions, 1.0);
        assert!(lists[0].is_empty());
        assert_eq!(lists[1].len(), 1);
    }
}
