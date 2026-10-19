// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Non-bonded topology of [`AmberForceField`]: 1-2 / 1-3 exclusions and 1-4
//! pairs derived from the bonded terms.
//!
//! The bond graph decides the treatment of a pair of atoms:
//!
//! * 1-2 (joined by a bond) and 1-3 (joined by two bonds, or the end atoms of
//!   an angle term) are **excluded** from the non-bonded interaction;
//! * 1-4 (separated by three bonds, or the end atoms of a dihedral term) interact
//!   with their energies scaled by [`AmberForceField::scaling_14`];
//! * when a pair is reachable both ways, as in a small ring, the shorter path
//!   wins: a 1-4 pair that is also 1-2 or 1-3 is excluded.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};

use crate::neighbor::ExclusionList;

use super::types::{AmberForceField, AngleTerm, Bond, DihedralTerm};

/// Exclusions and 1-4 pairs in effect for a force evaluation, borrowed from the
/// force field when its stored topology is up to date.
pub(super) type EffectiveTopology<'a> = (Cow<'a, ExclusionList>, Cow<'a, BTreeSet<(usize, usize)>>);

/// Exclusions and 1-4 pairs implied by a set of bonded terms.
#[derive(Debug, Clone, Default)]
pub(super) struct DerivedTopology {
    /// 1-2 and 1-3 pairs.
    pub exclusions: ExclusionList,
    /// 1-4 pairs `(i, j)` with `i < j` that are not excluded.
    pub pairs_14: BTreeSet<(usize, usize)>,
}

/// `(i, j)` ordered so that the smaller index comes first.
fn ordered(i: usize, j: usize) -> (usize, usize) {
    if i < j { (i, j) } else { (j, i) }
}

/// Bond graph as sorted, duplicate-free adjacency lists.
fn adjacency(bonds: &[Bond]) -> BTreeMap<usize, Vec<usize>> {
    let mut neighbours: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for bond in bonds.iter().filter(|bond| bond.i != bond.j) {
        neighbours.entry(bond.i).or_default().push(bond.j);
        neighbours.entry(bond.j).or_default().push(bond.i);
    }
    for list in neighbours.values_mut() {
        list.sort_unstable();
        list.dedup();
    }
    neighbours
}

/// Derive the exclusions and 1-4 pairs implied by the bonded terms.
///
/// Costs `O(Σ_bonds deg_j · deg_k)`, so it stays cheap for large molecules.
pub(super) fn derive_topology(
    bonds: &[Bond],
    angles: &[AngleTerm],
    dihedrals: &[DihedralTerm],
) -> DerivedTopology {
    let neighbours = adjacency(bonds);
    let mut exclusions = ExclusionList::new();
    let mut candidates_14: BTreeSet<(usize, usize)> = BTreeSet::new();

    for (&centre, list) in &neighbours {
        for (position, &first) in list.iter().enumerate() {
            // 1-2: every bond, recorded once from its lower-numbered atom.
            if centre < first {
                exclusions.add_12(centre, first);
            }
            // 1-3: two neighbours of a common atom.
            for &second in &list[position + 1..] {
                exclusions.add_12(first, second);
            }
            // 1-4: a neighbour of `centre` and a neighbour of `first` across
            // the bond centre-first.
            if centre < first {
                let far_ends = neighbours.get(&first).map(Vec::as_slice).unwrap_or(&[]);
                for &near in list.iter().filter(|&&atom| atom != first) {
                    for &far in far_ends.iter().filter(|&&atom| atom != centre) {
                        if near != far {
                            candidates_14.insert(ordered(near, far));
                        }
                    }
                }
            }
        }
    }
    // Explicit angle and dihedral terms define their end atoms as 1-3 / 1-4
    // neighbours even when the corresponding bonds were not listed.
    for angle in angles.iter().filter(|angle| angle.i != angle.k) {
        exclusions.add_12(angle.i, angle.k);
    }
    for dihedral in dihedrals.iter().filter(|dihedral| dihedral.i != dihedral.l) {
        candidates_14.insert(ordered(dihedral.i, dihedral.l));
    }

    let pairs_14 = candidates_14
        .into_iter()
        .filter(|&(i, j)| !exclusions.is_excluded(i, j))
        .collect();
    DerivedTopology {
        exclusions,
        pairs_14,
    }
}

impl AmberForceField {
    /// Numbers of bond, angle and dihedral terms.
    fn bonded_term_counts(&self) -> (usize, usize, usize) {
        (
            self.bonded.bonds.len(),
            self.angle.angles.len(),
            self.dihedral.dihedrals.len(),
        )
    }

    /// Derive the exclusions and 1-4 pairs from the bonded terms.
    ///
    /// 1-2 neighbours (bonds) and 1-3 neighbours (two bonds apart, or the end
    /// atoms of an angle term) are added to [`AmberForceField::exclusions`].
    /// 1-4 neighbours (three bonds apart, or the end atoms of a dihedral term)
    /// are added to [`AmberForceField::pairs_14`] unless they are excluded
    /// already, which happens for ring closures.  Entries that were added by
    /// hand are kept, so the call is idempotent and can be repeated after the
    /// bonded terms change.
    pub fn build_topology(&mut self) {
        let derived = self.derive_from_bonded_terms();
        self.exclusions.pairs.extend(derived.exclusions.pairs);
        let exclusions = &self.exclusions;
        self.pairs_14.extend(
            derived
                .pairs_14
                .into_iter()
                .filter(|&(i, j)| !exclusions.is_excluded(i, j)),
        );
        self.topology_term_counts = Some(self.bonded_term_counts());
    }

    /// Whether [`AmberForceField::build_topology`] has seen the current number
    /// of bonds, angles and dihedrals.
    ///
    /// Changing a term in place (rather than adding or removing one) is not
    /// detected; call [`AmberForceField::build_topology`] again after such an
    /// edit.
    pub fn topology_is_current(&self) -> bool {
        self.topology_term_counts == Some(self.bonded_term_counts())
    }

    /// Topology implied by the present bonded terms.
    fn derive_from_bonded_terms(&self) -> DerivedTopology {
        derive_topology(
            &self.bonded.bonds,
            &self.angle.angles,
            &self.dihedral.dihedrals,
        )
    }

    /// The exclusions and 1-4 pairs that the non-bonded interaction uses: the
    /// stored ones, extended by those derived from the bonded terms if these
    /// changed since the last [`AmberForceField::build_topology`] call.
    ///
    /// Borrows the stored lists in the usual, up-to-date case.
    pub(super) fn effective_topology(&self) -> EffectiveTopology<'_> {
        if self.topology_is_current() {
            return (
                Cow::Borrowed(&self.exclusions),
                Cow::Borrowed(&self.pairs_14),
            );
        }
        let derived = self.derive_from_bonded_terms();
        let mut exclusions = self.exclusions.clone();
        exclusions.pairs.extend(derived.exclusions.pairs);
        let mut pairs_14 = self.pairs_14.clone();
        pairs_14.extend(derived.pairs_14);
        (Cow::Owned(exclusions), Cow::Owned(pairs_14))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain_bonds(n_atoms: usize) -> Vec<Bond> {
        (0..n_atoms - 1)
            .map(|i| Bond {
                i,
                j: i + 1,
                k: 100.0,
                r0: 1.5,
            })
            .collect()
    }

    fn dihedral(i: usize, j: usize, k: usize, l: usize) -> DihedralTerm {
        DihedralTerm {
            i,
            j,
            k,
            l,
            v_n: 1.0,
            gamma: 0.0,
            n: 3,
        }
    }

    fn force_field_with_chain(n_atoms: usize) -> AmberForceField {
        let mut ff = AmberForceField::new();
        for bond in chain_bonds(n_atoms) {
            ff.bonded.add_bond(bond.i, bond.j, bond.k, bond.r0);
        }
        ff
    }

    #[test]
    fn butane_chain_has_five_exclusions_and_one_one_four_pair() {
        let ff_topology = derive_topology(&chain_bonds(4), &[], &[]);
        let mut excluded: Vec<(usize, usize)> =
            ff_topology.exclusions.pairs.iter().copied().collect();
        excluded.sort_unstable();
        assert_eq!(excluded, vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]);
        let pairs_14: Vec<(usize, usize)> = ff_topology.pairs_14.iter().copied().collect();
        assert_eq!(pairs_14, vec![(0, 3)]);
    }

    #[test]
    fn pentane_chain_classifies_every_pair_by_topological_distance() {
        let topology = derive_topology(&chain_bonds(5), &[], &[]);
        for i in 0..5_usize {
            for j in (i + 1)..5 {
                let separation = j - i;
                assert_eq!(
                    topology.exclusions.is_excluded(i, j),
                    separation <= 2,
                    "({i}, {j}) exclusion"
                );
                assert_eq!(
                    topology.pairs_14.contains(&(i, j)),
                    separation == 3,
                    "({i}, {j}) 1-4 membership"
                );
            }
        }
    }

    #[test]
    fn ring_closure_pairs_are_excluded_rather_than_scaled() {
        // Cyclobutane: the dihedral 0-1-2-3 has bonded end atoms (3-0).
        let mut bonds = chain_bonds(4);
        bonds.push(Bond {
            i: 3,
            j: 0,
            k: 100.0,
            r0: 1.5,
        });
        let topology = derive_topology(&bonds, &[], &[dihedral(0, 1, 2, 3)]);
        assert!(topology.pairs_14.is_empty(), "{:?}", topology.pairs_14);
        // All six pairs of the four ring atoms are 1-2 or 1-3 neighbours.
        assert_eq!(topology.exclusions.len(), 6);
    }

    #[test]
    fn branched_molecule_has_one_four_pairs_through_every_branch() {
        // Isobutane-like skeleton: centre 1 carries atoms 0, 2, 3; atom 3 carries 4.
        let bonds = vec![
            Bond {
                i: 0,
                j: 1,
                k: 1.0,
                r0: 1.0,
            },
            Bond {
                i: 1,
                j: 2,
                k: 1.0,
                r0: 1.0,
            },
            Bond {
                i: 1,
                j: 3,
                k: 1.0,
                r0: 1.0,
            },
            Bond {
                i: 3,
                j: 4,
                k: 1.0,
                r0: 1.0,
            },
        ];
        let topology = derive_topology(&bonds, &[], &[]);
        let pairs_14: Vec<(usize, usize)> = topology.pairs_14.iter().copied().collect();
        assert_eq!(pairs_14, vec![(0, 4), (2, 4)]);
        // 1-3: all pairs of the centre's neighbours plus 1-4 (centre, 4).
        for pair in [(0, 2), (0, 3), (2, 3), (1, 4)] {
            assert!(topology.exclusions.is_excluded(pair.0, pair.1), "{pair:?}");
        }
        assert_eq!(topology.exclusions.len(), 4 + 4);
    }

    #[test]
    fn angle_and_dihedral_terms_define_their_end_atoms_without_bonds() {
        let angle = AngleTerm {
            i: 5,
            j: 6,
            k: 7,
            k_theta: 10.0,
            theta0: 1.9,
        };
        let topology = derive_topology(&[], &[angle], &[dihedral(10, 11, 12, 13)]);
        assert!(topology.exclusions.is_excluded(5, 7));
        assert_eq!(topology.exclusions.len(), 1);
        let pairs_14: Vec<(usize, usize)> = topology.pairs_14.iter().copied().collect();
        assert_eq!(pairs_14, vec![(10, 13)]);
    }

    #[test]
    fn degenerate_and_duplicate_terms_are_harmless() {
        let bonds = vec![
            Bond {
                i: 2,
                j: 2,
                k: 1.0,
                r0: 1.0,
            },
            Bond {
                i: 0,
                j: 1,
                k: 1.0,
                r0: 1.0,
            },
            Bond {
                i: 1,
                j: 0,
                k: 1.0,
                r0: 1.0,
            },
        ];
        let topology = derive_topology(&bonds, &[], &[dihedral(4, 5, 6, 4)]);
        assert_eq!(topology.exclusions.len(), 1);
        assert!(topology.exclusions.is_excluded(0, 1));
        assert!(topology.pairs_14.is_empty());
    }

    #[test]
    fn build_topology_is_idempotent_and_keeps_manual_entries() {
        let mut ff = force_field_with_chain(5);
        ff.exclusions.add_12(10, 11);
        ff.pairs_14.insert((20, 21));
        assert!(!ff.topology_is_current());
        ff.build_topology();
        assert!(ff.topology_is_current());
        let exclusions_once = ff.exclusions.clone();
        let pairs_once = ff.pairs_14.clone();
        ff.build_topology();
        assert_eq!(ff.exclusions.pairs, exclusions_once.pairs);
        assert_eq!(ff.pairs_14, pairs_once);
        assert!(ff.exclusions.is_excluded(10, 11), "manual exclusion kept");
        assert!(ff.pairs_14.contains(&(20, 21)), "manual 1-4 pair kept");
        assert!(ff.pairs_14.contains(&(0, 3)) && ff.pairs_14.contains(&(1, 4)));
    }

    #[test]
    fn topology_becomes_stale_when_bonded_terms_are_added() {
        let mut ff = force_field_with_chain(4);
        ff.build_topology();
        assert!(ff.topology_is_current());
        ff.bonded.add_bond(3, 4, 100.0, 1.5);
        assert!(!ff.topology_is_current());
        let (exclusions, pairs_14) = ff.effective_topology();
        assert!(
            exclusions.is_excluded(3, 4),
            "new bond is excluded on the fly"
        );
        assert!(pairs_14.contains(&(1, 4)), "new 1-4 pair found on the fly");
        // The stored lists are untouched until build_topology runs again.
        assert!(!ff.exclusions.is_excluded(3, 4));
        ff.build_topology();
        assert!(ff.exclusions.is_excluded(3, 4));
        assert!(ff.topology_is_current());
    }
}
