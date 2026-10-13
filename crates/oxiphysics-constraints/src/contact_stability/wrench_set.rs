// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Clean-up helpers shared by the exact force-closure test and the grasp
//! quality metrics.
//!
//! Both consume a wrench set as a *point set* in `R^6` (force closure and the
//! inscribed-ball radius only depend on the convex hull / positive span of the
//! wrenches), so zero wrenches carry no information and repeated wrenches are
//! redundant.  Removing them keeps the linear program and the facet
//! enumeration small; a friction-less contact, for instance, yields one
//! identical wrench per cone edge.

use super::{Wrench, vec6_len};

/// Wrenches whose Euclidean norm is not larger than this are treated as zero.
pub(super) const MIN_WRENCH_NORM: f64 = 1e-12;

/// Returns `true` when every component of every wrench is finite.
pub(super) fn all_finite(wrenches: &[Wrench]) -> bool {
    wrenches.iter().all(|w| w.iter().all(|c| c.is_finite()))
}

/// Total lexicographic order on wrenches (used to find exact duplicates).
fn lexicographic(a: &Wrench, b: &Wrench) -> std::cmp::Ordering {
    for (x, y) in a.iter().zip(b.iter()) {
        let order = x.total_cmp(y);
        if order != std::cmp::Ordering::Equal {
            return order;
        }
    }
    std::cmp::Ordering::Equal
}

/// Sort a wrench list and drop exact duplicates.
fn sorted_distinct(mut wrenches: Vec<Wrench>) -> Vec<Wrench> {
    wrenches.sort_by(lexicographic);
    wrenches.dedup();
    wrenches
}

/// The distinct non-zero wrenches of `wrenches`.
///
/// Returns `None` if any component is not finite.  The order of the result is
/// unspecified (it is sorted lexicographically, which makes the downstream
/// algorithms deterministic).
pub(super) fn distinct_nonzero(wrenches: &[Wrench]) -> Option<Vec<Wrench>> {
    if !all_finite(wrenches) {
        return None;
    }
    let kept: Vec<Wrench> = wrenches
        .iter()
        .copied()
        .filter(|w| vec6_len(w) > MIN_WRENCH_NORM)
        .collect();
    Some(sorted_distinct(kept))
}

/// The distinct non-zero wrenches of `wrenches`, each scaled to unit length.
///
/// Positive scaling does not change the positive span of a wrench set, so the
/// result has the same force-closure property as the input.  Returns `None`
/// if any component is not finite.
pub(super) fn distinct_unit(wrenches: &[Wrench]) -> Option<Vec<Wrench>> {
    if !all_finite(wrenches) {
        return None;
    }
    let unit: Vec<Wrench> = wrenches
        .iter()
        .filter_map(|w| {
            let len = vec6_len(w);
            if len > MIN_WRENCH_NORM {
                let inv = 1.0 / len;
                Some([
                    w[0] * inv,
                    w[1] * inv,
                    w[2] * inv,
                    w[3] * inv,
                    w[4] * inv,
                    w[5] * inv,
                ])
            } else {
                None
            }
        })
        .collect();
    Some(sorted_distinct(unit))
}
