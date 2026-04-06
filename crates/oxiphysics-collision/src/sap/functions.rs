//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg, clippy::type_complexity)]
use std::collections::HashSet;

use super::types::{Aabb3, IncrementalSap, SapEndpointU32, SapStats};

/// Update the AABB for a body and propagate to all axis lists in one call.
#[allow(dead_code)]
pub fn propagate_aabb_update(sap: &mut IncrementalSap, body_id: u32, new_aabb: Aabb3) {
    sap.update(body_id, new_aabb);
}
/// Move a body by the given displacement and update the SAP.
#[allow(dead_code)]
pub fn translate_body(sap: &mut IncrementalSap, body_id: u32, delta: [f64; 3]) {
    if let Some(old_aabb) = sap.aabbs.get(&body_id).cloned() {
        let new_aabb = Aabb3 {
            min: [
                old_aabb.min[0] + delta[0],
                old_aabb.min[1] + delta[1],
                old_aabb.min[2] + delta[2],
            ],
            max: [
                old_aabb.max[0] + delta[0],
                old_aabb.max[1] + delta[1],
                old_aabb.max[2] + delta[2],
            ],
        };
        sap.update(body_id, new_aabb);
    }
}
/// Expand the AABB of a body by a margin on each side and update the SAP.
#[allow(dead_code)]
pub fn expand_aabb(sap: &mut IncrementalSap, body_id: u32, margin: f64) {
    if let Some(old_aabb) = sap.aabbs.get(&body_id).cloned() {
        let new_aabb = Aabb3 {
            min: [
                old_aabb.min[0] - margin,
                old_aabb.min[1] - margin,
                old_aabb.min[2] - margin,
            ],
            max: [
                old_aabb.max[0] + margin,
                old_aabb.max[1] + margin,
                old_aabb.max[2] + margin,
            ],
        };
        sap.update(body_id, new_aabb);
    }
}
/// Move a single endpoint left or right using bubble-sort steps until the list
/// is sorted again.  Efficient when the endpoint moved only a few positions.
///
/// Returns the number of swaps performed.
#[allow(dead_code)]
pub fn bubble_sort_endpoint(endpoints: &mut Vec<SapEndpointU32>, idx: usize) -> usize {
    let n = endpoints.len();
    if n == 0 {
        return 0;
    }
    let mut swaps = 0usize;
    let mut i = idx;
    while i > 0 && endpoints[i - 1].value > endpoints[i].value {
        endpoints.swap(i - 1, i);
        i -= 1;
        swaps += 1;
    }
    let start = if swaps > 0 { idx } else { i };
    let mut j = start;
    while j + 1 < n && endpoints[j].value > endpoints[j + 1].value {
        endpoints.swap(j, j + 1);
        j += 1;
        swaps += 1;
    }
    swaps
}
/// Find the index of the min-endpoint for `body_id` in a sorted endpoint list.
/// Returns `None` if not found.
#[allow(dead_code)]
pub fn find_min_endpoint(endpoints: &[SapEndpointU32], body_id: u32) -> Option<usize> {
    endpoints
        .iter()
        .position(|e| e.body_id == body_id && e.is_min)
}
/// Find the index of the max-endpoint for `body_id` in a sorted endpoint list.
/// Returns `None` if not found.
#[allow(dead_code)]
pub fn find_max_endpoint(endpoints: &[SapEndpointU32], body_id: u32) -> Option<usize> {
    endpoints
        .iter()
        .position(|e| e.body_id == body_id && !e.is_min)
}
/// Compute the union (merge) of two AABBs.
#[allow(dead_code)]
pub fn aabb_union(a: &Aabb3, b: &Aabb3) -> Aabb3 {
    Aabb3 {
        min: [
            a.min[0].min(b.min[0]),
            a.min[1].min(b.min[1]),
            a.min[2].min(b.min[2]),
        ],
        max: [
            a.max[0].max(b.max[0]),
            a.max[1].max(b.max[1]),
            a.max[2].max(b.max[2]),
        ],
    }
}
/// Compute the intersection of two AABBs.  Returns `None` if they do not overlap.
#[allow(dead_code)]
pub fn aabb_intersection(a: &Aabb3, b: &Aabb3) -> Option<Aabb3> {
    let min = [
        a.min[0].max(b.min[0]),
        a.min[1].max(b.min[1]),
        a.min[2].max(b.min[2]),
    ];
    let max = [
        a.max[0].min(b.max[0]),
        a.max[1].min(b.max[1]),
        a.max[2].min(b.max[2]),
    ];
    if min[0] <= max[0] && min[1] <= max[1] && min[2] <= max[2] {
        Some(Aabb3 { min, max })
    } else {
        None
    }
}
/// Surface area of an AABB (used for SAH cost heuristics).
#[allow(dead_code)]
pub fn aabb_surface_area(a: &Aabb3) -> f64 {
    let dx = a.max[0] - a.min[0];
    let dy = a.max[1] - a.min[1];
    let dz = a.max[2] - a.min[2];
    2.0 * (dx * dy + dy * dz + dz * dx)
}
/// Volume of an AABB.
#[allow(dead_code)]
pub fn aabb_volume(a: &Aabb3) -> f64 {
    let dx = (a.max[0] - a.min[0]).max(0.0);
    let dy = (a.max[1] - a.min[1]).max(0.0);
    let dz = (a.max[2] - a.min[2]).max(0.0);
    dx * dy * dz
}
/// Whether point `p` is inside AABB `a`.
#[allow(dead_code)]
pub fn aabb_contains_point(a: &Aabb3, p: [f64; 3]) -> bool {
    p[0] >= a.min[0]
        && p[0] <= a.max[0]
        && p[1] >= a.min[1]
        && p[1] <= a.max[1]
        && p[2] >= a.min[2]
        && p[2] <= a.max[2]
}
/// Whether AABB `a` fully contains AABB `b`.
#[allow(dead_code)]
pub fn aabb_contains_aabb(a: &Aabb3, b: &Aabb3) -> bool {
    a.min[0] <= b.min[0]
        && a.max[0] >= b.max[0]
        && a.min[1] <= b.min[1]
        && a.max[1] >= b.max[1]
        && a.min[2] <= b.min[2]
        && a.max[2] >= b.max[2]
}
/// Extend an AABB by `margin` on each side.
#[allow(dead_code)]
pub fn aabb_pad(a: &Aabb3, margin: f64) -> Aabb3 {
    Aabb3 {
        min: [a.min[0] - margin, a.min[1] - margin, a.min[2] - margin],
        max: [a.max[0] + margin, a.max[1] + margin, a.max[2] + margin],
    }
}
/// Centre of an AABB.
#[allow(dead_code)]
pub fn aabb_center(a: &Aabb3) -> [f64; 3] {
    [
        (a.min[0] + a.max[0]) * 0.5,
        (a.min[1] + a.max[1]) * 0.5,
        (a.min[2] + a.max[2]) * 0.5,
    ]
}
/// Half-extents of an AABB.
#[allow(dead_code)]
pub fn aabb_half_extents(a: &Aabb3) -> [f64; 3] {
    [
        (a.max[0] - a.min[0]) * 0.5,
        (a.max[1] - a.min[1]) * 0.5,
        (a.max[2] - a.min[2]) * 0.5,
    ]
}
/// Report the number of endpoints on each axis and the current pair count.
#[allow(dead_code)]
pub fn sap_endpoint_stats(sap: &IncrementalSap) -> SapStats {
    let endpoint_count = sap.endpoints_x.len();
    SapStats {
        pair_count: sap.active_pairs.len(),
        sweep_count: endpoint_count,
        body_count: sap.body_count(),
    }
}
/// Return `true` if the endpoint list is non-decreasingly sorted by value.
#[allow(dead_code)]
pub fn axis_is_sorted(endpoints: &[SapEndpointU32]) -> bool {
    endpoints.windows(2).all(|w| w[0].value <= w[1].value)
}
/// Count the number of open (active) bodies at a given sweep position.
///
/// Sweeps from left to right and counts bodies whose interval is currently
/// open at `pos`.
#[allow(dead_code)]
pub fn count_active_at(endpoints: &[SapEndpointU32], pos: f64) -> usize {
    let mut active = 0usize;
    for ep in endpoints {
        if ep.value > pos {
            break;
        }
        if ep.is_min {
            active += 1;
        } else {
            active = active.saturating_sub(1);
        }
    }
    active
}
/// Return the index range `[lo, hi)` of endpoints whose `value` falls in
/// `[range_min, range_max]` using binary search (assumes sorted input).
#[allow(dead_code)]
pub fn endpoint_range(
    endpoints: &[SapEndpointU32],
    range_min: f64,
    range_max: f64,
) -> (usize, usize) {
    let lo = endpoints.partition_point(|e| e.value < range_min);
    let hi = endpoints.partition_point(|e| e.value <= range_max);
    (lo, hi)
}
/// Bipartite 3-axis sweep-and-prune between two AABB sets A and B.
///
/// Builds per-axis endpoint lists containing *both* sets, sweeps each axis
/// independently, and intersects the three sets to return only cross-set pairs
/// `(a, b)` where `a ∈ set_a_ids` and `b ∈ set_b_ids`.
///
/// All AABBs are provided in the `aabbs` map (merged from both sets).
#[allow(dead_code)]
pub fn bipartite_sap_query(
    set_a_ids: &[u32],
    set_b_ids: &[u32],
    aabbs: &std::collections::HashMap<u32, Aabb3>,
) -> Vec<(u32, u32)> {
    let set_a_hs: std::collections::HashSet<u32> = set_a_ids.iter().copied().collect();
    let set_b_hs: std::collections::HashSet<u32> = set_b_ids.iter().copied().collect();
    let all_ids: Vec<u32> = set_a_ids.iter().chain(set_b_ids.iter()).copied().collect();
    let mut eps_x: Vec<SapEndpointU32> = Vec::new();
    let mut eps_y: Vec<SapEndpointU32> = Vec::new();
    let mut eps_z: Vec<SapEndpointU32> = Vec::new();
    for &id in &all_ids {
        if let Some(aabb) = aabbs.get(&id) {
            eps_x.push(SapEndpointU32 {
                value: aabb.min[0],
                body_id: id,
                is_min: true,
            });
            eps_x.push(SapEndpointU32 {
                value: aabb.max[0],
                body_id: id,
                is_min: false,
            });
            eps_y.push(SapEndpointU32 {
                value: aabb.min[1],
                body_id: id,
                is_min: true,
            });
            eps_y.push(SapEndpointU32 {
                value: aabb.max[1],
                body_id: id,
                is_min: false,
            });
            eps_z.push(SapEndpointU32 {
                value: aabb.min[2],
                body_id: id,
                is_min: true,
            });
            eps_z.push(SapEndpointU32 {
                value: aabb.max[2],
                body_id: id,
                is_min: false,
            });
        }
    }
    let pairs_x = IncrementalSap::sort_and_sweep_axis(&mut eps_x);
    let pairs_y = IncrementalSap::sort_and_sweep_axis(&mut eps_y);
    let pairs_z = IncrementalSap::sort_and_sweep_axis(&mut eps_z);
    let all_pairs: std::collections::HashSet<(u32, u32)> = pairs_x
        .intersection(&pairs_y)
        .copied()
        .collect::<std::collections::HashSet<_>>()
        .intersection(&pairs_z)
        .copied()
        .collect();
    let mut result: Vec<(u32, u32)> = all_pairs
        .into_iter()
        .filter(|&(a, b)| {
            (set_a_hs.contains(&a) && set_b_hs.contains(&b))
                || (set_a_hs.contains(&b) && set_b_hs.contains(&a))
        })
        .map(|(a, b)| if a < b { (a, b) } else { (b, a) })
        .collect();
    result.sort_unstable();
    result.dedup();
    result
}
/// Insert a single `SapEndpointU32` into an already-sorted slice using
/// binary-search for the insertion point, then shift elements right.
///
/// Returns the index at which the element was inserted.
#[allow(dead_code)]
pub fn sorted_insert(endpoints: &mut Vec<SapEndpointU32>, ep: SapEndpointU32) -> usize {
    let pos = endpoints
        .partition_point(|e| e.value < ep.value || (e.value == ep.value && e.is_min && !ep.is_min));
    endpoints.insert(pos, ep);
    pos
}
/// Return `true` if two `Aabb3` objects overlap on all three axes.
#[allow(dead_code)]
pub fn aabb3_overlaps(a: &Aabb3, b: &Aabb3) -> bool {
    for ax in 0..3 {
        if a.max[ax] < b.min[ax] || b.max[ax] < a.min[ax] {
            return false;
        }
    }
    true
}
/// Expand `a` to also cover `b` (in-place union).
#[allow(dead_code)]
pub fn aabb3_expand_to_include(a: &mut Aabb3, b: &Aabb3) {
    for ax in 0..3 {
        if b.min[ax] < a.min[ax] {
            a.min[ax] = b.min[ax];
        }
        if b.max[ax] > a.max[ax] {
            a.max[ax] = b.max[ax];
        }
    }
}
/// Compute the squared distance from point `p` to the surface of `a`
/// (0.0 if the point is inside).
#[allow(dead_code)]
pub fn aabb3_point_dist_sq(a: &Aabb3, p: [f64; 3]) -> f64 {
    let mut sq = 0.0_f64;
    for ax in 0..3 {
        let d = if p[ax] < a.min[ax] {
            a.min[ax] - p[ax]
        } else if p[ax] > a.max[ax] {
            p[ax] - a.max[ax]
        } else {
            0.0
        };
        sq += d * d;
    }
    sq
}
/// Return all body-ids whose X-axis interval overlaps the range `[lo, hi]`.
///
/// The endpoint list must be sorted before calling this.
#[allow(dead_code)]
pub fn sweep_window_query(endpoints: &[SapEndpointU32], lo: f64, hi: f64) -> Vec<u32> {
    let mut result = Vec::new();
    let mut active: Vec<u32> = Vec::new();
    for ep in endpoints {
        if ep.value > hi {
            break;
        }
        if ep.is_min {
            if ep.value <= hi {
                active.push(ep.body_id);
            }
        } else {
            if ep.value >= lo && active.contains(&ep.body_id) {
                result.push(ep.body_id);
            }
            active.retain(|&id| id != ep.body_id);
        }
    }
    for id in active {
        if !result.contains(&id) {
            result.push(id);
        }
    }
    result.sort_unstable();
    result.dedup();
    result
}
/// Compute the set of pairs that changed between two consecutive queries.
///
/// Returns `(new_pairs, removed_pairs)` by diffing `prev` and `current`.
#[allow(dead_code)]
pub fn pair_delta(
    prev: &HashSet<(u32, u32)>,
    current: &HashSet<(u32, u32)>,
) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
    let mut new_pairs: Vec<(u32, u32)> = current.difference(prev).copied().collect();
    let mut removed_pairs: Vec<(u32, u32)> = prev.difference(current).copied().collect();
    new_pairs.sort_unstable();
    removed_pairs.sort_unstable();
    (new_pairs, removed_pairs)
}
/// Clamp an AABB to fit within a world-space `world_min`/`world_max` boundary.
#[allow(dead_code)]
pub fn aabb3_clamp(a: &Aabb3, world_min: [f64; 3], world_max: [f64; 3]) -> Aabb3 {
    let mut out = a.clone();
    for ax in 0..3 {
        out.min[ax] = out.min[ax].max(world_min[ax]);
        out.max[ax] = out.max[ax].min(world_max[ax]);
        if out.min[ax] > out.max[ax] {
            out.max[ax] = out.min[ax];
        }
    }
    out
}
