// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hybrid GPU LBVH (Linear Bounding Volume Hierarchy) builder.
//!
//! This module ships the *hybrid* construction path: the expensive Morton-code
//! sort runs on the GPU radix sorter ([`crate::gpu_radix::radix_sort_pairs_gpu`]),
//! and the hierarchy is then assembled with the CPU reference builder
//! [`crate::bvh::compute_bvh_from_sorted`]. The result is correct *and*
//! GPU-accelerated for the sort phase, which is the dominant cost for large
//! scenes.
//!
//! # Bit-exact CPU parity
//!
//! The GPU pair-sort is **stable**: when two primitives share a Morton code it
//! preserves their input order. Because the payload carried alongside each key
//! is the primitive's *original index* (`0..n`), the emitted order is exactly
//! "sort by `(morton, original_index)`". Feeding that order into
//! [`compute_bvh_from_sorted`][crate::bvh::compute_bvh_from_sorted] therefore
//! yields a tree *bit-identical* to one produced by a CPU
//! `sort_by_key(|p| (p.morton, original_index))` followed by the same builder.
//! This holds whether the radix sort executed on the GPU or transparently fell
//! back to the CPU oracle (both are stable).
//!
//! The Morton codes themselves come from the canonical
//! [`LbvhPrimitive::new`][crate::bvh::LbvhPrimitive::new], so they match the CPU
//! BVH path exactly.
//!
//! A follow-on pass may replace the CPU hierarchy step with a fully on-GPU
//! Karras radix-tree kernel; this module intentionally ships the correct hybrid
//! first.

use crate::bvh::{Aabb, Bvh, BvhPrimitive, LbvhPrimitive, compute_bvh_from_sorted};
use crate::gpu_radix::radix_sort_pairs_gpu;

/// Build an LBVH whose Morton-code sort runs on the GPU radix sorter, then
/// assembles the hierarchy with the CPU [`compute_bvh_from_sorted`] reference.
///
/// The construction is split into two phases:
///
/// 1. [`gpu_morton_sort`] computes the scene AABB, derives canonical 30-bit
///    Morton codes via [`LbvhPrimitive::new`], and stably sorts the primitives
///    into Morton order on the GPU (carrying each primitive's original index as
///    the sort payload).
/// 2. [`compute_bvh_from_sorted`] builds the HLBVH hierarchy from the sorted
///    slice (`LEAF_SIZE == 4`, highest-differing-bit splits).
///
/// Because the GPU pair-sort is stable, equal-Morton primitives keep ascending
/// original-index order — i.e. the output is exactly the order of a CPU
/// `sort_by_key(|p| (p.morton, original_index))`. The resulting [`Bvh`] is
/// therefore bit-identical to the CPU reference path (see the module docs).
///
/// # Edge cases
///
/// * An empty input returns an empty [`Bvh`] (`root.is_none()`).
/// * A single primitive is routed through the same path with a degenerate scene
///   AABB equal to its own AABB.
///
/// Neither edge case panics.
pub fn gpu_lbvh_build(prims: Vec<BvhPrimitive>) -> Bvh {
    let sorted = gpu_morton_sort(&prims);
    compute_bvh_from_sorted(&sorted)
}

/// GPU-sort primitives into Morton order, stable by original index.
///
/// Returns the primitives as [`LbvhPrimitive`]s (carrying their canonical Morton
/// codes) ordered exactly as a CPU `sort_by_key(|p| (p.morton, original_index))`
/// would order them. Useful for parity tests against a CPU sort of the same
/// `(morton, index)` keys, and as the first phase of [`gpu_lbvh_build`].
///
/// The sort itself runs on [`radix_sort_pairs_gpu`]; if no GPU adapter is
/// available it transparently falls back to the CPU radix oracle (also stable),
/// so the returned order is identical either way.
pub fn gpu_morton_sort(prims: &[BvhPrimitive]) -> Vec<LbvhPrimitive> {
    let n = prims.len();
    if n == 0 {
        return Vec::new();
    }

    // Scene AABB: exact CPU fold over `Aabb::merge`. `LbvhPrimitive::new`
    // clamps each scene extent with 1e-10, so a single-point / zero-extent
    // scene is handled without division by zero.
    let mut scene = prims[0].aabb.clone();
    for p in &prims[1..] {
        scene = Aabb::merge(&scene, &p.aabb);
    }

    // Canonical Morton codes from the centroid normalised by the scene AABB.
    let lbvh_prims: Vec<LbvhPrimitive> = prims
        .iter()
        .map(|p| LbvhPrimitive::new(p.aabb.clone(), p.object_id, &scene))
        .collect();

    // When the primitive count exceeds `u32::MAX` the original index cannot be
    // carried as a `u32` payload; fall back to a stable CPU sort by
    // `(morton, original_index)`, which yields the same order as the GPU path.
    if n > u32::MAX as usize {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|&i| (lbvh_prims[i].morton, i));
        return order.into_iter().map(|i| lbvh_prims[i].clone()).collect();
    }

    let keys: Vec<u32> = lbvh_prims.iter().map(|lp| lp.morton).collect();
    let payload: Vec<u32> = (0..n as u32).collect();

    // Stable GPU radix sort: keys ascending, payload (original index) follows.
    let (_sorted_mortons, sorted_indices) = radix_sort_pairs_gpu(&keys, &payload);

    // Reorder the primitives by the sorted original indices. Each
    // `sorted_indices[k]` is a valid index into `lbvh_prims` (it is a stable
    // permutation of `0..n`), so the lookup never goes out of bounds.
    sorted_indices
        .iter()
        .map(|&i| lbvh_prims[i as usize].clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bvh::BvhNode;

    // ── Deterministic pseudo-random scene generation ─────────────────────────

    /// A tiny SplitMix64-style PRNG: deterministic, no external deps.
    struct DetRng {
        state: u64,
    }

    impl DetRng {
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

        /// Uniform `f32` in `[0, span)`.
        fn next_f32(&mut self, span: f32) -> f32 {
            let frac = (self.next_u64() >> 40) as f32 / (1u64 << 24) as f32;
            frac * span
        }
    }

    /// Build `n` small boxes at deterministic pseudo-random positions inside a
    /// `[0, extent]^3` scene. `object_id` is the original index.
    fn make_random_primitives(n: usize, extent: f32, seed: u64) -> Vec<BvhPrimitive> {
        let mut rng = DetRng::new(seed);
        (0..n)
            .map(|i| {
                let x = rng.next_f32(extent);
                let y = rng.next_f32(extent);
                let z = rng.next_f32(extent);
                let s = 0.25_f32;
                BvhPrimitive::new(Aabb::new([x, y, z], [x + s, y + s, z + s]), i)
            })
            .collect()
    }

    /// CPU reference order: indices sorted by `(morton, original_index)` using
    /// the canonical [`LbvhPrimitive::new`] Morton code.
    fn cpu_reference_sorted(prims: &[BvhPrimitive]) -> Vec<LbvhPrimitive> {
        if prims.is_empty() {
            return Vec::new();
        }
        let mut scene = prims[0].aabb.clone();
        for p in &prims[1..] {
            scene = Aabb::merge(&scene, &p.aabb);
        }
        let lbvh_prims: Vec<LbvhPrimitive> = prims
            .iter()
            .map(|p| LbvhPrimitive::new(p.aabb.clone(), p.object_id, &scene))
            .collect();
        let mut order: Vec<usize> = (0..lbvh_prims.len()).collect();
        order.sort_by_key(|&i| (lbvh_prims[i].morton, i));
        order.into_iter().map(|i| lbvh_prims[i].clone()).collect()
    }

    /// Collect the set of leaf `object_id`s, mapping leaf primitive indices
    /// through `bvh.primitives` (which `compute_bvh_from_sorted` stores in
    /// sorted order). Returned sorted for stable comparison.
    fn leaf_object_ids(bvh: &Bvh) -> Vec<usize> {
        fn walk(node: &BvhNode, bvh: &Bvh, out: &mut Vec<usize>) {
            if node.is_leaf() {
                for &idx in &node.primitives {
                    out.push(bvh.primitives[idx].object_id);
                }
            } else {
                if let Some(l) = &node.left {
                    walk(l, bvh, out);
                }
                if let Some(r) = &node.right {
                    walk(r, bvh, out);
                }
            }
        }
        let mut out = Vec::new();
        if let Some(root) = &bvh.root {
            walk(root, bvh, &mut out);
        }
        out.sort_unstable();
        out
    }

    // ── Edge cases ───────────────────────────────────────────────────────────

    #[test]
    fn empty_no_panic() {
        let bvh = gpu_lbvh_build(vec![]);
        assert!(bvh.root.is_none());
        assert_eq!(bvh.node_count(), 0);
        assert_eq!(bvh.depth(), 0);
    }

    #[test]
    fn single_no_panic() {
        let prims = vec![BvhPrimitive::new(
            Aabb::new([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]),
            42,
        )];
        let bvh = gpu_lbvh_build(prims);
        assert!(bvh.root.is_some());
        assert_eq!(bvh.primitives.len(), 1);
        assert_eq!(bvh.primitives[0].object_id, 42);
    }

    // ── GPU↔CPU Morton-sort parity ───────────────────────────────────────────

    #[test]
    fn gpu_morton_sort_matches_cpu_sort() {
        let prims = make_random_primitives(500, 100.0, 0xC001_CAFE);

        let gpu_sorted = gpu_morton_sort(&prims);
        let cpu_sorted = cpu_reference_sorted(&prims);

        assert_eq!(gpu_sorted.len(), cpu_sorted.len());
        assert_eq!(gpu_sorted.len(), 500);

        // The (morton, object_id) sequence must be identical. This is the
        // GPU↔CPU equivalence assertion: it holds whether the radix ran on the
        // GPU or fell back to the CPU oracle (both stable).
        for (k, (g, c)) in gpu_sorted.iter().zip(cpu_sorted.iter()).enumerate() {
            assert_eq!(
                (g.morton, g.object_id),
                (c.morton, c.object_id),
                "mismatch at sorted position {k}"
            );
        }

        // Mortons must be non-decreasing.
        for w in gpu_sorted.windows(2) {
            assert!(
                w[0].morton <= w[1].morton,
                "mortons not sorted at {} -> {}",
                w[0].morton,
                w[1].morton
            );
        }
    }

    // ── Hybrid topology parity vs the CPU reference builder ───────────────────

    #[test]
    fn hybrid_topology_matches_reference() {
        let prims = make_random_primitives(2000, 100.0, 0xBADC_0FFE);

        // Reference: CPU sort by (morton, index) then the same builder.
        let reference = compute_bvh_from_sorted(&cpu_reference_sorted(&prims));
        // GPU hybrid: GPU sort then the same builder.
        let gpu = gpu_lbvh_build(prims.clone());

        assert_eq!(
            reference.node_count(),
            gpu.node_count(),
            "node_count mismatch"
        );
        assert_eq!(reference.depth(), gpu.depth(), "depth mismatch");
        assert_eq!(reference.primitives.len(), gpu.primitives.len());

        // Identical leaf object_id sets.
        assert_eq!(
            leaf_object_ids(&reference),
            leaf_object_ids(&gpu),
            "leaf object_id sets differ"
        );

        // Identical hit sets along a sweep of rays.
        for i in 0..50 {
            let t = i as f32;
            let origin = [-1.0, t * 1.9, t * 0.7];
            let dir = [1.0, 0.0, 0.0];
            let mut ref_hits = reference.query_ray(origin, dir, 200.0);
            let mut gpu_hits = gpu.query_ray(origin, dir, 200.0);
            ref_hits.sort_unstable();
            gpu_hits.sort_unstable();
            assert_eq!(ref_hits, gpu_hits, "ray {i} hit sets differ");
        }
    }
}
