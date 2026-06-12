// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU/CPU parity integration tests for the hybrid GPU LBVH builder
//! (`oxiphysics_gpu::gpu_lbvh`).
//!
//! The hybrid builder GPU radix-sorts the `(morton, original_index)` pairs and
//! then assembles the hierarchy with the CPU reference [`compute_bvh_from_sorted`].
//! Because the GPU pair-sort is stable, the emitted order is exactly that of a
//! CPU `sort_by_key(|p| (p.morton, original_index))`, so the resulting tree is
//! bit-identical to the CPU reference path by construction. These tests assert
//! that bit-exact equality at the project's standard 10^5-leaf budget.
//!
//! Skip-not-fail on headless CI: if no adapter is available the GPU-gated tests
//! print "SKIPPED: no GPU adapter available" and return.

/// Always-available smoke test so the integration-test binary always links and
/// has at least one test even when compiled without the `wgpu-backend` feature.
#[test]
fn gpu_lbvh_link_smoke() {
    // Intentionally trivial — just proves the integration test binary links.
    assert_eq!(1 + 1, 2);
}

/// Halton sequence helper (base `base`, index `index`, 0-indexed).
///
/// Copied verbatim from `tests/wgpu_kernels.rs` so the generated scene matches
/// the project's conventions; duplicating this private helper avoids cross-file
/// test coupling.
fn halton(index: usize, base: u32) -> f32 {
    let mut n = index;
    let b = base as usize;
    let mut f = 1.0_f64;
    let mut r = 0.0_f64;
    while n > 0 {
        f /= b as f64;
        r += f * (n % b) as f64;
        n /= b;
    }
    r as f32
}

#[cfg(feature = "wgpu-backend")]
mod gpu_tests {
    use oxiphysics_gpu::bvh::{
        Aabb, Bvh, BvhNode, BvhPrimitive, LbvhPrimitive, compute_bvh_from_sorted,
    };
    use oxiphysics_gpu::compute::wgpu_backend::real::WgpuBackendReal;
    use oxiphysics_gpu::gpu_lbvh::{gpu_lbvh_build, gpu_morton_sort};

    use super::halton;

    /// Macro: obtain a real `WgpuBackendReal` or skip the test if no GPU adapter
    /// is available (headless CI).
    macro_rules! skip_if_no_gpu {
        () => {
            match WgpuBackendReal::try_new() {
                Ok(b) => b,
                Err(_) => {
                    eprintln!("SKIPPED: no GPU adapter available");
                    return;
                }
            }
        };
    }

    // ── Scene conventions (match `test_bvh_gpu_parity_10e5_leaves`) ────────────

    /// Edge length of the cubic scene primitives are scattered through.
    const SCENE_SIZE: f32 = 100.0;
    /// Half-extent of each leaf box (0.5-unit boxes that typically do not
    /// overlap their Halton-spaced neighbours in a 100-unit scene).
    const BOX_HALF: f32 = 0.25;

    /// Build `n` primitives at Halton(2,3,5)·SCENE_SIZE centres with BOX_HALF
    /// half-extents — the exact recipe used by `test_bvh_gpu_parity_10e5_leaves`.
    fn make_halton_primitives(n: usize) -> Vec<BvhPrimitive> {
        (0..n)
            .map(|i| {
                let x = halton(i, 2) * SCENE_SIZE;
                let y = halton(i, 3) * SCENE_SIZE;
                let z = halton(i, 5) * SCENE_SIZE;
                BvhPrimitive::new(
                    Aabb::new(
                        [x - BOX_HALF, y - BOX_HALF, z - BOX_HALF],
                        [x + BOX_HALF, y + BOX_HALF, z + BOX_HALF],
                    ),
                    i,
                )
            })
            .collect()
    }

    /// CPU reference: fold the scene AABB, derive canonical Morton codes via
    /// [`LbvhPrimitive::new`], sort the *indices* by `(morton, original_index)`,
    /// and materialise the sorted [`LbvhPrimitive`] slice. This is exactly what
    /// the GPU hybrid path reproduces (the GPU pair-sort is stable).
    fn cpu_reference_sorted(prims: &[BvhPrimitive]) -> Vec<LbvhPrimitive> {
        if prims.is_empty() {
            return Vec::new();
        }
        // Scene AABB = fold of all prim aabbs via `Aabb::merge`.
        let mut scene = prims[0].aabb.clone();
        for p in &prims[1..] {
            scene = Aabb::merge(&scene, &p.aabb);
        }
        let lbvh_prims: Vec<LbvhPrimitive> = prims
            .iter()
            .map(|p| LbvhPrimitive::new(p.aabb.clone(), p.object_id, &scene))
            .collect();
        // Sort INDICES by (morton, original_index), then materialise.
        let mut order: Vec<usize> = (0..lbvh_prims.len()).collect();
        order.sort_by_key(|&i| (lbvh_prims[i].morton, i));
        order.into_iter().map(|i| lbvh_prims[i].clone()).collect()
    }

    /// Recursively collect leaf `object_id`s, mapping each leaf primitive index
    /// through `bvh.primitives` (which `compute_bvh_from_sorted` stores in sorted
    /// order). Returned sorted for a stable, order-independent comparison.
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

    // ── A) Main 100k topology + ray parity test ───────────────────────────────

    /// The main parity test. At the project's standard 10^5-leaf budget, build a
    /// Halton-positioned scene and assert the GPU hybrid LBVH is bit-identical to
    /// the CPU reference: same node count, depth, primitive count, leaf set, and
    /// ray-traversal hit set across 10_000 deterministic rays.
    #[test]
    fn test_lbvh_100k_topology_and_ray_parity() {
        let _backend = skip_if_no_gpu!();

        const N_LEAVES: usize = 100_000;
        const N_RAYS: usize = 10_000;

        let prims = make_halton_primitives(N_LEAVES);

        // Reference: CPU sort by (morton, index) then the same builder.
        let sorted = cpu_reference_sorted(&prims);
        let reference = compute_bvh_from_sorted(&sorted);

        // GPU hybrid: GPU sort then the same builder.
        let gpu = gpu_lbvh_build(prims.clone());

        // Structural parity.
        assert_eq!(
            reference.node_count(),
            gpu.node_count(),
            "node_count mismatch: reference={} gpu={}",
            reference.node_count(),
            gpu.node_count()
        );
        assert_eq!(
            reference.depth(),
            gpu.depth(),
            "depth mismatch: reference={} gpu={}",
            reference.depth(),
            gpu.depth()
        );
        assert_eq!(
            reference.primitives.len(),
            N_LEAVES,
            "reference primitive count mismatch"
        );
        assert_eq!(
            gpu.primitives.len(),
            N_LEAVES,
            "gpu primitive count mismatch"
        );

        // Identical leaf set (the "identical leaf set" bar).
        let ref_leaves = leaf_object_ids(&reference);
        let gpu_leaves = leaf_object_ids(&gpu);
        assert_eq!(
            ref_leaves.len(),
            N_LEAVES,
            "reference leaf set should cover every primitive"
        );
        assert_eq!(
            ref_leaves, gpu_leaves,
            "leaf object_id sets differ between CPU reference and GPU hybrid"
        );

        // Ray-traversal hit parity: fire N_RAYS Halton(7,11)·SCENE_SIZE + 0.25
        // origins straight down +Z. The trees are bit-identical by construction,
        // so we assert EXACT hit-set equality per ray.
        let mut total_hits = 0usize;
        for i in 0..N_RAYS {
            let ox = halton(i, 7) * SCENE_SIZE + 0.25;
            let oy = halton(i, 11) * SCENE_SIZE + 0.25;
            let origin = [ox, oy, 0.0];
            let dir = [0.0, 0.0, 1.0];
            let max_t = SCENE_SIZE + 2.0;
            let mut ref_hits = reference.query_ray(origin, dir, max_t);
            let mut gpu_hits = gpu.query_ray(origin, dir, max_t);
            ref_hits.sort_unstable();
            gpu_hits.sort_unstable();
            total_hits += ref_hits.len();
            assert_eq!(
                ref_hits, gpu_hits,
                "ray {i} hit sets differ (origin={origin:?}): reference={ref_hits:?} gpu={gpu_hits:?}"
            );
        }
        eprintln!(
            "PARITY BAR USED: EXACT equality (node_count={}, depth={}, leaves={}, rays={}, total ray-hits={})",
            gpu.node_count(),
            gpu.depth(),
            N_LEAVES,
            N_RAYS,
            total_hits
        );
    }

    // ── B) Empty edge ─────────────────────────────────────────────────────────

    /// `gpu_lbvh_build(vec![])` must yield an empty tree without panicking.
    #[test]
    fn test_lbvh_empty_edge() {
        let _backend = skip_if_no_gpu!();

        let bvh = gpu_lbvh_build(vec![]);
        assert!(bvh.root.is_none(), "empty input must yield no root");
        assert_eq!(bvh.node_count(), 0, "empty input must yield zero nodes");
        assert!(
            bvh.primitives.is_empty(),
            "empty input must yield no primitives"
        );
    }

    // ── C) Single edge ────────────────────────────────────────────────────────

    /// A single primitive must produce a one-leaf tree without panicking.
    #[test]
    fn test_lbvh_single_edge() {
        let _backend = skip_if_no_gpu!();

        let prims = vec![BvhPrimitive::new(
            Aabb::new([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]),
            42,
        )];
        let bvh = gpu_lbvh_build(prims);
        assert!(bvh.root.is_some(), "single input must yield a root");
        assert_eq!(
            bvh.primitives.len(),
            1,
            "single input must keep one primitive"
        );
        assert_eq!(
            bvh.primitives[0].object_id, 42,
            "single primitive object_id must be preserved"
        );

        // The lone primitive must also be the only leaf object.
        let leaves = leaf_object_ids(&bvh);
        assert_eq!(leaves, vec![42usize], "single leaf object_id mismatch");
    }

    // ── D) Morton-sort parity at 100k ─────────────────────────────────────────

    /// For the same 100k Halton scene, assert that `gpu_morton_sort` produces the
    /// exact same `(morton, object_id)` sequence as the CPU `(morton, index)`
    /// sort — verifying the GPU sort underpinning the build at scale.
    #[test]
    fn test_gpu_morton_sort_parity_100k() {
        let _backend = skip_if_no_gpu!();

        const N_LEAVES: usize = 100_000;
        let prims = make_halton_primitives(N_LEAVES);

        let gpu_sorted = gpu_morton_sort(&prims);
        let cpu_sorted = cpu_reference_sorted(&prims);

        assert_eq!(gpu_sorted.len(), N_LEAVES, "gpu sort length mismatch");
        assert_eq!(cpu_sorted.len(), N_LEAVES, "cpu sort length mismatch");

        // The (morton, object_id) sequence must be identical position-by-position.
        for (k, (g, c)) in gpu_sorted.iter().zip(cpu_sorted.iter()).enumerate() {
            assert_eq!(
                (g.morton, g.object_id),
                (c.morton, c.object_id),
                "morton/object_id mismatch at sorted position {k}"
            );
        }

        // Mortons must be non-decreasing across the whole sequence.
        for w in gpu_sorted.windows(2) {
            assert!(
                w[0].morton <= w[1].morton,
                "mortons not sorted: {} > {}",
                w[0].morton,
                w[1].morton
            );
        }
    }
}
