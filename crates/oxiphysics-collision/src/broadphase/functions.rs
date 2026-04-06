//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::CollisionPair;
use oxiphysics_core::Aabb;
use oxiphysics_core::math::{Real, Vec3};

use super::types::{
    BroadphaseStats, BruteForceBroadPhase, BvhBroadphase, BvhNodeData, DynamicAabbTree, Frustum,
    ObjectType, PairCountHistogram, SweepAndPrune,
};

/// Trait for broad-phase collision detection.
pub trait BroadPhase {
    /// Return candidate collision pairs from a set of AABBs.
    fn find_pairs(&self, aabbs: &[Aabb]) -> Vec<CollisionPair>;
}
/// Perform a batch query: for each query AABB, find all overlapping objects.
///
/// Returns a vector of hit lists (one per query).
#[allow(dead_code)]
pub fn batch_query(bvh: &BvhBroadphase, queries: &[Aabb]) -> Vec<Vec<usize>> {
    queries.iter().map(|q| bvh.query(q)).collect()
}
/// Perform a batch ray query: for each ray, find all hit objects.
#[allow(dead_code)]
pub fn batch_ray_query(
    bvh: &BvhBroadphase,
    rays: &[(Vec3, Vec3)],
    max_toi: Real,
) -> Vec<Vec<usize>> {
    rays.iter()
        .map(|(origin, dir)| bvh.ray_query(origin, dir, max_toi))
        .collect()
}
/// Frustum culling: filter a BVH by a view frustum.
///
/// Returns indices of objects whose AABBs are potentially visible.
#[allow(dead_code)]
pub fn frustum_cull(_bvh: &BvhBroadphase, aabbs: &[Aabb], frustum: &Frustum) -> Vec<usize> {
    let all: Vec<usize> = (0..aabbs.len()).collect();
    all.into_iter()
        .filter(|&i| frustum.contains_aabb(&aabbs[i]))
        .collect()
}
/// Collect statistics from a brute-force query.
#[allow(dead_code)]
pub fn brute_force_with_stats(aabbs: &[Aabb]) -> (Vec<CollisionPair>, BroadphaseStats) {
    let bf = BruteForceBroadPhase;
    let pairs = bf.find_pairs(aabbs);
    let n = aabbs.len();
    let stats = BroadphaseStats {
        num_objects: n,
        num_pairs: pairs.len(),
        num_tests: n * (n.saturating_sub(1)) / 2,
    };
    (pairs, stats)
}
/// Collect statistics from a SAP query.
#[allow(dead_code)]
pub fn sap_with_stats(aabbs: &[Aabb], axis: usize) -> (Vec<CollisionPair>, BroadphaseStats) {
    let sap = SweepAndPrune::new(axis);
    let pairs = sap.find_pairs(aabbs);
    let stats = BroadphaseStats {
        num_objects: aabbs.len(),
        num_pairs: pairs.len(),
        num_tests: 0,
    };
    (pairs, stats)
}
/// Parallel Sweep-and-Prune: simulates data-parallel execution by
/// pre-sorting all objects along the chosen axis, then performing
/// the standard SAP sweep.  In a real multi-threaded implementation
/// each thread would own a contiguous slice of the sorted array.
///
/// This implementation is functionally identical to [`SweepAndPrune`];
/// the "parallel" aspect refers to the data-layout strategy.
#[allow(dead_code)]
pub fn parallel_sap(aabbs: &[Aabb], axis: usize) -> Vec<CollisionPair> {
    SweepAndPrune::new(axis).find_pairs(aabbs)
}
/// Parallel brute-force: chunks the pair index space.
/// Useful for small scenes where BVH overhead is not worth it.
#[allow(dead_code)]
pub fn parallel_brute_force(aabbs: &[Aabb]) -> Vec<CollisionPair> {
    let mut pairs = Vec::new();
    let n = aabbs.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if aabbs[i].intersects(&aabbs[j]) {
                pairs.push(CollisionPair::new(i, j));
            }
        }
    }
    pairs
}
/// Test whether `aabb` is potentially inside the frustum defined by `planes`.
pub(super) fn aabb_in_frustum(aabb: &Aabb, planes: &[(Vec3, Real); 6]) -> bool {
    for &(ref n, d) in planes {
        let px = if n.x >= 0.0 { aabb.max.x } else { aabb.min.x };
        let py = if n.y >= 0.0 { aabb.max.y } else { aabb.min.y };
        let pz = if n.z >= 0.0 { aabb.max.z } else { aabb.min.z };
        let dot = n.x * px + n.y * py + n.z * pz;
        if dot < d {
            return false;
        }
    }
    true
}
/// Expand `aabb` uniformly by `amount` on all sides.
pub(super) fn inflate_aabb(aabb: &Aabb, amount: Real) -> Aabb {
    let v = Vec3::new(amount, amount, amount);
    Aabb::new(aabb.min - v, aabb.max + v)
}
/// Check if an AABB fully contains another AABB.
///
/// This is an extension trait method that we implement as a free function
/// to avoid modifying the core crate.
#[allow(dead_code)]
pub(super) fn aabb_contains(outer: &Aabb, inner: &Aabb) -> bool {
    outer.min.x <= inner.min.x
        && outer.min.y <= inner.min.y
        && outer.min.z <= inner.min.z
        && outer.max.x >= inner.max.x
        && outer.max.y >= inner.max.y
        && outer.max.z >= inner.max.z
}
/// Extension trait for Aabb to check containment.
pub trait AabbContains {
    fn contains_aabb(&self, other: &Aabb) -> bool;
}
impl AabbContains for Aabb {
    fn contains_aabb(&self, other: &Aabb) -> bool {
        aabb_contains(self, other)
    }
}
/// Compute a histogram of collision pair counts broken down by object-type pair.
///
/// `aabbs` and `types` must have the same length.  `pairs` is the output of
/// any broadphase `find_pairs` call.
#[allow(dead_code)]
pub fn compute_pair_count_histogram(
    pairs: &[CollisionPair],
    types: &[ObjectType],
) -> PairCountHistogram {
    let mut hist = PairCountHistogram::default();
    for p in pairs {
        if p.a >= types.len() || p.b >= types.len() {
            continue;
        }
        let ta = types[p.a];
        let tb = types[p.b];
        let (lo, hi) = if (ta as u8) <= (tb as u8) {
            (ta, tb)
        } else {
            (tb, ta)
        };
        match (lo, hi) {
            (ObjectType::Static, ObjectType::Static) => hist.static_static += 1,
            (ObjectType::Static, ObjectType::Dynamic) => hist.static_dynamic += 1,
            (ObjectType::Static, ObjectType::Kinematic) => hist.static_kinematic += 1,
            (ObjectType::Dynamic, ObjectType::Dynamic) => hist.dynamic_dynamic += 1,
            (ObjectType::Dynamic, ObjectType::Kinematic) => hist.dynamic_kinematic += 1,
            (ObjectType::Kinematic, ObjectType::Kinematic) => {
                hist.kinematic_kinematic += 1;
            }
            _ => {}
        }
    }
    hist
}
/// Batch-update a [`DynamicAabbTree`] for a set of moved objects.
///
/// `updates` maps object index → new tight AABB.  Objects not present in
/// `updates` are left unchanged.  The tree is marked dirty and rebuilt lazily
/// on the next query.
#[allow(dead_code)]
pub fn update_batch(tree: &mut DynamicAabbTree, updates: &[(usize, Aabb)]) {
    for (idx, aabb) in updates {
        tree.update(*idx, aabb.clone());
    }
}
/// Bottom-up BVH refit: walk every internal node from the leaves upward and
/// recompute its AABB as the union of its two children.
///
/// This is cheaper than a full rebuild when leaf AABBs change only slightly.
/// The BVH topology (shape) is preserved; only node AABBs are updated.
///
/// Returns the number of internal nodes that were refitted.
#[allow(dead_code)]
pub fn refit_bottom_up(bvh: &mut BvhBroadphase) -> usize {
    let Some(root) = bvh.root else {
        return 0;
    };
    let mut order: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = vec![root];
    while let Some(idx) = stack.pop() {
        order.push(idx);
        match bvh.nodes[idx].data {
            BvhNodeData::Internal { left, right } => {
                stack.push(left);
                stack.push(right);
            }
            BvhNodeData::Leaf { .. } => {}
        }
    }
    let mut refitted = 0usize;
    for &idx in order.iter().rev() {
        if let BvhNodeData::Internal { left, right } = bvh.nodes[idx].data {
            let merged = bvh.nodes[left].aabb.merge(&bvh.nodes[right].aabb);
            bvh.nodes[idx].aabb = merged;
            refitted += 1;
        }
    }
    refitted
}
/// Slab method ray-AABB intersection test.
pub(super) fn ray_aabb_intersect(
    origin: &Vec3,
    direction: &Vec3,
    aabb: &Aabb,
    max_toi: Real,
) -> bool {
    let mut tmin = 0.0_f64;
    let mut tmax = max_toi;
    for i in 0..3 {
        if direction[i].abs() < 1e-12 {
            if origin[i] < aabb.min[i] || origin[i] > aabb.max[i] {
                return false;
            }
        } else {
            let inv_d = 1.0 / direction[i];
            let mut t1 = (aabb.min[i] - origin[i]) * inv_d;
            let mut t2 = (aabb.max[i] - origin[i]) * inv_d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            tmin = tmin.max(t1);
            tmax = tmax.min(t2);
            if tmin > tmax {
                return false;
            }
        }
    }
    true
}
