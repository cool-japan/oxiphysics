//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::collections::HashMap;

use super::types::{
    CcdBodyState, CcdHit, NarrowphaseBatch, ParAabb, ParBox, ParCollisionConfig,
    ParCollisionResult, ParContact, ParShapeKind, ParSphere, ParallelUnionFind, SoaAabbBuffer,
};

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Squared length of a 3-vector.
#[inline]
pub fn len_sq3(a: [f64; 3]) -> f64 {
    dot3(a, a)
}
/// Length of a 3-vector.
#[inline]
pub fn len3(a: [f64; 3]) -> f64 {
    len_sq3(a).sqrt()
}
/// Normalize a 3-vector; returns `[0,0,0]` for zero-length input.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}
/// Component-wise minimum.
#[inline]
pub fn min3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])]
}
/// Component-wise maximum.
#[inline]
pub fn max3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])]
}
/// Test one AABB against a slice of AABBs, returning indices of overlaps.
///
/// This function is structured to allow the compiler to auto-vectorize the
/// inner loop (all comparisons are simple scalar ops without early exit).
pub fn aabb_overlap_batch(query: &ParAabb, targets: &[ParAabb]) -> Vec<usize> {
    let mut result = Vec::new();
    let qminx = query.min[0];
    let qminy = query.min[1];
    let qminz = query.min[2];
    let qmaxx = query.max[0];
    let qmaxy = query.max[1];
    let qmaxz = query.max[2];
    for (i, t) in targets.iter().enumerate() {
        let overlap = (qminx <= t.max[0])
            & (qmaxx >= t.min[0])
            & (qminy <= t.max[1])
            & (qmaxy >= t.min[1])
            & (qminz <= t.max[2])
            & (qmaxz >= t.min[2]);
        if overlap {
            result.push(i);
        }
    }
    result
}
/// Test a slice of AABBs against another slice, collecting all overlapping pairs.
///
/// Returns pairs `(i, j)` with `i < j`.
pub fn aabb_all_pairs_overlap(aabbs: &[ParAabb]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..aabbs.len() {
        for j in (i + 1)..aabbs.len() {
            if aabbs[i].overlaps(&aabbs[j]) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}
/// Index into a node array; `u32::MAX` represents null.
pub type NodeIdx = u32;
/// Test sphere-sphere contact. Returns `Some(ParContact)` if penetrating.
pub fn par_sphere_sphere(
    a_idx: u32,
    b_idx: u32,
    a: &ParSphere,
    b: &ParSphere,
) -> Option<ParContact> {
    let delta = sub3(a.center, b.center);
    let dist_sq = len_sq3(delta);
    let sum_r = a.radius + b.radius;
    if dist_sq >= sum_r * sum_r {
        return None;
    }
    let dist = dist_sq.sqrt();
    let depth = sum_r - dist;
    let normal = if dist < 1e-14 {
        [0.0, 1.0, 0.0]
    } else {
        normalize3(delta)
    };
    let point = sub3(a.center, scale3(normal, a.radius - depth * 0.5));
    Some(ParContact {
        body_a: a_idx,
        body_b: b_idx,
        normal,
        depth,
        point,
    })
}
/// Test sphere-box (AABB) contact. Returns `Some(ParContact)` if penetrating.
pub fn par_sphere_box(
    a_idx: u32,
    b_idx: u32,
    sphere: &ParSphere,
    bx: &ParBox,
) -> Option<ParContact> {
    let mut closest = [0.0f64; 3];
    for i in 0..3 {
        let lo = bx.center[i] - bx.half_extents[i];
        let hi = bx.center[i] + bx.half_extents[i];
        closest[i] = sphere.center[i].max(lo).min(hi);
    }
    let delta = sub3(sphere.center, closest);
    let dist_sq = len_sq3(delta);
    if dist_sq >= sphere.radius * sphere.radius {
        return None;
    }
    let dist = dist_sq.sqrt();
    let depth = sphere.radius - dist;
    let normal = if dist < 1e-14 {
        [0.0, 1.0, 0.0]
    } else {
        normalize3(delta)
    };
    Some(ParContact {
        body_a: a_idx,
        body_b: b_idx,
        normal,
        depth,
        point: closest,
    })
}
/// Test box-box (AABB) contact using SAT on the three coordinate axes.
#[allow(clippy::too_many_arguments)]
pub fn par_box_box(a_idx: u32, b_idx: u32, a: &ParBox, b: &ParBox) -> Option<ParContact> {
    let mut min_depth = f64::MAX;
    let mut best_axis = 0usize;
    for i in 0..3 {
        let a_lo = a.center[i] - a.half_extents[i];
        let a_hi = a.center[i] + a.half_extents[i];
        let b_lo = b.center[i] - b.half_extents[i];
        let b_hi = b.center[i] + b.half_extents[i];
        let overlap = a_hi.min(b_hi) - a_lo.max(b_lo);
        if overlap <= 0.0 {
            return None;
        }
        if overlap < min_depth {
            min_depth = overlap;
            best_axis = i;
        }
    }
    let mut normal = [0.0f64; 3];
    let sign = if a.center[best_axis] < b.center[best_axis] {
        -1.0
    } else {
        1.0
    };
    normal[best_axis] = sign;
    let point = b.center;
    Some(ParContact {
        body_a: a_idx,
        body_b: b_idx,
        normal,
        depth: min_depth,
        point,
    })
}
/// Test two swept spheres for continuous collision.
///
/// Returns `Some(CcdHit)` with the earliest time of impact in `[0, 1]`.
pub fn ccd_sphere_sphere(a: &CcdBodyState, b: &CcdBodyState) -> Option<CcdHit> {
    let p = sub3(a.pos0, b.pos0);
    let d = sub3(sub3(a.pos1, a.pos0), sub3(b.pos1, b.pos0));
    let sum_r = a.radius + b.radius;
    let aa = dot3(d, d);
    let bb = 2.0 * dot3(p, d);
    let cc = dot3(p, p) - sum_r * sum_r;
    if aa < 1e-20 {
        if cc <= 0.0 {
            let normal = normalize3(p);
            return Some(CcdHit {
                toi: 0.0,
                body_a: a.index,
                body_b: b.index,
                normal,
            });
        }
        return None;
    }
    let disc = bb * bb - 4.0 * aa * cc;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let t0 = (-bb - sqrt_disc) / (2.0 * aa);
    let t1 = (-bb + sqrt_disc) / (2.0 * aa);
    let toi = if (0.0..=1.0).contains(&t0) {
        t0
    } else if (0.0..=1.0).contains(&t1) {
        t1
    } else {
        return None;
    };
    let pa = add3(a.pos0, scale3(sub3(a.pos1, a.pos0), toi));
    let pb = add3(b.pos0, scale3(sub3(b.pos1, b.pos0), toi));
    let normal = normalize3(sub3(pa, pb));
    Some(CcdHit {
        toi,
        body_a: a.index,
        body_b: b.index,
        normal,
    })
}
/// Parallel CCD sweep: test all pairs in a body list.
///
/// Returns a list of CCD hits sorted by `toi`.
pub fn par_ccd_sweep(bodies: &[CcdBodyState]) -> Vec<CcdHit> {
    let swept: Vec<ParAabb> = bodies.iter().map(|b| b.swept_aabb()).collect();
    let pairs = aabb_all_pairs_overlap(&swept);
    let mut hits = Vec::new();
    for (i, j) in pairs {
        if let Some(hit) = ccd_sphere_sphere(&bodies[i], &bodies[j]) {
            hits.push(hit);
        }
    }
    hits.sort_by(|a, b| {
        a.toi
            .partial_cmp(&b.toi)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    hits
}
/// Reduces a list of contacts to at most `max_contacts` per body pair.
///
/// Selects contacts that maximize coverage (deepest + most spread out).
pub fn par_contact_reduction(contacts: &[ParContact], max_per_pair: usize) -> Vec<ParContact> {
    if contacts.is_empty() {
        return Vec::new();
    }
    let mut groups: HashMap<(u32, u32), Vec<&ParContact>> = HashMap::new();
    for c in contacts {
        let key = (c.body_a.min(c.body_b), c.body_a.max(c.body_b));
        groups.entry(key).or_default().push(c);
    }
    let mut result = Vec::new();
    for (_key, mut group) in groups {
        group.sort_by(|a, b| {
            b.depth
                .partial_cmp(&a.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if group.len() <= max_per_pair {
            result.extend(group.iter().map(|&c| *c));
            continue;
        }
        let mut kept: Vec<&ParContact> = vec![group[0]];
        for candidate in group.iter().skip(1) {
            if kept.len() >= max_per_pair {
                break;
            }
            let min_dist_sq = kept
                .iter()
                .map(|k| len_sq3(sub3(candidate.point, k.point)))
                .fold(f64::MAX, f64::min);
            if min_dist_sq > 1e-6 {
                kept.push(candidate);
            }
        }
        result.extend(kept.into_iter().copied());
    }
    result
}
/// Build contact islands from a list of contact pairs.
///
/// Returns groups of body indices that form connected components.
pub fn build_islands(n_bodies: usize, contacts: &[ParContact]) -> Vec<Vec<usize>> {
    let uf = ParallelUnionFind::new(n_bodies);
    for c in contacts {
        let a = c.body_a as usize;
        let b = c.body_b as usize;
        if a < n_bodies && b < n_bodies {
            uf.union(a, b);
        }
    }
    uf.islands()
}
/// Run one frame of parallel collision detection.
///
/// Steps:
/// 1. Broadphase via `SoaAabbBuffer`.
/// 2. Narrowphase via `NarrowphaseBatch`.
/// 3. Contact reduction.
/// 4. Island detection.
/// 5. Optionally CCD sweep.
#[allow(clippy::too_many_arguments)]
pub fn run_parallel_collision(
    aabbs: &[ParAabb],
    shapes: &[ParShapeKind],
    ccd_bodies: &[CcdBodyState],
    config: &ParCollisionConfig,
) -> ParCollisionResult {
    let mut soa = SoaAabbBuffer::new();
    for (i, aabb) in aabbs.iter().enumerate() {
        soa.push(&aabb.expanded(config.aabb_margin), i as u32);
    }
    let pairs = aabb_all_pairs_overlap(aabbs);
    let mut batch = NarrowphaseBatch::new();
    for (i, j) in &pairs {
        if *i < shapes.len() && *j < shapes.len() {
            batch.push(*i as u32, shapes[*i], *j as u32, shapes[*j]);
        }
    }
    let raw_contacts = batch.process();
    let contacts = par_contact_reduction(&raw_contacts, config.max_contacts_per_pair);
    let n_bodies = aabbs.len();
    let islands = build_islands(n_bodies, &contacts);
    let ccd_hits = if config.enable_ccd {
        par_ccd_sweep(ccd_bodies)
    } else {
        Vec::new()
    };
    ParCollisionResult {
        contacts,
        islands,
        ccd_hits,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::parallel_collision::GpuAabbBufferDesc;
    use crate::parallel_collision::GpuBvhBuffer;
    use crate::parallel_collision::LockFreeContactBuffer;
    use crate::parallel_collision::ParallelAabbTree;
    use crate::parallel_collision::WorkStealingBroadphase;
    #[test]
    fn test_paraabb_overlaps_touching() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_paraabb_no_overlap_x() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_paraabb_no_overlap_y() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([0.0, 2.0, 0.0], [1.0, 3.0, 1.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_paraabb_no_overlap_z() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([0.0, 0.0, 2.0], [1.0, 1.0, 3.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_paraabb_center() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_paraabb_merge() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([0.5, 0.5, 0.5], [2.0, 2.0, 2.0]);
        let m = a.merge(&b);
        assert!((m.min[0]).abs() < 1e-10);
        assert!((m.max[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_paraabb_contains_point() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(a.contains_point([0.5, 0.5, 0.5]));
        assert!(!a.contains_point([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_paraabb_surface_area() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        assert!((a.surface_area() - 22.0).abs() < 1e-10);
    }
    #[test]
    fn test_paraabb_expanded() {
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let e = a.expanded(0.1);
        assert!((e.min[0] + 0.1).abs() < 1e-10);
        assert!((e.max[0] - 1.1).abs() < 1e-10);
    }
    #[test]
    fn test_soa_push_and_query() {
        let mut soa = SoaAabbBuffer::new();
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = ParAabb::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        soa.push(&a, 0);
        soa.push(&b, 1);
        let q = ParAabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        let hits = soa.query(&q);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0], 0);
    }
    #[test]
    fn test_soa_len_and_is_empty() {
        let mut soa = SoaAabbBuffer::new();
        assert!(soa.is_empty());
        soa.push(&ParAabb::new([0.0; 3], [1.0; 3]), 0);
        assert_eq!(soa.len(), 1);
    }
    #[test]
    fn test_soa_rebuild() {
        let mut soa = SoaAabbBuffer::new();
        soa.push(&ParAabb::new([0.0; 3], [1.0; 3]), 0);
        let entries = vec![(ParAabb::new([2.0; 3], [3.0; 3]), 1u32)];
        soa.rebuild(&entries);
        assert_eq!(soa.len(), 1);
    }
    #[test]
    fn test_parallel_tree_single_insert_query() {
        let tree = ParallelAabbTree::new();
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        tree.insert(a, 42);
        let result = tree.query(&a);
        assert!(result.contains(&42));
    }
    #[test]
    fn test_parallel_tree_no_overlap() {
        let tree = ParallelAabbTree::new();
        let a = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        tree.insert(a, 1);
        let q = ParAabb::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        assert!(tree.query(&q).is_empty());
    }
    #[test]
    fn test_parallel_tree_multiple_overlaps() {
        let tree = ParallelAabbTree::new();
        for i in 0..5u32 {
            let off = i as f64;
            tree.insert(ParAabb::new([off, 0.0, 0.0], [off + 1.0, 1.0, 1.0]), i);
        }
        let q = ParAabb::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        let res = tree.query(&q);
        assert!(res.len() >= 2);
    }
    #[test]
    fn test_work_stealing_no_pairs() {
        let aabbs = vec![
            ParAabb::new([0.0; 3], [1.0; 3]),
            ParAabb::new([5.0; 3], [6.0; 3]),
        ];
        let ws = WorkStealingBroadphase::new(aabbs, 8);
        assert!(ws.is_done());
    }
    #[test]
    fn test_work_stealing_has_pairs() {
        let aabbs = vec![
            ParAabb::new([0.0; 3], [2.0; 3]),
            ParAabb::new([1.0; 3], [3.0; 3]),
        ];
        let ws = WorkStealingBroadphase::new(aabbs, 8);
        assert!(!ws.is_done());
        let stolen = ws.steal();
        assert!(!stolen.is_empty());
    }
    #[test]
    fn test_work_stealing_chunk_size() {
        let mut aabbs = Vec::new();
        for i in 0..5usize {
            let off = i as f64;
            aabbs.push(ParAabb::new([off, 0.0, 0.0], [off + 1.5, 1.0, 1.0]));
        }
        let ws = WorkStealingBroadphase::new(aabbs, 2);
        let first = ws.steal();
        assert!(first.len() <= 2);
    }
    #[test]
    fn test_sphere_sphere_overlap() {
        let a = ParSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let b = ParSphere {
            center: [1.5, 0.0, 0.0],
            radius: 1.0,
        };
        let c = par_sphere_sphere(0, 1, &a, &b).unwrap();
        assert!(c.depth > 0.0);
        assert!((c.depth - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_sphere_sphere_separated() {
        let a = ParSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let b = ParSphere {
            center: [10.0, 0.0, 0.0],
            radius: 1.0,
        };
        assert!(par_sphere_sphere(0, 1, &a, &b).is_none());
    }
    #[test]
    fn test_sphere_box_overlap() {
        let s = ParSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let b = ParBox {
            center: [1.2, 0.0, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        };
        let c = par_sphere_box(0, 1, &s, &b).unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn test_sphere_box_separated() {
        let s = ParSphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.3,
        };
        let b = ParBox {
            center: [5.0, 0.0, 0.0],
            half_extents: [0.5, 0.5, 0.5],
        };
        assert!(par_sphere_box(0, 1, &s, &b).is_none());
    }
    #[test]
    fn test_box_box_overlap() {
        let a = ParBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let b = ParBox {
            center: [1.5, 0.0, 0.0],
            half_extents: [1.0, 1.0, 1.0],
        };
        let c = par_box_box(0, 1, &a, &b).unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn test_box_box_separated() {
        let a = ParBox {
            center: [0.0, 0.0, 0.0],
            half_extents: [0.4, 0.4, 0.4],
        };
        let b = ParBox {
            center: [2.0, 0.0, 0.0],
            half_extents: [0.4, 0.4, 0.4],
        };
        assert!(par_box_box(0, 1, &a, &b).is_none());
    }
    #[test]
    fn test_narrowphase_batch_mixed() {
        let mut batch = NarrowphaseBatch::new();
        let s1 = ParShapeKind::Sphere(ParSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        });
        let s2 = ParShapeKind::Sphere(ParSphere {
            center: [1.5, 0.0, 0.0],
            radius: 1.0,
        });
        batch.push(0, s1, 1, s2);
        let contacts = batch.process();
        assert_eq!(contacts.len(), 1);
    }
    #[test]
    fn test_lock_free_buffer_push_collect() {
        let buf = LockFreeContactBuffer::new(16);
        let c = ParContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            depth: 0.5,
            point: [0.0, 0.0, 0.0],
        };
        buf.push_contact(c).unwrap();
        assert_eq!(buf.len(), 1);
        let collected = buf.collect();
        assert_eq!(collected.len(), 1);
    }
    #[test]
    fn test_lock_free_buffer_overflow() {
        let buf = LockFreeContactBuffer::new(1);
        let c = ParContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            depth: 0.5,
            point: [0.0; 3],
        };
        assert!(buf.push_contact(c).is_some());
        assert!(buf.push_contact(c).is_none());
    }
    #[test]
    fn test_lock_free_buffer_reset() {
        let buf = LockFreeContactBuffer::new(8);
        let c = ParContact {
            body_a: 0,
            body_b: 1,
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
            point: [0.0; 3],
        };
        buf.push_contact(c).unwrap();
        buf.reset();
        assert_eq!(buf.len(), 0);
    }
    #[test]
    fn test_union_find_isolated() {
        let uf = ParallelUnionFind::new(4);
        assert_eq!(uf.num_components(), 4);
    }
    #[test]
    fn test_union_find_single_union() {
        let uf = ParallelUnionFind::new(4);
        uf.union(0, 1);
        assert_eq!(uf.num_components(), 3);
    }
    #[test]
    fn test_union_find_all_connected() {
        let uf = ParallelUnionFind::new(4);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        assert_eq!(uf.num_components(), 1);
    }
    #[test]
    fn test_union_find_islands() {
        let uf = ParallelUnionFind::new(6);
        uf.union(0, 1);
        uf.union(2, 3);
        let islands = uf.islands();
        assert_eq!(islands.len(), 4);
    }
    #[test]
    fn test_ccd_sphere_hit() {
        let a = CcdBodyState::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.5, 0);
        let b = CcdBodyState::new([3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1);
        let hit = ccd_sphere_sphere(&a, &b);
        assert!(hit.is_some());
        let h = hit.unwrap();
        assert!(h.toi >= 0.0 && h.toi <= 1.0);
    }
    #[test]
    fn test_ccd_sphere_no_hit() {
        let a = CcdBodyState::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3, 0);
        let b = CcdBodyState::new([10.0, 0.0, 0.0], [11.0, 0.0, 0.0], 0.3, 1);
        assert!(ccd_sphere_sphere(&a, &b).is_none());
    }
    #[test]
    fn test_ccd_swept_aabb() {
        let a = CcdBodyState::new([0.0, 0.0, 0.0], [5.0, 0.0, 0.0], 1.0, 0);
        let swept = a.swept_aabb();
        assert!(swept.max[0] >= 6.0);
        assert!(swept.min[0] <= -1.0);
    }
    #[test]
    fn test_par_ccd_sweep_finds_hit() {
        let bodies = vec![
            CcdBodyState::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], 0.5, 0),
            CcdBodyState::new([4.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1),
        ];
        let hits = par_ccd_sweep(&bodies);
        assert!(!hits.is_empty());
    }
    #[test]
    fn test_contact_reduction_empty() {
        let r = par_contact_reduction(&[], 4);
        assert!(r.is_empty());
    }
    #[test]
    fn test_contact_reduction_under_limit() {
        let contacts = vec![
            ParContact {
                body_a: 0,
                body_b: 1,
                normal: [0.0, 1.0, 0.0],
                depth: 0.1,
                point: [0.0, 0.0, 0.0],
            },
            ParContact {
                body_a: 0,
                body_b: 1,
                normal: [0.0, 1.0, 0.0],
                depth: 0.2,
                point: [1.0, 0.0, 0.0],
            },
        ];
        let r = par_contact_reduction(&contacts, 4);
        assert_eq!(r.len(), 2);
    }
    #[test]
    fn test_contact_reduction_over_limit() {
        let mut contacts = Vec::new();
        for i in 0..10u32 {
            contacts.push(ParContact {
                body_a: 0,
                body_b: 1,
                normal: [0.0, 1.0, 0.0],
                depth: i as f64 * 0.1,
                point: [i as f64, 0.0, 0.0],
            });
        }
        let r = par_contact_reduction(&contacts, 4);
        assert!(r.len() <= 4);
    }
    #[test]
    fn test_gpu_aabb_desc_total_bytes() {
        let desc = GpuAabbBufferDesc::new_f32_soa(100);
        assert_eq!(desc.total_bytes(), 100 * 4 * 6);
    }
    #[test]
    fn test_gpu_bvh_build_single() {
        let entries = vec![(ParAabb::new([0.0; 3], [1.0; 3]), 42u32)];
        let bvh = GpuBvhBuffer::build(&entries);
        assert_eq!(bvh.node_count(), 1);
        assert!(bvh.nodes[0].is_leaf());
        assert_eq!(bvh.nodes[0].right_or_handle, 42);
    }
    #[test]
    fn test_gpu_bvh_build_multiple() {
        let entries: Vec<_> = (0..4u32)
            .map(|i| {
                (
                    ParAabb::new([i as f64, 0.0, 0.0], [i as f64 + 1.0, 1.0, 1.0]),
                    i,
                )
            })
            .collect();
        let bvh = GpuBvhBuffer::build(&entries);
        assert_eq!(bvh.node_count(), 7);
    }
    #[test]
    fn test_gpu_bvh_to_bytes_size() {
        let entries = vec![(ParAabb::new([0.0; 3], [1.0; 3]), 0u32)];
        let bvh = GpuBvhBuffer::build(&entries);
        assert_eq!(bvh.to_bytes().len(), bvh.node_count() * 32);
    }
    #[test]
    fn test_soa_serialize_size() {
        let mut soa = SoaAabbBuffer::new();
        soa.push(&ParAabb::new([0.0; 3], [1.0; 3]), 0);
        soa.push(&ParAabb::new([2.0; 3], [3.0; 3]), 1);
        let bytes = GpuAabbBufferDesc::serialize_soa(&soa);
        assert_eq!(bytes.len(), 48);
    }
    #[test]
    fn test_run_parallel_collision_basic() {
        let aabbs = vec![
            ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            ParAabb::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
            ParAabb::new([10.0, 0.0, 0.0], [11.0, 1.0, 1.0]),
        ];
        let shapes = vec![
            ParShapeKind::Sphere(ParSphere {
                center: [0.5, 0.5, 0.5],
                radius: 0.5,
            }),
            ParShapeKind::Sphere(ParSphere {
                center: [1.0, 0.5, 0.5],
                radius: 0.5,
            }),
            ParShapeKind::Sphere(ParSphere {
                center: [10.5, 0.5, 0.5],
                radius: 0.5,
            }),
        ];
        let config = ParCollisionConfig::default();
        let result = run_parallel_collision(&aabbs, &shapes, &[], &config);
        assert!(!result.contacts.is_empty());
    }
    #[test]
    fn test_run_parallel_collision_ccd() {
        let aabbs = vec![
            ParAabb::new([0.0; 3], [1.0; 3]),
            ParAabb::new([3.0; 3], [4.0; 3]),
        ];
        let shapes = vec![
            ParShapeKind::Sphere(ParSphere {
                center: [0.5; 3],
                radius: 0.5,
            }),
            ParShapeKind::Sphere(ParSphere {
                center: [3.5; 3],
                radius: 0.5,
            }),
        ];
        let ccd_bodies = vec![
            CcdBodyState::new([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], 0.5, 0),
            CcdBodyState::new([4.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 1),
        ];
        let config = ParCollisionConfig {
            enable_ccd: true,
            ..Default::default()
        };
        let result = run_parallel_collision(&aabbs, &shapes, &ccd_bodies, &config);
        assert!(!result.ccd_hits.is_empty());
    }
    #[test]
    fn test_build_islands_basic() {
        let contacts = vec![
            ParContact {
                body_a: 0,
                body_b: 1,
                normal: [0.0, 1.0, 0.0],
                depth: 0.1,
                point: [0.0; 3],
            },
            ParContact {
                body_a: 2,
                body_b: 3,
                normal: [0.0, 1.0, 0.0],
                depth: 0.1,
                point: [0.0; 3],
            },
        ];
        let islands = build_islands(4, &contacts);
        assert_eq!(islands.len(), 2);
    }
    #[test]
    fn test_aabb_overlap_batch_multiple() {
        let query = ParAabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let targets = vec![
            ParAabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
            ParAabb::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
            ParAabb::new([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]),
        ];
        let hits = aabb_overlap_batch(&query, &targets);
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
        assert!(hits.contains(&2));
    }
    #[test]
    fn test_vec3_helpers() {
        assert!((len3([3.0, 4.0, 0.0]) - 5.0).abs() < 1e-10);
        let n = normalize3([1.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-10);
        let z = normalize3([0.0; 3]);
        assert_eq!(z, [0.0; 3]);
    }
}
