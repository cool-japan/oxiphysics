// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{AabbGpu, SupportPoint};

/// Expand a 10-bit integer into 30 bits by inserting two zero bits between
/// each original bit.  Used for 3-D Morton code generation.
pub(super) fn expand_bits(mut v: u32) -> u32 {
    v &= 0x0000_03FF;
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v << 8)) & 0x0300F00F;
    v = (v | (v << 4)) & 0x030C30C3;
    v = (v | (v << 2)) & 0x09249249;
    v
}
/// Compute a 30-bit Morton code for a normalised 3-D position `(x, y, z)`,
/// where each component is in `[0, 1]`.
pub fn morton_code(x: f32, y: f32, z: f32) -> u32 {
    let xi = (x.clamp(0.0, 1.0) * 1023.0) as u32;
    let yi = (y.clamp(0.0, 1.0) * 1023.0) as u32;
    let zi = (z.clamp(0.0, 1.0) * 1023.0) as u32;
    (expand_bits(xi) << 2) | (expand_bits(yi) << 1) | expand_bits(zi)
}
/// Dot product of two f32 3-vectors.
pub(super) fn dot3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Subtract two f32 3-vectors.
pub(super) fn sub3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two f32 3-vectors.
pub(super) fn add3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale an f32 3-vector.
pub(super) fn scale3f(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Squared length of an f32 3-vector.
pub(super) fn len_sq3f(a: [f32; 3]) -> f32 {
    dot3f(a, a)
}
/// Length of an f32 3-vector.
pub(super) fn len3f(a: [f32; 3]) -> f32 {
    len_sq3f(a).sqrt()
}
/// Normalise an f32 3-vector; returns zero vector if length < eps.
pub(super) fn norm3f(a: [f32; 3]) -> [f32; 3] {
    let l = len3f(a);
    if l < 1e-9 {
        [0.0; 3]
    } else {
        scale3f(a, 1.0 / l)
    }
}
/// Support function for an AABB along direction `d`.
pub(super) fn aabb_support(aabb: &AabbGpu, d: [f32; 3]) -> [f32; 3] {
    [
        if d[0] >= 0.0 {
            aabb.max[0]
        } else {
            aabb.min[0]
        },
        if d[1] >= 0.0 {
            aabb.max[1]
        } else {
            aabb.min[1]
        },
        if d[2] >= 0.0 {
            aabb.max[2]
        } else {
            aabb.min[2]
        },
    ]
}
/// Support in Minkowski difference A - B along `d`.
pub(super) fn minkowski_support(a: &AabbGpu, b: &AabbGpu, d: [f32; 3]) -> SupportPoint {
    let neg_d = [-d[0], -d[1], -d[2]];
    let pa = aabb_support(a, d);
    let pb = aabb_support(b, neg_d);
    SupportPoint::new(sub3f(pa, pb))
}
/// Update the GJK simplex and return `(intersecting, new_direction)`.
pub(super) fn update_simplex(simplex: &mut Vec<SupportPoint>) -> (bool, [f32; 3]) {
    match simplex.len() {
        2 => {
            let b = simplex[0].v;
            let a = simplex[1].v;
            let ab = sub3f(b, a);
            let ao = [-a[0], -a[1], -a[2]];
            if dot3f(ab, ao) > 0.0 {
                let t = dot3f(ao, ab) / dot3f(ab, ab);
                let closest = add3f(a, scale3f(ab, t));
                (false, [-closest[0], -closest[1], -closest[2]])
            } else {
                simplex.drain(0..1);
                (false, ao)
            }
        }
        3 => {
            let c = simplex[0].v;
            let b_pt = simplex[1].v;
            let a = simplex[2].v;
            let ab = sub3f(b_pt, a);
            let ac = sub3f(c, a);
            let ao = [-a[0], -a[1], -a[2]];
            let n = cross3f(ab, ac);
            if dot3f(n, ao) > 0.0 {
                (false, n)
            } else {
                (false, [-n[0], -n[1], -n[2]])
            }
        }
        4 => (true, [0.0; 3]),
        _ => {
            let a = simplex.last().expect("collection should not be empty").v;
            (false, [-a[0], -a[1], -a[2]])
        }
    }
}
/// Cross product for f32 3-vectors.
pub(super) fn cross3f(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::collision_gpu::BroadphaseGpuKernel;
    use crate::collision_gpu::BvhNodeGpu;
    use crate::collision_gpu::CollisionGpuPipeline;
    use crate::collision_gpu::CollisionKernelStats;
    use crate::collision_gpu::CollisionPair;
    use crate::collision_gpu::ContactResult;
    use crate::collision_gpu::GpuAabbTree;
    use crate::collision_gpu::GpuBroadphase;
    use crate::collision_gpu::GpuBvhBuilder;
    use crate::collision_gpu::GpuCollisionPipeline;
    use crate::collision_gpu::GpuContactCache;
    use crate::collision_gpu::GpuNarrowphase;
    use crate::collision_gpu::ManifoldPoint;
    use crate::collision_gpu::NarrowphaseGpuKernel;
    use crate::collision_gpu::PersistentManifoldGpu;
    fn unit_box() -> AabbGpu {
        AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    }
    #[test]
    fn test_aabb_new() {
        let a = unit_box();
        assert_eq!(a.min, [0.0, 0.0, 0.0]);
        assert_eq!(a.max, [1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_aabb_from_point() {
        let p = [3.0, 4.0, 5.0];
        let a = AabbGpu::from_point(p);
        assert_eq!(a.min, p);
        assert_eq!(a.max, p);
    }
    #[test]
    fn test_aabb_surface_area() {
        let a = unit_box();
        assert!((a.surface_area() - 6.0).abs() < 1e-5);
    }
    #[test]
    fn test_aabb_volume() {
        let a = unit_box();
        assert!((a.volume() - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_aabb_centre() {
        let a = unit_box();
        let c = a.centre();
        assert!((c[0] - 0.5).abs() < 1e-5);
        assert!((c[1] - 0.5).abs() < 1e-5);
        assert!((c[2] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_aabb_half_extents() {
        let a = unit_box();
        let h = a.half_extents();
        assert!((h[0] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_aabb_overlaps_true() {
        let a = unit_box();
        let b = AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_aabb_overlaps_touching() {
        let a = unit_box();
        let b = AabbGpu::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_aabb_overlaps_false() {
        let a = unit_box();
        let b = AabbGpu::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_aabb_contains_point_inside() {
        let a = unit_box();
        assert!(a.contains_point([0.5, 0.5, 0.5]));
    }
    #[test]
    fn test_aabb_contains_point_outside() {
        let a = unit_box();
        assert!(!a.contains_point([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_aabb_expanded() {
        let a = unit_box();
        let e = a.expanded(0.1);
        assert!((e.min[0] - (-0.1)).abs() < 1e-5);
        assert!((e.max[0] - 1.1).abs() < 1e-5);
    }
    #[test]
    fn test_aabb_merge() {
        let a = unit_box();
        let b = AabbGpu::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        let m = a.merge(&b);
        assert_eq!(m.min, [0.0, 0.0, 0.0]);
        assert_eq!(m.max, [3.0, 3.0, 3.0]);
    }
    #[test]
    fn test_morton_code_origin() {
        let code = morton_code(0.0, 0.0, 0.0);
        assert_eq!(code, 0);
    }
    #[test]
    fn test_morton_code_max() {
        let _code = morton_code(1.0, 1.0, 1.0);
    }
    #[test]
    fn test_morton_code_sorted_along_x() {
        let c1 = morton_code(0.1, 0.5, 0.5);
        let c2 = morton_code(0.9, 0.5, 0.5);
        assert_ne!(c1, c2);
    }
    #[test]
    fn test_bvh_node_leaf() {
        let node = BvhNodeGpu::leaf(unit_box(), 42);
        assert!(node.is_leaf);
        assert_eq!(node.primitive_index, 42);
    }
    #[test]
    fn test_bvh_node_internal() {
        let node = BvhNodeGpu::internal(unit_box(), 1, 2);
        assert!(!node.is_leaf);
        assert_eq!(node.left_child, 1);
        assert_eq!(node.right_child, 2);
    }
    #[test]
    fn test_bvh_builder_empty() {
        let mut builder = GpuBvhBuilder::new();
        let root = builder.build(&[]);
        assert_eq!(root, 0);
        assert!(builder.nodes().is_empty());
    }
    #[test]
    fn test_bvh_builder_single() {
        let mut builder = GpuBvhBuilder::new();
        let aabbs = vec![unit_box()];
        let root = builder.build(&aabbs);
        assert_eq!(builder.nodes()[root].primitive_index, 0);
        assert!(builder.nodes()[root].is_leaf);
    }
    #[test]
    fn test_bvh_builder_two() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let mut builder = GpuBvhBuilder::new();
        let root = builder.build(&aabbs);
        assert!(!builder.nodes()[root].is_leaf);
        assert_eq!(builder.nodes().len(), 3);
    }
    #[test]
    fn test_bvh_query_overlaps() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
            AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
        ];
        let mut builder = GpuBvhBuilder::new();
        let root = builder.build(&aabbs);
        let query = AabbGpu::new([0.2, 0.2, 0.2], [0.8, 0.8, 0.8]);
        let hits = builder.query_overlaps(root, &query);
        assert!(hits.contains(&0) || hits.contains(&2));
    }
    #[test]
    fn test_broadphase_no_overlap() {
        let kernel = BroadphaseGpuKernel::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let pairs = kernel.dispatch(&aabbs);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_broadphase_overlap() {
        let kernel = BroadphaseGpuKernel::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
        ];
        let pairs = kernel.dispatch(&aabbs);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].a, 0);
        assert_eq!(pairs[0].b, 1);
    }
    #[test]
    fn test_broadphase_margin_activates_pair() {
        let kernel = BroadphaseGpuKernel::new(1.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]),
        ];
        let pairs = kernel.dispatch(&aabbs);
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_broadphase_bvh_matches_brute() {
        let kernel = BroadphaseGpuKernel::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let brute = kernel.dispatch(&aabbs);
        let bvh = kernel.dispatch_bvh(&aabbs);
        assert_eq!(brute.len(), bvh.len());
    }
    #[test]
    fn test_narrowphase_no_contact() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let b = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        let result = kernel.test_aabb_pair(0, &a, 1, &b);
        assert!(result.is_none());
    }
    #[test]
    fn test_narrowphase_contact() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let b = AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        let result = kernel.test_aabb_pair(0, &a, 1, &b);
        assert!(result.is_some());
        let c = result.unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn test_narrowphase_dispatch_pairs() {
        let kernel = NarrowphaseGpuKernel::new();
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        let pairs = vec![CollisionPair::new(0, 1)];
        let contacts = kernel.dispatch(&pairs, &aabbs);
        assert_eq!(contacts.len(), 1);
    }
    #[test]
    fn test_gjk_intersect_overlapping() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let b = AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        let (intersect, _dist) = kernel.gjk_intersect(&a, &b);
        assert!(intersect);
    }
    #[test]
    fn test_gjk_separated() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let b = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        let (intersect, dist) = kernel.gjk_intersect(&a, &b);
        assert!(!intersect);
        assert!(dist > 0.0);
    }
    #[test]
    fn test_manifold_add_points() {
        let mut m = PersistentManifoldGpu::new(0, 1);
        let c = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.5, 0.0, 0.5],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        };
        m.add_contact(ManifoldPoint::from_contact(&c));
        assert_eq!(m.num_points, 1);
    }
    #[test]
    fn test_manifold_max_four_points() {
        let mut m = PersistentManifoldGpu::new(0, 1);
        for i in 0..6 {
            let c = ContactResult {
                prim_a: 0,
                prim_b: 1,
                contact_point: [i as f32, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                depth: 0.1 * (i as f32 + 1.0),
            };
            m.add_contact(ManifoldPoint::from_contact(&c));
        }
        assert!(m.num_points <= 4);
    }
    #[test]
    fn test_manifold_remove_stale() {
        let mut m = PersistentManifoldGpu::new(0, 1);
        let c = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.01,
        };
        m.add_contact(ManifoldPoint::from_contact(&c));
        m.remove_stale(0.05);
        assert_eq!(m.num_points, 0);
    }
    #[test]
    fn test_manifold_reset_warm_start() {
        let mut m = PersistentManifoldGpu::new(0, 1);
        let c = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
        };
        let mut pt = ManifoldPoint::from_contact(&c);
        pt.warm_impulse_normal = 5.0;
        m.add_contact(pt);
        m.reset_warm_start();
        for pt in m.active_points() {
            assert!((pt.warm_impulse_normal).abs() < 1e-6);
        }
    }
    #[test]
    fn test_pipeline_dispatch_no_contacts() {
        let mut pipeline = CollisionGpuPipeline::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([10.0, 0.0, 0.0], [11.0, 1.0, 1.0]),
        ];
        pipeline.dispatch(&aabbs);
        assert!(pipeline.contacts.is_empty());
    }
    #[test]
    fn test_pipeline_dispatch_with_contact() {
        let mut pipeline = CollisionGpuPipeline::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.dispatch(&aabbs);
        assert!(!pipeline.contacts.is_empty());
        assert!(pipeline.total_contact_points() > 0);
    }
    #[test]
    fn test_pipeline_reset() {
        let mut pipeline = CollisionGpuPipeline::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.dispatch(&aabbs);
        pipeline.reset();
        assert!(pipeline.contacts.is_empty());
        assert!(pipeline.manifolds.is_empty());
    }
    #[test]
    fn test_pipeline_bvh_dispatch() {
        let mut pipeline = CollisionGpuPipeline::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
            AabbGpu::new([10.0, 0.0, 0.0], [11.0, 1.0, 1.0]),
        ];
        pipeline.dispatch_bvh(&aabbs);
        assert!(!pipeline.contacts.is_empty());
    }
    #[test]
    fn test_pipeline_persistent_manifolds_accumulate() {
        let mut pipeline = CollisionGpuPipeline::new(0.0);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.dispatch(&aabbs);
        let first_count = pipeline.manifolds.len();
        pipeline.dispatch(&aabbs);
        assert!(pipeline.manifolds.len() >= first_count);
    }
    #[test]
    fn test_collision_pair_canonical() {
        let p1 = CollisionPair::new(3, 1);
        let p2 = CollisionPair::new(1, 3);
        assert_eq!(p1, p2);
        assert_eq!(p1.a, 1);
        assert_eq!(p1.b, 3);
    }
    #[test]
    fn test_stats_default_zero() {
        let s = CollisionKernelStats::new();
        assert_eq!(s.broadphase_pairs_tested, 0);
        assert_eq!(s.narrowphase_contacts, 0);
    }
    #[test]
    fn test_stats_broadphase_hit_rate() {
        let mut s = CollisionKernelStats::new();
        s.broadphase_pairs_tested = 10;
        s.broadphase_hits = 3;
        let rate = s.broadphase_hit_rate();
        assert!((rate - 0.3_f64).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_stats_broadphase_hit_rate_zero_tested() {
        let s = CollisionKernelStats::new();
        assert!((s.broadphase_hit_rate()).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_stats_narrowphase_contact_rate() {
        let mut s = CollisionKernelStats::new();
        s.narrowphase_queries = 8;
        s.narrowphase_contacts = 4;
        assert!((s.narrowphase_contact_rate() - 0.5_f64).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_stats_bandwidth() {
        let mut s = CollisionKernelStats::new();
        s.bytes_transferred = 1_000_000_000;
        s.elapsed_secs = 1.0_f64;
        assert!((s.bandwidth_gb_s() - 1.0_f64).abs() < 1.0e-9_f64);
    }
    #[test]
    fn test_stats_bandwidth_zero_time() {
        let s = CollisionKernelStats::new();
        assert!((s.bandwidth_gb_s()).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_stats_pair_throughput() {
        let mut s = CollisionKernelStats::new();
        s.broadphase_pairs_tested = 1000;
        s.elapsed_secs = 0.001_f64;
        assert!((s.pair_throughput() - 1_000_000.0_f64).abs() < 1.0_f64);
    }
    #[test]
    fn test_stats_accumulate() {
        let mut a = CollisionKernelStats::new();
        a.broadphase_pairs_tested = 5;
        let mut b = CollisionKernelStats::new();
        b.broadphase_pairs_tested = 3;
        a.accumulate(&b);
        assert_eq!(a.broadphase_pairs_tested, 8);
    }
    #[test]
    fn test_gpu_aabb_tree_empty() {
        let mut tree = GpuAabbTree::new();
        tree.build(&[]);
        assert_eq!(tree.leaf_count(), 0);
        assert_eq!(tree.depth(), 0);
    }
    #[test]
    fn test_gpu_aabb_tree_single() {
        let mut tree = GpuAabbTree::new();
        tree.build(&[unit_box()]);
        assert_eq!(tree.leaf_count(), 1);
        assert_eq!(tree.num_primitives, 1);
    }
    #[test]
    fn test_gpu_aabb_tree_two_nodes() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
        ];
        let mut tree = GpuAabbTree::new();
        tree.build(&aabbs);
        assert_eq!(tree.leaf_count(), 2);
        assert_eq!(tree.nodes.len(), 3);
    }
    #[test]
    fn test_gpu_aabb_tree_query_hit() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
        ];
        let mut tree = GpuAabbTree::new();
        tree.build(&aabbs);
        let query = AabbGpu::new([0.2, 0.2, 0.2], [0.8, 0.8, 0.8]);
        let hits = tree.query(&query);
        assert!(!hits.is_empty());
    }
    #[test]
    fn test_gpu_aabb_tree_query_miss() {
        let aabbs = vec![AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])];
        let mut tree = GpuAabbTree::new();
        tree.build(&aabbs);
        let query = AabbGpu::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0]);
        let hits = tree.query(&query);
        assert!(hits.is_empty());
    }
    #[test]
    fn test_gpu_aabb_tree_sah_cost_nonneg() {
        let aabbs = vec![
            unit_box(),
            AabbGpu::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
            AabbGpu::new([4.0, 0.0, 0.0], [5.0, 1.0, 1.0]),
        ];
        let mut tree = GpuAabbTree::new();
        tree.build(&aabbs);
        assert!(tree.sah_cost() >= 0.0_f32);
    }
    #[test]
    fn test_gpu_aabb_tree_depth_single() {
        let mut tree = GpuAabbTree::new();
        tree.build(&[unit_box()]);
        assert_eq!(tree.depth(), 1);
    }
    #[test]
    fn test_gpu_aabb_tree_morton_codes_computed() {
        let aabbs = vec![unit_box(), AabbGpu::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0])];
        let mut tree = GpuAabbTree::new();
        tree.build(&aabbs);
        assert_eq!(tree.morton_codes.len(), 2);
    }
    #[test]
    fn test_sap_broadphase_no_overlap() {
        let mut bp = GpuBroadphase::new(0, 0.0_f32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let pairs = bp.dispatch(&aabbs);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_sap_broadphase_overlap() {
        let mut bp = GpuBroadphase::new(0, 0.0_f32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        let pairs = bp.dispatch(&aabbs);
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_sap_broadphase_margin() {
        let mut bp = GpuBroadphase::new(0, 1.0_f32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]),
        ];
        let pairs = bp.dispatch(&aabbs);
        assert_eq!(pairs.len(), 1);
    }
    #[test]
    fn test_sap_broadphase_stats_updated() {
        let mut bp = GpuBroadphase::new(0, 0.0_f32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        bp.dispatch(&aabbs);
        assert!(bp.stats.broadphase_pairs_tested > 0);
        assert!(bp.stats.bytes_transferred > 0);
    }
    #[test]
    fn test_sap_choose_best_axis_x() {
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [100.0, 1.0, 1.0]),
            AabbGpu::new([200.0, 0.0, 0.0], [300.0, 1.0, 1.0]),
        ];
        let axis = GpuBroadphase::choose_best_axis(&aabbs);
        assert_eq!(axis, 0);
    }
    #[test]
    fn test_sap_choose_best_axis_empty() {
        let axis = GpuBroadphase::choose_best_axis(&[]);
        assert_eq!(axis, 0);
    }
    #[test]
    fn test_gjk_distance_intersecting() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        let res = np.gjk_distance(&a, &b);
        assert!(res.intersecting);
        assert!((res.distance).abs() < 1.0e-5_f32);
    }
    #[test]
    fn test_gjk_distance_separated() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);
        let res = np.gjk_distance(&a, &b);
        assert!(!res.intersecting);
        assert!(res.distance > 0.0_f32);
    }
    #[test]
    fn test_gjk_axis_is_unit() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]);
        let res = np.gjk_distance(&a, &b);
        let len = len3f(res.axis);
        assert!((len - 1.0_f32).abs() < 1.0e-5_f32);
    }
    #[test]
    fn test_sat_depth_overlap() {
        let np = GpuNarrowphase::new(16);
        let a = unit_box();
        let b = AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        let result = np.sat_depth(&a, &b);
        assert!(result.is_some());
        let (depth, _axis) = result.unwrap();
        assert!(depth > 0.0_f32);
    }
    #[test]
    fn test_sat_depth_no_overlap() {
        let np = GpuNarrowphase::new(16);
        let a = unit_box();
        let b = AabbGpu::new([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]);
        assert!(np.sat_depth(&a, &b).is_none());
    }
    #[test]
    fn test_gpu_narrowphase_dispatch() {
        let mut np = GpuNarrowphase::new(32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        let pairs = vec![CollisionPair::new(0, 1)];
        let contacts = np.dispatch(&pairs, &aabbs);
        assert_eq!(contacts.len(), 1);
        assert_eq!(np.stats.narrowphase_contacts, 1);
    }
    #[test]
    fn test_gpu_narrowphase_dispatch_no_contact() {
        let mut np = GpuNarrowphase::new(32);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let pairs = vec![CollisionPair::new(0, 1)];
        let contacts = np.dispatch(&pairs, &aabbs);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_contact_cache_empty() {
        let cache = GpuContactCache::new(3);
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
    }
    #[test]
    fn test_contact_cache_insert_and_find() {
        let mut cache = GpuContactCache::new(3);
        let pair = CollisionPair::new(0, 1);
        let contact = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.5_f32, 0.0_f32, 0.5_f32],
            normal: [0.0_f32, 1.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        cache.insert(&pair, &contact);
        assert_eq!(cache.len(), 1);
        assert!(cache.find(0, 1).is_some());
    }
    #[test]
    fn test_contact_cache_canonical_key() {
        let mut cache = GpuContactCache::new(3);
        let pair = CollisionPair::new(1, 0);
        let contact = ContactResult {
            prim_a: 1,
            prim_b: 0,
            contact_point: [0.0_f32; 3],
            normal: [0.0_f32, 1.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        cache.insert(&pair, &contact);
        assert!(cache.find(0, 1).is_some());
        assert!(cache.find(1, 0).is_some());
    }
    #[test]
    fn test_contact_cache_advance_frame_evicts() {
        let mut cache = GpuContactCache::new(2);
        let pair = CollisionPair::new(0, 1);
        let contact = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0_f32; 3],
            normal: [1.0_f32, 0.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        cache.insert(&pair, &contact);
        cache.advance_frame();
        cache.advance_frame();
        cache.advance_frame();
        assert!(cache.is_empty());
    }
    #[test]
    fn test_contact_cache_clear() {
        let mut cache = GpuContactCache::new(3);
        let pair = CollisionPair::new(0, 1);
        let contact = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0_f32; 3],
            normal: [1.0_f32, 0.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        cache.insert(&pair, &contact);
        cache.clear();
        assert!(cache.is_empty());
    }
    #[test]
    fn test_contact_cache_warmstart() {
        let mut cache = GpuContactCache::new(5);
        let pair = CollisionPair::new(0, 1);
        let contact = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0_f32; 3],
            normal: [1.0_f32, 0.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        cache.insert(&pair, &contact);
        if let Some(entry) = cache.find_mut(0, 1) {
            entry.update(5.0_f32, 1.0_f32, 2.0_f32);
        }
        let raw_contact = ContactResult {
            prim_a: 0,
            prim_b: 1,
            contact_point: [0.0_f32; 3],
            normal: [1.0_f32, 0.0_f32, 0.0_f32],
            depth: 0.1_f32,
        };
        let mut pt = ManifoldPoint::from_contact(&raw_contact);
        if let Some(entry) = cache.find(0, 1) {
            entry.apply_warm_start(&mut pt);
        }
        assert!((pt.warm_impulse_normal - 5.0_f32).abs() < 1.0e-6_f32);
    }
    #[test]
    fn test_gpu_pipeline_no_contacts() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([10.0, 10.0, 10.0], [11.0, 11.0, 11.0]),
        ];
        pipeline.step(&aabbs);
        assert!(pipeline.contacts.is_empty());
    }
    #[test]
    fn test_gpu_pipeline_contact_detected() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.step(&aabbs);
        assert!(!pipeline.contacts.is_empty());
    }
    #[test]
    fn test_gpu_pipeline_manifold_created() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.step(&aabbs);
        assert!(!pipeline.manifolds.is_empty());
    }
    #[test]
    fn test_gpu_pipeline_cache_populated() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.step(&aabbs);
        assert!(!pipeline.contact_cache.is_empty());
    }
    #[test]
    fn test_gpu_pipeline_reset() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.step(&aabbs);
        pipeline.reset();
        assert!(pipeline.contacts.is_empty());
        assert!(pipeline.manifolds.is_empty());
        assert!(pipeline.contact_cache.is_empty());
    }
    #[test]
    fn test_gpu_pipeline_cumulative_stats() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            unit_box(),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
            AabbGpu::new([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        pipeline.step(&aabbs);
        pipeline.step(&aabbs);
        assert!(pipeline.cumulative_stats.bytes_transferred > 0);
    }
    #[test]
    fn test_gpu_pipeline_total_contact_points() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]),
        ];
        pipeline.step(&aabbs);
        assert!(pipeline.total_contact_points() > 0);
    }
    #[test]
    fn test_gpu_pipeline_aabb_tree_rebuilt() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![unit_box()];
        pipeline.step(&aabbs);
        assert_eq!(pipeline.aabb_tree.num_primitives, 1);
    }
    #[test]
    fn test_gpu_pipeline_frame_stats_nonzero() {
        let mut pipeline = GpuCollisionPipeline::new(0.0_f32, 4);
        let aabbs = vec![unit_box(), AabbGpu::new([0.5, 0.0, 0.0], [1.5, 1.0, 1.0])];
        pipeline.step(&aabbs);
        let stats = pipeline.frame_stats();
        assert!(stats.narrowphase_queries > 0);
    }
}
