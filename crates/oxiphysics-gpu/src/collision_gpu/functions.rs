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
/// Negate an f32 3-vector.
pub(super) fn neg3f(a: [f32; 3]) -> [f32; 3] {
    [-a[0], -a[1], -a[2]]
}
/// Squared sine of the angle below which two directions count as parallel
/// (`sin ≈ 1e-5`, a few hundred f32 ulps).
const PARALLEL_SIN2_EPS: f32 = 1.0e-10;
/// Iteration budget of the boolean GJK query used by [`NarrowphaseGpuKernel`].
///
/// The Minkowski difference of two boxes is a box with eight vertices, so a
/// well-behaved query ends after a handful of iterations; the budget only
/// matters for degenerate (touching or numerically flat) configurations.
///
/// [`NarrowphaseGpuKernel`]: super::types::NarrowphaseGpuKernel
pub(super) const GJK_MAX_ITERATIONS: usize = 32;
/// Classification of a GJK search direction.
enum SearchDirection {
    /// The direction vanished: the origin lies on the current simplex.
    Zero,
    /// A unit-length search direction.
    Unit([f32; 3]),
    /// The direction is not a finite number (degenerate input).
    Invalid,
}
/// Classify and normalise a raw search direction.
///
/// The vector is rescaled by its largest component before the length is
/// taken, so tiny or huge scenes neither underflow nor overflow the squared
/// length.
fn classify_direction(d: [f32; 3]) -> SearchDirection {
    let m = d[0].abs().max(d[1].abs()).max(d[2].abs());
    if m == 0.0 {
        return SearchDirection::Zero;
    }
    if !m.is_finite() {
        return SearchDirection::Invalid;
    }
    let scaled = scale3f(d, 1.0 / m);
    SearchDirection::Unit(scale3f(scaled, 1.0 / len3f(scaled)))
}
/// Reduce a 2-point simplex `[B, A]` (`A` newest) to the Voronoi feature of the
/// segment that is nearest to the origin.
///
/// Returns `(true, _)` when the origin lies on the segment, otherwise
/// `(false, direction)` with `direction` pointing from the kept feature
/// towards the origin.
fn do_simplex_line(simplex: &mut Vec<SupportPoint>) -> (bool, [f32; 3]) {
    let &[b, a] = simplex.as_slice() else {
        return (false, [1.0, 0.0, 0.0]);
    };
    let ab = sub3f(b.v, a.v);
    let ao = neg3f(a.v);
    if dot3f(ab, ao) <= 0.0 {
        simplex.clear();
        simplex.push(a);
        return (false, ao);
    }
    let bo = neg3f(b.v);
    if dot3f(neg3f(ab), bo) <= 0.0 {
        simplex.clear();
        simplex.push(b);
        return (false, bo);
    }
    let normal = cross3f(ab, ao);
    if len_sq3f(normal) <= PARALLEL_SIN2_EPS * len_sq3f(ab) * len_sq3f(ao) {
        return (true, [0.0; 3]);
    }
    (false, cross3f(normal, ab))
}
/// Reduce a 3-point simplex `[C, B, A]` (`A` newest) to the Voronoi feature of
/// the triangle that is nearest to the origin.
///
/// Edges and vertices outside the face region drop the vertices that cannot
/// be nearest.  When the origin is on the far side of the face the winding is
/// flipped so that the face normal keeps pointing at the origin.
fn do_simplex_triangle(simplex: &mut Vec<SupportPoint>) -> (bool, [f32; 3]) {
    let &[c, b, a] = simplex.as_slice() else {
        return (false, [1.0, 0.0, 0.0]);
    };
    let ab = sub3f(b.v, a.v);
    let ac = sub3f(c.v, a.v);
    let ao = neg3f(a.v);
    let abc = cross3f(ab, ac);
    if len_sq3f(abc) <= PARALLEL_SIN2_EPS * len_sq3f(ab) * len_sq3f(ac) {
        let bc = sub3f(c.v, b.v);
        let (len_ab, len_ac, len_bc) = (len_sq3f(ab), len_sq3f(ac), len_sq3f(bc));
        let pair = if len_bc >= len_ab && len_bc >= len_ac {
            [b, c]
        } else if len_ac >= len_ab {
            [a, c]
        } else {
            [a, b]
        };
        simplex.clear();
        simplex.extend_from_slice(&pair);
        return do_simplex_line(simplex);
    }
    if dot3f(cross3f(abc, ac), ao) > 0.0 {
        if dot3f(ac, ao) > 0.0 {
            simplex.clear();
            simplex.extend_from_slice(&[c, a]);
            (false, cross3f(cross3f(ac, ao), ac))
        } else {
            simplex.clear();
            simplex.extend_from_slice(&[b, a]);
            do_simplex_line(simplex)
        }
    } else if dot3f(cross3f(ab, abc), ao) > 0.0 {
        simplex.clear();
        simplex.extend_from_slice(&[b, a]);
        do_simplex_line(simplex)
    } else if dot3f(abc, ao) > 0.0 {
        (false, abc)
    } else {
        simplex.clear();
        simplex.extend_from_slice(&[b, c, a]);
        (false, neg3f(abc))
    }
}
/// Reduce a 4-point simplex `[D, C, B, A]` (`A` newest) to the face that is
/// nearest to the origin, or report that the origin is enclosed.
///
/// The origin is enclosed only if it lies behind all three faces through the
/// newest vertex `A`; the fourth face (`BCD`) was already rejected when the
/// previous triangle was reduced.
fn do_simplex_tetrahedron(simplex: &mut Vec<SupportPoint>) -> (bool, [f32; 3]) {
    let &[d, c, b, a] = simplex.as_slice() else {
        return (false, [1.0, 0.0, 0.0]);
    };
    let ao = neg3f(a.v);
    let ab = sub3f(b.v, a.v);
    let ac = sub3f(c.v, a.v);
    let ad = sub3f(d.v, a.v);
    let abc = cross3f(ab, ac);
    let acd = cross3f(ac, ad);
    let adb = cross3f(ad, ab);
    if dot3f(abc, ao) > 0.0 {
        simplex.clear();
        simplex.extend_from_slice(&[c, b, a]);
        do_simplex_triangle(simplex)
    } else if dot3f(acd, ao) > 0.0 {
        simplex.clear();
        simplex.extend_from_slice(&[d, c, a]);
        do_simplex_triangle(simplex)
    } else if dot3f(adb, ao) > 0.0 {
        simplex.clear();
        simplex.extend_from_slice(&[b, d, a]);
        do_simplex_triangle(simplex)
    } else {
        (true, [0.0; 3])
    }
}
/// Refine the GJK simplex and return `(contains_origin, new_direction)`.
///
/// `simplex` is ordered oldest first, so its last element is the support point
/// added most recently.  Every case performs the Voronoi-region tests of the
/// matching sub-simplex and *drops the vertices that cannot be nearest to the
/// origin*; the returned direction points from the kept feature towards the
/// origin.  `contains_origin` is `true` only if the origin lies on a segment of
/// the simplex or behind all faces of a tetrahedron.
pub(super) fn update_simplex(simplex: &mut Vec<SupportPoint>) -> (bool, [f32; 3]) {
    match simplex.len() {
        0 => (false, [1.0, 0.0, 0.0]),
        1 => {
            let a = simplex[0].v;
            (false, neg3f(a))
        }
        2 => do_simplex_line(simplex),
        3 => do_simplex_triangle(simplex),
        4 => do_simplex_tetrahedron(simplex),
        len => {
            simplex.drain(..len - 4);
            do_simplex_tetrahedron(simplex)
        }
    }
}
/// Result of the iterative GJK boolean query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum GjkVerdict {
    /// The simplex encloses (or touches) the origin.
    Intersecting,
    /// A plane through the origin separates the shapes.
    Separated,
    /// The iteration budget ran out or the simplex degenerated.
    Inconclusive,
}
/// Run GJK on the Minkowski difference `a − b` of two boxes.
///
/// Each iteration asks for the support point of `a − b` along the current
/// direction.  If that point does not reach the origin, the plane through it
/// separates the shapes.  Otherwise it joins the simplex, which is then
/// reduced by [`update_simplex`] to the feature nearest to the origin.
pub(super) fn gjk_verdict(a: &AabbGpu, b: &AabbGpu, max_iterations: usize) -> GjkVerdict {
    let mut d = match classify_direction(sub3f(a.centre(), b.centre())) {
        SearchDirection::Unit(u) => u,
        SearchDirection::Zero => [1.0, 0.0, 0.0],
        SearchDirection::Invalid => return GjkVerdict::Inconclusive,
    };
    let first = minkowski_support(a, b, d);
    let mut simplex = vec![first];
    d = match classify_direction(neg3f(first.v)) {
        SearchDirection::Unit(u) => u,
        SearchDirection::Zero => return GjkVerdict::Intersecting,
        SearchDirection::Invalid => return GjkVerdict::Inconclusive,
    };
    for _ in 0..max_iterations {
        let support = minkowski_support(a, b, d);
        if dot3f(support.v, d) < 0.0 {
            return GjkVerdict::Separated;
        }
        if simplex.iter().any(|p| p.v == support.v) {
            return GjkVerdict::Inconclusive;
        }
        simplex.push(support);
        let (contains_origin, next) = update_simplex(&mut simplex);
        if contains_origin {
            return GjkVerdict::Intersecting;
        }
        d = match classify_direction(next) {
            SearchDirection::Unit(u) => u,
            SearchDirection::Zero => return GjkVerdict::Intersecting,
            SearchDirection::Invalid => return GjkVerdict::Inconclusive,
        };
    }
    GjkVerdict::Inconclusive
}
/// Closest pair of points of two boxes, computed per axis.
///
/// Along an axis where the slabs are disjoint the points are the facing
/// bounds; where they overlap both points take the middle of the shared
/// interval, so the pair realises the exact box-to-box distance.
pub(super) fn aabb_closest_points(a: &AabbGpu, b: &AabbGpu) -> ([f32; 3], [f32; 3]) {
    let mut closest_a = [0.0_f32; 3];
    let mut closest_b = [0.0_f32; 3];
    for k in 0..3 {
        if a.max[k] < b.min[k] {
            closest_a[k] = a.max[k];
            closest_b[k] = b.min[k];
        } else if b.max[k] < a.min[k] {
            closest_a[k] = a.min[k];
            closest_b[k] = b.max[k];
        } else {
            let lo = a.min[k].max(b.min[k]);
            let hi = a.max[k].min(b.max[k]);
            let mid = 0.5 * (lo + hi);
            closest_a[k] = mid;
            closest_b[k] = mid;
        }
    }
    (closest_a, closest_b)
}
/// Result of [`aabb_gjk_query`].
#[derive(Debug, Clone, Copy)]
pub(super) struct AabbGjkOutcome {
    /// Whether the boxes intersect (touching counts as intersecting).
    pub(super) intersecting: bool,
    /// Exact distance between the boxes (`0.0` if they intersect).
    pub(super) distance: f32,
    /// Closest point on box `a` (box centre if the boxes intersect).
    pub(super) closest_a: [f32; 3],
    /// Closest point on box `b` (box centre if the boxes intersect).
    pub(super) closest_b: [f32; 3],
}
/// GJK intersection and distance query for a pair of boxes.
///
/// The intersection verdict comes from [`gjk_verdict`]; if the iteration is
/// inconclusive (budget exhausted, degenerate simplex) the exact box overlap
/// test decides.  A separation certificate is trusted only if the exact gap
/// between the boxes is positive: boxes that merely touch intersect, as in
/// [`AabbGpu::overlaps`].  The distance of separated boxes is the exact norm of
/// the per-axis gap `max(0, b.min − a.max, a.min − b.max)`.
pub(super) fn aabb_gjk_query(a: &AabbGpu, b: &AabbGpu, max_iterations: usize) -> AabbGjkOutcome {
    let verdict = gjk_verdict(a, b, max_iterations);
    let intersecting = match verdict {
        GjkVerdict::Intersecting => true,
        GjkVerdict::Separated => false,
        GjkVerdict::Inconclusive => a.overlaps(b),
    };
    if !intersecting {
        let (closest_a, closest_b) = aabb_closest_points(a, b);
        let distance = len3f(sub3f(closest_a, closest_b));
        if distance != 0.0 {
            return AabbGjkOutcome {
                intersecting: false,
                distance,
                closest_a,
                closest_b,
            };
        }
    }
    AabbGjkOutcome {
        intersecting: true,
        distance: 0.0,
        closest_a: a.centre(),
        closest_b: b.centre(),
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
        // The gap is 4 along every axis: distance 4 * sqrt(3).  A search
        // direction that is not normalised used to report 6 * sqrt(3).
        let expected = 4.0_f32 * 3.0_f32.sqrt();
        assert!(
            (dist - expected).abs() < 1.0e-5,
            "distance {dist} != {expected}"
        );
    }
    #[test]
    fn test_gjk_single_axis_gap_is_reported_exactly() {
        // A = [0, 1]^3 and B = [1.1, 2.1] x [0, 1]^2 are 0.1 apart along x.  The
        // Minkowski-difference search used to report 2.1 here.
        let kernel = NarrowphaseGpuKernel::new();
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([1.1, 0.0, 0.0], [2.1, 1.0, 1.0]);
        let (intersect, dist) = kernel.gjk_intersect(&a, &b);
        assert!(!intersect);
        assert!((dist - 0.1).abs() < 1.0e-6, "distance {dist}");
        let res = np.gjk_distance(&a, &b);
        assert!(!res.intersecting);
        assert!(
            (res.distance - 0.1).abs() < 1.0e-6,
            "distance {}",
            res.distance
        );
        // The same query with the shapes swapped.
        let (intersect_ba, dist_ba) = kernel.gjk_intersect(&b, &a);
        assert!(!intersect_ba);
        assert!((dist_ba - 0.1).abs() < 1.0e-6, "distance {dist_ba}");
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
        let expected = 4.0_f32 * 3.0_f32.sqrt();
        assert!(
            (res.distance - expected).abs() < 1.0e-5_f32,
            "distance {} != {expected}",
            res.distance
        );
        for k in 0..3 {
            assert_eq!(res.closest_a[k], 1.0);
            assert_eq!(res.closest_b[k], 5.0);
        }
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
    // ── GJK simplex handling and distance (F1865) ──────────────────────────

    /// Deterministic generator of `f32` values in `[0, 1)`.
    struct Lcg(u64);
    impl Lcg {
        fn next_unit(&mut self) -> f32 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 40) as f32) / ((1u64 << 24) as f32)
        }
        fn random_box(&mut self) -> AabbGpu {
            let mut min = [0.0_f32; 3];
            let mut max = [0.0_f32; 3];
            for k in 0..3 {
                min[k] = self.next_unit() * 4.0 - 2.0;
                max[k] = min[k] + 0.2 + self.next_unit() * 2.0;
            }
            AabbGpu::new(min, max)
        }
    }
    /// Exact box-to-box distance in double precision.
    fn exact_gap(a: &AabbGpu, b: &AabbGpu) -> f64 {
        let mut sq = 0.0_f64;
        for k in 0..3 {
            let g = (f64::from(b.min[k]) - f64::from(a.max[k]))
                .max(f64::from(a.min[k]) - f64::from(b.max[k]))
                .max(0.0);
            sq += g * g;
        }
        sq.sqrt()
    }
    /// Smallest per-axis overlap depth (negative if the boxes are separated).
    fn exact_depth(a: &AabbGpu, b: &AabbGpu) -> f64 {
        (0..3)
            .map(|k| f64::from(a.max[k].min(b.max[k])) - f64::from(a.min[k].max(b.min[k])))
            .fold(f64::INFINITY, f64::min)
    }
    /// Order `points` as a GJK simplex (oldest first, newest last) whose
    /// tetrahedron has outward-facing face normals, the winding the
    /// tetrahedron case relies on.
    fn tetra_simplex(points: [[f32; 3]; 4]) -> Vec<SupportPoint> {
        let [p0, p1, p2, p3] = points;
        let a = p3;
        let (mut b, mut c, d) = (p0, p1, p2);
        let volume = dot3f(cross3f(sub3f(b, a), sub3f(c, a)), sub3f(d, a));
        if volume > 0.0 {
            std::mem::swap(&mut b, &mut c);
        }
        vec![
            SupportPoint::new(d),
            SupportPoint::new(c),
            SupportPoint::new(b),
            SupportPoint::new(a),
        ]
    }
    #[test]
    fn test_gjk_separated_cubes_report_exact_distance() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let b = AabbGpu::new([1.2, 1.2, 1.2], [2.0, 2.0, 2.0]);
        let (intersecting, distance) = kernel.gjk_intersect(&a, &b);
        assert!(!intersecting);
        let expected = 0.2_f32 * 3.0_f32.sqrt();
        assert!(
            (distance - expected).abs() < 1.0e-5,
            "distance {distance} != {expected}"
        );
    }
    #[test]
    fn test_gjk_distance_separated_cubes_closest_points_and_axis() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([1.2, 1.2, 1.2], [2.0, 2.0, 2.0]);
        let res = np.gjk_distance(&a, &b);
        assert!(!res.intersecting);
        assert!((res.distance - 0.2 * 3.0_f32.sqrt()).abs() < 1.0e-5);
        for k in 0..3 {
            assert!((res.closest_a[k] - 1.0).abs() < 1.0e-6);
            assert!((res.closest_b[k] - 1.2).abs() < 1.0e-6);
        }
        // The axis points from B towards A.
        let expected = -1.0 / 3.0_f32.sqrt();
        for k in 0..3 {
            assert!(
                (res.axis[k] - expected).abs() < 1.0e-5,
                "axis {:?}",
                res.axis
            );
        }
    }
    #[test]
    fn test_gjk_diagonal_separation_distance_and_shared_axis_midpoint() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([1.5, 1.5, 0.5], [2.5, 2.5, 1.5]);
        let res = np.gjk_distance(&a, &b);
        assert!(!res.intersecting);
        assert!((res.distance - 0.5_f32.sqrt()).abs() < 1.0e-5);
        // z-slabs overlap on [0.5, 1.0]: both witnesses take the middle.
        assert!((res.closest_a[2] - 0.75).abs() < 1.0e-6);
        assert!((res.closest_b[2] - 0.75).abs() < 1.0e-6);
    }
    #[test]
    fn test_gjk_touching_boxes_intersect_with_zero_distance() {
        let kernel = NarrowphaseGpuKernel::new();
        let a = unit_box();
        let face = AabbGpu::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        let edge = AabbGpu::new([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]);
        let vertex = AabbGpu::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        for b in [&face, &edge, &vertex] {
            let (intersecting, distance) = kernel.gjk_intersect(&a, b);
            assert!(intersecting, "touching boxes {b:?} must intersect");
            assert_eq!(distance, 0.0);
        }
    }
    #[test]
    fn test_gjk_nested_boxes_intersect() {
        let kernel = NarrowphaseGpuKernel::new();
        let outer = AabbGpu::new([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]);
        let inner = AabbGpu::new([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]);
        assert!(kernel.gjk_intersect(&outer, &inner).0);
        assert!(kernel.gjk_intersect(&inner, &outer).0);
        // Concentric boxes: the centre difference is zero.
        assert!(kernel.gjk_intersect(&inner, &inner).0);
    }
    #[test]
    fn test_gjk_distance_is_symmetric_up_to_axis_sign() {
        let np = GpuNarrowphase::new(32);
        let a = unit_box();
        let b = AabbGpu::new([3.0, -2.0, 0.25], [4.0, -1.0, 0.75]);
        let ab = np.gjk_distance(&a, &b);
        let ba = np.gjk_distance(&b, &a);
        assert!(!ab.intersecting && !ba.intersecting);
        assert!((ab.distance - ba.distance).abs() < 1.0e-6);
        for k in 0..3 {
            assert!((ab.axis[k] + ba.axis[k]).abs() < 1.0e-6);
        }
    }
    #[test]
    fn test_raw_gjk_verdict_agrees_with_exact_overlap_on_random_pairs() {
        let mut rng = Lcg(0x5EED_1234_ABCD_0001);
        let mut inconclusive = 0usize;
        let mut decided = 0usize;
        for _ in 0..4000 {
            let a = rng.random_box();
            let b = rng.random_box();
            let gap = exact_gap(&a, &b);
            let depth = exact_depth(&a, &b);
            // Skip configurations within f32 rounding of touching.
            if gap < 1.0e-3 && depth < 1.0e-3 {
                continue;
            }
            match gjk_verdict(&a, &b, GJK_MAX_ITERATIONS) {
                GjkVerdict::Intersecting => {
                    assert!(a.overlaps(&b), "false intersection: {a:?} vs {b:?}");
                    decided += 1;
                }
                GjkVerdict::Separated => {
                    assert!(!a.overlaps(&b), "false separation: {a:?} vs {b:?}");
                    decided += 1;
                }
                GjkVerdict::Inconclusive => inconclusive += 1,
            }
        }
        assert!(decided > 3000, "only {decided} decided queries");
        assert!(
            inconclusive * 100 <= decided,
            "{inconclusive} inconclusive vs {decided} decided"
        );
    }
    #[test]
    fn test_gjk_queries_match_exact_box_distance_on_random_pairs() {
        let kernel = NarrowphaseGpuKernel::new();
        let np = GpuNarrowphase::new(32);
        let mut rng = Lcg(0xFEED_0000_7777_0002);
        let mut separated = 0usize;
        let mut overlapping = 0usize;
        for _ in 0..4000 {
            let a = rng.random_box();
            let b = rng.random_box();
            let gap = exact_gap(&a, &b);
            let depth = exact_depth(&a, &b);
            if gap < 1.0e-3 && depth < 1.0e-3 {
                continue;
            }
            let (intersecting, distance) = kernel.gjk_intersect(&a, &b);
            let res = np.gjk_distance(&a, &b);
            assert_eq!(intersecting, a.overlaps(&b));
            assert_eq!(res.intersecting, intersecting);
            assert!(
                (f64::from(distance) - gap).abs() < 1.0e-4,
                "distance {distance} vs exact {gap}"
            );
            assert!((f64::from(res.distance) - gap).abs() < 1.0e-4);
            if intersecting {
                overlapping += 1;
            } else {
                separated += 1;
                let witness = len3f(sub3f(res.closest_a, res.closest_b));
                assert!((f64::from(witness) - gap).abs() < 1.0e-4);
            }
        }
        assert!(
            separated > 500 && overlapping > 200,
            "{separated} separated / {overlapping} overlapping pairs"
        );
    }
    #[test]
    fn test_gjk_grid_aligned_boxes_match_exact_overlap() {
        // Exhaustive lattice of integer offsets: face, edge and vertex contacts
        // all occur exactly (gap == 0).
        let np = GpuNarrowphase::new(32);
        let a = AabbGpu::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let mut touching = 0usize;
        for size in [1.0_f32, 2.0] {
            for ix in -3..=4 {
                for iy in -3..=4 {
                    for iz in -3..=4 {
                        let min = [ix as f32, iy as f32, iz as f32];
                        let b = AabbGpu::new(min, [min[0] + size, min[1] + size, min[2] + size]);
                        let gap = exact_gap(&a, &b);
                        let res = np.gjk_distance(&a, &b);
                        assert_eq!(res.intersecting, a.overlaps(&b), "verdict for {b:?}");
                        assert!(
                            (f64::from(res.distance) - gap).abs() < 1.0e-6,
                            "distance for {b:?}: {} vs {gap}",
                            res.distance
                        );
                        if gap == 0.0 && exact_depth(&a, &b) == 0.0 {
                            touching += 1;
                        }
                        // A conclusive raw verdict never contradicts the exact test.
                        match gjk_verdict(&a, &b, GJK_MAX_ITERATIONS) {
                            GjkVerdict::Intersecting => assert!(a.overlaps(&b), "{b:?}"),
                            GjkVerdict::Separated => assert!(!a.overlaps(&b), "{b:?}"),
                            GjkVerdict::Inconclusive => {}
                        }
                    }
                }
            }
        }
        assert!(
            touching > 20,
            "lattice produced only {touching} touching cases"
        );
    }
    #[test]
    fn test_gjk_exhausted_iteration_budget_falls_back_to_exact_test() {
        // The smallest allowed budget cannot finish every query; the verdict
        // must still equal the exact box test and never default to "intersecting".
        let np = GpuNarrowphase::new(1);
        assert_eq!(np.max_iterations, 4);
        let mut rng = Lcg(0x0DD_BA11_0003);
        for _ in 0..1500 {
            let a = rng.random_box();
            let b = rng.random_box();
            let res = np.gjk_distance(&a, &b);
            assert_eq!(res.intersecting, a.overlaps(&b));
            assert!((f64::from(res.distance) - exact_gap(&a, &b)).abs() < 1.0e-4);
        }
    }
    #[test]
    fn test_update_simplex_tetrahedron_enclosing_origin() {
        let mut simplex = tetra_simplex([
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ]);
        let (contains, _) = update_simplex(&mut simplex);
        assert!(contains);
    }
    #[test]
    fn test_update_simplex_tetrahedron_not_enclosing_origin() {
        // Every vertex has x >= 1: the origin lies beyond the face x = 1, which
        // is spanned by the vertices at index 0, 2 and 3.  GJK guarantees that
        // the face opposite the newest vertex was rejected earlier, so the apex
        // (index 1) is never the newest vertex of a real iteration.
        for newest in [0usize, 2, 3] {
            let mut points = [
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 1.0],
            ];
            points.swap(newest, 3);
            let mut simplex = tetra_simplex(points);
            let (contains, direction) = update_simplex(&mut simplex);
            assert!(!contains, "newest vertex {newest}: origin is outside");
            assert!(
                simplex.len() <= 3,
                "a face must replace the tetrahedron, got {} vertices",
                simplex.len()
            );
            for vertex in &simplex {
                assert!(
                    dot3f(direction, vertex.v) < 0.0,
                    "search direction {direction:?} must point from {:?} towards the origin",
                    vertex.v
                );
            }
        }
    }
    #[test]
    fn test_update_simplex_triangle_vertex_region_keeps_only_that_vertex() {
        // Origin behind vertex A = (1, 0, 0): the nearest feature is A alone.
        let mut simplex = vec![
            SupportPoint::new([2.0, -1.0, 0.0]),
            SupportPoint::new([2.0, 1.0, 0.0]),
            SupportPoint::new([1.0, 0.0, 0.0]),
        ];
        let (contains, direction) = update_simplex(&mut simplex);
        assert!(!contains);
        assert_eq!(simplex.len(), 1);
        assert_eq!(simplex[0].v, [1.0, 0.0, 0.0]);
        assert_eq!(direction, [-1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_update_simplex_triangle_edge_region_keeps_the_edge() {
        // Triangle in the plane x = 1; the origin projects across edge AB.
        let c = [1.0, 3.0, 0.0];
        let b = [1.0, 1.0, 1.0];
        let a = [1.0, 1.0, -1.0];
        let mut simplex = vec![
            SupportPoint::new(c),
            SupportPoint::new(b),
            SupportPoint::new(a),
        ];
        let (contains, direction) = update_simplex(&mut simplex);
        assert!(!contains);
        assert_eq!(simplex.len(), 2, "the far vertex C must be dropped");
        assert_eq!(simplex[0].v, b);
        assert_eq!(simplex[1].v, a);
        // Perpendicular to AB and towards the origin: parallel to (-1, -1, 0).
        assert_eq!(direction[2], 0.0);
        assert!(direction[0] < 0.0 && (direction[0] - direction[1]).abs() < 1.0e-6);
    }
    #[test]
    fn test_update_simplex_line_projection_inside_segment() {
        let mut simplex = vec![
            SupportPoint::new([1.0, 1.0, 0.0]),
            SupportPoint::new([1.0, -1.0, 0.0]),
        ];
        let (contains, direction) = update_simplex(&mut simplex);
        assert!(!contains);
        assert_eq!(simplex.len(), 2);
        assert!(direction[0] < 0.0);
        assert_eq!(direction[1], 0.0);
        assert_eq!(direction[2], 0.0);
    }
    #[test]
    fn test_update_simplex_line_origin_on_segment() {
        let mut simplex = vec![
            SupportPoint::new([-1.0, 0.0, 0.0]),
            SupportPoint::new([2.0, 0.0, 0.0]),
        ];
        assert!(update_simplex(&mut simplex).0);
    }
    #[test]
    fn test_update_simplex_line_origin_beyond_segment_end() {
        let mut simplex = vec![
            SupportPoint::new([2.0, 0.0, 0.0]),
            SupportPoint::new([1.0, 0.0, 0.0]),
        ];
        let (contains, direction) = update_simplex(&mut simplex);
        assert!(!contains);
        assert_eq!(simplex.len(), 1);
        assert_eq!(simplex[0].v, [1.0, 0.0, 0.0]);
        assert_eq!(direction, [-1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_update_simplex_single_point_points_to_origin() {
        let mut simplex = vec![SupportPoint::new([0.0, 3.0, 0.0])];
        let (contains, direction) = update_simplex(&mut simplex);
        assert!(!contains);
        assert_eq!(direction, [0.0, -3.0, 0.0]);
    }
    // ── Contact cache keys (F1868) ─────────────────────────────────────────

    #[test]
    fn test_gpu_pipeline_cache_keys_follow_surviving_contacts() {
        // Box 1 lies 0.05 beyond box 0: inside the 0.1 broadphase margin, so the
        // pair (0, 1) reaches the narrowphase, which rejects it.  The surviving
        // contacts (0, 2) and (1, 2) are therefore *not* index-aligned with the
        // broadphase pair list.
        let mut pipeline = GpuCollisionPipeline::new(0.1_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([1.05, 0.0, 0.0], [2.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
        ];
        let pairs = pipeline.broadphase.dispatch(&aabbs);
        assert!(
            pairs.contains(&CollisionPair::new(0, 1)),
            "margin must let the near-miss pair through the broadphase"
        );
        pipeline.step(&aabbs);
        assert_eq!(pipeline.contacts.len(), 2);
        assert_eq!(pipeline.contact_cache.len(), pipeline.contacts.len());
        for entry in &pipeline.contact_cache.entries {
            let (a, b) = entry.key;
            assert!(
                aabbs[a as usize].overlaps(&aabbs[b as usize]),
                "cache entry {:?} does not belong to an overlapping pair",
                entry.key
            );
        }
        assert!(pipeline.contact_cache.find(0, 1).is_none());
        assert!(pipeline.contact_cache.find(0, 2).is_some());
        assert!(pipeline.contact_cache.find(1, 2).is_some());
    }
    #[test]
    fn test_gpu_pipeline_warm_start_lands_on_the_right_manifold() {
        let mut pipeline = GpuCollisionPipeline::new(0.1_f32, 4);
        let aabbs = vec![
            AabbGpu::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            AabbGpu::new([1.05, 0.0, 0.0], [2.0, 1.0, 1.0]),
            AabbGpu::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
        ];
        pipeline.step(&aabbs);
        let entry = pipeline
            .contact_cache
            .find_mut(1, 2)
            .expect("pair (1, 2) is in contact");
        entry.update(7.0_f32, 0.0_f32, 0.0_f32);
        pipeline.manifolds.clear();
        pipeline.step(&aabbs);
        for manifold in &pipeline.manifolds {
            let key = if manifold.body_a <= manifold.body_b {
                (manifold.body_a, manifold.body_b)
            } else {
                (manifold.body_b, manifold.body_a)
            };
            let expected = if key == (1, 2) { 7.0_f32 } else { 0.0_f32 };
            for point in manifold.active_points() {
                assert!(
                    (point.warm_impulse_normal - expected).abs() < 1.0e-6,
                    "manifold {key:?} carries warm impulse {}",
                    point.warm_impulse_normal
                );
            }
        }
        assert!(pipeline.manifolds.iter().any(|m| m.num_points > 0));
    }
}
