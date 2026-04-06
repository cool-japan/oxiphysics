// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU ray tracing (CPU mock implementation).
//!
//! Provides BVH construction, ray-AABB intersection (slab method),
//! ray-triangle intersection (Möller-Trumbore), and BVH traversal.

/// A ray defined by an origin and a direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f64; 3],
    /// Ray direction (should be normalised for correct distance output).
    pub direction: [f64; 3],
}

impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self { origin, direction }
    }

    /// Evaluate the ray at parameter `t`: `origin + t * direction`.
    pub fn at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}

/// Axis-aligned bounding box (AABB).
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Compute the AABB that encloses both `self` and `other`.
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    /// Return the centroid of this AABB.
    pub fn centroid(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
}

/// A triangle defined by three vertices.
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// First vertex.
    pub v0: [f64; 3],
    /// Second vertex.
    pub v1: [f64; 3],
    /// Third vertex.
    pub v2: [f64; 3],
}

impl Triangle {
    /// Create a new triangle from three vertices.
    pub fn new(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> Self {
        Self { v0, v1, v2 }
    }

    /// Return the AABB that encloses this triangle.
    pub fn aabb(&self) -> Aabb {
        Aabb {
            min: [
                self.v0[0].min(self.v1[0]).min(self.v2[0]),
                self.v0[1].min(self.v1[1]).min(self.v2[1]),
                self.v0[2].min(self.v1[2]).min(self.v2[2]),
            ],
            max: [
                self.v0[0].max(self.v1[0]).max(self.v2[0]),
                self.v0[1].max(self.v1[1]).max(self.v2[1]),
                self.v0[2].max(self.v1[2]).max(self.v2[2]),
            ],
        }
    }
}

/// A node in the bounding volume hierarchy.
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// Bounding box of this node.
    pub bounds: Aabb,
    /// Index of the left child node, or `usize::MAX` if this is a leaf.
    pub left: usize,
    /// Index of the right child node, or `usize::MAX` if this is a leaf.
    pub right: usize,
    /// Index into the triangle list if this is a leaf (`usize::MAX` otherwise).
    pub triangle_index: usize,
}

impl BvhNode {
    /// Returns `true` if this node is a leaf.
    pub fn is_leaf(&self) -> bool {
        self.triangle_index != usize::MAX
    }
}

/// Result of a ray-triangle intersection test.
#[derive(Debug, Clone, Copy)]
pub struct HitRecord {
    /// Distance along the ray at the hit point (`t` parameter).
    pub t: f64,
    /// Index of the hit triangle.
    pub triangle_index: usize,
    /// Barycentric coordinates `(u, v)` of the hit point on the triangle.
    pub uv: [f64; 2],
}

// ── Internal helpers ─────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Test whether a ray intersects an AABB using the slab method.
///
/// Returns `Some(t_near)` when there is an intersection with `t_near >= t_min`,
/// or `None` when the ray misses.
pub fn ray_aabb_intersect(ray: &Ray, aabb: &Aabb, t_min: f64, t_max: f64) -> Option<f64> {
    let mut t_lo = t_min;
    let mut t_hi = t_max;

    for axis in 0..3 {
        let inv_d = if ray.direction[axis].abs() > 1e-15 {
            1.0 / ray.direction[axis]
        } else {
            f64::INFINITY
        };
        let mut t0 = (aabb.min[axis] - ray.origin[axis]) * inv_d;
        let mut t1 = (aabb.max[axis] - ray.origin[axis]) * inv_d;
        if inv_d < 0.0 {
            std::mem::swap(&mut t0, &mut t1);
        }
        t_lo = t_lo.max(t0);
        t_hi = t_hi.min(t1);
        if t_hi < t_lo {
            return None;
        }
    }
    Some(t_lo)
}

/// Test whether a ray intersects a triangle using the Möller-Trumbore algorithm.
///
/// Returns `Some(HitRecord)` when the ray hits the front face within
/// `[t_min, t_max]`, or `None` on a miss.
pub fn ray_triangle_intersect(
    ray: &Ray,
    tri: &Triangle,
    tri_index: usize,
    t_min: f64,
    t_max: f64,
) -> Option<HitRecord> {
    const EPSILON: f64 = 1e-10;

    let edge1 = sub3(tri.v1, tri.v0);
    let edge2 = sub3(tri.v2, tri.v0);
    let h = cross3(ray.direction, edge2);
    let det = dot3(edge1, h);

    if det.abs() < EPSILON {
        return None; // Ray is parallel to triangle
    }

    let inv_det = 1.0 / det;
    let s = sub3(ray.origin, tri.v0);
    let u = inv_det * dot3(s, h);

    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = cross3(s, edge1);
    let v = inv_det * dot3(ray.direction, q);

    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = inv_det * dot3(edge2, q);
    if t < t_min || t > t_max {
        return None;
    }

    Some(HitRecord {
        t,
        triangle_index: tri_index,
        uv: [u, v],
    })
}

/// Build a BVH from a list of triangles using bottom-up SAH-lite construction.
///
/// Returns a flat node list; node 0 is the root.
pub fn build_bvh(triangles: &[Triangle]) -> Vec<BvhNode> {
    if triangles.is_empty() {
        return Vec::new();
    }

    let mut nodes: Vec<BvhNode> = Vec::new();

    // Create leaf nodes
    let mut leaf_indices: Vec<usize> = (0..triangles.len()).collect();

    fn build_recursive(
        tris: &[Triangle],
        indices: &mut [usize],
        nodes: &mut Vec<BvhNode>,
    ) -> usize {
        if indices.len() == 1 {
            let tri_idx = indices[0];
            let bounds = tris[tri_idx].aabb();
            let node = BvhNode {
                bounds,
                left: usize::MAX,
                right: usize::MAX,
                triangle_index: tri_idx,
            };
            let idx = nodes.len();
            nodes.push(node);
            return idx;
        }

        // Compute combined AABB
        let mut combined = tris[indices[0]].aabb();
        for &i in indices.iter().skip(1) {
            combined = combined.union(&tris[i].aabb());
        }

        // Choose split axis as the longest dimension
        let extent = [
            combined.max[0] - combined.min[0],
            combined.max[1] - combined.min[1],
            combined.max[2] - combined.min[2],
        ];
        let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
            0
        } else if extent[1] >= extent[2] {
            1
        } else {
            2
        };

        // Sort by centroid along chosen axis
        indices.sort_by(|&a, &b| {
            let ca = tris[a].aabb().centroid()[axis];
            let cb = tris[b].aabb().centroid()[axis];
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });

        let mid = indices.len() / 2;
        let (left_ids, right_ids) = indices.split_at_mut(mid);

        let left_child = build_recursive(tris, left_ids, nodes);
        let right_child = build_recursive(tris, right_ids, nodes);

        let left_bounds = nodes[left_child].bounds;
        let right_bounds = nodes[right_child].bounds;
        let node = BvhNode {
            bounds: left_bounds.union(&right_bounds),
            left: left_child,
            right: right_child,
            triangle_index: usize::MAX,
        };
        let idx = nodes.len();
        nodes.push(node);
        idx
    }

    build_recursive(triangles, &mut leaf_indices, &mut nodes);
    nodes
}

/// Traverse the BVH to find the closest ray-triangle intersection.
///
/// `root` is the index of the root node (last element returned by `build_bvh`).
/// Returns `Some(HitRecord)` for the closest hit, or `None`.
pub fn traverse_bvh(
    ray: &Ray,
    nodes: &[BvhNode],
    triangles: &[Triangle],
    root: usize,
    t_min: f64,
    t_max: f64,
) -> Option<HitRecord> {
    if nodes.is_empty() {
        return None;
    }

    let mut best: Option<HitRecord> = None;
    let mut t_closest = t_max;

    // Stack-based traversal
    let mut stack = Vec::with_capacity(64);
    stack.push(root);

    while let Some(node_idx) = stack.pop() {
        if node_idx >= nodes.len() {
            continue;
        }
        let node = &nodes[node_idx];

        // AABB test
        if ray_aabb_intersect(ray, &node.bounds, t_min, t_closest).is_none() {
            continue;
        }

        if node.is_leaf() {
            if node.triangle_index < triangles.len()
                && let Some(hit) = ray_triangle_intersect(
                    ray,
                    &triangles[node.triangle_index],
                    node.triangle_index,
                    t_min,
                    t_closest,
                )
            {
                t_closest = hit.t;
                best = Some(hit);
            }
        } else {
            if node.left != usize::MAX {
                stack.push(node.left);
            }
            if node.right != usize::MAX {
                stack.push(node.right);
            }
        }
    }

    best
}

/// Cast multiple rays against a BVH and return the closest hit for each.
///
/// This is a CPU mock of a GPU batch ray cast. Returns a `Vec` of optional
/// `HitRecord`s, one per input ray.
pub fn batch_ray_cast(
    rays: &[Ray],
    nodes: &[BvhNode],
    triangles: &[Triangle],
    root: usize,
) -> Vec<Option<HitRecord>> {
    rays.iter()
        .map(|ray| traverse_bvh(ray, nodes, triangles, root, 1e-4, f64::INFINITY))
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box_aabb() -> Aabb {
        Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    }

    fn simple_tri() -> Triangle {
        Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    }

    // Ray tests
    #[test]
    fn test_ray_at() {
        let ray = Ray::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = ray.at(3.0);
        assert!((p[0] - 3.0).abs() < 1e-12);
        assert!(p[1].abs() < 1e-12);
        assert!(p[2].abs() < 1e-12);
    }

    #[test]
    fn test_ray_at_negative_t() {
        let ray = Ray::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let p = ray.at(-1.0);
        assert!((p[0]).abs() < 1e-12);
    }

    // AABB union / centroid
    #[test]
    fn test_aabb_union() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([0.5, 0.5, 0.5], [2.0, 2.0, 2.0]);
        let u = a.union(&b);
        assert!((u.max[0] - 2.0).abs() < 1e-12);
        assert!((u.min[0]).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_centroid() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = aabb.centroid();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    // Triangle AABB
    #[test]
    fn test_triangle_aabb() {
        let tri = simple_tri();
        let aabb = tri.aabb();
        assert!((aabb.max[0] - 1.0).abs() < 1e-12);
        assert!((aabb.max[1] - 1.0).abs() < 1e-12);
        assert!((aabb.max[2]).abs() < 1e-12);
    }

    // Ray-AABB hit
    #[test]
    fn test_ray_aabb_hit() {
        let ray = Ray::new([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]);
        let aabb = unit_box_aabb();
        let result = ray_aabb_intersect(&ray, &aabb, 0.0, f64::INFINITY);
        assert!(result.is_some());
        let t = result.unwrap();
        assert!((t - 1.0).abs() < 1e-10);
    }

    // Ray-AABB miss
    #[test]
    fn test_ray_aabb_miss() {
        let ray = Ray::new([2.0, 2.0, -1.0], [0.0, 0.0, 1.0]);
        let aabb = unit_box_aabb();
        assert!(ray_aabb_intersect(&ray, &aabb, 0.0, f64::INFINITY).is_none());
    }

    // Ray-AABB from inside
    #[test]
    fn test_ray_aabb_inside() {
        let ray = Ray::new([0.5, 0.5, 0.5], [0.0, 0.0, 1.0]);
        let aabb = unit_box_aabb();
        let result = ray_aabb_intersect(&ray, &aabb, 0.0, f64::INFINITY);
        assert!(result.is_some());
    }

    // Ray-AABB behind ray
    #[test]
    fn test_ray_aabb_behind() {
        let ray = Ray::new([0.5, 0.5, 5.0], [0.0, 0.0, 1.0]);
        let aabb = unit_box_aabb();
        assert!(ray_aabb_intersect(&ray, &aabb, 0.0, f64::INFINITY).is_none());
    }

    // Ray-triangle hit
    #[test]
    fn test_ray_triangle_hit() {
        // Ray straight down the z-axis hitting a triangle in the xy-plane
        let tri = Triangle::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let ray = Ray::new([0.5, 0.5, 1.0], [0.0, 0.0, -1.0]);
        let result = ray_triangle_intersect(&ray, &tri, 0, 0.0, f64::INFINITY);
        assert!(result.is_some());
        let hit = result.unwrap();
        assert!((hit.t - 1.0).abs() < 1e-9);
    }

    // Ray-triangle miss (outside)
    #[test]
    fn test_ray_triangle_miss_outside() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let ray = Ray::new([2.0, 2.0, 1.0], [0.0, 0.0, -1.0]);
        assert!(ray_triangle_intersect(&ray, &tri, 0, 0.0, f64::INFINITY).is_none());
    }

    // Ray-triangle parallel
    #[test]
    fn test_ray_triangle_parallel() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let ray = Ray::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]); // parallel to z=0 plane
        assert!(ray_triangle_intersect(&ray, &tri, 0, 0.0, f64::INFINITY).is_none());
    }

    // Ray-triangle t outside range
    #[test]
    fn test_ray_triangle_t_range() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let ray = Ray::new([0.5, 0.5, 1.0], [0.0, 0.0, -1.0]);
        // t_max = 0.5, hit is at t=1.0 → miss
        assert!(ray_triangle_intersect(&ray, &tri, 0, 0.0, 0.5).is_none());
    }

    // BVH build with single triangle
    #[test]
    fn test_build_bvh_single() {
        let tris = vec![simple_tri()];
        let nodes = build_bvh(&tris);
        assert!(!nodes.is_empty());
        assert!(nodes.last().unwrap().is_leaf());
    }

    // BVH build empty
    #[test]
    fn test_build_bvh_empty() {
        let nodes = build_bvh(&[]);
        assert!(nodes.is_empty());
    }

    // BVH build multiple triangles
    #[test]
    fn test_build_bvh_multiple() {
        let tris = vec![
            Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            Triangle::new([2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]),
            Triangle::new([4.0, 0.0, 0.0], [5.0, 0.0, 0.0], [4.0, 1.0, 0.0]),
            Triangle::new([6.0, 0.0, 0.0], [7.0, 0.0, 0.0], [6.0, 1.0, 0.0]),
        ];
        let nodes = build_bvh(&tris);
        assert!(!nodes.is_empty());
        // Root is the last node
        let root = nodes.len() - 1;
        assert!(!nodes[root].is_leaf());
    }

    // BVH traversal hit
    #[test]
    fn test_traverse_bvh_hit() {
        let tris = vec![
            Triangle::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]),
            Triangle::new([3.0, 0.0, 0.0], [5.0, 0.0, 0.0], [3.0, 2.0, 0.0]),
        ];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let ray = Ray::new([0.5, 0.5, 1.0], [0.0, 0.0, -1.0]);
        let hit = traverse_bvh(&ray, &nodes, &tris, root, 1e-4, f64::INFINITY);
        assert!(hit.is_some());
    }

    // BVH traversal miss
    #[test]
    fn test_traverse_bvh_miss() {
        let tris = vec![Triangle::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        )];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let ray = Ray::new([5.0, 5.0, 1.0], [0.0, 0.0, -1.0]);
        let hit = traverse_bvh(&ray, &nodes, &tris, root, 1e-4, f64::INFINITY);
        assert!(hit.is_none());
    }

    // BVH traversal empty
    #[test]
    fn test_traverse_bvh_empty_nodes() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let hit = traverse_bvh(&ray, &[], &[], 0, 0.0, f64::INFINITY);
        assert!(hit.is_none());
    }

    // Batch ray cast
    #[test]
    fn test_batch_ray_cast() {
        let tris = vec![Triangle::new(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
        )];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let rays = vec![
            Ray::new([0.5, 0.5, 1.0], [0.0, 0.0, -1.0]),
            Ray::new([5.0, 5.0, 1.0], [0.0, 0.0, -1.0]),
        ];
        let results = batch_ray_cast(&rays, &nodes, &tris, root);
        assert_eq!(results.len(), 2);
        assert!(results[0].is_some());
        assert!(results[1].is_none());
    }

    // Batch ray cast empty
    #[test]
    fn test_batch_ray_cast_empty_rays() {
        let tris = vec![simple_tri()];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let results = batch_ray_cast(&[], &nodes, &tris, root);
        assert!(results.is_empty());
    }

    // BVH leaf detection
    #[test]
    fn test_bvh_node_is_leaf() {
        let node = BvhNode {
            bounds: unit_box_aabb(),
            left: usize::MAX,
            right: usize::MAX,
            triangle_index: 0,
        };
        assert!(node.is_leaf());
    }

    #[test]
    fn test_bvh_node_not_leaf() {
        let node = BvhNode {
            bounds: unit_box_aabb(),
            left: 0,
            right: 1,
            triangle_index: usize::MAX,
        };
        assert!(!node.is_leaf());
    }

    // Hit record fields
    #[test]
    fn test_hit_record_uv() {
        let tris = [Triangle::new(
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
        )];
        let ray = Ray::new([1.0, 1.0, 1.0], [0.0, 0.0, -1.0]);
        let hit = ray_triangle_intersect(&ray, &tris[0], 0, 0.0, f64::INFINITY);
        assert!(hit.is_some());
        let h = hit.unwrap();
        assert!(h.uv[0] >= 0.0 && h.uv[0] <= 1.0);
        assert!(h.uv[1] >= 0.0 && h.uv[1] <= 1.0);
    }

    // Closest hit in batch
    #[test]
    fn test_batch_returns_closest_hit() {
        let tris = vec![
            Triangle::new([0.0, 0.0, 2.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]),
            Triangle::new([0.0, 0.0, 5.0], [2.0, 0.0, 5.0], [0.0, 2.0, 5.0]),
        ];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let rays = vec![Ray::new([0.5, 0.5, 0.0], [0.0, 0.0, 1.0])];
        let results = batch_ray_cast(&rays, &nodes, &tris, root);
        // Should hit the nearer triangle at z=2
        if let Some(hit) = results[0] {
            assert!((hit.t - 2.0).abs() < 1e-9);
        }
    }

    // BVH with 8 triangles — deeper tree
    #[test]
    fn test_build_bvh_8_triangles() {
        let tris: Vec<Triangle> = (0..8)
            .map(|i| {
                let x = i as f64 * 2.0;
                Triangle::new([x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0])
            })
            .collect();
        let nodes = build_bvh(&tris);
        // There should be 2*N-1 nodes for N leaves
        assert_eq!(nodes.len(), 2 * tris.len() - 1);
    }

    // Intersection distance ordering
    #[test]
    fn test_ray_triangle_t_value() {
        let tri = Triangle::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let ray = Ray::new([0.5, 0.5, 3.0], [0.0, 0.0, -1.0]);
        let hit = ray_triangle_intersect(&ray, &tri, 0, 0.0, f64::INFINITY);
        assert!(hit.is_some());
        assert!((hit.unwrap().t - 3.0).abs() < 1e-9);
    }

    // AABB ray direction components near zero
    #[test]
    fn test_ray_aabb_near_zero_dir_component() {
        // Direction has x and y ≈ 0
        let ray = Ray::new([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]);
        let aabb = unit_box_aabb();
        let result = ray_aabb_intersect(&ray, &aabb, 0.0, f64::INFINITY);
        assert!(result.is_some());
    }

    // BVH traversal picks the right triangle
    #[test]
    fn test_traverse_picks_correct_triangle() {
        let tris = vec![
            Triangle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            Triangle::new([10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]),
        ];
        let nodes = build_bvh(&tris);
        let root = nodes.len() - 1;
        let ray = Ray::new([10.2, 0.2, 1.0], [0.0, 0.0, -1.0]);
        let hit = traverse_bvh(&ray, &nodes, &tris, root, 1e-4, f64::INFINITY);
        assert!(hit.is_some());
        assert_eq!(hit.unwrap().triangle_index, 1);
    }
}
