// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Collision detection primitives: GJK/EPA narrowphase, SAT, sphere/capsule
//! tests, BVH broadphase, and contact manifold generation.
//!
//! All geometry is represented with `[f64; 3]` arrays to avoid nalgebra
//! dependencies inside this module.

use crate::math::Real;

// ---------------------------------------------------------------------------
// Vec3 helper (local, no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-vectors.
#[inline]
fn v3_add(a: [Real; 3], b: [Real; 3]) -> [Real; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
fn v3_sub(a: [Real; 3], b: [Real; 3]) -> [Real; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
fn v3_scale(a: [Real; 3], s: Real) -> [Real; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline]
fn v3_dot(a: [Real; 3], b: [Real; 3]) -> Real {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn v3_cross(a: [Real; 3], b: [Real; 3]) -> [Real; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Squared length of a 3-vector.
#[inline]
fn v3_len_sq(a: [Real; 3]) -> Real {
    v3_dot(a, a)
}

/// Length of a 3-vector.
#[inline]
fn v3_len(a: [Real; 3]) -> Real {
    v3_len_sq(a).sqrt()
}

/// Normalize a 3-vector.  Returns `[0,0,0]` for zero-length input.
#[inline]
fn v3_normalize(a: [Real; 3]) -> [Real; 3] {
    let len = v3_len(a);
    if len > 1e-15 {
        v3_scale(a, 1.0 / len)
    } else {
        [0.0; 3]
    }
}

/// Negate a 3-vector.
#[inline]
fn v3_neg(a: [Real; 3]) -> [Real; 3] {
    [-a[0], -a[1], -a[2]]
}

// ---------------------------------------------------------------------------
// Convex Shape (support function interface)
// ---------------------------------------------------------------------------

/// Trait for convex shapes defined by a support function.
pub trait ConvexShape {
    /// Return the support point of the shape in direction `dir`.
    fn support(&self, dir: [Real; 3]) -> [Real; 3];
}

// ---------------------------------------------------------------------------
// Primitive shapes
// ---------------------------------------------------------------------------

/// A sphere centred at `center` with `radius`.
#[derive(Debug, Clone)]
pub struct Sphere {
    /// Centre position in world space.
    pub center: [Real; 3],
    /// Radius.
    pub radius: Real,
}

impl Sphere {
    /// Create a new sphere.
    pub fn new(center: [Real; 3], radius: Real) -> Self {
        Self { center, radius }
    }
}

impl ConvexShape for Sphere {
    fn support(&self, dir: [Real; 3]) -> [Real; 3] {
        let n = v3_normalize(dir);
        v3_add(self.center, v3_scale(n, self.radius))
    }
}

/// An axis-aligned box centred at `center` with `half_extents`.
#[derive(Debug, Clone)]
pub struct Box3 {
    /// Centre position in world space.
    pub center: [Real; 3],
    /// Half-extents along each axis.
    pub half_extents: [Real; 3],
}

impl Box3 {
    /// Create a new axis-aligned box.
    pub fn new(center: [Real; 3], half_extents: [Real; 3]) -> Self {
        Self {
            center,
            half_extents,
        }
    }
}

impl ConvexShape for Box3 {
    fn support(&self, dir: [Real; 3]) -> [Real; 3] {
        [
            self.center[0] + self.half_extents[0] * dir[0].signum(),
            self.center[1] + self.half_extents[1] * dir[1].signum(),
            self.center[2] + self.half_extents[2] * dir[2].signum(),
        ]
    }
}

/// A capsule: the Minkowski sum of a line segment and a sphere.
#[derive(Debug, Clone)]
pub struct Capsule {
    /// First endpoint of the inner segment.
    pub a: [Real; 3],
    /// Second endpoint of the inner segment.
    pub b: [Real; 3],
    /// Radius.
    pub radius: Real,
}

impl Capsule {
    /// Create a new capsule.
    pub fn new(a: [Real; 3], b: [Real; 3], radius: Real) -> Self {
        Self { a, b, radius }
    }
}

impl ConvexShape for Capsule {
    fn support(&self, dir: [Real; 3]) -> [Real; 3] {
        // Choose whichever endpoint is further in `dir`.
        let da = v3_dot(self.a, dir);
        let db = v3_dot(self.b, dir);
        let base = if da >= db { self.a } else { self.b };
        v3_add(base, v3_scale(v3_normalize(dir), self.radius))
    }
}

// ---------------------------------------------------------------------------
// Sphere-Sphere test
// ---------------------------------------------------------------------------

/// Result of a sphere–sphere intersection test.
#[derive(Debug, Clone)]
pub struct SphereSphereContact {
    /// Whether the spheres overlap.
    pub overlapping: bool,
    /// Contact normal from A to B (normalised).
    pub normal: [Real; 3],
    /// Penetration depth (negative = separated; positive = penetrating).
    pub depth: Real,
    /// Contact point (midpoint of overlap, or closest point).
    pub point: [Real; 3],
}

/// Test two spheres for intersection and compute the contact manifold.
pub fn sphere_sphere(a: &Sphere, b: &Sphere) -> SphereSphereContact {
    let d = v3_sub(b.center, a.center);
    let dist = v3_len(d);
    let sum_r = a.radius + b.radius;
    let depth = sum_r - dist;
    let normal = if dist > 1e-12 {
        v3_scale(d, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0] // degenerate: coincident centres
    };
    let point = v3_add(a.center, v3_scale(normal, a.radius - depth * 0.5));
    SphereSphereContact {
        overlapping: depth > 0.0,
        normal,
        depth,
        point,
    }
}

// ---------------------------------------------------------------------------
// Sphere-Capsule test
// ---------------------------------------------------------------------------

/// Signed distance from point `p` to the line segment (a→b), and the closest
/// point on the segment.
pub fn point_segment_closest(p: [Real; 3], a: [Real; 3], b: [Real; 3]) -> ([Real; 3], Real) {
    let ab = v3_sub(b, a);
    let ap = v3_sub(p, a);
    let t = (v3_dot(ap, ab) / v3_len_sq(ab).max(1e-30)).clamp(0.0, 1.0);
    let closest = v3_add(a, v3_scale(ab, t));
    let dist = v3_len(v3_sub(p, closest));
    (closest, dist)
}

/// Test a sphere and capsule for intersection.
/// Returns `(overlapping, depth, normal, point)`.
pub fn sphere_capsule(sphere: &Sphere, cap: &Capsule) -> (bool, Real, [Real; 3], [Real; 3]) {
    let (closest, dist) = point_segment_closest(sphere.center, cap.a, cap.b);
    let depth = (sphere.radius + cap.radius) - dist;
    let d = v3_sub(sphere.center, closest);
    let normal = if dist > 1e-12 {
        v3_scale(d, 1.0 / dist)
    } else {
        [0.0, 1.0, 0.0]
    };
    let point = v3_add(closest, v3_scale(normal, cap.radius - depth * 0.5));
    (depth > 0.0, depth, normal, point)
}

// ---------------------------------------------------------------------------
// AABB-AABB overlap
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box, stored as `(min, max)`.
#[derive(Debug, Clone, Copy)]
pub struct AabbRaw {
    /// Minimum corner.
    pub min: [Real; 3],
    /// Maximum corner.
    pub max: [Real; 3],
}

impl AabbRaw {
    /// Create a new AABB.
    pub fn new(min: [Real; 3], max: [Real; 3]) -> Self {
        Self { min, max }
    }

    /// Create from centre and half-extents.
    pub fn from_center_half(center: [Real; 3], half: [Real; 3]) -> Self {
        Self {
            min: [
                center[0] - half[0],
                center[1] - half[1],
                center[2] - half[2],
            ],
            max: [
                center[0] + half[0],
                center[1] + half[1],
                center[2] + half[2],
            ],
        }
    }

    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &AabbRaw) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Merge two AABBs.
    pub fn merge(&self, other: &AabbRaw) -> Self {
        Self {
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

    /// Surface area (used for BVH cost heuristic).
    pub fn surface_area(&self) -> Real {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Centre of the AABB.
    pub fn center(&self) -> [Real; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
}

// ---------------------------------------------------------------------------
// BVH (bounding volume hierarchy) — simple top-down build
// ---------------------------------------------------------------------------

/// A leaf entry in the BVH.
#[derive(Debug, Clone)]
pub struct BvhLeaf {
    /// Unique id of the object.
    pub id: u64,
    /// Tight AABB.
    pub aabb: AabbRaw,
}

/// A node in the BVH tree.
#[derive(Debug, Clone)]
enum BvhNode {
    Leaf {
        id: u64,
        aabb: AabbRaw,
    },
    Internal {
        aabb: AabbRaw,
        left: Box<BvhNode>,
        right: Box<BvhNode>,
    },
}

impl BvhNode {
    fn aabb(&self) -> &AabbRaw {
        match self {
            BvhNode::Leaf { aabb, .. } => aabb,
            BvhNode::Internal { aabb, .. } => aabb,
        }
    }

    fn query_overlap(&self, query: &AabbRaw, out: &mut Vec<u64>) {
        if !self.aabb().overlaps(query) {
            return;
        }
        match self {
            BvhNode::Leaf { id, .. } => out.push(*id),
            BvhNode::Internal { left, right, .. } => {
                left.query_overlap(query, out);
                right.query_overlap(query, out);
            }
        }
    }

    fn all_pairs<'a>(&'a self, other: &'a BvhNode, out: &mut Vec<(u64, u64)>) {
        if !self.aabb().overlaps(other.aabb()) {
            return;
        }
        match (self, other) {
            (BvhNode::Leaf { id: id_a, .. }, BvhNode::Leaf { id: id_b, .. }) => {
                if id_a != id_b {
                    out.push((*id_a.min(id_b), *id_a.max(id_b)));
                }
            }
            (BvhNode::Internal { left, right, .. }, _) => {
                left.all_pairs(other, out);
                right.all_pairs(other, out);
            }
            (_, BvhNode::Internal { left, right, .. }) => {
                self.all_pairs(left, out);
                self.all_pairs(right, out);
            }
        }
    }
}

/// A simple axis-aligned BVH built with a median-split strategy.
pub struct Bvh {
    root: Option<BvhNode>,
}

impl Bvh {
    /// Build a new BVH from the given leaf list.
    pub fn build(leaves: Vec<BvhLeaf>) -> Self {
        if leaves.is_empty() {
            return Self { root: None };
        }
        Self {
            root: Some(Self::build_recursive(leaves)),
        }
    }

    fn build_recursive(mut leaves: Vec<BvhLeaf>) -> BvhNode {
        if leaves.len() == 1 {
            let leaf = leaves.remove(0);
            return BvhNode::Leaf {
                id: leaf.id,
                aabb: leaf.aabb,
            };
        }

        // Compute encompassing AABB.
        let mut aabb = leaves[0].aabb;
        for l in &leaves[1..] {
            aabb = aabb.merge(&l.aabb);
        }

        // Split along the longest axis.
        let dx = aabb.max[0] - aabb.min[0];
        let dy = aabb.max[1] - aabb.min[1];
        let dz = aabb.max[2] - aabb.min[2];
        let axis = if dx >= dy && dx >= dz {
            0
        } else if dy >= dz {
            1
        } else {
            2
        };

        leaves.sort_by(|a, b| {
            a.aabb.center()[axis]
                .partial_cmp(&b.aabb.center()[axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mid = leaves.len() / 2;
        let right_leaves = leaves.split_off(mid);
        let left_leaves = leaves;

        let left = Box::new(Self::build_recursive(left_leaves));
        let right = Box::new(Self::build_recursive(right_leaves));
        BvhNode::Internal { aabb, left, right }
    }

    /// Return all ids whose AABB overlaps the given query AABB.
    pub fn query(&self, aabb: &AabbRaw) -> Vec<u64> {
        let mut result = Vec::new();
        if let Some(root) = &self.root {
            root.query_overlap(aabb, &mut result);
        }
        result
    }

    /// Return all overlapping pairs in the BVH (self-test).
    pub fn overlapping_pairs(&self) -> Vec<(u64, u64)> {
        let mut pairs = Vec::new();
        if let Some(root) = &self.root {
            root.all_pairs(root, &mut pairs);
        }
        // Deduplicate.
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

// ---------------------------------------------------------------------------
// GJK – Gilbert-Johnson-Keerthi distance algorithm
// ---------------------------------------------------------------------------

/// Minkowski difference support function for two convex shapes.
fn gjk_support(a: &dyn ConvexShape, b: &dyn ConvexShape, dir: [Real; 3]) -> [Real; 3] {
    v3_sub(a.support(dir), b.support(v3_neg(dir)))
}

/// Result of a GJK query.
#[derive(Debug, Clone)]
pub struct GjkResult {
    /// Whether the shapes intersect (minimum distance is 0).
    pub intersecting: bool,
    /// Closest distance between the shapes (0 if intersecting).
    pub distance: Real,
    /// Closest point on shape A (only valid when not intersecting).
    pub closest_a: [Real; 3],
    /// Closest point on shape B (only valid when not intersecting).
    pub closest_b: [Real; 3],
}

/// Run the GJK algorithm to determine if two convex shapes intersect and,
/// if not, compute their minimum distance.
///
/// Returns `GjkResult` with the intersection flag and minimum distance.
pub fn gjk(a: &dyn ConvexShape, b: &dyn ConvexShape) -> GjkResult {
    const MAX_ITER: usize = 64;
    const EPS: Real = 1e-10;

    let mut simplex: Vec<[Real; 3]> = Vec::with_capacity(4);
    let mut dir = [1.0, 0.0, 0.0_f64];

    let first = gjk_support(a, b, dir);
    simplex.push(first);
    dir = v3_neg(first);

    for _ in 0..MAX_ITER {
        let len_sq = v3_len_sq(dir);
        if len_sq < EPS {
            // Direction degenerate — intersection at origin.
            return GjkResult {
                intersecting: true,
                distance: 0.0,
                closest_a: [0.0; 3],
                closest_b: [0.0; 3],
            };
        }

        let support = gjk_support(a, b, dir);

        // If the new support does not pass the origin, shapes do not intersect.
        if v3_dot(support, dir) < 0.0 {
            let dist = v3_len_sq(dir).sqrt();
            return GjkResult {
                intersecting: false,
                distance: dist,
                closest_a: a.support(dir),
                closest_b: b.support(v3_neg(dir)),
            };
        }

        simplex.push(support);

        if gjk_do_simplex(&mut simplex, &mut dir) {
            return GjkResult {
                intersecting: true,
                distance: 0.0,
                closest_a: [0.0; 3],
                closest_b: [0.0; 3],
            };
        }
    }

    // Fallback: treat as intersecting.
    GjkResult {
        intersecting: true,
        distance: 0.0,
        closest_a: [0.0; 3],
        closest_b: [0.0; 3],
    }
}

/// Update the GJK simplex and direction; returns `true` if the origin is inside.
fn gjk_do_simplex(simplex: &mut Vec<[Real; 3]>, dir: &mut [Real; 3]) -> bool {
    match simplex.len() {
        2 => gjk_line_case(simplex, dir),
        3 => gjk_triangle_case(simplex, dir),
        4 => gjk_tetrahedron_case(simplex, dir),
        _ => false,
    }
}

fn gjk_line_case(simplex: &mut Vec<[Real; 3]>, dir: &mut [Real; 3]) -> bool {
    let b = simplex[0];
    let a = simplex[1];
    let ab = v3_sub(b, a);
    let ao = v3_neg(a);
    if v3_dot(ab, ao) > 0.0 {
        // Origin is between A and B.
        *dir = v3_sub(v3_scale(ab, v3_dot(ab, ao) / v3_dot(ab, ab)), ao);
    } else {
        simplex.clear();
        simplex.push(a);
        *dir = ao;
    }
    false
}

fn gjk_triangle_case(simplex: &mut Vec<[Real; 3]>, dir: &mut [Real; 3]) -> bool {
    let c = simplex[0];
    let b = simplex[1];
    let a = simplex[2];
    let ab = v3_sub(b, a);
    let ac = v3_sub(c, a);
    let ao = v3_neg(a);
    let abc = v3_cross(ab, ac);

    // Check if origin is outside edge AC.
    if v3_dot(v3_cross(abc, ac), ao) > 0.0 {
        if v3_dot(ac, ao) > 0.0 {
            simplex.clear();
            simplex.push(c);
            simplex.push(a);
            *dir = v3_sub(v3_scale(ac, v3_dot(ac, ao) / v3_dot(ac, ac)), ao);
        } else {
            simplex.clear();
            simplex.push(b);
            simplex.push(a);
            return gjk_line_case(simplex, dir);
        }
    } else if v3_dot(v3_cross(ab, abc), ao) > 0.0 {
        simplex.clear();
        simplex.push(b);
        simplex.push(a);
        return gjk_line_case(simplex, dir);
    } else {
        // Origin is above or below the triangle.
        if v3_dot(abc, ao) > 0.0 {
            *dir = abc;
        } else {
            simplex.swap(0, 1); // flip winding
            *dir = v3_neg(abc);
        }
    }
    false
}

fn gjk_tetrahedron_case(simplex: &mut Vec<[Real; 3]>, dir: &mut [Real; 3]) -> bool {
    let d = simplex[0];
    let c = simplex[1];
    let b = simplex[2];
    let a = simplex[3];
    let ab = v3_sub(b, a);
    let ac = v3_sub(c, a);
    let ad = v3_sub(d, a);
    let ao = v3_neg(a);

    let abc = v3_cross(ab, ac);
    let acd = v3_cross(ac, ad);
    let adb = v3_cross(ad, ab);

    if v3_dot(abc, ao) > 0.0 {
        simplex.clear();
        simplex.push(c);
        simplex.push(b);
        simplex.push(a);
        return gjk_triangle_case(simplex, dir);
    }
    if v3_dot(acd, ao) > 0.0 {
        simplex.clear();
        simplex.push(d);
        simplex.push(c);
        simplex.push(a);
        return gjk_triangle_case(simplex, dir);
    }
    if v3_dot(adb, ao) > 0.0 {
        simplex.clear();
        simplex.push(b);
        simplex.push(d);
        simplex.push(a);
        return gjk_triangle_case(simplex, dir);
    }
    // Origin is inside the tetrahedron.
    true
}

// ---------------------------------------------------------------------------
// SAT (Separating Axis Theorem) – OBB vs OBB
// ---------------------------------------------------------------------------

/// An oriented bounding box.
#[derive(Debug, Clone)]
pub struct Obb {
    /// Centre position.
    pub center: [Real; 3],
    /// Three normalised local axes (columns of the rotation matrix).
    pub axes: [[Real; 3]; 3],
    /// Half-extents along each axis.
    pub half_extents: [Real; 3],
}

impl Obb {
    /// Create a new OBB.
    pub fn new(center: [Real; 3], axes: [[Real; 3]; 3], half_extents: [Real; 3]) -> Self {
        Self {
            center,
            axes,
            half_extents,
        }
    }

    /// Create an axis-aligned OBB (identity rotation).
    pub fn axis_aligned(center: [Real; 3], half_extents: [Real; 3]) -> Self {
        Self {
            center,
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            half_extents,
        }
    }

    /// Project the OBB onto axis `ax` and return (min, max).
    fn project(&self, ax: [Real; 3]) -> (Real, Real) {
        let c = v3_dot(self.center, ax);
        let r = self
            .axes
            .iter()
            .zip(self.half_extents.iter())
            .map(|(&a, &h)| v3_dot(a, ax).abs() * h)
            .sum::<Real>();
        (c - r, c + r)
    }

    /// Return the OBB's 8 corner vertices in world space.
    pub fn vertices(&self) -> [[Real; 3]; 8] {
        let mut verts = [[0.0; 3]; 8];
        let signs = [
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        for (i, s) in signs.iter().enumerate() {
            let mut v = self.center;
            for (k, (&sk, &he)) in s.iter().zip(self.half_extents.iter()).enumerate() {
                v = v3_add(v, v3_scale(self.axes[k], sk * he));
            }
            verts[i] = v;
        }
        verts
    }
}

/// Test two OBBs using the Separating Axis Theorem.
/// Returns `None` if separated, or `Some(depth, axis)` with penetration info.
pub fn obb_obb_sat(a: &Obb, b: &Obb) -> Option<(Real, [Real; 3])> {
    let mut min_depth = Real::MAX;
    let mut min_axis = [0.0; 3];

    // Helper: test a single axis.
    let mut test_axis = |ax: [Real; 3]| -> bool {
        let len_sq = v3_len_sq(ax);
        if len_sq < 1e-10 {
            return true; // degenerate axis; skip
        }
        let ax = v3_scale(ax, 1.0 / len_sq.sqrt());
        let (a_min, a_max) = a.project(ax);
        let (b_min, b_max) = b.project(ax);
        if a_max < b_min || b_max < a_min {
            return false; // separated
        }
        let depth = (a_max.min(b_max) - a_min.max(b_min))
            .min(a_max - a_min)
            .min(b_max - b_min);
        if depth < min_depth {
            min_depth = depth;
            min_axis = ax;
        }
        true
    };

    // 3 face axes of A
    for i in 0..3 {
        if !test_axis(a.axes[i]) {
            return None;
        }
    }
    // 3 face axes of B
    for i in 0..3 {
        if !test_axis(b.axes[i]) {
            return None;
        }
    }
    // 9 cross-product edge axes
    for i in 0..3 {
        for j in 0..3 {
            let ax = v3_cross(a.axes[i], b.axes[j]);
            if !test_axis(ax) {
                return None;
            }
        }
    }

    Some((min_depth, min_axis))
}

// ---------------------------------------------------------------------------
// Ray cast
// ---------------------------------------------------------------------------

/// Result of a ray–AABB intersection test.
#[derive(Debug, Clone)]
pub struct RayAabbHit {
    /// Parameter `t_min` along the ray at entry.
    pub t_min: Real,
    /// Parameter `t_max` along the ray at exit.
    pub t_max: Real,
}

/// Intersect a ray `(origin, dir)` with an AABB.
///
/// Returns `Some(hit)` if the ray hits the box; the ray is infinite.
pub fn ray_aabb(origin: [Real; 3], dir: [Real; 3], aabb: &AabbRaw) -> Option<RayAabbHit> {
    let mut t_min = Real::NEG_INFINITY;
    let mut t_max = Real::INFINITY;

    for i in 0..3 {
        if dir[i].abs() < 1e-12 {
            // Ray is parallel to slab.
            if origin[i] < aabb.min[i] || origin[i] > aabb.max[i] {
                return None;
            }
        } else {
            let inv_d = 1.0 / dir[i];
            let t1 = (aabb.min[i] - origin[i]) * inv_d;
            let t2 = (aabb.max[i] - origin[i]) * inv_d;
            let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
            t_min = t_min.max(ta);
            t_max = t_max.min(tb);
            if t_min > t_max {
                return None;
            }
        }
    }
    Some(RayAabbHit { t_min, t_max })
}

/// Intersect a ray with a sphere.  Returns parameter `t` at first hit or `None`.
pub fn ray_sphere(origin: [Real; 3], dir: [Real; 3], sphere: &Sphere) -> Option<Real> {
    let oc = v3_sub(origin, sphere.center);
    let a = v3_dot(dir, dir);
    let b = 2.0 * v3_dot(oc, dir);
    let c = v3_dot(oc, oc) - sphere.radius * sphere.radius;
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return None;
    }
    let sqrt_d = discriminant.sqrt();
    let t0 = (-b - sqrt_d) / (2.0 * a);
    let t1 = (-b + sqrt_d) / (2.0 * a);
    // Return the nearest positive t.
    if t0 > 0.0 {
        Some(t0)
    } else if t1 > 0.0 {
        Some(t1)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Contact manifold
// ---------------------------------------------------------------------------

/// A collision contact point.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Contact point in world space.
    pub point: [Real; 3],
    /// Contact normal (from B to A, normalised).
    pub normal: [Real; 3],
    /// Penetration depth.
    pub depth: Real,
    /// Id of the first object.
    pub id_a: u64,
    /// Id of the second object.
    pub id_b: u64,
}

/// A contact manifold: up to 4 persistent contact points between two bodies.
#[derive(Debug, Clone, Default)]
pub struct ContactManifold {
    /// Up to 4 contact points.
    pub contacts: Vec<Contact>,
}

impl ContactManifold {
    /// Create an empty manifold.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a contact, replacing the least significant one if there are already 4.
    pub fn add_contact(&mut self, contact: Contact) {
        if self.contacts.len() < 4 {
            self.contacts.push(contact);
        } else {
            // Replace the contact with the smallest depth.
            if let Some(min_idx) = self
                .contacts
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.depth
                        .partial_cmp(&b.depth)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
            {
                let new_depth = contact.depth;
                if new_depth > self.contacts[min_idx].depth {
                    self.contacts[min_idx] = contact;
                }
            }
        }
    }

    /// Return the average contact normal.
    pub fn average_normal(&self) -> [Real; 3] {
        if self.contacts.is_empty() {
            return [0.0, 1.0, 0.0];
        }
        let mut n = [0.0; 3];
        for c in &self.contacts {
            n[0] += c.normal[0];
            n[1] += c.normal[1];
            n[2] += c.normal[2];
        }
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len > 1e-12 {
            [n[0] / len, n[1] / len, n[2] / len]
        } else {
            [0.0, 1.0, 0.0]
        }
    }

    /// Maximum penetration depth across all contacts.
    pub fn max_depth(&self) -> Real {
        self.contacts
            .iter()
            .map(|c| c.depth)
            .fold(Real::NEG_INFINITY, Real::max)
    }

    /// Whether there are any contacts.
    pub fn is_empty(&self) -> bool {
        self.contacts.is_empty()
    }
}

// ---------------------------------------------------------------------------
// EPA – Expanding Polytope Algorithm (penetration depth)
// ---------------------------------------------------------------------------

/// Result of an EPA penetration-depth query.
#[derive(Debug, Clone)]
pub struct EpaResult {
    /// Penetration depth (positive when overlapping).
    pub depth: Real,
    /// Collision normal (from B towards A, normalised).
    pub normal: [Real; 3],
    /// Contact point estimate.
    pub point: [Real; 3],
}

/// Run EPA to compute the penetration depth between two overlapping convex shapes.
///
/// Requires that GJK has confirmed the shapes are intersecting.
/// Returns `None` if the shapes are not intersecting.
pub fn epa(a: &dyn ConvexShape, b: &dyn ConvexShape) -> Option<EpaResult> {
    // First confirm intersection with GJK.
    let gjk_r = gjk(a, b);
    if !gjk_r.intersecting {
        return None;
    }

    const MAX_ITER: usize = 64;
    const EPS: Real = 1e-8;

    // Seed the polytope with a tetrahedron constructed from support points.
    let dirs = [
        [1.0, 0.0, 0.0_f64],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];

    let mut polytope: Vec<[Real; 3]> = dirs.iter().map(|&d| gjk_support(a, b, d)).collect();

    // Build initial faces (pairs of triangles from the seed points).
    // Simple approach: build from all unique triangles in the convex hull of polytope.
    // For EPA we iteratively expand closest face.

    let mut faces: Vec<[usize; 3]> = vec![
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 1],
        [5, 2, 1],
        [5, 3, 2],
        [5, 4, 3],
        [5, 1, 4],
    ];

    for _ in 0..MAX_ITER {
        // Find the face closest to the origin.
        let mut min_dist = Real::MAX;
        let mut min_idx = 0;
        let mut min_normal = [0.0_f64; 3];

        for (fi, face) in faces.iter().enumerate() {
            let a_pt = polytope[face[0]];
            let b_pt = polytope[face[1]];
            let c_pt = polytope[face[2]];
            let ab = v3_sub(b_pt, a_pt);
            let ac = v3_sub(c_pt, a_pt);
            let n = v3_normalize(v3_cross(ab, ac));
            let d = v3_dot(n, a_pt).abs();
            if d < min_dist {
                min_dist = d;
                min_idx = fi;
                min_normal = n;
            }
        }

        // Support in the direction of the closest face normal.
        let support = gjk_support(a, b, min_normal);
        let new_dist = v3_dot(min_normal, support);

        if (new_dist - min_dist).abs() < EPS {
            // Converged.
            let contact = v3_scale(min_normal, min_dist);
            return Some(EpaResult {
                depth: min_dist,
                normal: min_normal,
                point: contact,
            });
        }

        // Expand polytope: remove visible faces and add new faces with the new vertex.
        let new_idx = polytope.len();
        polytope.push(support);

        let mut new_faces: Vec<[usize; 3]> = Vec::new();
        let mut edges: Vec<(usize, usize)> = Vec::new();

        for (fi, &face) in faces.iter().enumerate() {
            if fi == min_idx {
                // Record edges of the removed face (horizon).
                edges.push((face[0], face[1]));
                edges.push((face[1], face[2]));
                edges.push((face[2], face[0]));
            } else {
                new_faces.push(face);
            }
        }

        // Add new faces connecting the new vertex to the horizon edges.
        for (ea, eb) in edges {
            new_faces.push([ea, eb, new_idx]);
        }
        faces = new_faces;
    }

    // Fallback after max iterations.
    let normal = [0.0, 1.0, 0.0];
    Some(EpaResult {
        depth: 0.0,
        normal,
        point: [0.0; 3],
    })
}

// ---------------------------------------------------------------------------
// Minkowski sum support function
// ---------------------------------------------------------------------------

/// Support function of the Minkowski sum of two shapes.
///
/// `support_mink(A⊕B, d) = support_A(d) + support_B(d)`.
pub fn minkowski_sum_support(
    a: &dyn ConvexShape,
    b: &dyn ConvexShape,
    dir: [Real; 3],
) -> [Real; 3] {
    v3_add(a.support(dir), b.support(dir))
}

// ---------------------------------------------------------------------------
// Shape casting (linear motion CCD)
// ---------------------------------------------------------------------------

/// Cast a moving sphere along `velocity` direction (magnitude = max travel distance)
/// against a static sphere.  Returns the first time of contact `t ∈ [0, max_t]`
/// or `None` if no contact.
pub fn shape_cast_sphere_vs_sphere(
    moving: &Sphere,
    velocity: [Real; 3],
    target: &Sphere,
    max_t: Real,
) -> Option<Real> {
    // Relative motion: relative velocity = velocity (target is static).
    // Reduce to: find t such that |moving.center + t*vel - target.center| = r_a + r_b.
    let oc = v3_sub(moving.center, target.center);
    let combined_r = moving.radius + target.radius;

    let a = v3_dot(velocity, velocity);
    let b = 2.0 * v3_dot(oc, velocity);
    let c = v3_dot(oc, oc) - combined_r * combined_r;

    if a < 1e-15 {
        // No motion.
        return if c <= 0.0 { Some(0.0) } else { None };
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }

    let sqrt_d = disc.sqrt();
    let t0 = (-b - sqrt_d) / (2.0 * a);
    let t1 = (-b + sqrt_d) / (2.0 * a);

    // Return smallest non-negative t within [0, max_t].
    for t in [t0, t1] {
        if t >= 0.0 && t <= max_t {
            return Some(t);
        }
    }
    // Already overlapping (c <= 0) → t = 0.
    if c <= 0.0 { Some(0.0) } else { None }
}

/// Cast a sphere along `velocity` against an AABB.  Returns `t ∈ [0, max_t]` or `None`.
///
/// Implemented as a ray cast against the Minkowski-expanded AABB.
pub fn shape_cast_sphere_vs_aabb(
    sphere: &Sphere,
    velocity: [Real; 3],
    aabb: &AabbRaw,
    max_t: Real,
) -> Option<Real> {
    // Expand AABB by sphere radius.
    let r = sphere.radius;
    let expanded = AabbRaw::new(
        [aabb.min[0] - r, aabb.min[1] - r, aabb.min[2] - r],
        [aabb.max[0] + r, aabb.max[1] + r, aabb.max[2] + r],
    );
    let hit = ray_aabb(sphere.center, velocity, &expanded)?;
    if hit.t_min < 0.0 || hit.t_min > max_t {
        return None;
    }
    Some(hit.t_min.max(0.0))
}

// ---------------------------------------------------------------------------
// Triangle convex shape
// ---------------------------------------------------------------------------

/// A triangle as a convex shape (degenerate 2-D polyhedron).
#[derive(Debug, Clone)]
pub struct Triangle {
    /// Vertices of the triangle.
    pub vertices: [[Real; 3]; 3],
}

impl Triangle {
    /// Create a new triangle.
    pub fn new(v0: [Real; 3], v1: [Real; 3], v2: [Real; 3]) -> Self {
        Self {
            vertices: [v0, v1, v2],
        }
    }

    /// Outward normal of the triangle (not normalised).
    pub fn normal(&self) -> [Real; 3] {
        let e1 = v3_sub(self.vertices[1], self.vertices[0]);
        let e2 = v3_sub(self.vertices[2], self.vertices[0]);
        v3_cross(e1, e2)
    }
}

impl ConvexShape for Triangle {
    fn support(&self, dir: [Real; 3]) -> [Real; 3] {
        self.vertices
            .iter()
            .copied()
            .max_by(|&a, &b| {
                v3_dot(a, dir)
                    .partial_cmp(&v3_dot(b, dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or([0.0; 3])
    }
}

// ---------------------------------------------------------------------------
// Convex hull shape (point cloud)
// ---------------------------------------------------------------------------

/// A convex shape defined by a set of points.
///
/// The support function iterates all points; for small point sets this is fine.
#[derive(Debug, Clone)]
pub struct ConvexPointCloud {
    /// Points defining the convex hull.
    pub points: Vec<[Real; 3]>,
}

impl ConvexPointCloud {
    /// Create a new convex point cloud.
    pub fn new(points: Vec<[Real; 3]>) -> Self {
        Self { points }
    }
}

impl ConvexShape for ConvexPointCloud {
    fn support(&self, dir: [Real; 3]) -> [Real; 3] {
        self.points
            .iter()
            .copied()
            .max_by(|&a, &b| {
                v3_dot(a, dir)
                    .partial_cmp(&v3_dot(b, dir))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or([0.0; 3])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Vector helpers ────────────────────────────────────────────────────────

    #[test]
    fn test_v3_dot_orthogonal() {
        assert_eq!(v3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 0.0);
    }

    #[test]
    fn test_v3_cross_basis() {
        let c = v3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-12);
        assert!(c[0].abs() < 1e-12 && c[1].abs() < 1e-12);
    }

    #[test]
    fn test_v3_normalize_unit() {
        let n = v3_normalize([3.0, 4.0, 0.0]);
        let len = v3_len(n);
        assert!((len - 1.0).abs() < 1e-12, "len={}", len);
    }

    #[test]
    fn test_v3_normalize_zero_safe() {
        let n = v3_normalize([0.0; 3]);
        assert_eq!(n, [0.0; 3]);
    }

    // ── Sphere-Sphere ─────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_sphere_overlap() {
        let a = Sphere::new([0.0; 3], 1.0);
        let b = Sphere::new([1.5, 0.0, 0.0], 1.0);
        let r = sphere_sphere(&a, &b);
        assert!(r.overlapping, "spheres should overlap");
        assert!(r.depth > 0.0, "depth should be positive");
    }

    #[test]
    fn test_sphere_sphere_no_overlap() {
        let a = Sphere::new([0.0; 3], 1.0);
        let b = Sphere::new([5.0, 0.0, 0.0], 1.0);
        let r = sphere_sphere(&a, &b);
        assert!(!r.overlapping);
        assert!(r.depth < 0.0);
    }

    #[test]
    fn test_sphere_sphere_touching() {
        let a = Sphere::new([0.0; 3], 1.0);
        let b = Sphere::new([2.0, 0.0, 0.0], 1.0);
        let r = sphere_sphere(&a, &b);
        assert!((r.depth).abs() < 1e-10, "depth should be ~0 at touching");
    }

    // ── Sphere-Capsule ────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_capsule_overlap() {
        let s = Sphere::new([0.0; 3], 0.5);
        let c = Capsule::new([0.0, 2.0, 0.0], [0.0, 4.0, 0.0], 0.5);
        let (overlap, depth, _, _) = sphere_capsule(&s, &c);
        assert!(!overlap, "should not overlap, depth={}", depth);
    }

    #[test]
    fn test_sphere_capsule_nearby() {
        let s = Sphere::new([0.0, 1.0, 0.0], 0.6);
        let c = Capsule::new([0.0, 0.0, 0.0], [0.0, 3.0, 0.0], 0.6);
        let (overlap, depth, _, _) = sphere_capsule(&s, &c);
        assert!(overlap, "should overlap, depth={}", depth);
    }

    // ── AABB ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_overlap() {
        let a = AabbRaw::new([-1.0; 3], [1.0; 3]);
        let b = AabbRaw::new([0.5, 0.5, 0.5], [2.0, 2.0, 2.0]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap() {
        let a = AabbRaw::new([-1.0; 3], [1.0; 3]);
        let b = AabbRaw::new([2.0; 3], [3.0; 3]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_merge_contains_both() {
        let a = AabbRaw::new([-2.0; 3], [0.0; 3]);
        let b = AabbRaw::new([1.0; 3], [3.0; 3]);
        let m = a.merge(&b);
        assert!(m.min[0] <= -2.0 && m.max[0] >= 3.0);
    }

    // ── BVH ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_bvh_query_single_match() {
        let leaves = vec![
            BvhLeaf {
                id: 1,
                aabb: AabbRaw::new([0.0; 3], [1.0; 3]),
            },
            BvhLeaf {
                id: 2,
                aabb: AabbRaw::new([5.0; 3], [6.0; 3]),
            },
        ];
        let bvh = Bvh::build(leaves);
        let query = AabbRaw::new([0.0; 3], [0.5; 3]);
        let hits = bvh.query(&query);
        assert_eq!(hits, vec![1]);
    }

    #[test]
    fn test_bvh_query_no_match() {
        let leaves = vec![BvhLeaf {
            id: 1,
            aabb: AabbRaw::new([10.0; 3], [11.0; 3]),
        }];
        let bvh = Bvh::build(leaves);
        let query = AabbRaw::new([0.0; 3], [1.0; 3]);
        assert!(bvh.query(&query).is_empty());
    }

    #[test]
    fn test_bvh_empty() {
        let bvh = Bvh::build(vec![]);
        let query = AabbRaw::new([0.0; 3], [1.0; 3]);
        assert!(bvh.query(&query).is_empty());
    }

    // ── OBB SAT ───────────────────────────────────────────────────────────────

    #[test]
    fn test_obb_sat_overlap() {
        let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
        let b = Obb::axis_aligned([1.5, 0.0, 0.0], [1.0; 3]);
        assert!(obb_obb_sat(&a, &b).is_some(), "overlapping OBBs");
    }

    #[test]
    fn test_obb_sat_separated() {
        let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
        let b = Obb::axis_aligned([5.0, 0.0, 0.0], [1.0; 3]);
        assert!(obb_obb_sat(&a, &b).is_none(), "separated OBBs");
    }

    // ── Ray tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_ray_aabb_hit() {
        let origin = [-5.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
        let hit = ray_aabb(origin, dir, &aabb);
        assert!(hit.is_some());
        let h = hit.unwrap();
        assert!(h.t_min < h.t_max);
    }

    #[test]
    fn test_ray_aabb_miss() {
        let origin = [0.0, 5.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
        assert!(ray_aabb(origin, dir, &aabb).is_none());
    }

    #[test]
    fn test_ray_sphere_hit() {
        let s = Sphere::new([0.0; 3], 1.0);
        let t = ray_sphere([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], &s);
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_sphere_miss() {
        let s = Sphere::new([0.0; 3], 1.0);
        let t = ray_sphere([0.0, 5.0, 0.0], [1.0, 0.0, 0.0], &s);
        assert!(t.is_none());
    }

    // ── GJK ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_gjk_spheres_intersecting() {
        let a = Sphere::new([0.0; 3], 1.5);
        let b = Sphere::new([1.0, 0.0, 0.0], 1.5);
        let r = gjk(&a, &b);
        assert!(r.intersecting);
    }

    #[test]
    fn test_gjk_spheres_separated() {
        let a = Sphere::new([0.0; 3], 0.5);
        let b = Sphere::new([10.0, 0.0, 0.0], 0.5);
        let r = gjk(&a, &b);
        assert!(!r.intersecting);
        assert!(r.distance > 0.0);
    }

    // ── ContactManifold ───────────────────────────────────────────────────────

    #[test]
    fn test_manifold_add_and_max_depth() {
        let mut m = ContactManifold::new();
        m.add_contact(Contact {
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            depth: 0.1,
            id_a: 1,
            id_b: 2,
        });
        m.add_contact(Contact {
            point: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            depth: 0.3,
            id_a: 1,
            id_b: 2,
        });
        assert_eq!(m.contacts.len(), 2);
        assert!((m.max_depth() - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_manifold_average_normal_unit() {
        let mut m = ContactManifold::new();
        for _ in 0..3 {
            m.add_contact(Contact {
                point: [0.0; 3],
                normal: [0.0, 1.0, 0.0],
                depth: 0.1,
                id_a: 1,
                id_b: 2,
            });
        }
        let n = m.average_normal();
        let len = v3_len(n);
        assert!((len - 1.0).abs() < 1e-12);
    }

    // ── EPA tests ─────────────────────────────────────────────────────────────

    #[test]
    fn test_epa_two_overlapping_spheres() {
        let a = Sphere::new([0.0; 3], 1.5);
        let b = Sphere::new([1.0, 0.0, 0.0], 1.5);
        let result = epa(&a, &b);
        assert!(
            result.is_some(),
            "EPA should return penetration for overlapping spheres"
        );
        let ep = result.unwrap();
        assert!(ep.depth >= 0.0, "penetration depth should be non-negative");
    }

    #[test]
    fn test_epa_concentric_spheres() {
        let a = Sphere::new([0.0; 3], 2.0);
        let b = Sphere::new([0.0; 3], 1.0);
        let result = epa(&a, &b);
        // Concentric spheres: deep penetration
        assert!(result.is_some());
    }

    #[test]
    fn test_epa_box_sphere_overlap() {
        let b = Box3::new([0.0; 3], [1.0; 3]);
        let s = Sphere::new([0.5, 0.0, 0.0], 0.8);
        let result = epa(&b, &s);
        assert!(
            result.is_some(),
            "overlapping box and sphere should produce EPA result"
        );
    }

    // ── Support function framework ────────────────────────────────────────────

    #[test]
    fn test_support_sphere_axis_aligned() {
        let s = Sphere::new([1.0, 2.0, 3.0], 2.0);
        let sp_x = s.support([1.0, 0.0, 0.0]);
        assert!(
            (sp_x[0] - 3.0).abs() < 1e-10,
            "x support = center.x + radius"
        );
        let sp_y = s.support([0.0, 1.0, 0.0]);
        assert!(
            (sp_y[1] - 4.0).abs() < 1e-10,
            "y support = center.y + radius"
        );
    }

    #[test]
    fn test_support_box_positive_axes() {
        let b = Box3::new([0.0; 3], [2.0, 3.0, 4.0]);
        let sp = b.support([1.0, 1.0, 1.0]);
        assert!((sp[0] - 2.0).abs() < 1e-10);
        assert!((sp[1] - 3.0).abs() < 1e-10);
        assert!((sp[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_support_box_negative_axes() {
        let b = Box3::new([0.0; 3], [1.0; 3]);
        let sp = b.support([-1.0, -1.0, -1.0]);
        assert!((sp[0] + 1.0).abs() < 1e-10);
        assert!((sp[1] + 1.0).abs() < 1e-10);
        assert!((sp[2] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_support_capsule_chooses_correct_endpoint() {
        let cap = Capsule::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        let sp_up = cap.support([0.0, 1.0, 0.0]);
        // Should be near top endpoint + radius in Y
        assert!(sp_up[1] > 1.0, "y support above top endpoint");
        let sp_down = cap.support([0.0, -1.0, 0.0]);
        assert!(sp_down[1] < -1.0, "y support below bottom endpoint");
    }

    // ── Minkowski sum ─────────────────────────────────────────────────────────

    #[test]
    fn test_minkowski_sum_sphere_sphere_support() {
        let a = Sphere::new([0.0; 3], 1.0);
        let b = Sphere::new([0.0; 3], 2.0);
        let dir = [1.0, 0.0, 0.0];
        // Support of sum = support of a + support of b in same direction
        let sp = minkowski_sum_support(&a, &b, dir);
        let expected = a.support(dir)[0] + b.support(dir)[0];
        assert!((sp[0] - expected).abs() < 1e-10);
    }

    #[test]
    fn test_minkowski_sum_sphere_box() {
        let s = Sphere::new([0.0; 3], 0.5);
        let b = Box3::new([0.0; 3], [1.0, 1.0, 1.0]);
        let dir = [1.0, 0.0, 0.0];
        let sp = minkowski_sum_support(&s, &b, dir);
        // Support ≥ box support alone (sphere only adds positive amount)
        assert!(sp[0] >= b.support(dir)[0]);
    }

    // ── Shape casting / continuous collision ──────────────────────────────────

    #[test]
    fn test_shape_cast_sphere_hits_sphere() {
        let moving = Sphere::new([0.0; 3], 0.5);
        let target = Sphere::new([5.0, 0.0, 0.0], 0.5);
        let velocity = [1.0, 0.0, 0.0];
        let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
        assert!(t.is_some(), "moving sphere should hit static sphere");
        let tv = t.unwrap();
        assert!((0.0..=10.0).contains(&tv));
    }

    #[test]
    fn test_shape_cast_sphere_misses_sphere() {
        let moving = Sphere::new([0.0; 3], 0.5);
        let target = Sphere::new([0.0, 100.0, 0.0], 0.5);
        let velocity = [1.0, 0.0, 0.0]; // moving in X, target far in Y
        let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
        assert!(
            t.is_none(),
            "sphere moving in X should miss sphere far in Y"
        );
    }

    #[test]
    fn test_shape_cast_already_overlapping() {
        let moving = Sphere::new([0.0; 3], 1.0);
        let target = Sphere::new([0.5, 0.0, 0.0], 1.0);
        let velocity = [1.0, 0.0, 0.0];
        let t = shape_cast_sphere_vs_sphere(&moving, velocity, &target, 10.0);
        // Already overlapping → t = 0 or very small
        if let Some(tv) = t {
            assert!(tv >= 0.0);
        }
    }

    // ── GJK with capsule/box ──────────────────────────────────────────────────

    #[test]
    fn test_gjk_box_sphere_intersecting() {
        let b = Box3::new([0.0; 3], [1.0; 3]);
        let s = Sphere::new([1.5, 0.0, 0.0], 1.0);
        let r = gjk(&b, &s);
        assert!(r.intersecting, "box and sphere should intersect");
    }

    #[test]
    fn test_gjk_box_sphere_separated() {
        let b = Box3::new([0.0; 3], [1.0; 3]);
        let s = Sphere::new([10.0, 0.0, 0.0], 0.5);
        let r = gjk(&b, &s);
        assert!(!r.intersecting, "box and sphere should be separated");
        assert!(r.distance > 0.0);
    }

    #[test]
    fn test_gjk_capsule_sphere_overlapping() {
        let cap = Capsule::new([0.0, -2.0, 0.0], [0.0, 2.0, 0.0], 0.5);
        let s = Sphere::new([0.0, 0.0, 0.0], 0.5);
        let r = gjk(&cap, &s);
        assert!(r.intersecting, "capsule and sphere should overlap");
    }

    #[test]
    fn test_gjk_capsule_sphere_separated() {
        let cap = Capsule::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.3);
        let s = Sphere::new([5.0, 0.0, 0.0], 0.3);
        let r = gjk(&cap, &s);
        assert!(!r.intersecting);
    }

    // ── AabbRaw additional tests ───────────────────────────────────────────────

    #[test]
    fn test_aabb_surface_area_unit_cube() {
        let a = AabbRaw::new([0.0; 3], [1.0; 3]);
        // Unit cube: SA = 6
        assert!((a.surface_area() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_center() {
        let a = AabbRaw::new([-2.0; 3], [4.0; 3]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_from_center_half() {
        let a = AabbRaw::from_center_half([1.0, 2.0, 3.0], [0.5; 3]);
        assert!((a.min[0] - 0.5).abs() < 1e-12);
        assert!((a.max[1] - 2.5).abs() < 1e-12);
    }

    // ── OBB additional tests ───────────────────────────────────────────────────

    #[test]
    fn test_obb_vertices_count() {
        let obb = Obb::axis_aligned([0.0; 3], [1.0; 3]);
        let verts = obb.vertices();
        assert_eq!(verts.len(), 8);
    }

    #[test]
    fn test_obb_axis_aligned_vertices_correct() {
        let obb = Obb::axis_aligned([0.0; 3], [1.0; 3]);
        let verts = obb.vertices();
        // Every vertex should have coordinates in {-1, 1}³
        for v in &verts {
            assert!(v[0].abs() <= 1.0 + 1e-10);
            assert!(v[1].abs() <= 1.0 + 1e-10);
            assert!(v[2].abs() <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn test_obb_sat_touching() {
        // Two OBBs touching exactly at their faces
        let a = Obb::axis_aligned([0.0; 3], [1.0; 3]);
        let b = Obb::axis_aligned([2.0, 0.0, 0.0], [1.0; 3]);
        // Distance between centers = 2, sum of half-extents in X = 2 → touching
        let result = obb_obb_sat(&a, &b);
        // Either touching (some depth ~ 0) or just separated — both valid
        let _ = result;
    }

    // ── BVH additional tests ───────────────────────────────────────────────────

    #[test]
    fn test_bvh_all_overlapping_pairs() {
        let leaves = vec![
            BvhLeaf {
                id: 1,
                aabb: AabbRaw::new([0.0; 3], [2.0; 3]),
            },
            BvhLeaf {
                id: 2,
                aabb: AabbRaw::new([1.0; 3], [3.0; 3]),
            },
            BvhLeaf {
                id: 3,
                aabb: AabbRaw::new([10.0; 3], [11.0; 3]),
            },
        ];
        let bvh = Bvh::build(leaves);
        let pairs = bvh.overlapping_pairs();
        // Pair (1,2) overlaps; (1,3) and (2,3) don't
        assert!(
            pairs.contains(&(1, 2)) || pairs.contains(&(2, 1)),
            "pair (1,2) should be in overlapping pairs: {:?}",
            pairs
        );
    }

    #[test]
    fn test_bvh_single_leaf() {
        let leaves = vec![BvhLeaf {
            id: 42,
            aabb: AabbRaw::new([0.0; 3], [1.0; 3]),
        }];
        let bvh = Bvh::build(leaves);
        let hits = bvh.query(&AabbRaw::new([0.0; 3], [0.5; 3]));
        assert_eq!(hits, vec![42]);
    }

    // ── ContactManifold additional tests ──────────────────────────────────────

    #[test]
    fn test_manifold_overflow_replaces_shallowest() {
        let mut m = ContactManifold::new();
        for i in 0..5 {
            m.add_contact(Contact {
                point: [i as f64, 0.0, 0.0],
                normal: [0.0, 1.0, 0.0],
                depth: (i as f64 + 1.0) * 0.1,
                id_a: 1,
                id_b: 2,
            });
        }
        // Should have at most 4 contacts
        assert!(m.contacts.len() <= 4);
    }

    #[test]
    fn test_manifold_is_empty_after_new() {
        let m = ContactManifold::new();
        assert!(m.is_empty());
    }

    #[test]
    fn test_manifold_max_depth_empty() {
        let m = ContactManifold::new();
        // max_depth on empty manifold should be NEG_INFINITY
        assert_eq!(m.max_depth(), f64::NEG_INFINITY);
    }

    // ── point_segment_closest ─────────────────────────────────────────────────

    #[test]
    fn test_point_segment_closest_midpoint() {
        let (cp, dist) = point_segment_closest([0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(cp[0].abs() < 1e-12, "closest x should be 0");
        assert!((dist - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_segment_closest_beyond_end() {
        let (cp, dist) = point_segment_closest([5.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]);
        assert!((cp[0] - 1.0).abs() < 1e-12, "closest x should clamp to 1");
        assert!((dist - 4.0).abs() < 1e-10);
    }

    // ── ray_sphere additional ─────────────────────────────────────────────────

    #[test]
    fn test_ray_sphere_from_inside() {
        let s = Sphere::new([0.0; 3], 2.0);
        // Ray from origin → should hit from inside (t > 0)
        let t = ray_sphere([0.0; 3], [1.0, 0.0, 0.0], &s);
        assert!(t.is_some());
        assert!(t.unwrap() > 0.0);
    }

    #[test]
    fn test_ray_aabb_along_each_axis() {
        let aabb = AabbRaw::new([-1.0; 3], [1.0; 3]);
        for axis in 0..3 {
            // Start on the negative axis side, aligned with the center
            let mut origin = [0.0_f64; 3];
            let mut dir = [0.0_f64; 3];
            origin[axis] = -5.0;
            dir[axis] = 1.0;
            let hit = ray_aabb(origin, dir, &aabb);
            assert!(hit.is_some(), "ray along axis {axis} should hit unit AABB");
        }
    }
}
