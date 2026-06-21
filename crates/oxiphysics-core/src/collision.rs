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
    gjk_simplex(a, b).0
}

/// GJK core that also returns the terminating simplex when the shapes overlap.
///
/// On intersection the second element is the enclosing tetrahedron (the simplex
/// that contained the origin), which EPA uses to seed a polytope that genuinely
/// encloses the origin.  When the shapes are separated it is `None`.
fn gjk_simplex(a: &dyn ConvexShape, b: &dyn ConvexShape) -> (GjkResult, Option<Vec<[Real; 3]>>) {
    const MAX_ITER: usize = 64;
    const EPS: Real = 1e-10;

    let mut simplex: Vec<[Real; 3]> = Vec::with_capacity(4);
    let mut dir = [1.0, 0.0, 0.0_f64];

    let first = gjk_support(a, b, dir);
    simplex.push(first);
    dir = v3_neg(first);

    // Best separation lower bound found so far, together with the witness
    // direction that produced it.  `v3_dot(support, dir_unit)` is a valid lower
    // bound on the true minimum distance along that direction, so the largest
    // such value seen is our best honest distance estimate if the algorithm
    // never converges.
    let mut best_dist = Real::MAX; // distance of the closest simplex feature
    let mut best_dir = dir; // direction toward the origin from that feature

    for _ in 0..MAX_ITER {
        let len_sq = v3_len_sq(dir);
        if len_sq < EPS {
            // Direction degenerate — the closest simplex feature is the origin,
            // i.e. the Minkowski difference contains the origin: intersection.
            let s = (simplex.len() == 4).then(|| simplex.clone());
            return (
                GjkResult {
                    intersecting: true,
                    distance: 0.0,
                    closest_a: [0.0; 3],
                    closest_b: [0.0; 3],
                },
                s,
            );
        }

        // Track the closest-feature distance: |dir| is the distance from the
        // origin to the closest point on the current simplex (do_simplex sets
        // `dir` to the origin-ward vector from that closest feature).
        let cur_dist = len_sq.sqrt();
        if cur_dist < best_dist {
            best_dist = cur_dist;
            best_dir = dir;
        }

        let support = gjk_support(a, b, dir);

        // If the new support does not pass the origin, shapes do not intersect.
        if v3_dot(support, dir) < 0.0 {
            let dist = cur_dist;
            return (
                GjkResult {
                    intersecting: false,
                    distance: dist,
                    closest_a: a.support(dir),
                    closest_b: b.support(v3_neg(dir)),
                },
                None,
            );
        }

        simplex.push(support);

        if gjk_do_simplex(&mut simplex, &mut dir) {
            // Origin enclosed: `simplex` holds the terminating tetrahedron.
            let s = (simplex.len() == 4).then(|| simplex.clone());
            return (
                GjkResult {
                    intersecting: true,
                    distance: 0.0,
                    closest_a: [0.0; 3],
                    closest_b: [0.0; 3],
                },
                s,
            );
        }
    }

    // Non-convergence fallback (slow / near-degenerate geometry): report the
    // best estimate actually computed rather than a fabricated constant.  If the
    // closest simplex feature is essentially at the origin we honestly call it an
    // intersection; otherwise we report the best measured separation and witness
    // support points instead of pretending `intersecting = true` with zeros.
    if best_dist <= EPS.sqrt() {
        let s = (simplex.len() == 4).then(|| simplex.clone());
        (
            GjkResult {
                intersecting: true,
                distance: 0.0,
                closest_a: [0.0; 3],
                closest_b: [0.0; 3],
            },
            s,
        )
    } else {
        (
            GjkResult {
                intersecting: false,
                distance: best_dist,
                closest_a: a.support(best_dir),
                closest_b: b.support(v3_neg(best_dir)),
            },
            None,
        )
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

/// Outward-oriented normal and origin-distance of a polytope face.
struct EpaFace {
    /// Vertex indices into the polytope (consistent winding).
    verts: [usize; 3],
    /// Unit outward normal (points away from the origin / polytope interior).
    normal: [Real; 3],
    /// Signed distance from the origin to the face plane (≥ 0 once oriented).
    dist: Real,
}

/// Build an [`EpaFace`] from three polytope vertices in the order given,
/// computing the outward normal from that winding (CCW as seen from outside).
///
/// The polytope is kept with globally consistent CCW winding, so the cross
/// product `(b−a)×(c−a)` already points outward; the face distance is the
/// (possibly negative) projection of the plane onto that normal.  Returns
/// `None` for a degenerate (zero-area) triangle.
fn epa_make_face(polytope: &[[Real; 3]], i: usize, j: usize, k: usize) -> Option<EpaFace> {
    let a = polytope[i];
    let b = polytope[j];
    let c = polytope[k];
    let n = v3_cross(v3_sub(b, a), v3_sub(c, a));
    let len = v3_len(n);
    if len < 1e-12 {
        return None;
    }
    let normal = v3_scale(n, 1.0 / len);
    let dist = v3_dot(normal, a);
    Some(EpaFace {
        verts: [i, j, k],
        normal,
        dist,
    })
}

/// Build a seed-tetrahedron face `(i, j, k)` with winding chosen so its normal
/// points *away* from the opposite vertex `l` (i.e. outward).  This establishes
/// the globally consistent CCW orientation that the horizon-stitching expansion
/// then preserves.
fn epa_seed_face(
    polytope: &[[Real; 3]],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Option<EpaFace> {
    let f = epa_make_face(polytope, i, j, k)?;
    // If the opposite vertex is on the positive side of this face's plane, the
    // winding is inward — swap two vertices to flip it outward.
    if v3_dot(f.normal, v3_sub(polytope[l], polytope[i])) > 0.0 {
        epa_make_face(polytope, i, k, j)
    } else {
        Some(f)
    }
}

/// Check that a 4-vertex GJK simplex is a non-degenerate tetrahedron (has real
/// volume), so it can seed EPA.  Returns `false` for flat / collinear simplices.
fn epa_simplex_is_valid(s: &[[Real; 3]]) -> bool {
    if s.len() != 4 {
        return false;
    }
    let e1 = v3_sub(s[1], s[0]);
    let e2 = v3_sub(s[2], s[0]);
    let e3 = v3_sub(s[3], s[0]);
    // Six times the signed tetrahedron volume.
    v3_dot(v3_cross(e1, e2), e3).abs() > 1e-12
}

/// Run EPA to compute the penetration depth between two overlapping convex shapes.
///
/// Requires that GJK has confirmed the shapes are intersecting.
/// Returns `None` if the shapes are not intersecting.
///
/// Implements the genuine Expanding Polytope Algorithm: starting from a
/// non-degenerate tetrahedron enclosing the origin in Minkowski-difference
/// space, it repeatedly (1) selects the face closest to the origin, (2) queries
/// the support point along that face's outward normal, and (3) if that support
/// lies beyond the face, expands the polytope by deleting *all* faces visible
/// from the new point and re-triangulating the resulting horizon with consistent
/// outward winding.  The iteration converges when no face can be pushed further
/// out; the penetration depth is then the distance to the closest face and the
/// collision normal is that face's outward normal.
pub fn epa(a: &dyn ConvexShape, b: &dyn ConvexShape) -> Option<EpaResult> {
    // First confirm intersection with GJK and recover its terminating simplex.
    let (gjk_r, gjk_simplex) = gjk_simplex(a, b);
    if !gjk_r.intersecting {
        return None;
    }

    const MAX_ITER: usize = 128;
    const EPS: Real = 1e-8;

    // ── Seed the polytope with a tetrahedron that encloses the origin ────────
    // The GJK terminating simplex genuinely contains the origin; prefer it.
    // If it is unavailable or degenerate, fall back to an explicit construction.
    let mut polytope: Vec<[Real; 3]> = match gjk_simplex {
        Some(s) if epa_simplex_is_valid(&s) && epa_origin_in_tetra(&s) => s,
        _ => epa_seed_tetrahedron(a, b)?,
    };

    // Seed faces with globally consistent outward (CCW) winding: each face is
    // oriented away from the tetrahedron's opposite vertex.
    let mut faces: Vec<EpaFace> = Vec::new();
    for &[i, j, k, l] in &[[0usize, 1, 2, 3], [0, 1, 3, 2], [0, 2, 3, 1], [1, 2, 3, 0]] {
        if let Some(f) = epa_seed_face(&polytope, i, j, k, l) {
            faces.push(f);
        }
    }
    if faces.len() < 4 {
        return None;
    }

    // Running closest-face estimate.  As the polytope is refined, the distance
    // to its closest face increases monotonically toward the true penetration
    // depth, so the LAST (most-refined) estimate is the best one — this is what
    // the honest non-convergence fallback returns (never a fabricated constant).
    let mut last_dist = 0.0_f64;
    let mut last_normal = [0.0_f64; 3];

    for _ in 0..MAX_ITER {
        // Closest face to the origin (only outward-oriented, dist ≥ 0 faces).
        let mut min_idx = 0usize;
        let mut min_dist = Real::MAX;
        for (fi, f) in faces.iter().enumerate() {
            if f.dist < min_dist {
                min_dist = f.dist;
                min_idx = fi;
            }
        }
        let min_normal = faces[min_idx].normal;
        if v3_len_sq(min_normal) > 0.5 {
            last_dist = min_dist;
            last_normal = min_normal;
        }

        // Support along the outward normal of the closest face.
        let support = gjk_support(a, b, min_normal);
        let support_dist = v3_dot(min_normal, support);

        // Converged when the support cannot push the closest face outward any
        // further (absolute or relative to the current depth, so smooth shapes
        // such as spheres terminate in finitely many steps).
        let tol = EPS + 1e-6 * min_dist.abs();
        if support_dist - min_dist < tol {
            let contact = v3_scale(min_normal, min_dist);
            return Some(EpaResult {
                depth: min_dist,
                normal: min_normal,
                point: contact,
            });
        }

        // ── Expand: remove every face visible from `support`, build horizon ──
        let new_idx = polytope.len();
        polytope.push(support);

        // A face is visible if the new vertex is in front of its plane.
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        let mut kept: Vec<EpaFace> = Vec::with_capacity(faces.len());
        for f in faces.drain(..) {
            let visible = v3_dot(f.normal, v3_sub(support, polytope[f.verts[0]])) > EPS;
            if visible {
                // Add the three directed edges to the horizon set, cancelling
                // edges shared by two visible faces (interior edges).
                for &(p, q) in &[
                    (f.verts[0], f.verts[1]),
                    (f.verts[1], f.verts[2]),
                    (f.verts[2], f.verts[0]),
                ] {
                    epa_add_horizon_edge(&mut horizon, p, q);
                }
            } else {
                kept.push(f);
            }
        }

        faces = kept;
        // Stitch the horizon to the new vertex.
        for (p, q) in horizon {
            if let Some(f) = epa_make_face(&polytope, p, q, new_idx) {
                faces.push(f);
            }
        }
        if faces.is_empty() {
            break;
        }
    }

    // Non-convergence fallback (e.g. a smooth shape refined to the iteration
    // cap): return the LAST closest-face estimate actually computed — a real
    // depth and normal, never a hardcoded 0-depth axis-aligned constant.
    if v3_len_sq(last_normal) > 0.5 {
        let contact = v3_scale(last_normal, last_dist);
        Some(EpaResult {
            depth: last_dist,
            normal: last_normal,
            point: contact,
        })
    } else {
        None
    }
}

/// Insert a directed edge into the horizon list, cancelling its reverse if
/// already present (an interior edge shared by two visible faces).
fn epa_add_horizon_edge(horizon: &mut Vec<(usize, usize)>, p: usize, q: usize) {
    if let Some(pos) = horizon.iter().position(|&(a, b)| a == q && b == p) {
        horizon.swap_remove(pos);
    } else {
        horizon.push((p, q));
    }
}

/// Construct a non-degenerate seed tetrahedron of the Minkowski difference
/// `A ⊖ B` that **encloses the origin**.
///
/// Builds the simplex incrementally (point → segment spanning the origin →
/// triangle → tetrahedron straddling the origin) using support queries along
/// progressively constrained directions, exactly as a from-scratch GJK seed
/// would.  The fourth vertex is taken on the far side of the triangle from the
/// origin so the resulting tetrahedron contains it.  Returns `None` only when
/// the difference is genuinely lower-dimensional (no enclosing volume exists,
/// e.g. a boundary touch).
fn epa_seed_tetrahedron(a: &dyn ConvexShape, b: &dyn ConvexShape) -> Option<Vec<[Real; 3]>> {
    // Candidate spanning directions covering all octants.
    const AXES: [[Real; 3]; 6] = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];

    // Vertex 0: any support.
    let v0 = gjk_support(a, b, [1.0, 0.0, 0.0]);

    // Vertex 1: search the opposite direction; require a real edge length.
    let mut dir1 = v3_neg(v0);
    if v3_len_sq(dir1) < 1e-9 {
        dir1 = [-1.0, 0.0, 0.0];
    }
    let v1 = gjk_support(a, b, dir1);
    let e01 = v3_sub(v1, v0);
    if v3_len_sq(e01) < 1e-12 {
        return None;
    }

    // Vertex 2: maximise the area perpendicular to the segment v0→v1.  Try every
    // axis cross-product and keep the largest off-line support.
    let mut v2 = v0;
    let mut best_area = 0.0;
    for axis in AXES.iter() {
        let d = v3_cross(e01, *axis);
        if v3_len_sq(d) < 1e-12 {
            continue;
        }
        for cand_dir in [d, v3_neg(d)] {
            let cand = gjk_support(a, b, cand_dir);
            let area = v3_len_sq(v3_cross(v3_sub(cand, v0), e01));
            if area > best_area {
                best_area = area;
                v2 = cand;
            }
        }
    }
    if best_area < 1e-18 {
        return None;
    }

    // Vertex 3: off the triangle plane, on the side *away* from the origin so the
    // tetrahedron straddles it.
    let tri_n = v3_cross(v3_sub(v1, v0), v3_sub(v2, v0));
    if v3_len_sq(tri_n) < 1e-18 {
        return None;
    }
    // The origin sits at signed height -dot(tri_n, v0) relative to the plane;
    // pick the search direction that points away from the origin.
    let toward_origin = v3_dot(tri_n, v0); // >0 ⇒ origin is below the plane
    let search = if toward_origin > 0.0 {
        v3_neg(tri_n)
    } else {
        tri_n
    };
    let mut v3 = gjk_support(a, b, search);
    if v3_dot(tri_n, v3_sub(v3, v0)).abs() < 1e-9 {
        // Degenerate on this side — try the other.
        v3 = gjk_support(a, b, v3_neg(search));
        if v3_dot(tri_n, v3_sub(v3, v0)).abs() < 1e-9 {
            return None;
        }
    }

    let tetra = vec![v0, v1, v2, v3];
    if epa_simplex_is_valid(&tetra) && epa_origin_in_tetra(&tetra) {
        Some(tetra)
    } else {
        None
    }
}

/// Test whether the origin lies inside (or on) the tetrahedron `t` via the
/// same-sign-of-signed-volumes criterion.
fn epa_origin_in_tetra(t: &[[Real; 3]]) -> bool {
    if t.len() != 4 {
        return false;
    }
    // Signed volume (×6) of the tetra formed by replacing vertex `omit` with the
    // origin must share the sign of the full tetra for the origin to be inside.
    let signed = |p: [Real; 3], q: [Real; 3], r: [Real; 3], s: [Real; 3]| -> Real {
        v3_dot(v3_cross(v3_sub(q, p), v3_sub(r, p)), v3_sub(s, p))
    };
    let full = signed(t[0], t[1], t[2], t[3]);
    if full.abs() < 1e-15 {
        return false;
    }
    let o = [0.0_f64; 3];
    let d0 = signed(o, t[1], t[2], t[3]);
    let d1 = signed(t[0], o, t[2], t[3]);
    let d2 = signed(t[0], t[1], o, t[3]);
    let d3 = signed(t[0], t[1], t[2], o);
    let sgn = full.signum();
    // Allow a tiny tolerance so boundary cases (origin on a face) still seed.
    let tol = -1e-12 * full.abs();
    d0 * sgn >= tol && d1 * sgn >= tol && d2 * sgn >= tol && d3 * sgn >= tol
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
#[path = "collision_tests.rs"]
mod tests;
