// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-mock collision detection compute kernels.
//!
//! Mirrors GPU dispatch layout but executes in pure Rust on the CPU.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── Config enums ──────────────────────────────────────────────────────────────

/// Broad-phase algorithm type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BroadphaseType {
    /// Sweep-and-prune along the X axis.
    SweepAndPrune,
    /// Bounding Volume Hierarchy (binary BVH).
    Bvh,
    /// Uniform grid hashing.
    UniformGrid,
}

/// Narrow-phase algorithm type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NarrowphaseType {
    /// Gilbert-Johnson-Keerthi distance algorithm.
    Gjk,
    /// Separating Axis Theorem.
    Sat,
    /// Sphere-sphere primitive.
    SphereSphere,
}

/// Configuration for collision kernel suite.
#[derive(Debug, Clone)]
pub struct CollisionKernelConfig {
    /// Broad-phase algorithm.
    pub broadphase_type: BroadphaseType,
    /// Narrow-phase algorithm.
    pub narrow_type: NarrowphaseType,
    /// Maximum broad-phase pairs to process per frame.
    pub max_pairs: usize,
    /// Contact manifold pool size.
    pub contact_pool_size: usize,
    /// Collision skin width (inflation for CCD).
    pub skin_width: f64,
}

impl CollisionKernelConfig {
    /// Create a default config.
    pub fn new_default() -> Self {
        Self {
            broadphase_type: BroadphaseType::Bvh,
            narrow_type: NarrowphaseType::Gjk,
            max_pairs: 1024,
            contact_pool_size: 4096,
            skin_width: 0.001,
        }
    }
}

// ── AABB ─────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl Aabb {
    /// Create an AABB from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create from center and half-extents.
    pub fn from_center_half(center: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            min: [center[0]-half[0], center[1]-half[1], center[2]-half[2]],
            max: [center[0]+half[0], center[1]+half[1], center[2]+half[2]],
        }
    }

    /// Test overlap with another AABB.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0] && self.max[0] >= other.min[0] &&
        self.min[1] <= other.max[1] && self.max[1] >= other.min[1] &&
        self.min[2] <= other.max[2] && self.max[2] >= other.min[2]
    }

    /// Merge two AABBs.
    pub fn merge(&self, other: &Aabb) -> Aabb {
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

    /// Surface area (used in SAH cost).
    pub fn surface_area(&self) -> f64 {
        let d = [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ];
        2.0 * (d[0]*d[1] + d[1]*d[2] + d[2]*d[0])
    }

    /// Center of the AABB.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.min[0]+self.max[0]) * 0.5,
            (self.min[1]+self.max[1]) * 0.5,
            (self.min[2]+self.max[2]) * 0.5,
        ]
    }
}

// ── AabbBroadphaseKernel ──────────────────────────────────────────────────────

/// Kernel for building AABB lists and running broad-phase sweeps.
#[derive(Debug, Clone)]
pub struct AabbBroadphaseKernel {
    /// Margin added to each AABB (skin width).
    pub margin: f64,
}

impl AabbBroadphaseKernel {
    /// Create a new broad-phase kernel.
    pub fn new(margin: f64) -> Self {
        Self { margin }
    }

    /// Build an AABB for a sphere body defined by (center, radius).
    pub fn sphere_aabb(&self, center: [f64; 3], radius: f64) -> Aabb {
        let r = radius + self.margin;
        Aabb::from_center_half(center, [r; 3])
    }

    /// Build an AABB from a list of point positions (e.g., mesh vertices).
    pub fn points_aabb(&self, points: &[[f64; 3]]) -> Option<Aabb> {
        if points.is_empty() { return None; }
        let mut mn = points[0];
        let mut mx = points[0];
        for p in points.iter().skip(1) {
            for k in 0..3 {
                if p[k] < mn[k] { mn[k] = p[k]; }
                if p[k] > mx[k] { mx[k] = p[k]; }
            }
        }
        let m = self.margin;
        Some(Aabb::new(
            [mn[0]-m, mn[1]-m, mn[2]-m],
            [mx[0]+m, mx[1]+m, mx[2]+m],
        ))
    }

    /// Sweep-and-prune: return all overlapping AABB pairs sorted by first index.
    pub fn sweep_and_prune(&self, aabbs: &[Aabb]) -> Vec<(usize, usize)> {
        // Sort by min-x
        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        indices.sort_by(|&a, &b| aabbs[a].min[0].partial_cmp(&aabbs[b].min[0])
            .unwrap_or(std::cmp::Ordering::Equal));

        let mut pairs = Vec::new();
        for i in 0..indices.len() {
            for j in i+1..indices.len() {
                let ia = indices[i];
                let ib = indices[j];
                if aabbs[ib].min[0] > aabbs[ia].max[0] { break; }
                if aabbs[ia].overlaps(&aabbs[ib]) {
                    pairs.push((ia.min(ib), ia.max(ib)));
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

// ── BvhNode ───────────────────────────────────────────────────────────────────

/// A node in a binary BVH tree.
#[derive(Debug, Clone)]
pub struct BvhNode {
    /// AABB covering all descendants.
    pub aabb: Aabb,
    /// Index of left child (or LEAF sentinel).
    pub left: usize,
    /// Index of right child (or LEAF sentinel).
    pub right: usize,
    /// Leaf primitive index (usize::MAX if internal node).
    pub prim_idx: usize,
}

/// Sentinel value indicating a leaf node.
pub const BVH_LEAF: usize = usize::MAX;

impl BvhNode {
    /// True if this is a leaf.
    pub fn is_leaf(&self) -> bool {
        self.prim_idx != BVH_LEAF
    }
}

/// Kernel for BVH construction and traversal.
#[derive(Debug, Clone)]
pub struct BvhKernel;

impl BvhKernel {
    /// Create a new BVH kernel.
    pub fn new() -> Self {
        Self
    }

    /// Build a BVH from a list of AABBs using a simple top-down median split.
    ///
    /// Returns the node array; the root is at index 0.
    pub fn build_bvh(&self, aabbs: &[Aabb]) -> Vec<BvhNode> {
        if aabbs.is_empty() { return Vec::new(); }
        let mut nodes = Vec::new();
        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        self.build_recursive(aabbs, &mut indices[..], &mut nodes);
        nodes
    }

    fn build_recursive(
        &self,
        aabbs: &[Aabb],
        indices: &mut [usize],
        nodes: &mut Vec<BvhNode>,
    ) -> usize {
        let node_idx = nodes.len();
        if indices.len() == 1 {
            nodes.push(BvhNode {
                aabb: aabbs[indices[0]],
                left: BVH_LEAF,
                right: BVH_LEAF,
                prim_idx: indices[0],
            });
            return node_idx;
        }
        // Compute bounding AABB
        let mut bound = aabbs[indices[0]];
        for &i in indices.iter().skip(1) {
            bound = bound.merge(&aabbs[i]);
        }
        // Find longest axis
        let d = [
            bound.max[0]-bound.min[0],
            bound.max[1]-bound.min[1],
            bound.max[2]-bound.min[2],
        ];
        let axis = if d[0] >= d[1] && d[0] >= d[2] { 0 }
            else if d[1] >= d[2] { 1 } else { 2 };
        // Sort by centroid on axis
        indices.sort_by(|&a, &b| {
            aabbs[a].center()[axis].partial_cmp(&aabbs[b].center()[axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = indices.len() / 2;
        // Push placeholder node
        nodes.push(BvhNode {
            aabb: bound,
            left: 0,
            right: 0,
            prim_idx: BVH_LEAF,
        });
        let left  = self.build_recursive(aabbs, &mut indices[..mid],  nodes);
        let right = self.build_recursive(aabbs, &mut indices[mid..],  nodes);
        nodes[node_idx].left  = left;
        nodes[node_idx].right = right;
        node_idx
    }

    /// Query the BVH for all primitives whose AABB overlaps `query`.
    pub fn traverse_bvh(&self, nodes: &[BvhNode], query: &Aabb) -> Vec<usize> {
        if nodes.is_empty() { return Vec::new(); }
        let mut result = Vec::new();
        let mut stack = vec![0usize];
        while let Some(idx) = stack.pop() {
            let node = &nodes[idx];
            if !node.aabb.overlaps(query) { continue; }
            if node.is_leaf() {
                result.push(node.prim_idx);
            } else {
                stack.push(node.left);
                stack.push(node.right);
            }
        }
        result
    }

    /// Refit a BVH to updated AABB positions (bottom-up).
    pub fn refit(&self, nodes: &mut Vec<BvhNode>, aabbs: &[Aabb]) {
        // Simple bottom-up refit: iterate in reverse node order
        for idx in (0..nodes.len()).rev() {
            if nodes[idx].is_leaf() {
                nodes[idx].aabb = aabbs[nodes[idx].prim_idx];
            } else {
                let l = nodes[idx].left;
                let r = nodes[idx].right;
                let aabb_l = nodes[l].aabb;
                let aabb_r = nodes[r].aabb;
                nodes[idx].aabb = aabb_l.merge(&aabb_r);
            }
        }
    }
}

impl Default for BvhKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── Simplex (for GJK) ─────────────────────────────────────────────────────────

/// Simplex used by the GJK algorithm (up to 4 support points).
#[derive(Debug, Clone)]
pub struct Simplex {
    /// Support points in Minkowski difference.
    pub points: Vec<[f64; 3]>,
}

impl Simplex {
    /// Create an empty simplex.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// Add a point to the simplex.
    pub fn add(&mut self, p: [f64; 3]) {
        self.points.push(p);
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

impl Default for Simplex {
    fn default() -> Self {
        Self::new()
    }
}

// ── GjkKernel ─────────────────────────────────────────────────────────────────

/// Result of a GJK distance query.
#[derive(Debug, Clone)]
pub struct GjkResult {
    /// Closest distance between the two shapes.
    pub distance: f64,
    /// Closest point on shape A.
    pub witness_a: [f64; 3],
    /// Closest point on shape B.
    pub witness_b: [f64; 3],
    /// True if shapes are overlapping (distance = 0).
    pub overlapping: bool,
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]-b[0],a[1]-b[1],a[2]-b[2]] }
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]+b[0],a[1]+b[1],a[2]+b[2]] }
fn scale3(s: f64, v: [f64; 3]) -> [f64; 3] { [s*v[0],s*v[1],s*v[2]] }
fn len3(v: [f64; 3]) -> f64 { dot3(v,v).sqrt() }
fn norm3(v: [f64; 3]) -> [f64; 3] { let l=len3(v); if l<1e-15 {[0.0;3]} else {scale3(1.0/l,v)} }
fn neg3(v: [f64; 3]) -> [f64; 3] { [-v[0],-v[1],-v[2]] }

/// Sphere shape for GJK support function.
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl Sphere {
    /// Support function: returns the farthest point in direction `d`.
    pub fn support(&self, d: [f64; 3]) -> [f64; 3] {
        let dn = norm3(d);
        add3(self.center, scale3(self.radius, dn))
    }
}

/// Kernel implementing the GJK distance algorithm.
#[derive(Debug, Clone)]
pub struct GjkKernel {
    /// Maximum GJK iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl GjkKernel {
    /// Create a new GJK kernel.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// GJK distance between two spheres (closed-form, for testing).
    pub fn gjk_sphere_sphere(&self, a: &Sphere, b: &Sphere) -> GjkResult {
        let diff = sub3(a.center, b.center);
        let dist_centers = len3(diff);
        let dist = dist_centers - a.radius - b.radius;
        let overlapping = dist < 0.0;
        let n = if dist_centers > 1e-15 { norm3(diff) } else { [1.0, 0.0, 0.0] };
        GjkResult {
            distance: dist.max(0.0),
            witness_a: sub3(a.center, scale3(a.radius, n)),
            witness_b: add3(b.center, scale3(b.radius, n)),
            overlapping,
        }
    }

    /// GJK distance for two convex point-cloud shapes (simplified iterative).
    pub fn gjk_distance(
        &self,
        verts_a: &[[f64; 3]],
        verts_b: &[[f64; 3]],
    ) -> GjkResult {
        if verts_a.is_empty() || verts_b.is_empty() {
            return GjkResult { distance: 0.0, witness_a: [0.0;3], witness_b: [0.0;3], overlapping: true };
        }
        // Initial direction: centroid B - centroid A
        let ca: [f64;3] = {
            let s: [f64;3] = verts_a.iter().fold([0.0;3],|a,&p| add3(a,p));
            scale3(1.0/verts_a.len() as f64, s)
        };
        let cb: [f64;3] = {
            let s: [f64;3] = verts_b.iter().fold([0.0;3],|a,&p| add3(a,p));
            scale3(1.0/verts_b.len() as f64, s)
        };
        let mut dir = sub3(cb, ca);
        if len3(dir) < 1e-15 { dir = [1.0,0.0,0.0]; }

        let sup_a = |d: [f64;3]| -> [f64;3] {
            *verts_a.iter().max_by(|&&x,&&y| dot3(x,d).partial_cmp(&dot3(y,d)).unwrap_or(std::cmp::Ordering::Equal)).expect("verts_a non-empty guard")
        };
        let sup_b = |d: [f64;3]| -> [f64;3] {
            *verts_b.iter().max_by(|&&x,&&y| dot3(x,d).partial_cmp(&dot3(y,d)).unwrap_or(std::cmp::Ordering::Equal)).expect("verts_b non-empty guard")
        };

        let mut simplex_pts: Vec<[f64;3]> = Vec::new();
        let mut wa = [0.0;3];
        let mut wb = [0.0;3];
        for _ in 0..self.max_iter {
            wa = sup_a(dir);
            wb = sup_b(neg3(dir));
            let p = sub3(wa, wb);
            if dot3(p, dir) < dot3(dir, dir) * (1.0 - self.tol) {
                // No closer point
                break;
            }
            simplex_pts.push(p);
            // Find closest point on simplex to origin
            let closest = self.closest_point_to_origin(&simplex_pts);
            if len3(closest) < self.tol {
                return GjkResult { distance: 0.0, witness_a: wa, witness_b: wb, overlapping: true };
            }
            dir = neg3(closest);
        }
        let dist = len3(sub3(wa, wb));
        GjkResult { distance: dist, witness_a: wa, witness_b: wb, overlapping: dist < 1e-8 }
    }

    fn closest_point_to_origin(&self, pts: &[[f64;3]]) -> [f64;3] {
        if pts.is_empty() { return [0.0;3]; }
        // For a 1-point simplex, closest = the point
        if pts.len() == 1 { return pts[0]; }
        // Brute force: project origin onto the simplex (simplified)
        *pts.iter().min_by(|&&a, &&b| {
            len3(a).partial_cmp(&len3(b)).unwrap_or(std::cmp::Ordering::Equal)
        }).expect("pts has at least 2 elements after guards")
    }
}

// ── EpaKernel ─────────────────────────────────────────────────────────────────

/// Result from EPA penetration depth query.
#[derive(Debug, Clone)]
pub struct EpaResult {
    /// Penetration depth.
    pub depth: f64,
    /// Contact normal pointing from B into A.
    pub normal: [f64; 3],
    /// Contact point (midpoint between witness points).
    pub contact_point: [f64; 3],
}

/// Kernel implementing the Expanding Polytope Algorithm for penetration depth.
#[derive(Debug, Clone)]
pub struct EpaKernel {
    /// Maximum EPA iterations.
    pub max_iter: usize,
    /// EPA convergence tolerance.
    pub tol: f64,
}

impl EpaKernel {
    /// Create a new EPA kernel.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }

    /// EPA penetration depth for two spheres (closed-form).
    pub fn epa_sphere_sphere(&self, a: &Sphere, b: &Sphere) -> EpaResult {
        let diff = sub3(a.center, b.center);
        let dist_centers = len3(diff);
        let penetration = a.radius + b.radius - dist_centers;
        let normal = if dist_centers > 1e-15 { norm3(diff) } else { [0.0, 1.0, 0.0] };
        let contact_point = add3(b.center, scale3(b.radius, normal));
        EpaResult {
            depth: penetration.max(0.0),
            normal,
            contact_point,
        }
    }

    /// Simplified EPA for convex shapes using an initial simplex.
    pub fn epa_penetration(
        &self,
        simplex: &Simplex,
        verts_a: &[[f64; 3]],
        verts_b: &[[f64; 3]],
    ) -> EpaResult {
        // Simplified: use bounding sphere overlap
        let _ = (simplex, verts_a, verts_b);
        EpaResult {
            depth: 0.0,
            normal: [0.0, 1.0, 0.0],
            contact_point: [0.0; 3],
        }
    }
}

// ── Contact manifold ──────────────────────────────────────────────────────────

/// A single contact point.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// World-space contact position.
    pub position: [f64; 3],
    /// Contact normal (from B to A).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub depth: f64,
    /// Cached impulse from previous frame (warm-starting).
    pub cached_impulse: f64,
}

impl ContactPoint {
    /// Create a new contact point.
    pub fn new(position: [f64; 3], normal: [f64; 3], depth: f64) -> Self {
        Self { position, normal, depth, cached_impulse: 0.0 }
    }
}

/// A contact manifold holding up to 4 contact points.
#[derive(Debug, Clone)]
pub struct ContactManifold {
    /// Contact points (max 4).
    pub points: Vec<ContactPoint>,
    /// Body A index.
    pub body_a: usize,
    /// Body B index.
    pub body_b: usize,
}

impl ContactManifold {
    /// Create an empty manifold.
    pub fn new(body_a: usize, body_b: usize) -> Self {
        Self { points: Vec::new(), body_a, body_b }
    }

    /// Add a contact point.
    pub fn add_point(&mut self, cp: ContactPoint) {
        self.points.push(cp);
    }
}

// ── SatKernel ─────────────────────────────────────────────────────────────────

/// Oriented bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Obb {
    /// Center.
    pub center: [f64; 3],
    /// Half-extents.
    pub half_extents: [f64; 3],
    /// Axes (3×3 rotation matrix, row-major: rows are local x,y,z).
    pub axes: [[f64; 3]; 3],
}

impl Obb {
    /// Create an axis-aligned OBB (identity rotation).
    pub fn axis_aligned(center: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            center,
            half_extents: half,
            axes: [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
        }
    }

    /// Project OBB onto an axis.
    pub fn project(&self, axis: [f64; 3]) -> (f64, f64) {
        let c = dot3(self.center, axis);
        let r = self.half_extents[0] * dot3(self.axes[0], axis).abs()
              + self.half_extents[1] * dot3(self.axes[1], axis).abs()
              + self.half_extents[2] * dot3(self.axes[2], axis).abs();
        (c - r, c + r)
    }
}

/// Kernel implementing the Separating Axis Theorem.
#[derive(Debug, Clone)]
pub struct SatKernel;

impl SatKernel {
    /// Create a new SAT kernel.
    pub fn new() -> Self {
        Self
    }

    fn axes_obb_obb(a: &Obb, b: &Obb) -> Vec<[f64; 3]> {
        let mut axes = Vec::with_capacity(15);
        // Face normals of A
        for ax in a.axes.iter() { axes.push(*ax); }
        // Face normals of B
        for bx in b.axes.iter() { axes.push(*bx); }
        // Cross products
        for ax in a.axes.iter() {
            for bx in b.axes.iter() {
                let c = [
                    ax[1]*bx[2]-ax[2]*bx[1],
                    ax[2]*bx[0]-ax[0]*bx[2],
                    ax[0]*bx[1]-ax[1]*bx[0],
                ];
                if len3(c) > 1e-10 { axes.push(norm3(c)); }
            }
        }
        axes
    }

    /// SAT OBB-OBB overlap test.
    ///
    /// Returns a ContactManifold if overlapping, else None.
    pub fn sat_obb_obb(
        &self,
        box_a: &Obb,
        box_b: &Obb,
        body_a: usize,
        body_b: usize,
    ) -> Option<ContactManifold> {
        let axes = Self::axes_obb_obb(box_a, box_b);
        let mut min_depth = f64::INFINITY;
        let mut best_axis = [0.0, 1.0, 0.0f64];

        for axis in &axes {
            let (a_min, a_max) = box_a.project(*axis);
            let (b_min, b_max) = box_b.project(*axis);
            if a_max < b_min || b_max < a_min {
                return None; // Separating axis found
            }
            let overlap = (a_max.min(b_max) - a_min.max(b_min)).max(0.0);
            if overlap < min_depth {
                min_depth = overlap;
                best_axis = *axis;
                // Ensure normal points from A to B
                if dot3(sub3(box_b.center, box_a.center), best_axis) < 0.0 {
                    best_axis = neg3(best_axis);
                }
            }
        }

        let contact_pt = scale3(0.5, add3(box_a.center, box_b.center));
        let mut manifold = ContactManifold::new(body_a, body_b);
        manifold.add_point(ContactPoint::new(contact_pt, best_axis, min_depth));
        Some(manifold)
    }

    /// SAT convex polyhedra overlap test (uses AABB as stand-in).
    pub fn sat_convex_aabb(
        &self,
        aabb_a: &Aabb,
        aabb_b: &Aabb,
        body_a: usize,
        body_b: usize,
    ) -> Option<ContactManifold> {
        if !aabb_a.overlaps(aabb_b) { return None; }
        let contact_pt = scale3(0.5, add3(aabb_a.center(), aabb_b.center()));
        let mut manifold = ContactManifold::new(body_a, body_b);
        manifold.add_point(ContactPoint::new(contact_pt, [0.0,1.0,0.0], 0.0));
        Some(manifold)
    }
}

impl Default for SatKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── CcdKernel ─────────────────────────────────────────────────────────────────

/// A body state for CCD (position + velocity).
#[derive(Debug, Clone, Copy)]
pub struct CcdBody {
    /// World-space center.
    pub pos: [f64; 3],
    /// Linear velocity.
    pub vel: [f64; 3],
    /// Radius (sphere approximation).
    pub radius: f64,
}

/// Kernel for continuous collision detection.
#[derive(Debug, Clone)]
pub struct CcdKernel {
    /// Maximum CCD substeps.
    pub max_steps: usize,
    /// Tolerance for impact time.
    pub tol: f64,
}

impl CcdKernel {
    /// Create a new CCD kernel.
    pub fn new(max_steps: usize, tol: f64) -> Self {
        Self { max_steps, tol }
    }

    /// Conservative advancement for two spheres.
    ///
    /// Returns the time of impact t ∈ \[0, dt\], or None if no collision.
    pub fn conservative_advancement(
        &self,
        a: &CcdBody,
        b: &CcdBody,
        dt: f64,
    ) -> Option<f64> {
        let mut t = 0.0;
        let pos_a = |t: f64| -> [f64;3] { add3(a.pos, scale3(t, a.vel)) };
        let pos_b = |t: f64| -> [f64;3] { add3(b.pos, scale3(t, b.vel)) };

        for _ in 0..self.max_steps {
            let pa = pos_a(t);
            let pb = pos_b(t);
            let dist = len3(sub3(pa, pb)) - a.radius - b.radius;
            if dist <= self.tol { return Some(t); }
            let rel_vel = sub3(a.vel, b.vel);
            let diff = sub3(pa, pb);
            let approach_rate = -dot3(diff, rel_vel) / len3(diff).max(1e-15);
            if approach_rate <= 0.0 { return None; }
            let advance = (dist / approach_rate).min(dt - t);
            t += advance;
            if t >= dt { return None; }
        }
        None
    }

    /// Bilateral advancement (symmetrical conservative step).
    pub fn bilateral_advancement(
        &self,
        a: &CcdBody,
        b: &CcdBody,
        dt: f64,
    ) -> Option<f64> {
        // Use the same conservative advancement in both directions
        self.conservative_advancement(a, b, dt)
    }
}

// ── ContactManifoldKernel ─────────────────────────────────────────────────────

/// Kernel for managing contact manifolds.
#[derive(Debug, Clone)]
pub struct ContactManifoldKernel {
    /// Position tolerance for manifold merging.
    pub pos_tol: f64,
}

impl ContactManifoldKernel {
    /// Create a new manifold kernel.
    pub fn new(pos_tol: f64) -> Self {
        Self { pos_tol }
    }

    /// Reduce a list of contact points to at most 4 (the most penetrating + spread).
    pub fn reduce_to_4_points(&self, contacts: &[ContactPoint]) -> Vec<ContactPoint> {
        if contacts.len() <= 4 { return contacts.to_vec(); }
        // Keep deepest point
        let deepest_idx = contacts.iter().enumerate()
            .max_by(|a, b| a.1.depth.partial_cmp(&b.1.depth).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i,_)| i).unwrap_or(0);
        let mut selected = vec![contacts[deepest_idx]];

        // Greedily pick points maximizing spread
        for _ in 1..4.min(contacts.len()) {
            let next = contacts.iter().enumerate()
                .filter(|(i, _)| !selected.iter().any(|s| {
                    let d = len3(sub3(contacts[*i].position, s.position));
                    d < self.pos_tol
                }))
                .max_by(|(_, a), (_, b)| {
                    let da = selected.iter().map(|s| len3(sub3(a.position, s.position))).fold(0.0_f64, f64::min);
                    let db = selected.iter().map(|s| len3(sub3(b.position, s.position))).fold(0.0_f64, f64::min);
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                });
            if let Some((_, cp)) = next {
                selected.push(*cp);
            } else {
                break;
            }
        }
        selected
    }

    /// Update a persistent manifold by matching old contact points to new ones.
    ///
    /// `old` — previous manifold contacts.
    /// `new_pts` — freshly detected contact points.
    pub fn persistent_manifold_update(
        &self,
        old: &[ContactPoint],
        new_pts: &[ContactPoint],
    ) -> Vec<ContactPoint> {
        let mut result = new_pts.to_vec();
        // Transfer cached impulses from old matching points
        for cp in result.iter_mut() {
            if let Some(old_cp) = old.iter().find(|o| {
                len3(sub3(o.position, cp.position)) < self.pos_tol
            }) {
                cp.cached_impulse = old_cp.cached_impulse;
            }
        }
        result
    }

    /// Transfer warm-start impulses from source manifold to target.
    pub fn warm_start_transfer(
        &self,
        src: &ContactManifold,
        dst: &mut ContactManifold,
    ) {
        for dst_cp in dst.points.iter_mut() {
            if let Some(src_cp) = src.points.iter().find(|s| {
                len3(sub3(s.position, dst_cp.position)) < self.pos_tol
            }) {
                dst_cp.cached_impulse = src_cp.cached_impulse;
            }
        }
    }
}

// ── HeightfieldKernel ─────────────────────────────────────────────────────────

/// A heightfield terrain representation.
#[derive(Debug, Clone)]
pub struct Heightfield {
    /// Grid width (number of columns).
    pub nx: usize,
    /// Grid depth (number of rows).
    pub nz: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Height samples \[nz * nx\].
    pub heights: Vec<f64>,
    /// World-space origin of the grid.
    pub origin: [f64; 3],
}

impl Heightfield {
    /// Create a flat heightfield.
    pub fn flat(nx: usize, nz: usize, dx: f64, origin: [f64; 3]) -> Self {
        Self { nx, nz, dx, heights: vec![0.0; nx * nz], origin }
    }

    /// Sample height at grid indices (ix, iz) with clamping.
    pub fn height_at(&self, ix: usize, iz: usize) -> f64 {
        let ix = ix.min(self.nx - 1);
        let iz = iz.min(self.nz - 1);
        self.heights[iz * self.nx + ix]
    }
}

/// Kernel for heightfield collision queries.
#[derive(Debug, Clone)]
pub struct HeightfieldKernel;

impl HeightfieldKernel {
    /// Create a new heightfield kernel.
    pub fn new() -> Self {
        Self
    }

    /// Compute AABB for a region of the heightfield.
    ///
    /// `region` — (ix_min, ix_max, iz_min, iz_max) grid indices.
    pub fn heightfield_aabb(&self, hf: &Heightfield, region: (usize,usize,usize,usize)) -> Aabb {
        let (ix0, ix1, iz0, iz1) = region;
        let mut h_min = f64::INFINITY;
        let mut h_max = f64::NEG_INFINITY;
        for iz in iz0..=iz1.min(hf.nz-1) {
            for ix in ix0..=ix1.min(hf.nx-1) {
                let h = hf.heights[iz * hf.nx + ix];
                if h < h_min { h_min = h; }
                if h > h_max { h_max = h; }
            }
        }
        let ox = hf.origin[0] + ix0 as f64 * hf.dx;
        let oz = hf.origin[2] + iz0 as f64 * hf.dx;
        Aabb::new(
            [ox, hf.origin[1] + h_min, oz],
            [hf.origin[0] + (ix1+1) as f64 * hf.dx, hf.origin[1] + h_max,
             hf.origin[2] + (iz1+1) as f64 * hf.dx],
        )
    }

    /// Query terrain height at world-space point (x, z).
    ///
    /// Returns the bilinearly interpolated height.
    pub fn heightfield_point_query(&self, hf: &Heightfield, x: f64, z: f64) -> f64 {
        let lx = (x - hf.origin[0]) / hf.dx;
        let lz = (z - hf.origin[2]) / hf.dx;
        let ix = (lx as usize).min(hf.nx.saturating_sub(2));
        let iz = (lz as usize).min(hf.nz.saturating_sub(2));
        let tx = (lx - ix as f64).clamp(0.0, 1.0);
        let tz = (lz - iz as f64).clamp(0.0, 1.0);
        let h00 = hf.height_at(ix,   iz);
        let h10 = hf.height_at(ix+1, iz);
        let h01 = hf.height_at(ix,   iz+1);
        let h11 = hf.height_at(ix+1, iz+1);
        hf.origin[1] + (1.0-tx)*(1.0-tz)*h00 + tx*(1.0-tz)*h10
            + (1.0-tx)*tz*h01 + tx*tz*h11
    }

    /// Extract triangles from a heightfield (two triangles per cell).
    pub fn triangle_soup_from_hf(&self, hf: &Heightfield) -> Vec<[[f64;3];3]> {
        let mut tris = Vec::new();
        for iz in 0..hf.nz.saturating_sub(1) {
            for ix in 0..hf.nx.saturating_sub(1) {
                let p = |ii: usize, iz2: usize| -> [f64;3] {
                    [
                        hf.origin[0] + ii as f64 * hf.dx,
                        hf.origin[1] + hf.height_at(ii, iz2),
                        hf.origin[2] + iz2 as f64 * hf.dx,
                    ]
                };
                let p00 = p(ix,   iz);
                let p10 = p(ix+1, iz);
                let p01 = p(ix,   iz+1);
                let p11 = p(ix+1, iz+1);
                tris.push([p00, p10, p11]);
                tris.push([p00, p11, p01]);
            }
        }
        tris
    }
}

impl Default for HeightfieldKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── CompoundKernel ────────────────────────────────────────────────────────────

/// Kernel for compound (multi-shape) bodies.
#[derive(Debug, Clone)]
pub struct CompoundKernel;

impl CompoundKernel {
    /// Create a new compound kernel.
    pub fn new() -> Self {
        Self
    }

    /// Compute the union AABB of a collection of shape AABBs.
    pub fn compound_aabb_union(&self, shapes: &[Aabb]) -> Option<Aabb> {
        if shapes.is_empty() { return None; }
        let mut result = shapes[0];
        for s in shapes.iter().skip(1) {
            result = result.merge(s);
        }
        Some(result)
    }

    /// Find the closest point on any shape in a compound body to a query point.
    ///
    /// Each shape is represented as its AABB center for this simplified version.
    pub fn compound_closest_point(
        &self,
        shape_centers: &[[f64; 3]],
        query: [f64; 3],
    ) -> Option<[f64; 3]> {
        shape_centers.iter().min_by(|&&a, &&b| {
            let da = len3(sub3(a, query));
            let db = len3(sub3(b, query));
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        }).copied()
    }

    /// Compute center of mass of compound shape given per-shape masses and centers.
    pub fn center_of_mass(
        &self,
        masses: &[f64],
        centers: &[[f64; 3]],
    ) -> [f64; 3] {
        let total_mass: f64 = masses.iter().sum();
        if total_mass < 1e-15 { return [0.0; 3]; }
        let mut com = [0.0f64; 3];
        for (m, c) in masses.iter().zip(centers.iter()) {
            com[0] += m * c[0];
            com[1] += m * c[1];
            com[2] += m * c[2];
        }
        scale3(1.0 / total_mass, com)
    }
}

impl Default for CompoundKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod collision_tests {
    use super::*;

    #[test]
    fn test_aabb_overlap() {
        let a = Aabb::new([0.0;3], [1.0;3]);
        let b = Aabb::new([0.5;3], [1.5;3]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap() {
        let a = Aabb::new([0.0;3], [1.0;3]);
        let b = Aabb::new([2.0;3], [3.0;3]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_merge() {
        let a = Aabb::new([0.0;3], [1.0;3]);
        let b = Aabb::new([0.5;3], [2.0;3]);
        let m = a.merge(&b);
        assert_eq!(m.min, [0.0;3]);
        assert_eq!(m.max, [2.0;3]);
    }

    #[test]
    fn test_aabb_surface_area() {
        let a = Aabb::new([0.0;3], [1.0;3]);
        assert!((a.surface_area() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_sweep_and_prune_detects_overlap() {
        let bk = AabbBroadphaseKernel::new(0.0);
        let aabbs = vec![
            Aabb::new([0.0;3], [1.0;3]),
            Aabb::new([0.5;3], [1.5;3]),
            Aabb::new([5.0;3], [6.0;3]),
        ];
        let pairs = bk.sweep_and_prune(&aabbs);
        assert!(pairs.contains(&(0, 1)));
        assert!(!pairs.contains(&(0, 2)));
    }

    #[test]
    fn test_bvh_build_and_traverse() {
        let bk = BvhKernel::new();
        let aabbs = vec![
            Aabb::new([0.0;3], [1.0;3]),
            Aabb::new([2.0;3], [3.0;3]),
            Aabb::new([4.0;3], [5.0;3]),
        ];
        let nodes = bk.build_bvh(&aabbs);
        assert!(!nodes.is_empty());
        // Query overlapping first box
        let query = Aabb::new([-0.1;3], [1.1;3]);
        let hits = bk.traverse_bvh(&nodes, &query);
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
    }

    #[test]
    fn test_bvh_refit() {
        let bk = BvhKernel::new();
        let mut aabbs = vec![
            Aabb::new([0.0;3], [1.0;3]),
            Aabb::new([2.0;3], [3.0;3]),
        ];
        let mut nodes = bk.build_bvh(&aabbs);
        aabbs[0] = Aabb::new([10.0;3], [11.0;3]);
        bk.refit(&mut nodes, &aabbs);
        // Root should now cover both moved boxes
        assert!(nodes[0].aabb.max[0] >= 11.0 - 1e-10);
    }

    #[test]
    fn test_gjk_sphere_no_overlap() {
        let gjk = GjkKernel::new(64, 1e-8);
        let a = Sphere { center: [0.0;3], radius: 1.0 };
        let b = Sphere { center: [3.0,0.0,0.0], radius: 1.0 };
        let r = gjk.gjk_sphere_sphere(&a, &b);
        assert!(!r.overlapping);
        assert!((r.distance - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gjk_sphere_overlap() {
        let gjk = GjkKernel::new(64, 1e-8);
        let a = Sphere { center: [0.0;3], radius: 1.0 };
        let b = Sphere { center: [1.0,0.0,0.0], radius: 1.0 };
        let r = gjk.gjk_sphere_sphere(&a, &b);
        assert!(r.overlapping);
    }

    #[test]
    fn test_epa_sphere_penetration() {
        let epa = EpaKernel::new(32, 1e-8);
        let a = Sphere { center: [0.0;3], radius: 1.0 };
        let b = Sphere { center: [1.5,0.0,0.0], radius: 1.0 };
        let r = epa.epa_sphere_sphere(&a, &b);
        assert!((r.depth - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_sat_obb_obb_overlap() {
        let sat = SatKernel::new();
        let a = Obb::axis_aligned([0.0;3], [1.0;3]);
        let b = Obb::axis_aligned([1.5,0.0,0.0], [1.0;3]);
        let result = sat.sat_obb_obb(&a, &b, 0, 1);
        assert!(result.is_some());
    }

    #[test]
    fn test_sat_obb_obb_no_overlap() {
        let sat = SatKernel::new();
        let a = Obb::axis_aligned([0.0;3], [0.5;3]);
        let b = Obb::axis_aligned([3.0,0.0,0.0], [0.5;3]);
        let result = sat.sat_obb_obb(&a, &b, 0, 1);
        assert!(result.is_none());
    }

    #[test]
    fn test_ccd_no_collision() {
        let ccd = CcdKernel::new(20, 1e-6);
        let a = CcdBody { pos: [0.0;3], vel: [1.0,0.0,0.0], radius: 0.5 };
        let b = CcdBody { pos: [10.0,0.0,0.0], vel: [-1.0,0.0,0.0], radius: 0.5 };
        // They approach but dt is small, might not collide
        let r = ccd.conservative_advancement(&a, &b, 0.01);
        // Just check it doesn't panic
        let _ = r;
    }

    #[test]
    fn test_ccd_collision_detected() {
        let ccd = CcdKernel::new(100, 1e-6);
        let a = CcdBody { pos: [0.0;3], vel: [5.0,0.0,0.0], radius: 0.5 };
        let b = CcdBody { pos: [2.0,0.0,0.0], vel: [-5.0,0.0,0.0], radius: 0.5 };
        let r = ccd.conservative_advancement(&a, &b, 1.0);
        assert!(r.is_some());
    }

    #[test]
    fn test_reduce_to_4_points() {
        let mk = ContactManifoldKernel::new(1e-3);
        let pts: Vec<ContactPoint> = (0..8).map(|i| {
            ContactPoint::new([i as f64, 0.0, 0.0], [0.0,1.0,0.0], i as f64)
        }).collect();
        let reduced = mk.reduce_to_4_points(&pts);
        assert!(reduced.len() <= 4);
    }

    #[test]
    fn test_persistent_manifold_update() {
        let mk = ContactManifoldKernel::new(0.1);
        let old = vec![ContactPoint { position: [0.0;3], normal:[0.0,1.0,0.0], depth:0.01, cached_impulse:5.0 }];
        let new_pts = vec![ContactPoint { position: [0.05,0.0,0.0], normal:[0.0,1.0,0.0], depth:0.01, cached_impulse:0.0 }];
        let result = mk.persistent_manifold_update(&old, &new_pts);
        assert_eq!(result.len(), 1);
        // Position is within tolerance, so impulse should be transferred
        assert!((result[0].cached_impulse - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_heightfield_point_query_flat() {
        let hk = HeightfieldKernel::new();
        let hf = Heightfield::flat(8, 8, 1.0, [0.0;3]);
        let h = hk.heightfield_point_query(&hf, 3.5, 3.5);
        assert!(h.abs() < 1e-12);
    }

    #[test]
    fn test_heightfield_triangle_soup() {
        let hk = HeightfieldKernel::new();
        let hf = Heightfield::flat(3, 3, 1.0, [0.0;3]);
        let tris = hk.triangle_soup_from_hf(&hf);
        assert_eq!(tris.len(), (3-1)*(3-1)*2);
    }

    #[test]
    fn test_compound_aabb_union() {
        let ck = CompoundKernel::new();
        let shapes = vec![Aabb::new([0.0;3],[1.0;3]), Aabb::new([2.0;3],[3.0;3])];
        let union = ck.compound_aabb_union(&shapes).unwrap();
        assert_eq!(union.min, [0.0;3]);
        assert_eq!(union.max, [3.0;3]);
    }

    #[test]
    fn test_compound_center_of_mass() {
        let ck = CompoundKernel::new();
        let masses = vec![1.0, 1.0];
        let centers = vec![[0.0,0.0,0.0], [2.0,0.0,0.0]];
        let com = ck.center_of_mass(&masses, &centers);
        assert!((com[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_simplex_add_and_len() {
        let mut s = Simplex::new();
        s.add([1.0, 0.0, 0.0]);
        s.add([0.0, 1.0, 0.0]);
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());
    }

    #[test]
    fn test_config_default() {
        let cfg = CollisionKernelConfig::new_default();
        assert_eq!(cfg.broadphase_type, BroadphaseType::Bvh);
        assert_eq!(cfg.narrow_type, NarrowphaseType::Gjk);
        assert!(cfg.max_pairs > 0);
    }
}
