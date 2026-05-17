// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core BVH types: AABB, primitives, nodes, flat nodes, and ray-hit results.

// ============================================================================
// Aabb
// ============================================================================

/// An axis-aligned bounding box stored as two `[f32; 3]` corners.
#[derive(Debug, Clone, PartialEq)]
pub struct Aabb {
    /// Minimum corner (component-wise).
    pub min: [f32; 3],
    /// Maximum corner (component-wise).
    pub max: [f32; 3],
}

impl Aabb {
    /// Construct an `Aabb` from explicit min/max corners.
    pub fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self { min, max }
    }

    /// Construct a degenerate `Aabb` that contains exactly one point.
    pub fn point(p: [f32; 3]) -> Self {
        Self { min: p, max: p }
    }

    /// Return the smallest `Aabb` that contains both `a` and `b`.
    pub fn merge(a: &Aabb, b: &Aabb) -> Aabb {
        Aabb {
            min: [
                a.min[0].min(b.min[0]),
                a.min[1].min(b.min[1]),
                a.min[2].min(b.min[2]),
            ],
            max: [
                a.max[0].max(b.max[0]),
                a.max[1].max(b.max[1]),
                a.max[2].max(b.max[2]),
            ],
        }
    }

    /// Returns `true` if this box overlaps `other` (touching counts).
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }

    /// Returns `true` if point `p` lies inside or on the surface of this box.
    pub fn contains(&self, p: [f32; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }

    /// Surface area of the box (sum of all six face areas).
    pub fn surface_area(&self) -> f32 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dz * dx)
    }

    /// Geometric centre of the box.
    pub fn center(&self) -> [f32; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Return a copy of this box expanded uniformly by `margin` on every side.
    pub fn expand(&self, margin: f32) -> Aabb {
        Aabb {
            min: [
                self.min[0] - margin,
                self.min[1] - margin,
                self.min[2] - margin,
            ],
            max: [
                self.max[0] + margin,
                self.max[1] + margin,
                self.max[2] + margin,
            ],
        }
    }
}

// ============================================================================
// BvhPrimitive
// ============================================================================

/// A leaf primitive: an AABB together with the logical object it belongs to.
#[derive(Debug, Clone)]
pub struct BvhPrimitive {
    /// Bounding box of the primitive.
    pub aabb: Aabb,
    /// Caller-defined identifier (returned by queries).
    pub object_id: usize,
}

impl BvhPrimitive {
    /// Construct a `BvhPrimitive`.
    pub fn new(aabb: Aabb, object_id: usize) -> Self {
        Self { aabb, object_id }
    }
}

// ============================================================================
// BvhNode
// ============================================================================

/// A node in the recursive BVH tree.
#[derive(Debug)]
pub struct BvhNode {
    /// Bounding box that contains all children / primitives.
    pub aabb: Aabb,
    /// Left subtree (internal nodes only).
    pub left: Option<Box<BvhNode>>,
    /// Right subtree (internal nodes only).
    pub right: Option<Box<BvhNode>>,
    /// Indices into `Bvh::primitives` (leaf nodes only).
    pub primitives: Vec<usize>,
}

impl BvhNode {
    /// Returns `true` if this node is a leaf (holds primitives directly).
    pub fn is_leaf(&self) -> bool {
        self.left.is_none() && self.right.is_none()
    }
}

// ============================================================================
// Flat BVH node
// ============================================================================

/// A single node in the linearised (flat) BVH representation.
///
/// Layout:
/// * If `count == 0` this is an **internal** node.
///   - The **left** child is always at index `node_idx + 1` (i.e. stored
///     immediately after the parent in DFS pre-order).
///   - `left_first` holds the index of the **right** child.
/// * If `count > 0` this is a **leaf**; `left_first` is the start index into
///   the accompanying primitive-index slice and `count` is the number of
///   entries.
#[derive(Debug, Clone)]
pub struct FlatBvhNode {
    /// Bounding box of this node.
    pub aabb: Aabb,
    /// Right-child index (internal) or first-primitive index (leaf).
    pub left_first: u32,
    /// 0 for internal nodes; number of primitives for leaf nodes.
    pub count: u32,
}

// ============================================================================
// RayHit
// ============================================================================

/// Result of a closest-hit ray traversal.
#[derive(Debug, Clone)]
pub struct RayHit {
    /// Object ID of the closest hit primitive.
    pub object_id: usize,
    /// Ray parameter at which the hit occurred.
    pub t: f32,
}

// ============================================================================
// GpuRay
// ============================================================================

/// A ray for GPU BVH traversal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuRay {
    /// Ray origin.
    pub origin: [f32; 3],
    /// Ray direction (does not need to be normalised).
    pub direction: [f32; 3],
    /// Maximum ray parameter (hit farther than this is ignored).
    pub max_t: f32,
}

impl GpuRay {
    /// Construct a new ray.
    pub fn new(origin: [f32; 3], direction: [f32; 3], max_t: f32) -> Self {
        Self {
            origin,
            direction,
            max_t,
        }
    }
}

// ============================================================================
// BVH statistics types
// ============================================================================

/// Runtime statistics about a BVH tree.
#[derive(Debug, Clone)]
pub struct BvhStats {
    /// Total number of nodes (internal + leaf).
    pub node_count: usize,
    /// Number of leaf nodes.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
    /// Maximum tree depth.
    pub max_depth: usize,
    /// Total number of primitives stored across all leaves.
    pub total_primitives: usize,
    /// Average primitives per leaf.
    pub avg_primitives_per_leaf: f32,
}

/// Extended BVH tree statistics including average fan-out.
#[derive(Debug, Clone)]
pub struct BvhTreeStatistics {
    /// Total node count.
    pub node_count: usize,
    /// Number of leaf nodes.
    pub leaf_count: usize,
    /// Number of internal nodes.
    pub internal_count: usize,
    /// Maximum depth from root (1-indexed).
    pub max_depth: usize,
    /// Total number of primitives across all leaves.
    pub total_primitives: usize,
    /// Average number of children per internal node (fan-out).
    /// For a binary tree this is at most 2.
    pub avg_fanout: f32,
    /// Total surface area of all leaf AABBs.
    pub total_leaf_surface_area: f32,
}

// ============================================================================
// Morton / LBVH types
// ============================================================================

/// LBVH primitive: AABB + Morton code.
#[derive(Debug, Clone)]
pub struct LbvhPrimitive {
    /// Bounding box.
    pub aabb: Aabb,
    /// Caller-defined object ID.
    pub object_id: usize,
    /// 30-bit Morton code computed from the AABB centroid.
    pub morton: u32,
}

/// A cluster of Morton-coded primitives, with a pre-computed bounding radius.
#[derive(Debug, Clone)]
pub struct MortonCluster {
    /// Indices of the primitives in this cluster (into the parent slice).
    pub indices: Vec<usize>,
    /// Axis-aligned bounding box of the cluster.
    pub aabb: Aabb,
    /// Bounding sphere radius (centred at the AABB centre).
    pub radius: f32,
}
