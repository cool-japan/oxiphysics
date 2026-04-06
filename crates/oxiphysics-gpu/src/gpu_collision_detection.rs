// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated collision detection (CPU mock backend via Rayon).
//!
//! Implements a broadphase AABB sweep-and-prune followed by a narrowphase
//! sphere-sphere test and a stub GJK dispatcher.  The "GPU" dispatch is mocked
//! using Rayon parallel iterators so the module compiles and runs on CPU.

use rayon::prelude::*;

// ---------------------------------------------------------------------------
// Aabb
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum corner (x, y, z).
    pub min: [f64; 3],
    /// Maximum corner (x, y, z).
    pub max: [f64; 3],
}

impl Aabb {
    /// Create a new AABB from min/max corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Create an AABB centred at `centre` with half-extents `half`.
    pub fn from_center_half(centre: [f64; 3], half: [f64; 3]) -> Self {
        Self {
            min: [
                centre[0] - half[0],
                centre[1] - half[1],
                centre[2] - half[2],
            ],
            max: [
                centre[0] + half[0],
                centre[1] + half[1],
                centre[2] + half[2],
            ],
        }
    }

    /// Centre of the AABB.
    pub fn centre(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }

    /// Returns `true` if this AABB overlaps `other`.
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
}

// ---------------------------------------------------------------------------
// GpuCollisionBuffer
// ---------------------------------------------------------------------------

/// Buffer holding the AABB list and broadphase potential pairs.
#[derive(Debug, Clone)]
pub struct GpuCollisionBuffer {
    /// Axis-aligned bounding boxes for each object.
    pub aabbs: Vec<Aabb>,
    /// Potential collision pairs `(i, j)` with `i < j` from the broadphase.
    pub potential_pairs: Vec<(usize, usize)>,
}

impl GpuCollisionBuffer {
    /// Create an empty collision buffer.
    pub fn new() -> Self {
        Self {
            aabbs: Vec::new(),
            potential_pairs: Vec::new(),
        }
    }

    /// Add an AABB to the buffer, returning its index.
    pub fn add_aabb(&mut self, aabb: Aabb) -> usize {
        let idx = self.aabbs.len();
        self.aabbs.push(aabb);
        idx
    }

    /// Number of objects in the buffer.
    pub fn n_objects(&self) -> usize {
        self.aabbs.len()
    }
}

impl Default for GpuCollisionBuffer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// gpu_broadphase_sort
// ---------------------------------------------------------------------------

/// Sort AABBs by their minimum x-coordinate (mock GPU sort).
///
/// Returns a permutation index array such that `perm[0]` is the index of the
/// AABB with the smallest `min.x`.  The buffer itself is not reordered.
///
/// # Arguments
/// * `buffer` - The collision buffer to sort.
pub fn gpu_broadphase_sort(buffer: &GpuCollisionBuffer) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..buffer.aabbs.len()).collect();
    indices.sort_unstable_by(|&a, &b| {
        buffer.aabbs[a].min[0]
            .partial_cmp(&buffer.aabbs[b].min[0])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
}

// ---------------------------------------------------------------------------
// gpu_aabb_overlap
// ---------------------------------------------------------------------------

/// AABB-AABB overlap test kernel (mock GPU dispatch).
///
/// Tests all pairs `(i, j)` with `i < j` in parallel and returns those that
/// overlap.
///
/// # Arguments
/// * `buffer` - Collision buffer with the AABB list.
pub fn gpu_aabb_overlap(buffer: &GpuCollisionBuffer) -> Vec<(usize, usize)> {
    let n = buffer.aabbs.len();
    if n < 2 {
        return Vec::new();
    }

    // Generate all pairs and test in parallel
    let pairs: Vec<(usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let aabb_i = &buffer.aabbs[i];
            (i + 1..n)
                .filter(|&j| aabb_i.overlaps(&buffer.aabbs[j]))
                .map(|j| (i, j))
                .collect::<Vec<_>>()
        })
        .collect();
    pairs
}

// ---------------------------------------------------------------------------
// gpu_sphere_collision
// ---------------------------------------------------------------------------

/// Sphere-sphere narrowphase collision result.
#[derive(Debug, Clone, Copy)]
pub struct SphereContact {
    /// Index of the first sphere.
    pub i: usize,
    /// Index of the second sphere.
    pub j: usize,
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// Contact normal pointing from sphere `i` to sphere `j`.
    pub normal: [f64; 3],
}

/// Sphere descriptor: centre + radius.
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    /// Centre position (x, y, z).
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
}

impl Sphere {
    /// Create a sphere at `centre` with the given `radius`.
    pub fn new(centre: [f64; 3], radius: f64) -> Self {
        Self { centre, radius }
    }
}

/// Test sphere-sphere collisions for the given pairs (mock GPU narrowphase).
///
/// For each pair `(i, j)` in `pairs`, tests whether spheres `i` and `j`
/// overlap and, if so, computes a [`SphereContact`].
///
/// # Arguments
/// * `spheres` - Sphere descriptors, indexed by object id.
/// * `pairs` - Candidate pairs from the broadphase.
pub fn gpu_sphere_collision(spheres: &[Sphere], pairs: &[(usize, usize)]) -> Vec<SphereContact> {
    pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            if i >= spheres.len() || j >= spheres.len() {
                return None;
            }
            let a = &spheres[i];
            let b = &spheres[j];
            let dx = b.centre[0] - a.centre[0];
            let dy = b.centre[1] - a.centre[1];
            let dz = b.centre[2] - a.centre[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            let r_sum = a.radius + b.radius;
            if dist_sq >= r_sum * r_sum {
                return None;
            }
            let dist = dist_sq.sqrt().max(1e-15);
            let depth = r_sum - dist;
            let normal = [dx / dist, dy / dist, dz / dist];
            Some(SphereContact {
                i,
                j,
                depth,
                normal,
            })
        })
        .collect()
}

// ---------------------------------------------------------------------------
// build_collision_pairs
// ---------------------------------------------------------------------------

/// Build the full broadphase collision pair list from a sorted AABB list.
///
/// Uses a sweep-and-prune approach along the x-axis: after sorting by `min.x`,
/// only pairs where the x-extents overlap need to be checked further.
///
/// Updates `buffer.potential_pairs` in place and returns the count of pairs.
///
/// # Arguments
/// * `buffer` - Collision buffer (sorted order applied internally).
pub fn build_collision_pairs(buffer: &mut GpuCollisionBuffer) -> usize {
    let sorted = gpu_broadphase_sort(buffer);
    let mut pairs = Vec::new();

    for (si, &i) in sorted.iter().enumerate() {
        let aabb_i = &buffer.aabbs[i];
        // Only check objects whose min.x <= max.x of i (sweep prune)
        for &j in sorted[si + 1..].iter() {
            let aabb_j = &buffer.aabbs[j];
            if aabb_j.min[0] > aabb_i.max[0] {
                break;
            }
            if aabb_i.overlaps(aabb_j) {
                let pair = if i < j { (i, j) } else { (j, i) };
                pairs.push(pair);
            }
        }
    }

    pairs.sort_unstable();
    pairs.dedup();
    buffer.potential_pairs = pairs;
    buffer.potential_pairs.len()
}

// ---------------------------------------------------------------------------
// GjkResult
// ---------------------------------------------------------------------------

/// Stub result from the GJK narrowphase.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GjkResult {
    /// The two shapes are separated; `distance` is the gap.
    Separated {
        /// Minimum distance between the two shapes.
        distance: f64,
    },
    /// The two shapes intersect.
    Intersecting,
}

/// Dispatch GJK for a batch of collision pairs (stub implementation).
///
/// In a real GPU pipeline, this would launch one thread per pair.  Here we
/// use a heuristic based on AABB centre distance as a placeholder.
///
/// # Arguments
/// * `buffer` - Collision buffer with the AABB list and pair list.
pub fn batch_gjk_dispatch(buffer: &GpuCollisionBuffer) -> Vec<GjkResult> {
    buffer
        .potential_pairs
        .par_iter()
        .map(|&(i, j)| {
            let ci = buffer.aabbs[i].centre();
            let cj = buffer.aabbs[j].centre();
            let dx = cj[0] - ci[0];
            let dy = cj[1] - ci[1];
            let dz = cj[2] - ci[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            // Approximate "radius" as half the diagonal of the AABB
            let diag_i = {
                let d = [
                    buffer.aabbs[i].max[0] - buffer.aabbs[i].min[0],
                    buffer.aabbs[i].max[1] - buffer.aabbs[i].min[1],
                    buffer.aabbs[i].max[2] - buffer.aabbs[i].min[2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() * 0.5
            };
            let diag_j = {
                let d = [
                    buffer.aabbs[j].max[0] - buffer.aabbs[j].min[0],
                    buffer.aabbs[j].max[1] - buffer.aabbs[j].min[1],
                    buffer.aabbs[j].max[2] - buffer.aabbs[j].min[2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() * 0.5
            };
            if dist < diag_i + diag_j {
                GjkResult::Intersecting
            } else {
                GjkResult::Separated {
                    distance: dist - diag_i - diag_j,
                }
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Aabb ──────────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_new() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(a.min, [0.0; 3]);
        assert_eq!(a.max, [1.0; 3]);
    }

    #[test]
    fn test_aabb_from_center_half() {
        let a = Aabb::from_center_half([1.0, 2.0, 3.0], [0.5, 0.5, 0.5]);
        assert!((a.min[0] - 0.5).abs() < 1e-12);
        assert!((a.max[0] - 1.5).abs() < 1e-12);
        assert!((a.min[1] - 1.5).abs() < 1e-12);
        assert!((a.max[1] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_centre() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = a.centre();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_aabb_overlaps_true() {
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([0.5; 3], [1.5; 3]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_overlaps_false_separated_x() {
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }

    #[test]
    fn test_aabb_overlaps_touching_edge() {
        // Touching at x=1.0 counts as overlap (<=)
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
        assert!(a.overlaps(&b));
    }

    #[test]
    fn test_aabb_no_overlap_y() {
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([0.0, 2.0, 0.0], [1.0, 3.0, 1.0]);
        assert!(!a.overlaps(&b));
    }

    // ── GpuCollisionBuffer ────────────────────────────────────────────────

    #[test]
    fn test_buffer_new_empty() {
        let buf = GpuCollisionBuffer::new();
        assert_eq!(buf.n_objects(), 0);
        assert!(buf.potential_pairs.is_empty());
    }

    #[test]
    fn test_buffer_add_aabb() {
        let mut buf = GpuCollisionBuffer::new();
        let idx = buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        assert_eq!(idx, 0);
        assert_eq!(buf.n_objects(), 1);
    }

    #[test]
    fn test_buffer_default() {
        let buf = GpuCollisionBuffer::default();
        assert_eq!(buf.n_objects(), 0);
    }

    // ── gpu_broadphase_sort ───────────────────────────────────────────────

    #[test]
    fn test_broadphase_sort_order() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]));
        buf.add_aabb(Aabb::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]));
        buf.add_aabb(Aabb::new([0.0, 0.0, 0.0], [0.5, 1.0, 1.0]));
        let perm = gpu_broadphase_sort(&buf);
        assert_eq!(perm, vec![2, 1, 0]);
    }

    #[test]
    fn test_broadphase_sort_empty() {
        let buf = GpuCollisionBuffer::new();
        let perm = gpu_broadphase_sort(&buf);
        assert!(perm.is_empty());
    }

    // ── gpu_aabb_overlap ──────────────────────────────────────────────────

    #[test]
    fn test_aabb_overlap_two_overlapping() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        buf.add_aabb(Aabb::new([0.5; 3], [1.5; 3]));
        let pairs = gpu_aabb_overlap(&buf);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 1));
    }

    #[test]
    fn test_aabb_overlap_two_separated() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        buf.add_aabb(Aabb::new([2.0; 3], [3.0; 3]));
        let pairs = gpu_aabb_overlap(&buf);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_aabb_overlap_three_objects() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [2.0; 3]));
        buf.add_aabb(Aabb::new([1.0; 3], [3.0; 3]));
        buf.add_aabb(Aabb::new([5.0; 3], [6.0; 3]));
        let pairs = gpu_aabb_overlap(&buf);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 1));
    }

    #[test]
    fn test_aabb_overlap_single_object_no_pairs() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        let pairs = gpu_aabb_overlap(&buf);
        assert!(pairs.is_empty());
    }

    // ── gpu_sphere_collision ──────────────────────────────────────────────

    #[test]
    fn test_sphere_collision_overlapping() {
        let spheres = vec![
            Sphere::new([0.0, 0.0, 0.0], 1.0),
            Sphere::new([1.5, 0.0, 0.0], 1.0),
        ];
        let pairs = vec![(0usize, 1usize)];
        let contacts = gpu_sphere_collision(&spheres, &pairs);
        assert_eq!(contacts.len(), 1);
        assert!((contacts[0].depth - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_collision_separated() {
        let spheres = vec![
            Sphere::new([0.0, 0.0, 0.0], 1.0),
            Sphere::new([5.0, 0.0, 0.0], 1.0),
        ];
        let pairs = vec![(0usize, 1usize)];
        let contacts = gpu_sphere_collision(&spheres, &pairs);
        assert!(contacts.is_empty());
    }

    #[test]
    fn test_sphere_collision_normal_direction() {
        let spheres = vec![
            Sphere::new([0.0, 0.0, 0.0], 1.0),
            Sphere::new([1.0, 0.0, 0.0], 1.0),
        ];
        let pairs = vec![(0usize, 1usize)];
        let contacts = gpu_sphere_collision(&spheres, &pairs);
        assert_eq!(contacts.len(), 1);
        assert!((contacts[0].normal[0] - 1.0).abs() < 1e-10);
        assert!(contacts[0].normal[1].abs() < 1e-10);
        assert!(contacts[0].normal[2].abs() < 1e-10);
    }

    #[test]
    fn test_sphere_collision_out_of_bounds_index() {
        let spheres = vec![Sphere::new([0.0; 3], 1.0)];
        let pairs = vec![(0usize, 99usize)];
        let contacts = gpu_sphere_collision(&spheres, &pairs);
        assert!(contacts.is_empty());
    }

    // ── build_collision_pairs ─────────────────────────────────────────────

    #[test]
    fn test_build_pairs_two_overlapping() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        buf.add_aabb(Aabb::new([0.5; 3], [1.5; 3]));
        let n = build_collision_pairs(&mut buf);
        assert_eq!(n, 1);
        assert_eq!(buf.potential_pairs[0], (0, 1));
    }

    #[test]
    fn test_build_pairs_no_overlap() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        buf.add_aabb(Aabb::new([10.0; 3], [11.0; 3]));
        let n = build_collision_pairs(&mut buf);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_build_pairs_multiple_objects() {
        let mut buf = GpuCollisionBuffer::new();
        // All three overlap with each other
        buf.add_aabb(Aabb::new([0.0; 3], [2.0; 3]));
        buf.add_aabb(Aabb::new([1.0; 3], [3.0; 3]));
        buf.add_aabb(Aabb::new([1.5; 3], [3.5; 3]));
        let n = build_collision_pairs(&mut buf);
        assert!(n >= 2); // At minimum (0,1) and (1,2) should be found
    }

    // ── batch_gjk_dispatch ─────────────────────────────────────────────────

    #[test]
    fn test_gjk_dispatch_intersecting_pair() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [1.0; 3]));
        buf.add_aabb(Aabb::new([0.5; 3], [1.5; 3]));
        build_collision_pairs(&mut buf);
        let results = batch_gjk_dispatch(&buf);
        assert!(!results.is_empty());
        assert_eq!(results[0], GjkResult::Intersecting);
    }

    #[test]
    fn test_gjk_dispatch_empty_pairs() {
        let buf = GpuCollisionBuffer::new();
        let results = batch_gjk_dispatch(&buf);
        assert!(results.is_empty());
    }

    #[test]
    fn test_gjk_dispatch_separated_pair() {
        let mut buf = GpuCollisionBuffer::new();
        buf.add_aabb(Aabb::new([0.0; 3], [0.1; 3]));
        buf.add_aabb(Aabb::new([10.0; 3], [10.1; 3]));
        // Manually add a pair (broadphase won't find this since they don't overlap)
        buf.potential_pairs = vec![(0, 1)];
        let results = batch_gjk_dispatch(&buf);
        assert_eq!(results.len(), 1);
        matches!(results[0], GjkResult::Separated { .. });
    }
}
