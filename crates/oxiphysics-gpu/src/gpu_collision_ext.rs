// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU collision detection extensions (CPU mock backend).
//!
//! Provides additional broadphase and narrowphase utilities that complement
//! [`gpu_collision_detection`](super::gpu_collision_detection):
//! - [`GpuBroadphaseGrid`]: uniform-grid broadphase
//! - [`GpuContactManifold`]: contact point storage
//! - Free functions: AABB overlap, sphere overlap, point containment,
//!   ray-AABB / ray-sphere intersection, GJK distance estimate.

// ── GpuBroadphaseGrid ────────────────────────────────────────────────────────

/// A uniform-grid broadphase accelerator.
///
/// The 3-D domain is partitioned into `nx × ny × nz` axis-aligned cells.
/// Particles are inserted by position and can be queried by neighbourhood.
#[derive(Debug, Clone)]
pub struct GpuBroadphaseGrid {
    /// Size of one cell along each axis.
    pub cell_size: f64,
    /// Number of cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z.
    pub nz: usize,
    /// Flattened cell storage: `cells[flat_index]` is the list of particle IDs.
    pub cells: Vec<Vec<usize>>,
}

impl GpuBroadphaseGrid {
    /// Create a new broadphase grid.
    ///
    /// # Arguments
    /// * `cell_size` – edge length of each cubic cell.
    /// * `nx`, `ny`, `nz` – cell counts along each axis.
    pub fn new(cell_size: f64, nx: usize, ny: usize, nz: usize) -> Self {
        let total = nx * ny * nz;
        Self {
            cell_size,
            nx,
            ny,
            nz,
            cells: vec![Vec::new(); total],
        }
    }

    /// Insert `particle_id` at world-space position `pos`.
    ///
    /// Positions that fall outside the grid are clamped to the nearest cell.
    pub fn insert(&mut self, particle_id: usize, pos: [f64; 3]) {
        let (ix, iy, iz) = self.cell_index(pos);
        let flat = iz * self.nx * self.ny + iy * self.nx + ix;
        if flat < self.cells.len() {
            self.cells[flat].push(particle_id);
        }
    }

    /// Map a world-space position to grid-cell indices `(ix, iy, iz)`.
    ///
    /// Result is clamped so it always lies within `[0, nx-1] × [0, ny-1] × [0, nz-1]`.
    pub fn cell_index(&self, pos: [f64; 3]) -> (usize, usize, usize) {
        let clamp = |v: f64, n: usize| {
            let i = (v / self.cell_size).floor() as isize;
            i.max(0).min(n as isize - 1) as usize
        };
        let ix = clamp(pos[0], self.nx.max(1));
        let iy = clamp(pos[1], self.ny.max(1));
        let iz = clamp(pos[2], self.nz.max(1));
        (ix, iy, iz)
    }

    /// Return all particle IDs in the cell containing `pos` and its 26 neighbours.
    pub fn query_neighbors(&self, pos: [f64; 3]) -> Vec<usize> {
        let (cx, cy, cz) = self.cell_index(pos);
        let mut result = Vec::new();
        for dz in -1isize..=1 {
            for dy in -1isize..=1 {
                for dx in -1isize..=1 {
                    let ix = cx as isize + dx;
                    let iy = cy as isize + dy;
                    let iz = cz as isize + dz;
                    if ix < 0
                        || iy < 0
                        || iz < 0
                        || ix >= self.nx as isize
                        || iy >= self.ny as isize
                        || iz >= self.nz as isize
                    {
                        continue;
                    }
                    let flat =
                        iz as usize * self.nx * self.ny + iy as usize * self.nx + ix as usize;
                    result.extend_from_slice(&self.cells[flat]);
                }
            }
        }
        result
    }

    /// Total number of cells in the grid.
    pub fn cell_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Remove all particle IDs from every cell.
    pub fn clear(&mut self) {
        for cell in &mut self.cells {
            cell.clear();
        }
    }
}

// ── Overlap tests ─────────────────────────────────────────────────────────────

/// Test whether two axis-aligned bounding boxes overlap.
///
/// Returns `true` when the boxes touch or interpenetrate on all three axes.
pub fn gpu_aabb_overlap(
    a_min: [f32; 3],
    a_max: [f32; 3],
    b_min: [f32; 3],
    b_max: [f32; 3],
) -> bool {
    a_min[0] <= b_max[0]
        && a_max[0] >= b_min[0]
        && a_min[1] <= b_max[1]
        && a_max[1] >= b_min[1]
        && a_min[2] <= b_max[2]
        && a_max[2] >= b_min[2]
}

/// Test whether two spheres overlap (or touch).
///
/// Returns `true` when the distance between centres is ≤ `r1 + r2`.
pub fn gpu_sphere_sphere_overlap(c1: [f32; 3], r1: f32, c2: [f32; 3], r2: f32) -> bool {
    let dx = c1[0] - c2[0];
    let dy = c1[1] - c2[1];
    let dz = c1[2] - c2[2];
    let dist2 = dx * dx + dy * dy + dz * dz;
    let rsum = r1 + r2;
    dist2 <= rsum * rsum
}

/// Test whether point `p` is inside (or on the boundary of) the AABB.
pub fn gpu_point_in_aabb(p: [f32; 3], min: [f32; 3], max: [f32; 3]) -> bool {
    p[0] >= min[0]
        && p[0] <= max[0]
        && p[1] >= min[1]
        && p[1] <= max[1]
        && p[2] >= min[2]
        && p[2] <= max[2]
}

// ── Ray intersection ──────────────────────────────────────────────────────────

/// Intersect a ray with an AABB using the slab method.
///
/// Returns `Some(t)` (the entry distance along the ray) if the ray hits the
/// box, or `None` on a miss.  The ray direction need not be normalised.
pub fn gpu_ray_aabb_intersect(
    ray_origin: [f32; 3],
    ray_dir: [f32; 3],
    min: [f32; 3],
    max: [f32; 3],
) -> Option<f32> {
    let mut t_min = f32::NEG_INFINITY;
    let mut t_max = f32::INFINITY;
    for i in 0..3 {
        let inv_d = if ray_dir[i].abs() > 1e-12 {
            1.0 / ray_dir[i]
        } else {
            // Ray parallel to slab – check if origin is inside
            if ray_origin[i] < min[i] || ray_origin[i] > max[i] {
                return None;
            }
            continue;
        };
        let t1 = (min[i] - ray_origin[i]) * inv_d;
        let t2 = (max[i] - ray_origin[i]) * inv_d;
        let (ta, tb) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
        t_min = t_min.max(ta);
        t_max = t_max.min(tb);
        if t_min > t_max {
            return None;
        }
    }
    if t_max < 0.0 {
        return None;
    }
    Some(if t_min >= 0.0 { t_min } else { t_max })
}

/// Intersect a ray with a sphere.
///
/// Returns `Some(t)` for the nearest positive intersection, or `None` on miss.
pub fn gpu_ray_sphere_intersect(
    ray_origin: [f32; 3],
    ray_dir: [f32; 3],
    center: [f32; 3],
    radius: f32,
) -> Option<f32> {
    let ocx = ray_origin[0] - center[0];
    let ocy = ray_origin[1] - center[1];
    let ocz = ray_origin[2] - center[2];
    let a = ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1] + ray_dir[2] * ray_dir[2];
    let b = 2.0 * (ray_dir[0] * ocx + ray_dir[1] * ocy + ray_dir[2] * ocz);
    let c = ocx * ocx + ocy * ocy + ocz * ocz - radius * radius;
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return None;
    }
    let sqrt_d = discriminant.sqrt();
    let t1 = (-b - sqrt_d) / (2.0 * a);
    let t2 = (-b + sqrt_d) / (2.0 * a);
    if t1 >= 0.0 {
        Some(t1)
    } else if t2 >= 0.0 {
        Some(t2)
    } else {
        None
    }
}

// ── GpuContactManifold ────────────────────────────────────────────────────────

/// Contact manifold between two rigid bodies.
///
/// Stores the list of contact points, surface normals, and penetration depths
/// generated during narrowphase collision resolution.
#[derive(Debug, Clone)]
pub struct GpuContactManifold {
    /// Index of the first body.
    pub body_a: usize,
    /// Index of the second body.
    pub body_b: usize,
    /// World-space contact points.
    pub contact_points: Vec<[f32; 3]>,
    /// Outward contact normals (from A towards B).
    pub normals: Vec<[f32; 3]>,
    /// Penetration depths (positive = overlapping).
    pub depths: Vec<f32>,
}

impl GpuContactManifold {
    /// Create an empty manifold for the given body pair.
    pub fn new(body_a: usize, body_b: usize) -> Self {
        Self {
            body_a,
            body_b,
            contact_points: Vec::new(),
            normals: Vec::new(),
            depths: Vec::new(),
        }
    }

    /// Add a single contact to the manifold.
    ///
    /// * `point`  – world-space contact point.
    /// * `normal` – surface normal pointing from A to B.
    /// * `depth`  – penetration depth (positive value).
    pub fn add_contact(&mut self, point: [f32; 3], normal: [f32; 3], depth: f32) {
        self.contact_points.push(point);
        self.normals.push(normal);
        self.depths.push(depth);
    }

    /// Number of contact points in this manifold.
    pub fn contact_count(&self) -> usize {
        self.contact_points.len()
    }

    /// Maximum penetration depth across all contacts, or `0.0` if empty.
    pub fn max_penetration(&self) -> f32 {
        self.depths.iter().copied().fold(0.0_f32, f32::max)
    }
}

// ── GJK distance (simplified) ────────────────────────────────────────────────

/// Estimate the minimum distance between two convex hulls using a simplified
/// GJK approach (brute-force vertex-pair minimum for the mock GPU backend).
///
/// Returns `0.0` if either set is empty.
pub fn gpu_gjk_distance(vertices_a: &[[f32; 3]], vertices_b: &[[f32; 3]]) -> f32 {
    if vertices_a.is_empty() || vertices_b.is_empty() {
        return 0.0;
    }
    let mut min_dist2 = f32::MAX;
    for &va in vertices_a {
        for &vb in vertices_b {
            let dx = va[0] - vb[0];
            let dy = va[1] - vb[1];
            let dz = va[2] - vb[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < min_dist2 {
                min_dist2 = d2;
            }
        }
    }
    min_dist2.sqrt()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── GpuBroadphaseGrid ────────────────────────────────────────────────

    #[test]
    fn broadphase_cell_count_basic() {
        let g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        assert_eq!(g.cell_count(), 64);
    }

    #[test]
    fn broadphase_cell_count_non_uniform() {
        let g = GpuBroadphaseGrid::new(0.5, 2, 3, 5);
        assert_eq!(g.cell_count(), 30);
    }

    #[test]
    fn broadphase_insert_and_query() {
        let mut g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        g.insert(0, [0.5, 0.5, 0.5]);
        let neighbors = g.query_neighbors([0.5, 0.5, 0.5]);
        assert!(!neighbors.is_empty());
        assert!(neighbors.contains(&0));
    }

    #[test]
    fn broadphase_query_nearby_cell() {
        let mut g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        g.insert(42, [1.5, 1.5, 1.5]);
        // Query from adjacent cell should still find it
        let neighbors = g.query_neighbors([0.5, 0.5, 0.5]);
        assert!(neighbors.contains(&42));
    }

    #[test]
    fn broadphase_clear_removes_particles() {
        let mut g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        g.insert(0, [0.5, 0.5, 0.5]);
        g.clear();
        let neighbors = g.query_neighbors([0.5, 0.5, 0.5]);
        assert!(neighbors.is_empty());
    }

    #[test]
    fn broadphase_multiple_particles_same_cell() {
        let mut g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        g.insert(1, [0.1, 0.1, 0.1]);
        g.insert(2, [0.9, 0.9, 0.9]);
        let neighbors = g.query_neighbors([0.5, 0.5, 0.5]);
        assert!(neighbors.contains(&1));
        assert!(neighbors.contains(&2));
    }

    #[test]
    fn broadphase_cell_index_origin() {
        let g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        let (ix, iy, iz) = g.cell_index([0.0, 0.0, 0.0]);
        assert_eq!((ix, iy, iz), (0, 0, 0));
    }

    #[test]
    fn broadphase_cell_index_clamped() {
        let g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        // position outside grid is clamped
        let (ix, _iy, _iz) = g.cell_index([100.0, 0.0, 0.0]);
        assert_eq!(ix, 3);
    }

    #[test]
    fn broadphase_query_empty_grid() {
        let g = GpuBroadphaseGrid::new(1.0, 4, 4, 4);
        let neighbors = g.query_neighbors([0.5, 0.5, 0.5]);
        assert!(neighbors.is_empty());
    }

    #[test]
    fn broadphase_insert_many() {
        let mut g = GpuBroadphaseGrid::new(1.0, 8, 8, 8);
        for i in 0..20 {
            g.insert(i, [(i % 8) as f64 * 0.9, 0.5, 0.5]);
        }
        // Grid should still have cell_count = 512
        assert_eq!(g.cell_count(), 512);
    }

    // ── gpu_aabb_overlap ─────────────────────────────────────────────────

    #[test]
    fn aabb_overlap_basic() {
        let a_min = [0.0f32; 3];
        let a_max = [1.0f32; 3];
        let b_min = [0.5f32; 3];
        let b_max = [1.5f32; 3];
        assert!(gpu_aabb_overlap(a_min, a_max, b_min, b_max));
    }

    #[test]
    fn aabb_overlap_no_overlap() {
        let a_min = [0.0f32; 3];
        let a_max = [1.0f32; 3];
        let b_min = [2.0f32; 3];
        let b_max = [3.0f32; 3];
        assert!(!gpu_aabb_overlap(a_min, a_max, b_min, b_max));
    }

    #[test]
    fn aabb_overlap_touching_face() {
        // Boxes touching at x=1
        let a_min = [0.0f32, 0.0, 0.0];
        let a_max = [1.0f32, 1.0, 1.0];
        let b_min = [1.0f32, 0.0, 0.0];
        let b_max = [2.0f32, 1.0, 1.0];
        assert!(gpu_aabb_overlap(a_min, a_max, b_min, b_max));
    }

    #[test]
    fn aabb_overlap_separated_in_y() {
        let a_min = [0.0f32, 0.0, 0.0];
        let a_max = [1.0f32, 1.0, 1.0];
        let b_min = [0.0f32, 2.0, 0.0];
        let b_max = [1.0f32, 3.0, 1.0];
        assert!(!gpu_aabb_overlap(a_min, a_max, b_min, b_max));
    }

    #[test]
    fn aabb_overlap_contained() {
        let outer_min = [0.0f32; 3];
        let outer_max = [10.0f32; 3];
        let inner_min = [2.0f32; 3];
        let inner_max = [3.0f32; 3];
        assert!(gpu_aabb_overlap(outer_min, outer_max, inner_min, inner_max));
    }

    // ── gpu_sphere_sphere_overlap ────────────────────────────────────────

    #[test]
    fn sphere_overlap_clearly_overlapping() {
        let c1 = [0.0f32; 3];
        let c2 = [1.0f32, 0.0, 0.0];
        assert!(gpu_sphere_sphere_overlap(c1, 1.0, c2, 1.0));
    }

    #[test]
    fn sphere_overlap_just_touching() {
        // Centers 2.0 apart, radii 1.0 each → exactly touching
        let c1 = [0.0f32; 3];
        let c2 = [2.0f32, 0.0, 0.0];
        assert!(gpu_sphere_sphere_overlap(c1, 1.0, c2, 1.0));
    }

    #[test]
    fn sphere_overlap_separated() {
        let c1 = [0.0f32; 3];
        let c2 = [3.0f32, 0.0, 0.0];
        assert!(!gpu_sphere_sphere_overlap(c1, 1.0, c2, 1.0));
    }

    #[test]
    fn sphere_overlap_same_center() {
        let c = [0.0f32; 3];
        assert!(gpu_sphere_sphere_overlap(c, 1.0, c, 1.0));
    }

    // ── gpu_point_in_aabb ────────────────────────────────────────────────

    #[test]
    fn point_in_aabb_inside() {
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(gpu_point_in_aabb([0.5, 0.5, 0.5], min, max));
    }

    #[test]
    fn point_in_aabb_outside() {
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(!gpu_point_in_aabb([2.0, 0.5, 0.5], min, max));
    }

    #[test]
    fn point_in_aabb_corner() {
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(gpu_point_in_aabb([0.0, 0.0, 0.0], min, max));
        assert!(gpu_point_in_aabb([1.0, 1.0, 1.0], min, max));
    }

    #[test]
    fn point_in_aabb_just_outside() {
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(!gpu_point_in_aabb([1.0001, 0.5, 0.5], min, max));
    }

    // ── gpu_ray_aabb_intersect ────────────────────────────────────────────

    #[test]
    fn ray_aabb_hit() {
        let origin = [-2.0f32, 0.5, 0.5];
        let dir = [1.0f32, 0.0, 0.0];
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        let t = gpu_ray_aabb_intersect(origin, dir, min, max);
        assert!(t.is_some());
        assert!((t.unwrap() - 2.0).abs() < 1e-5);
    }

    #[test]
    fn ray_aabb_miss() {
        let origin = [-2.0f32, 5.0, 0.5];
        let dir = [1.0f32, 0.0, 0.0];
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(gpu_ray_aabb_intersect(origin, dir, min, max).is_none());
    }

    #[test]
    fn ray_aabb_origin_inside() {
        let origin = [0.5f32; 3];
        let dir = [1.0f32, 0.0, 0.0];
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        // Origin inside box — should return positive t
        let t = gpu_ray_aabb_intersect(origin, dir, min, max);
        assert!(t.is_some());
        assert!(t.unwrap() >= 0.0);
    }

    #[test]
    fn ray_aabb_opposite_direction() {
        // Ray pointing away from box
        let origin = [5.0f32, 0.5, 0.5];
        let dir = [1.0f32, 0.0, 0.0];
        let min = [0.0f32; 3];
        let max = [1.0f32; 3];
        assert!(gpu_ray_aabb_intersect(origin, dir, min, max).is_none());
    }

    // ── gpu_ray_sphere_intersect ──────────────────────────────────────────

    #[test]
    fn ray_sphere_hit() {
        let origin = [-5.0f32, 0.0, 0.0];
        let dir = [1.0f32, 0.0, 0.0];
        let center = [0.0f32; 3];
        let t = gpu_ray_sphere_intersect(origin, dir, center, 1.0);
        assert!(t.is_some());
        // Entry at x = -1, so t = 5 - 1 = 4
        assert!((t.unwrap() - 4.0).abs() < 1e-4);
    }

    #[test]
    fn ray_sphere_miss() {
        let origin = [0.0f32, 5.0, 0.0];
        let dir = [1.0f32, 0.0, 0.0];
        let center = [0.0f32; 3];
        assert!(gpu_ray_sphere_intersect(origin, dir, center, 1.0).is_none());
    }

    #[test]
    fn ray_sphere_origin_inside() {
        let origin = [0.0f32; 3];
        let dir = [1.0f32, 0.0, 0.0];
        let center = [0.0f32; 3];
        let t = gpu_ray_sphere_intersect(origin, dir, center, 2.0);
        // Should return the exit distance
        assert!(t.is_some());
        assert!(t.unwrap() > 0.0);
    }

    // ── GpuContactManifold ───────────────────────────────────────────────

    #[test]
    fn manifold_new_empty() {
        let m = GpuContactManifold::new(0, 1);
        assert_eq!(m.body_a, 0);
        assert_eq!(m.body_b, 1);
        assert_eq!(m.contact_count(), 0);
    }

    #[test]
    fn manifold_add_contact_increases_count() {
        let mut m = GpuContactManifold::new(0, 1);
        m.add_contact([0.0; 3], [0.0, 1.0, 0.0], 0.1);
        assert_eq!(m.contact_count(), 1);
        m.add_contact([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.2);
        assert_eq!(m.contact_count(), 2);
    }

    #[test]
    fn manifold_max_penetration_empty() {
        let m = GpuContactManifold::new(0, 1);
        assert!((m.max_penetration() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn manifold_max_penetration_multiple() {
        let mut m = GpuContactManifold::new(0, 1);
        m.add_contact([0.0; 3], [0.0, 1.0, 0.0], 0.1);
        m.add_contact([0.0; 3], [0.0, 1.0, 0.0], 0.5);
        m.add_contact([0.0; 3], [0.0, 1.0, 0.0], 0.3);
        assert!((m.max_penetration() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn manifold_contact_points_stored() {
        let mut m = GpuContactManifold::new(2, 3);
        let pt = [1.0, 2.0, 3.0];
        m.add_contact(pt, [0.0, 1.0, 0.0], 0.05);
        assert!((m.contact_points[0][0] - 1.0).abs() < 1e-6);
    }

    // ── gpu_gjk_distance ─────────────────────────────────────────────────

    #[test]
    fn gjk_distance_touching() {
        let verts_a: Vec<[f32; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let verts_b: Vec<[f32; 3]> = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let d = gpu_gjk_distance(&verts_a, &verts_b);
        assert!(d < 1e-5, "expected ~0, got {d}");
    }

    #[test]
    fn gjk_distance_separated() {
        let verts_a: Vec<[f32; 3]> = vec![[0.0; 3]];
        let verts_b: Vec<[f32; 3]> = vec![[3.0, 0.0, 0.0]];
        let d = gpu_gjk_distance(&verts_a, &verts_b);
        assert!((d - 3.0).abs() < 1e-5);
    }

    #[test]
    fn gjk_distance_empty_a() {
        let d = gpu_gjk_distance(&[], &[[1.0, 0.0, 0.0]]);
        assert!((d - 0.0).abs() < 1e-6);
    }

    #[test]
    fn gjk_distance_empty_b() {
        let d = gpu_gjk_distance(&[[1.0, 0.0, 0.0]], &[]);
        assert!((d - 0.0).abs() < 1e-6);
    }

    #[test]
    fn gjk_distance_same_point() {
        let v = [[0.0f32; 3]];
        let d = gpu_gjk_distance(&v, &v);
        assert!(d < 1e-6);
    }
}
