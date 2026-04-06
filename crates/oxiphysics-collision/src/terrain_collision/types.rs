//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Terrain sampling utilities: elevation profiles, slope maps.
pub struct TerrainSampling;
impl TerrainSampling {
    /// Sample terrain height along a path of world-space (x, z) points.
    pub fn elevation_profile(
        hf: &HeightfieldShape,
        path: &[[f64; 2]],
        samples_per_segment: usize,
    ) -> Vec<ElevationSample> {
        let mut result = Vec::new();
        let mut cumulative_dist = 0.0;
        for seg in path.windows(2) {
            let a = seg[0];
            let b = seg[1];
            let dx = b[0] - a[0];
            let dz = b[1] - a[1];
            let seg_len = (dx * dx + dz * dz).sqrt();
            for i in 0..samples_per_segment {
                let t = i as f64 / samples_per_segment as f64;
                let x = a[0] + dx * t;
                let z = a[1] + dz * t;
                let h = hf.height_at(x, z);
                let n = hf.normal_at(x, z);
                let slope = n[1].clamp(-1.0, 1.0).acos();
                let dist_along = t * seg_len;
                result.push(ElevationSample {
                    position: [x, h, z],
                    distance: cumulative_dist + dist_along,
                    height: h,
                    slope,
                });
            }
            cumulative_dist += seg_len;
        }
        result
    }
    /// Compute a slope map for the entire heightfield.
    /// Returns a Vec of slope angles (radians) in row-major order.
    pub fn slope_map(hf: &HeightfieldShape) -> Vec<f64> {
        let mut slopes = Vec::with_capacity(hf.nx * hf.nz);
        for iz in 0..hf.nz {
            for ix in 0..hf.nx {
                let x = ix as f64 * hf.scale_x;
                let z = iz as f64 * hf.scale_z;
                let n = hf.normal_at(x, z);
                let slope = n[1].clamp(-1.0, 1.0).acos();
                slopes.push(slope);
            }
        }
        slopes
    }
    /// Find the highest point within a rectangle of the heightfield.
    pub fn find_peak(
        hf: &HeightfieldShape,
        x_min: f64,
        z_min: f64,
        x_max: f64,
        z_max: f64,
    ) -> [f64; 3] {
        let ix_min = (x_min / hf.scale_x) as usize;
        let iz_min = (z_min / hf.scale_z) as usize;
        let ix_max = ((x_max / hf.scale_x) as usize).min(hf.nx - 1);
        let iz_max = ((z_max / hf.scale_z) as usize).min(hf.nz - 1);
        let mut best_h = f64::MIN;
        let mut best_pos = [0.0; 3];
        for iz in iz_min..=iz_max {
            for ix in ix_min..=ix_max {
                let h = hf.height_at_cell(ix, iz) + hf.y_offset;
                if h > best_h {
                    best_h = h;
                    best_pos = [ix as f64 * hf.scale_x, h, iz as f64 * hf.scale_z];
                }
            }
        }
        best_pos
    }
}
/// Slope category classification.
#[derive(Debug, Clone, PartialEq)]
pub enum SlopeCategory {
    /// Slope < 5 degrees.
    Flat,
    /// Slope 5–15 degrees.
    Gentle,
    /// Slope 15–30 degrees.
    Moderate,
    /// Slope 30–45 degrees.
    Steep,
    /// Slope > 45 degrees.
    Cliff,
}
/// A BVH node covering a rectangular patch of the heightfield.
#[derive(Debug, Clone)]
pub struct HfBvhNode {
    /// AABB of this patch.
    pub aabb_min: [f64; 3],
    /// AABB of this patch.
    pub aabb_max: [f64; 3],
    /// Start cell index in X.
    pub x0: usize,
    /// Start cell index in Z.
    pub z0: usize,
    /// End cell index in X (exclusive).
    pub x1: usize,
    /// End cell index in Z (exclusive).
    pub z1: usize,
    /// Children indices or None for leaf.
    pub children: Option<[usize; 2]>,
}
/// A cached terrain contact result with invalidation tracking.
#[derive(Debug, Clone)]
pub struct CachedTerrainContact {
    /// World-space query point used to generate this cache entry.
    pub query_point: [f64; 3],
    /// Cached contact result (None if no contact).
    pub contact: Option<TerrainContact>,
    /// Whether this cache entry is still valid.
    pub valid: bool,
    /// Frame/tick when this cache entry was created.
    pub frame: u64,
}
/// Terrain shadow map computed by casting rays from a light direction.
///
/// For each terrain cell, determines whether it is shadowed by tracing a ray
/// toward the light source and checking for terrain intersections.
pub struct TerrainShadowMap {
    /// Shadow values: 1.0 = lit, 0.0 = shadowed.
    pub shadow: Vec<f64>,
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub nz: usize,
}
impl TerrainShadowMap {
    /// Compute the shadow map for a given heightfield and light direction.
    ///
    /// `light_dir` should point from terrain toward the light (un-normalized ok).
    pub fn compute(hf: &HeightfieldShape, light_dir: [f64; 3]) -> Self {
        let ld = normalize3(light_dir);
        let nx = hf.nx;
        let nz = hf.nz;
        let mut shadow = vec![1.0f64; nx * nz];
        for iz in 0..nz {
            for ix in 0..nx {
                let x0 = ix as f64 * hf.scale_x;
                let z0 = iz as f64 * hf.scale_z;
                let y0 = hf.height_at_cell(ix, iz) + hf.y_offset + 0.01;
                let max_dist = (hf.nx as f64 * hf.scale_x + hf.nz as f64 * hf.scale_z) * 0.5;
                let n_steps = 20usize;
                let mut occluded = false;
                for step in 1..=n_steps {
                    let t = (step as f64 / n_steps as f64) * max_dist;
                    let rx = x0 + ld[0] * t;
                    let rz = z0 + ld[2] * t;
                    let ry = y0 + ld[1] * t;
                    if rx < 0.0 || rz < 0.0 {
                        break;
                    }
                    if rx > (nx - 1) as f64 * hf.scale_x || rz > (nz - 1) as f64 * hf.scale_z {
                        break;
                    }
                    let terrain_h = hf.height_at(rx, rz);
                    if ry < terrain_h {
                        occluded = true;
                        break;
                    }
                }
                shadow[iz * nx + ix] = if occluded { 0.0 } else { 1.0 };
            }
        }
        TerrainShadowMap { shadow, nx, nz }
    }
    /// Get shadow value at cell (ix, iz): 1.0 = lit, 0.0 = shadowed.
    pub fn get(&self, ix: usize, iz: usize) -> f64 {
        self.shadow[iz * self.nx + ix]
    }
    /// Fraction of terrain cells that are lit.
    pub fn lit_fraction(&self) -> f64 {
        let total = self.nx * self.nz;
        if total == 0 {
            return 0.0;
        }
        let lit: f64 = self.shadow.iter().sum();
        lit / total as f64
    }
}
/// Contact between a shape and the terrain.
#[derive(Debug, Clone)]
pub struct TerrainContact {
    /// Contact point on terrain surface.
    pub point: [f64; 3],
    /// Terrain surface normal at contact.
    pub normal: [f64; 3],
    /// Penetration depth (positive = penetrating).
    pub depth: f64,
}
/// A single elevation sample along a path.
#[derive(Debug, Clone)]
pub struct ElevationSample {
    /// World-space position along the path.
    pub position: [f64; 3],
    /// Accumulated distance from path start.
    pub distance: f64,
    /// Terrain height at this position.
    pub height: f64,
    /// Slope angle in radians at this position.
    pub slope: f64,
}
/// LOD-aware terrain collision with adaptive resolution.
///
/// Selects the appropriate collision resolution based on the distance
/// between the query point and the terrain reference point.
pub struct LodTerrainCollider<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
    /// World-space reference point (e.g., camera or physics object).
    pub reference_point: [f64; 3],
}
impl<'a> LodTerrainCollider<'a> {
    /// Construct with a heightfield and reference point.
    pub fn new(hf: &'a HeightfieldShape, reference_point: [f64; 3]) -> Self {
        LodTerrainCollider {
            hf,
            reference_point,
        }
    }
    /// Determine LOD level for a given world-space query point.
    pub fn lod_for(&self, query: [f64; 3]) -> LodLevel {
        let dx = query[0] - self.reference_point[0];
        let dz = query[2] - self.reference_point[2];
        let dist = (dx * dx + dz * dz).sqrt();
        LodLevel::from_distance(dist)
    }
    /// Sphere-terrain contact using LOD-adaptive sampling.
    ///
    /// Uses coarser height queries when object is far from the reference point.
    pub fn sphere_contact_lod(&self, center: [f64; 3], radius: f64) -> Option<TerrainContact> {
        let lod = self.lod_for(center);
        let stride = lod.stride() as f64;
        let snap_x = (center[0] / (self.hf.scale_x * stride)).round() * self.hf.scale_x * stride;
        let snap_z = (center[2] / (self.hf.scale_z * stride)).round() * self.hf.scale_z * stride;
        let terrain_h = self.hf.height_at(snap_x, snap_z);
        let depth = terrain_h + radius - center[1];
        if depth > 0.0 {
            let normal = self.hf.normal_at(snap_x, snap_z);
            Some(TerrainContact {
                point: [snap_x, terrain_h, snap_z],
                normal,
                depth,
            })
        } else {
            None
        }
    }
    /// Build a coarse height grid at the given LOD level.
    ///
    /// Returns a Vec of (world_x, world_z, height) samples.
    pub fn coarse_height_grid(&self, lod: LodLevel) -> Vec<(f64, f64, f64)> {
        let stride = lod.stride();
        let mut samples = Vec::new();
        let mut iz = 0;
        while iz < self.hf.nz {
            let mut ix = 0;
            while ix < self.hf.nx {
                let x = ix as f64 * self.hf.scale_x;
                let z = iz as f64 * self.hf.scale_z;
                let h = self.hf.height_at_cell(ix, iz) + self.hf.y_offset;
                samples.push((x, z, h));
                ix += stride;
            }
            iz += stride;
        }
        samples
    }
}
/// Hierarchical BVH over terrain patches for fast ray traversal.
pub struct HeightfieldBvh {
    /// BVH nodes.
    pub nodes: Vec<HfBvhNode>,
    /// Reference to the heightfield.
    pub hf: HeightfieldShape,
}
impl HeightfieldBvh {
    /// Construct a flat (1-level) BVH from the heightfield.
    pub fn build(hf: HeightfieldShape) -> Self {
        let nx = hf.nx;
        let nz = hf.nz;
        let mut nodes = Vec::new();
        let (h_min, h_max) = hf.height_range();
        let root = HfBvhNode {
            aabb_min: [0.0, h_min, 0.0],
            aabb_max: [
                (nx - 1) as f64 * hf.scale_x,
                h_max,
                (nz - 1) as f64 * hf.scale_z,
            ],
            x0: 0,
            z0: 0,
            x1: nx - 1,
            z1: nz - 1,
            children: None,
        };
        nodes.push(root);
        HeightfieldBvh { nodes, hf }
    }
    /// Test a ray against the BVH root AABB.
    pub fn ray_vs_root(&self, origin: [f64; 3], dir: [f64; 3]) -> bool {
        let node = &self.nodes[0];
        ray_aabb(origin, dir, node.aabb_min, node.aabb_max).is_some()
    }
}
/// Precomputed terrain contact cache for frequently queried locations.
///
/// Caches terrain contact results to avoid redundant height queries in
/// tight physics loops, with spatial tolerance-based invalidation.
pub struct TerrainContactCache {
    /// Cached entries.
    pub entries: Vec<CachedTerrainContact>,
    /// Spatial tolerance: if query moves less than this, cache is valid.
    pub position_tolerance: f64,
    /// Frame age tolerance: invalidate entries older than this many frames.
    pub max_age_frames: u64,
    /// Current frame counter.
    pub current_frame: u64,
}
impl TerrainContactCache {
    /// Create a new contact cache.
    pub fn new(position_tolerance: f64, max_age_frames: u64) -> Self {
        TerrainContactCache {
            entries: Vec::new(),
            position_tolerance,
            max_age_frames,
            current_frame: 0,
        }
    }
    /// Advance to the next frame (call once per physics tick).
    pub fn next_frame(&mut self) {
        self.current_frame += 1;
    }
    /// Look up a cached contact for a given query point.
    ///
    /// Returns `Some(&CachedTerrainContact)` if a valid, recent, nearby cache entry exists.
    pub fn lookup(&self, query: [f64; 3]) -> Option<&CachedTerrainContact> {
        let tol_sq = self.position_tolerance * self.position_tolerance;
        for entry in &self.entries {
            if !entry.valid {
                continue;
            }
            if self.current_frame - entry.frame > self.max_age_frames {
                continue;
            }
            let dx = entry.query_point[0] - query[0];
            let dy = entry.query_point[1] - query[1];
            let dz = entry.query_point[2] - query[2];
            if dx * dx + dy * dy + dz * dz < tol_sq {
                return Some(entry);
            }
        }
        None
    }
    /// Insert or update a cache entry for a query point.
    pub fn insert(&mut self, query: [f64; 3], contact: Option<TerrainContact>) {
        let tol_sq = self.position_tolerance * self.position_tolerance;
        for entry in &mut self.entries {
            let dx = entry.query_point[0] - query[0];
            let dy = entry.query_point[1] - query[1];
            let dz = entry.query_point[2] - query[2];
            if dx * dx + dy * dy + dz * dz < tol_sq {
                entry.query_point = query;
                entry.contact = contact;
                entry.valid = true;
                entry.frame = self.current_frame;
                return;
            }
        }
        self.entries.push(CachedTerrainContact {
            query_point: query,
            contact,
            valid: true,
            frame: self.current_frame,
        });
    }
    /// Invalidate all cache entries older than `max_age_frames`.
    pub fn invalidate_old(&mut self) {
        for entry in &mut self.entries {
            if self.current_frame - entry.frame > self.max_age_frames {
                entry.valid = false;
            }
        }
    }
    /// Clear all cache entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Number of valid cache entries.
    pub fn num_valid(&self) -> usize {
        self.entries.iter().filter(|e| e.valid).count()
    }
}
/// Detects underwater terrain and computes buoyancy corrections.
///
/// Combines a water surface height with the terrain heightfield to determine
/// which terrain cells are submerged, and by how much.
pub struct UnderwaterTerrainDetector<'a> {
    /// Terrain heightfield.
    pub hf: &'a HeightfieldShape,
    /// Global water surface level (flat water assumption).
    pub water_level: f64,
}
impl<'a> UnderwaterTerrainDetector<'a> {
    /// Construct with a heightfield and water level.
    pub fn new(hf: &'a HeightfieldShape, water_level: f64) -> Self {
        UnderwaterTerrainDetector { hf, water_level }
    }
    /// Check whether world-space point (x, y, z) is underwater.
    pub fn is_underwater(&self, x: f64, y: f64, z: f64) -> bool {
        let terrain_h = self.hf.height_at(x, z);
        y > terrain_h && y < self.water_level
    }
    /// Compute water depth above terrain at (x, z). Returns 0 if not submerged.
    pub fn water_depth_at(&self, x: f64, z: f64) -> f64 {
        let terrain_h = self.hf.height_at(x, z);
        (self.water_level - terrain_h).max(0.0)
    }
    /// Buoyancy force magnitude for a sphere of `radius` at (x, y, z).
    ///
    /// Uses Archimedes' principle: F = rho_water * g * V_submerged.
    pub fn buoyancy_force(
        &self,
        _x: f64,
        y: f64,
        _z: f64,
        radius: f64,
        rho_water: f64,
        gravity: f64,
    ) -> f64 {
        let depth = self.water_level - (y - radius);
        if depth <= 0.0 {
            return 0.0;
        }
        let submerged_fraction = (depth / (2.0 * radius)).clamp(0.0, 1.0);
        let v_full = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        let v_sub = v_full * submerged_fraction;
        rho_water * gravity * v_sub
    }
    /// Count the number of fully submerged terrain cells.
    pub fn count_submerged_cells(&self) -> usize {
        let mut count = 0;
        for iz in 0..self.hf.nz {
            for ix in 0..self.hf.nx {
                let h = self.hf.height_at_cell(ix, iz) + self.hf.y_offset;
                if h < self.water_level {
                    count += 1;
                }
            }
        }
        count
    }
}
/// Ray-terrain intersection result.
#[derive(Debug, Clone)]
pub struct TerrainRayHit {
    /// World-space hit point.
    pub point: [f64; 3],
    /// Surface normal at hit point.
    pub normal: [f64; 3],
    /// Ray parameter t.
    pub t: f64,
    /// Cell index in X.
    pub cell_x: usize,
    /// Cell index in Z.
    pub cell_z: usize,
}
/// Voxel grid terrain with marching cubes isosurface extraction and AABB queries.
pub struct VoxelTerrain {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Grid dimensions.
    pub nz: usize,
    /// Voxel size.
    pub voxel_size: f64,
    /// World-space origin.
    pub origin: [f64; 3],
    /// Density values (positive = solid, negative = empty).
    pub densities: Vec<f64>,
}
impl VoxelTerrain {
    /// Construct an empty voxel terrain.
    pub fn new(nx: usize, ny: usize, nz: usize, voxel_size: f64, origin: [f64; 3]) -> Self {
        VoxelTerrain {
            nx,
            ny,
            nz,
            voxel_size,
            origin,
            densities: vec![0.0; nx * ny * nz],
        }
    }
    fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }
    /// Set density at voxel (ix, iy, iz).
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, density: f64) {
        let idx = self.index(ix, iy, iz);
        self.densities[idx] = density;
    }
    /// Get density at voxel (ix, iy, iz).
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.densities[self.index(ix, iy, iz)]
    }
    /// Test if a voxel is solid (density >= 0).
    pub fn is_solid(&self, ix: usize, iy: usize, iz: usize) -> bool {
        self.get(ix, iy, iz) >= 0.0
    }
    /// AABB vs voxel query: find all solid voxels overlapping the AABB.
    pub fn aabb_query(&self, aabb_min: [f64; 3], aabb_max: [f64; 3]) -> Vec<(usize, usize, usize)> {
        let vs = self.voxel_size;
        let ix0 = (((aabb_min[0] - self.origin[0]) / vs).floor() as isize).max(0) as usize;
        let iy0 = (((aabb_min[1] - self.origin[1]) / vs).floor() as isize).max(0) as usize;
        let iz0 = (((aabb_min[2] - self.origin[2]) / vs).floor() as isize).max(0) as usize;
        let ix1 = ((((aabb_max[0] - self.origin[0]) / vs).ceil() as isize).max(0) as usize)
            .min(self.nx - 1);
        let iy1 = ((((aabb_max[1] - self.origin[1]) / vs).ceil() as isize).max(0) as usize)
            .min(self.ny - 1);
        let iz1 = ((((aabb_max[2] - self.origin[2]) / vs).ceil() as isize).max(0) as usize)
            .min(self.nz - 1);
        let mut result = Vec::new();
        for ix in ix0..=ix1 {
            for iy in iy0..=iy1 {
                for iz in iz0..=iz1 {
                    if self.is_solid(ix, iy, iz) {
                        result.push((ix, iy, iz));
                    }
                }
            }
        }
        result
    }
    /// Extract the isosurface at density = 0 using a simplified marching cubes
    /// (returns triangles as vertex triples).
    pub fn extract_isosurface(&self) -> Vec<[[f64; 3]; 3]> {
        let mut tris = Vec::new();
        let vs = self.voxel_size;
        for iz in 0..self.nz.saturating_sub(1) {
            for iy in 0..self.ny.saturating_sub(1) {
                for ix in 0..self.nx.saturating_sub(1) {
                    let d000 = self.get(ix, iy, iz);
                    let d100 = self.get(ix + 1, iy, iz);
                    let d010 = self.get(ix, iy + 1, iz);
                    let has_positive = d000 > 0.0 || d100 > 0.0 || d010 > 0.0;
                    let has_negative = d000 < 0.0 || d100 < 0.0 || d010 < 0.0;
                    if has_positive && has_negative {
                        let base = voxel_to_world(ix, iy, iz, vs, self.origin);
                        let v0 = base;
                        let v1 = add3(base, [vs, 0.0, 0.0]);
                        let v2 = add3(base, [0.0, vs, 0.0]);
                        tris.push([v0, v1, v2]);
                    }
                }
            }
        }
        tris
    }
}
/// Cliff detection based on gradient thresholding.
///
/// A cell is considered a cliff if the local height gradient magnitude
/// exceeds a configurable threshold angle.
pub struct CliffDetector<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
    /// Gradient angle threshold in radians above which a cell is a cliff.
    pub threshold_angle: f64,
}
impl<'a> CliffDetector<'a> {
    /// Construct a cliff detector.
    pub fn new(hf: &'a HeightfieldShape, threshold_angle: f64) -> Self {
        CliffDetector {
            hf,
            threshold_angle,
        }
    }
    /// Check whether the terrain at (x, z) is a cliff.
    pub fn is_cliff(&self, x: f64, z: f64) -> bool {
        let n = self.hf.normal_at(x, z);
        let angle = n[1].clamp(-1.0, 1.0).acos();
        angle >= self.threshold_angle
    }
    /// Build a cliff mask for the entire heightfield.
    ///
    /// Returns a Vec`bool` in row-major order (true = cliff cell).
    pub fn cliff_mask(&self) -> Vec<bool> {
        let mut mask = Vec::with_capacity(self.hf.nx * self.hf.nz);
        for iz in 0..self.hf.nz {
            for ix in 0..self.hf.nx {
                let x = ix as f64 * self.hf.scale_x;
                let z = iz as f64 * self.hf.scale_z;
                mask.push(self.is_cliff(x, z));
            }
        }
        mask
    }
    /// Count cliff cells in the heightfield.
    pub fn num_cliffs(&self) -> usize {
        self.cliff_mask().iter().filter(|&&c| c).count()
    }
    /// Find the steepest cell in a rectangular region.
    pub fn steepest_cell(
        &self,
        x_min: f64,
        z_min: f64,
        x_max: f64,
        z_max: f64,
    ) -> (usize, usize, f64) {
        let ix_min = (x_min / self.hf.scale_x) as usize;
        let iz_min = (z_min / self.hf.scale_z) as usize;
        let ix_max = ((x_max / self.hf.scale_x) as usize).min(self.hf.nx - 1);
        let iz_max = ((z_max / self.hf.scale_z) as usize).min(self.hf.nz - 1);
        let mut max_angle = 0.0_f64;
        let mut best_ix = ix_min;
        let mut best_iz = iz_min;
        for iz in iz_min..=iz_max {
            for ix in ix_min..=ix_max {
                let x = ix as f64 * self.hf.scale_x;
                let z = iz as f64 * self.hf.scale_z;
                let angle = self.hf.normal_at(x, z)[1].clamp(-1.0, 1.0).acos();
                if angle > max_angle {
                    max_angle = angle;
                    best_ix = ix;
                    best_iz = iz;
                }
            }
        }
        (best_ix, best_iz, max_angle)
    }
}
/// Local heightfield modification: craters, tire tracks, SDF blending.
pub struct TerrainDeformation {
    /// Target heightfield.
    pub hf: HeightfieldShape,
}
impl TerrainDeformation {
    /// Construct from a heightfield.
    pub fn new(hf: HeightfieldShape) -> Self {
        TerrainDeformation { hf }
    }
    /// Create a circular crater at world position (cx, cz) with given radius
    /// and depth.
    pub fn apply_crater(&mut self, cx: f64, cz: f64, radius: f64, depth: f64) {
        let inv_sx = 1.0 / self.hf.scale_x;
        let inv_sz = 1.0 / self.hf.scale_z;
        let cell_r = (radius * inv_sx.max(inv_sz)).ceil() as usize;
        let ix_center = (cx * inv_sx) as usize;
        let iz_center = (cz * inv_sz) as usize;
        let nx = self.hf.nx;
        let nz = self.hf.nz;
        for dix in 0..=(2 * cell_r) {
            for diz in 0..=(2 * cell_r) {
                let ix = ix_center.saturating_sub(cell_r) + dix;
                let iz = iz_center.saturating_sub(cell_r) + diz;
                if ix >= nx || iz >= nz {
                    continue;
                }
                let wx = ix as f64 * self.hf.scale_x;
                let wz = iz as f64 * self.hf.scale_z;
                let dx = wx - cx;
                let dz_ = wz - cz;
                let dist = (dx * dx + dz_ * dz_).sqrt();
                if dist < radius {
                    let factor = 1.0 - (dist / radius).powi(2);
                    self.hf.heights[iz * nx + ix] -= depth * factor;
                }
            }
        }
    }
    /// Apply a linear track (tire track) along a path defined by a start/end
    /// point, width, and depth.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_track(&mut self, start: [f64; 2], end: [f64; 2], width: f64, depth: f64) {
        let d = [end[0] - start[0], end[1] - start[1]];
        let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
        if len < 1e-14 {
            return;
        }
        let nx = self.hf.nx;
        let nz = self.hf.nz;
        for ix in 0..nx {
            for iz in 0..nz {
                let wx = ix as f64 * self.hf.scale_x;
                let wz = iz as f64 * self.hf.scale_z;
                let px = wx - start[0];
                let pz = wz - start[1];
                let t = ((px * d[0] + pz * d[1]) / (len * len)).clamp(0.0, 1.0);
                let proj_x = start[0] + t * d[0];
                let proj_z = start[1] + t * d[1];
                let dist = ((wx - proj_x).powi(2) + (wz - proj_z).powi(2)).sqrt();
                if dist < width * 0.5 {
                    let factor = 1.0 - 2.0 * dist / width;
                    self.hf.heights[iz * nx + ix] -= depth * factor;
                }
            }
        }
    }
    /// Blend in a raised mound at (cx, cz) with given radius and height.
    pub fn apply_mound(&mut self, cx: f64, cz: f64, radius: f64, height: f64) {
        let inv_sx = 1.0 / self.hf.scale_x;
        let inv_sz = 1.0 / self.hf.scale_z;
        let cell_r = (radius * inv_sx.max(inv_sz)).ceil() as usize;
        let ix_center = (cx * inv_sx) as usize;
        let iz_center = (cz * inv_sz) as usize;
        let nx = self.hf.nx;
        let nz = self.hf.nz;
        for dix in 0..=(2 * cell_r) {
            for diz in 0..=(2 * cell_r) {
                let ix = ix_center.saturating_sub(cell_r) + dix;
                let iz = iz_center.saturating_sub(cell_r) + diz;
                if ix >= nx || iz >= nz {
                    continue;
                }
                let wx = ix as f64 * self.hf.scale_x;
                let wz = iz as f64 * self.hf.scale_z;
                let dx = wx - cx;
                let dz_ = wz - cz;
                let dist = (dx * dx + dz_ * dz_).sqrt();
                if dist < radius {
                    let factor = 1.0 - dist / radius;
                    self.hf.heights[iz * nx + ix] += height * factor;
                }
            }
        }
    }
}
/// A heightfield terrain shape on an nx × nz grid with uniform cell spacing.
pub struct HeightfieldShape {
    /// Number of cells in the X direction.
    pub nx: usize,
    /// Number of cells in the Z direction.
    pub nz: usize,
    /// World-space scale per cell in X.
    pub scale_x: f64,
    /// World-space scale per cell in Z.
    pub scale_z: f64,
    /// Heights stored in row-major order (z-row, x-column).
    pub heights: Vec<f64>,
    /// World-space Y offset.
    pub y_offset: f64,
}
impl HeightfieldShape {
    /// Construct a heightfield.
    pub fn new(nx: usize, nz: usize, scale_x: f64, scale_z: f64, heights: Vec<f64>) -> Self {
        assert_eq!(heights.len(), nx * nz, "heights length must be nx*nz");
        HeightfieldShape {
            nx,
            nz,
            scale_x,
            scale_z,
            heights,
            y_offset: 0.0,
        }
    }
    /// Get height at grid cell (ix, iz).
    pub fn height_at_cell(&self, ix: usize, iz: usize) -> f64 {
        self.heights[iz * self.nx + ix]
    }
    /// Interpolated height at world-space (x, z).
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let cx = (x / self.scale_x).max(0.0);
        let cz = (z / self.scale_z).max(0.0);
        let ix = (cx as usize).min(self.nx - 2);
        let iz = (cz as usize).min(self.nz - 2);
        let tx = cx - ix as f64;
        let tz = cz - iz as f64;
        let h00 = self.height_at_cell(ix, iz);
        let h10 = self.height_at_cell(ix + 1, iz);
        let h01 = self.height_at_cell(ix, iz + 1);
        let h11 = self.height_at_cell(ix + 1, iz + 1);
        bilinear_height_interp(h00, h10, h01, h11, tx, tz) + self.y_offset
    }
    /// Surface normal at world-space (x, z) using central differences.
    pub fn normal_at(&self, x: f64, z: f64) -> [f64; 3] {
        let h = 0.5_f64.max(self.scale_x.min(self.scale_z));
        let h_nx = self.height_at(x - h, z);
        let h_px = self.height_at(x + h, z);
        let h_nz = self.height_at(x, z - h);
        let h_pz = self.height_at(x, z + h);
        terrain_normal_from_heights(h_nx, h_px, h_nz, h_pz, h, h)
    }
    /// Minimum and maximum height in the field.
    pub fn height_range(&self) -> (f64, f64) {
        let mut mn = f64::MAX;
        let mut mx = f64::MIN;
        for &h in &self.heights {
            mn = mn.min(h);
            mx = mx.max(h);
        }
        (mn + self.y_offset, mx + self.y_offset)
    }
    /// World-space extent of the heightfield.
    pub fn world_extent(&self) -> ([f64; 3], [f64; 3]) {
        let (h_min, h_max) = self.height_range();
        (
            [0.0, h_min, 0.0],
            [
                (self.nx - 1) as f64 * self.scale_x,
                h_max,
                (self.nz - 1) as f64 * self.scale_z,
            ],
        )
    }
}
/// Slope and aspect analysis for terrain navigation and physics.
///
/// Slope is the steepest angle of descent; aspect is the compass direction
/// of steepest descent. Used for friction modulation, vegetation placement,
/// and AI pathfinding.
pub struct TerrainSlopeAspect<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
}
impl<'a> TerrainSlopeAspect<'a> {
    /// Construct for a given heightfield.
    pub fn new(hf: &'a HeightfieldShape) -> Self {
        TerrainSlopeAspect { hf }
    }
    /// Compute slope angle (radians) at world-space (x, z).
    ///
    /// Returns the angle between the terrain normal and the up vector \[0,1,0\].
    pub fn slope_angle(&self, x: f64, z: f64) -> f64 {
        let n = self.hf.normal_at(x, z);
        n[1].clamp(-1.0, 1.0).acos()
    }
    /// Compute aspect angle (radians from north = -z axis) at world-space (x, z).
    ///
    /// Returns angle in \[0, 2π\] measured clockwise from north.
    pub fn aspect_angle(&self, x: f64, z: f64) -> f64 {
        let h = self.hf.scale_x.min(self.hf.scale_z) * 0.5;
        let dhdx = (self.hf.height_at(x + h, z) - self.hf.height_at(x - h, z)) / (2.0 * h);
        let dhdz = (self.hf.height_at(x, z + h) - self.hf.height_at(x, z - h)) / (2.0 * h);
        let angle = dhdz.atan2(dhdx);
        if angle < 0.0 {
            angle + 2.0 * std::f64::consts::PI
        } else {
            angle
        }
    }
    /// Classify slope category.
    ///
    /// Returns one of: Flat (< 5°), Gentle (5–15°), Moderate (15–30°),
    /// Steep (30–45°), Cliff (> 45°).
    pub fn classify_slope(&self, x: f64, z: f64) -> SlopeCategory {
        let deg = self.slope_angle(x, z).to_degrees();
        if deg < 5.0 {
            SlopeCategory::Flat
        } else if deg < 15.0 {
            SlopeCategory::Gentle
        } else if deg < 30.0 {
            SlopeCategory::Moderate
        } else if deg < 45.0 {
            SlopeCategory::Steep
        } else {
            SlopeCategory::Cliff
        }
    }
    /// Compute full slope map with categories.
    pub fn slope_category_map(&self) -> Vec<SlopeCategory> {
        let mut cats = Vec::with_capacity(self.hf.nx * self.hf.nz);
        for iz in 0..self.hf.nz {
            for ix in 0..self.hf.nx {
                let x = ix as f64 * self.hf.scale_x;
                let z = iz as f64 * self.hf.scale_z;
                cats.push(self.classify_slope(x, z));
            }
        }
        cats
    }
}
/// A fully hierarchical BVH over a heightfield, built by recursive subdivision.
///
/// Splits the terrain into quadtree patches until a minimum patch size is reached,
/// storing AABB bounds at each node for fast ray traversal and broadphase queries.
pub struct HeightfieldBvhFull {
    /// All BVH nodes.
    pub nodes: Vec<HfBvhNode>,
    /// The underlying heightfield.
    pub hf: HeightfieldShape,
    /// Maximum cells per leaf patch.
    pub leaf_size: usize,
}
impl HeightfieldBvhFull {
    /// Build an iterative BVH over the heightfield.
    pub fn build(hf: HeightfieldShape, leaf_size: usize) -> Self {
        let ls = leaf_size.max(1);
        let mut nodes = Vec::new();
        let (h_min, h_max) = hf.height_range();
        let nx = hf.nx;
        let nz = hf.nz;
        let root = HfBvhNode {
            aabb_min: [0.0, h_min, 0.0],
            aabb_max: [
                (nx - 1) as f64 * hf.scale_x,
                h_max,
                (nz - 1) as f64 * hf.scale_z,
            ],
            x0: 0,
            z0: 0,
            x1: nx - 1,
            z1: nz - 1,
            children: None,
        };
        nodes.push(root);
        let mut work_stack: Vec<usize> = vec![0];
        while let Some(idx) = work_stack.pop() {
            let (x0, x1, z0, z1) = {
                let n = &nodes[idx];
                (n.x0, n.x1, n.z0, n.z1)
            };
            let x_span = x1.saturating_sub(x0);
            let z_span = z1.saturating_sub(z0);
            if x_span <= ls && z_span <= ls {
                continue;
            }
            let (ca, cb) = if x_span >= z_span {
                let xm = x0 + x_span / 2;
                (
                    Self::make_node(&hf, x0, xm, z0, z1),
                    Self::make_node(&hf, xm, x1, z0, z1),
                )
            } else {
                let zm = z0 + z_span / 2;
                (
                    Self::make_node(&hf, x0, x1, z0, zm),
                    Self::make_node(&hf, x0, x1, zm, z1),
                )
            };
            let ca_x = ca.x1.saturating_sub(ca.x0);
            let ca_z = ca.z1.saturating_sub(ca.z0);
            let cb_x = cb.x1.saturating_sub(cb.x0);
            let cb_z = cb.z1.saturating_sub(cb.z0);
            let ca_ok = ca_x < x_span || ca_z < z_span;
            let cb_ok = cb_x < x_span || cb_z < z_span;
            if !ca_ok || !cb_ok {
                continue;
            }
            let ca_idx = nodes.len();
            nodes.push(ca);
            let cb_idx = nodes.len();
            nodes.push(cb);
            nodes[idx].children = Some([ca_idx, cb_idx]);
            work_stack.push(ca_idx);
            work_stack.push(cb_idx);
        }
        HeightfieldBvhFull {
            nodes,
            hf,
            leaf_size: ls,
        }
    }
    fn make_node(hf: &HeightfieldShape, x0: usize, x1: usize, z0: usize, z1: usize) -> HfBvhNode {
        let mut h_min = f64::MAX;
        let mut h_max = f64::MIN;
        for iz in z0..=z1.min(hf.nz - 1) {
            for ix in x0..=x1.min(hf.nx - 1) {
                let h = hf.height_at_cell(ix, iz);
                h_min = h_min.min(h);
                h_max = h_max.max(h);
            }
        }
        HfBvhNode {
            aabb_min: [
                x0 as f64 * hf.scale_x,
                h_min + hf.y_offset,
                z0 as f64 * hf.scale_z,
            ],
            aabb_max: [
                x1 as f64 * hf.scale_x,
                h_max + hf.y_offset,
                z1 as f64 * hf.scale_z,
            ],
            x0,
            z0,
            x1,
            z1,
            children: None,
        }
    }
    /// Traverse the BVH to collect leaf nodes overlapping a query AABB.
    pub fn query_aabb(&self, qmin: [f64; 3], qmax: [f64; 3]) -> Vec<usize> {
        let mut result = Vec::new();
        self.traverse_aabb(0, qmin, qmax, &mut result);
        result
    }
    fn traverse_aabb(&self, idx: usize, qmin: [f64; 3], qmax: [f64; 3], result: &mut Vec<usize>) {
        let node = &self.nodes[idx];
        for k in 0..3 {
            if qmax[k] < node.aabb_min[k] || qmin[k] > node.aabb_max[k] {
                return;
            }
        }
        if let Some([ca, cb]) = node.children {
            self.traverse_aabb(ca, qmin, qmax, result);
            self.traverse_aabb(cb, qmin, qmax, result);
        } else {
            result.push(idx);
        }
    }
    /// Total number of BVH nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
}
/// CCD result for a fast-moving object vs terrain.
#[derive(Debug, Clone)]
pub struct TerrainCcdResult {
    /// Time of impact in \[0, 1\].
    pub toi: f64,
    /// World-space contact point.
    pub point: [f64; 3],
    /// Normal at contact.
    pub normal: [f64; 3],
}
/// Animated water heightfield using superposition of sine waves.
pub struct WaterSurface {
    /// Wave components: (amplitude, frequency_x, frequency_z, phase, speed).
    pub waves: Vec<(f64, f64, f64, f64, f64)>,
    /// Current simulation time.
    pub time: f64,
    /// Base water level.
    pub base_level: f64,
}
impl WaterSurface {
    /// Construct with a base water level.
    pub fn new(base_level: f64) -> Self {
        WaterSurface {
            waves: Vec::new(),
            time: 0.0,
            base_level,
        }
    }
    /// Add a wave component.
    pub fn add_wave(&mut self, amplitude: f64, freq_x: f64, freq_z: f64, phase: f64, speed: f64) {
        self.waves.push((amplitude, freq_x, freq_z, phase, speed));
    }
    /// Advance time.
    pub fn step(&mut self, dt: f64) {
        self.time += dt;
    }
    /// Water surface height at (x, z) at the current time.
    pub fn height_at(&self, x: f64, z: f64) -> f64 {
        let mut h = self.base_level;
        for &(amp, fx, fz, phase, speed) in &self.waves {
            h += amp * (fx * x + fz * z + phase + speed * self.time).sin();
        }
        h
    }
    /// Test if a point is floating (above water), submerged, or at surface.
    /// Returns (is_submerged, depth_below_surface).
    pub fn float_test(&self, point: [f64; 3]) -> (bool, f64) {
        let wh = self.height_at(point[0], point[2]);
        let depth = wh - point[1];
        (depth > 0.0, depth.max(0.0))
    }
}
/// Computes terrain normals using bilinear interpolation of corner normals.
///
/// More accurate than central-differences for highly curved terrain patches,
/// by blending analytic normals from the four cell corners.
pub struct TerrainNormalBilinear<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
}
impl<'a> TerrainNormalBilinear<'a> {
    /// Construct for a given heightfield.
    pub fn new(hf: &'a HeightfieldShape) -> Self {
        TerrainNormalBilinear { hf }
    }
    /// Compute the corner normal for cell vertex (ix, iz).
    pub fn corner_normal(&self, ix: usize, iz: usize) -> [f64; 3] {
        let x = ix as f64 * self.hf.scale_x;
        let z = iz as f64 * self.hf.scale_z;
        self.hf.normal_at(x, z)
    }
    /// Bilinear-interpolated normal at world-space (x, z).
    ///
    /// Blends the four corner normals of the containing cell.
    pub fn normal_at(&self, x: f64, z: f64) -> [f64; 3] {
        let cx = (x / self.hf.scale_x).max(0.0);
        let cz = (z / self.hf.scale_z).max(0.0);
        let ix = (cx as usize).min(self.hf.nx.saturating_sub(2));
        let iz = (cz as usize).min(self.hf.nz.saturating_sub(2));
        let tx = cx - ix as f64;
        let tz = cz - iz as f64;
        let n00 = self.corner_normal(ix, iz);
        let n10 = self.corner_normal((ix + 1).min(self.hf.nx - 1), iz);
        let n01 = self.corner_normal(ix, (iz + 1).min(self.hf.nz - 1));
        let n11 = self.corner_normal((ix + 1).min(self.hf.nx - 1), (iz + 1).min(self.hf.nz - 1));
        let nx = n00[0] * (1.0 - tx) * (1.0 - tz)
            + n10[0] * tx * (1.0 - tz)
            + n01[0] * (1.0 - tx) * tz
            + n11[0] * tx * tz;
        let ny = n00[1] * (1.0 - tx) * (1.0 - tz)
            + n10[1] * tx * (1.0 - tz)
            + n01[1] * (1.0 - tx) * tz
            + n11[1] * tx * tz;
        let nz = n00[2] * (1.0 - tx) * (1.0 - tz)
            + n10[2] * tx * (1.0 - tz)
            + n01[2] * (1.0 - tx) * tz
            + n11[2] * tx * tz;
        normalize3([nx, ny, nz])
    }
}
/// Terrain collider for spheres, capsules, and boxes against a heightfield.
pub struct TerrainCollider<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
}
impl<'a> TerrainCollider<'a> {
    /// Construct for a given heightfield.
    pub fn new(hf: &'a HeightfieldShape) -> Self {
        TerrainCollider { hf }
    }
    /// Sphere vs heightfield contact.
    pub fn sphere_vs_terrain(&self, center: [f64; 3], radius: f64) -> Option<TerrainContact> {
        let terrain_h = self.hf.height_at(center[0], center[2]);
        let depth = terrain_h + radius - center[1];
        if depth > 0.0 {
            let normal = self.hf.normal_at(center[0], center[2]);
            let point = [center[0], terrain_h, center[2]];
            Some(TerrainContact {
                point,
                normal,
                depth,
            })
        } else {
            None
        }
    }
    /// Capsule vs heightfield contact (tests both endpoints and midpoint).
    #[allow(clippy::too_many_arguments)]
    pub fn capsule_vs_terrain(
        &self,
        cap_a: [f64; 3],
        cap_b: [f64; 3],
        radius: f64,
    ) -> Vec<TerrainContact> {
        let mut contacts = Vec::new();
        for frac in [0.0, 0.5, 1.0] {
            let center = add3(cap_a, scale3(sub3(cap_b, cap_a), frac));
            if let Some(c) = self.sphere_vs_terrain(center, radius) {
                contacts.push(c);
            }
        }
        contacts
    }
    /// Box vs heightfield contact (tests all 8 corners).
    pub fn box_vs_terrain(&self, box_min: [f64; 3], box_max: [f64; 3]) -> Vec<TerrainContact> {
        let mut contacts = Vec::new();
        for &xi in &[box_min[0], box_max[0]] {
            for &zi in &[box_min[2], box_max[2]] {
                let terrain_h = self.hf.height_at(xi, zi);
                if box_min[1] < terrain_h {
                    let depth = terrain_h - box_min[1];
                    let normal = self.hf.normal_at(xi, zi);
                    contacts.push(TerrainContact {
                        point: [xi, terrain_h, zi],
                        normal,
                        depth,
                    });
                }
            }
        }
        contacts
    }
}
/// Material layer identifiers.
#[derive(Debug, Clone, PartialEq)]
pub enum TerrainLayer {
    /// Bedrock (granite, basalt, etc.). Highest hardness.
    Rock,
    /// Soil (loam, clay). Medium hardness, deformable.
    Soil,
    /// Vegetation (grass, moss). Low hardness, high energy absorption.
    Vegetation,
    /// Sand / gravel. Low hardness, high rolling friction.
    Sand,
}
/// Procedural crater impact deformation with realistic ejecta rim.
///
/// Models the parabolic crater bowl and elevated ejecta rim that form
/// from a ballistic impact. Supports multiple overlapping craters.
pub struct ProceduralCraterImpact;
impl ProceduralCraterImpact {
    /// Apply a single impact crater with an ejecta rim.
    ///
    /// * `hf` – mutable heightfield to deform
    /// * `cx`, `cz` – impact center in world space
    /// * `radius` – crater radius
    /// * `depth` – maximum crater depth at center
    /// * `rim_height` – height of the ejecta rim
    /// * `rim_width_fraction` – rim width as a fraction of radius
    pub fn apply_crater_with_rim(
        hf: &mut HeightfieldShape,
        cx: f64,
        cz: f64,
        radius: f64,
        depth: f64,
        rim_height: f64,
        rim_width_fraction: f64,
    ) {
        let inv_sx = 1.0 / hf.scale_x;
        let inv_sz = 1.0 / hf.scale_z;
        let cell_r = ((radius * 1.5) * inv_sx.max(inv_sz)).ceil() as usize;
        let ix_center = (cx * inv_sx) as usize;
        let iz_center = (cz * inv_sz) as usize;
        let nx = hf.nx;
        let nz = hf.nz;
        for dix in 0..=(2 * cell_r) {
            for diz in 0..=(2 * cell_r) {
                let ix = ix_center.saturating_sub(cell_r) + dix;
                let iz = iz_center.saturating_sub(cell_r) + diz;
                if ix >= nx || iz >= nz {
                    continue;
                }
                let wx = ix as f64 * hf.scale_x;
                let wz = iz as f64 * hf.scale_z;
                let dist = ((wx - cx).powi(2) + (wz - cz).powi(2)).sqrt();
                let r_norm = dist / radius;
                let delta = if r_norm < 1.0 {
                    -depth * (1.0 - r_norm * r_norm)
                } else if r_norm < 1.0 + rim_width_fraction {
                    let rim_r = (r_norm - 1.0) / rim_width_fraction;
                    rim_height * (-2.0 * rim_r * rim_r).exp()
                } else {
                    0.0
                };
                hf.heights[iz * nx + ix] += delta;
            }
        }
    }
    /// Apply multiple overlapping craters from an impactor list.
    ///
    /// Each impactor is `(cx, cz, radius, depth)`.
    pub fn apply_multiple(
        hf: &mut HeightfieldShape,
        impactors: &[(f64, f64, f64, f64)],
        rim_height_fraction: f64,
    ) {
        for &(cx, cz, radius, depth) in impactors {
            Self::apply_crater_with_rim(
                hf,
                cx,
                cz,
                radius,
                depth,
                depth * rim_height_fraction,
                0.3,
            );
        }
    }
}
/// Efficient ray traversal over a heightfield using a DDA (Digital
/// Differential Analyzer) grid walk.
pub struct TerrainRaycast<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
}
impl<'a> TerrainRaycast<'a> {
    /// Construct for a given heightfield.
    pub fn new(hf: &'a HeightfieldShape) -> Self {
        TerrainRaycast { hf }
    }
    /// Cast a ray and return the first terrain hit using DDA traversal.
    pub fn cast(&self, origin: [f64; 3], dir: [f64; 3]) -> Option<TerrainRayHit> {
        let inv_sx = 1.0 / self.hf.scale_x;
        let inv_sz = 1.0 / self.hf.scale_z;
        let mut cx = ((origin[0] * inv_sx) as isize).max(0) as usize;
        let mut cz = ((origin[2] * inv_sz) as isize).max(0) as usize;
        let dx = dir[0] * inv_sx;
        let dz = dir[2] * inv_sz;
        let step_x: isize = if dx >= 0.0 { 1 } else { -1 };
        let step_z: isize = if dz >= 0.0 { 1 } else { -1 };
        let t_delta_x = if dx.abs() > 1e-14 {
            1.0 / dx.abs()
        } else {
            f64::MAX
        };
        let t_delta_z = if dz.abs() > 1e-14 {
            1.0 / dz.abs()
        } else {
            f64::MAX
        };
        let mut t_max_x = if dx.abs() > 1e-14 {
            let boundary = if dx >= 0.0 {
                ((cx + 1) as f64 - origin[0] * inv_sx) / dx
            } else {
                (cx as f64 - origin[0] * inv_sx) / dx
            };
            boundary.abs()
        } else {
            f64::MAX
        };
        let mut t_max_z = if dz.abs() > 1e-14 {
            let boundary = if dz >= 0.0 {
                ((cz + 1) as f64 - origin[2] * inv_sz) / dz
            } else {
                (cz as f64 - origin[2] * inv_sz) / dz
            };
            boundary.abs()
        } else {
            f64::MAX
        };
        let max_steps = (self.hf.nx + self.hf.nz) * 2;
        for _ in 0..max_steps {
            if cx >= self.hf.nx - 1 || cz >= self.hf.nz - 1 {
                break;
            }
            let h = self.hf.height_at(
                (cx as f64 + 0.5) * self.hf.scale_x,
                (cz as f64 + 0.5) * self.hf.scale_z,
            );
            let t = if dir[1].abs() > 1e-14 {
                (h - origin[1]) / dir[1]
            } else {
                continue;
            };
            if t >= 0.0 {
                let hit_x = origin[0] + dir[0] * t;
                let hit_z = origin[2] + dir[2] * t;
                let cell_x = (hit_x * inv_sx) as usize;
                let cell_z = (hit_z * inv_sz) as usize;
                if cell_x == cx && cell_z == cz {
                    let point = [hit_x, h, hit_z];
                    let normal = self.hf.normal_at(hit_x, hit_z);
                    return Some(TerrainRayHit {
                        point,
                        normal,
                        t,
                        cell_x: cx,
                        cell_z: cz,
                    });
                }
            }
            if t_max_x < t_max_z {
                t_max_x += t_delta_x;
                if step_x > 0 {
                    cx += 1;
                } else if cx > 0 {
                    cx -= 1;
                } else {
                    break;
                }
            } else {
                t_max_z += t_delta_z;
                if step_z > 0 {
                    cz += 1;
                } else if cz > 0 {
                    cz -= 1;
                } else {
                    break;
                }
            }
        }
        None
    }
}
/// LOD level: determines cell sampling resolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LodLevel {
    /// Full resolution.
    High,
    /// 2x coarser.
    Medium,
    /// 4x coarser.
    Low,
    /// 8x coarser.
    VeryLow,
}
impl LodLevel {
    /// Stride (step size) in cells for this LOD level.
    pub fn stride(&self) -> usize {
        match self {
            LodLevel::High => 1,
            LodLevel::Medium => 2,
            LodLevel::Low => 4,
            LodLevel::VeryLow => 8,
        }
    }
    /// Select LOD level based on distance from camera/object.
    pub fn from_distance(dist: f64) -> Self {
        if dist < 20.0 {
            LodLevel::High
        } else if dist < 60.0 {
            LodLevel::Medium
        } else if dist < 150.0 {
            LodLevel::Low
        } else {
            LodLevel::VeryLow
        }
    }
}
/// Decompose a heightfield cell into two triangles (lower-left and upper-right)
/// and select the correct one for a given fractional coordinate.
pub struct TerrainTriangle;
impl TerrainTriangle {
    /// Get the two triangles for cell (ix, iz) in world coordinates.
    pub fn cell_triangles(
        hf: &HeightfieldShape,
        ix: usize,
        iz: usize,
    ) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
        let x0 = ix as f64 * hf.scale_x;
        let z0 = iz as f64 * hf.scale_z;
        let x1 = (ix + 1) as f64 * hf.scale_x;
        let z1 = (iz + 1) as f64 * hf.scale_z;
        let y00 = hf.height_at_cell(ix, iz) + hf.y_offset;
        let y10 = hf.height_at_cell(ix + 1, iz) + hf.y_offset;
        let y01 = hf.height_at_cell(ix, iz + 1) + hf.y_offset;
        let y11 = hf.height_at_cell(ix + 1, iz + 1) + hf.y_offset;
        let t0 = [[x0, y00, z0], [x1, y10, z0], [x0, y01, z1]];
        let t1 = [[x1, y10, z0], [x1, y11, z1], [x0, y01, z1]];
        (t0, t1)
    }
    /// Select which triangle a point (tx, tz) ∈ \[0,1\]² falls in.
    /// Returns the lower-left triangle if tx + tz <= 1.
    pub fn select_triangle(
        hf: &HeightfieldShape,
        ix: usize,
        iz: usize,
        tx: f64,
        tz: f64,
    ) -> [[f64; 3]; 3] {
        let (t0, t1) = Self::cell_triangles(hf, ix, iz);
        if tx + tz <= 1.0 { t0 } else { t1 }
    }
}
/// A multi-layer terrain model with per-cell material assignment.
///
/// Supports up to 4 material layers that can be queried for per-layer
/// physical properties such as friction, restitution, and deformability.
pub struct MultiLayerTerrain {
    /// Underlying heightfield.
    pub hf: HeightfieldShape,
    /// Per-cell top-layer material assignment (row-major).
    pub layers: Vec<TerrainLayer>,
}
impl MultiLayerTerrain {
    /// Construct from a heightfield, defaulting all cells to Rock.
    pub fn new(hf: HeightfieldShape) -> Self {
        let n = hf.nx * hf.nz;
        let layers = vec![TerrainLayer::Rock; n];
        MultiLayerTerrain { hf, layers }
    }
    fn cell_index(&self, ix: usize, iz: usize) -> usize {
        iz * self.hf.nx + ix
    }
    /// Set the material layer for cell (ix, iz).
    pub fn set_layer(&mut self, ix: usize, iz: usize, layer: TerrainLayer) {
        let i = self.cell_index(ix, iz);
        self.layers[i] = layer;
    }
    /// Get the material layer at cell (ix, iz).
    pub fn layer_at_cell(&self, ix: usize, iz: usize) -> &TerrainLayer {
        &self.layers[self.cell_index(ix, iz)]
    }
    /// Get the material layer at world-space (x, z) using nearest-cell lookup.
    pub fn layer_at(&self, x: f64, z: f64) -> &TerrainLayer {
        let ix = ((x / self.hf.scale_x).round() as usize).min(self.hf.nx - 1);
        let iz = ((z / self.hf.scale_z).round() as usize).min(self.hf.nz - 1);
        self.layer_at_cell(ix, iz)
    }
    /// Friction coefficient for a given layer.
    pub fn layer_friction(layer: &TerrainLayer) -> f64 {
        match layer {
            TerrainLayer::Rock => 0.7,
            TerrainLayer::Soil => 0.5,
            TerrainLayer::Vegetation => 0.4,
            TerrainLayer::Sand => 0.35,
        }
    }
    /// Restitution coefficient for a given layer.
    pub fn layer_restitution(layer: &TerrainLayer) -> f64 {
        match layer {
            TerrainLayer::Rock => 0.3,
            TerrainLayer::Soil => 0.1,
            TerrainLayer::Vegetation => 0.05,
            TerrainLayer::Sand => 0.05,
        }
    }
    /// Effective friction at world-space (x, z).
    pub fn friction_at(&self, x: f64, z: f64) -> f64 {
        Self::layer_friction(self.layer_at(x, z))
    }
    /// Count cells with a given layer type.
    pub fn count_layer(&self, target: &TerrainLayer) -> usize {
        self.layers.iter().filter(|l| *l == target).count()
    }
}
/// Continuous collision detection for a fast-moving sphere vs terrain.
pub struct TerrainCcd<'a> {
    /// Reference heightfield.
    pub hf: &'a HeightfieldShape,
    /// Number of substeps for linear search.
    pub steps: usize,
}
impl<'a> TerrainCcd<'a> {
    /// Construct with a heightfield and number of CCD substeps.
    pub fn new(hf: &'a HeightfieldShape, steps: usize) -> Self {
        TerrainCcd { hf, steps }
    }
    /// CCD for a sphere moving linearly from `c0` to `c1`.
    pub fn sphere_ccd(&self, c0: [f64; 3], c1: [f64; 3], radius: f64) -> Option<TerrainCcdResult> {
        let dc = sub3(c1, c0);
        for i in 0..=self.steps {
            let t = i as f64 / self.steps as f64;
            let c = add3(c0, scale3(dc, t));
            let h = self.hf.height_at(c[0], c[2]);
            if c[1] - radius < h {
                let normal = self.hf.normal_at(c[0], c[2]);
                return Some(TerrainCcdResult {
                    toi: t,
                    point: [c[0], h, c[2]],
                    normal,
                });
            }
        }
        None
    }
}
/// Erosion-aware collision with variable friction based on terrain material.
///
/// Models eroded vs. hard rock surfaces, wet soil, and gravel beds,
/// each with distinct friction coefficients and collision response.
pub struct ErosionFriction {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub nz: usize,
    /// Per-cell hardness in \[0, 1\]: 1.0 = solid rock, 0.0 = soft soil.
    pub hardness: Vec<f64>,
    /// Per-cell moisture in \[0, 1\]: 1.0 = saturated (low friction).
    pub moisture: Vec<f64>,
}
impl ErosionFriction {
    /// Create an erosion friction map.
    pub fn new(nx: usize, nz: usize) -> Self {
        ErosionFriction {
            nx,
            nz,
            hardness: vec![1.0; nx * nz],
            moisture: vec![0.0; nx * nz],
        }
    }
    fn index(&self, ix: usize, iz: usize) -> usize {
        iz * self.nx + ix
    }
    /// Set hardness at cell (ix, iz).
    pub fn set_hardness(&mut self, ix: usize, iz: usize, h: f64) {
        let i = self.index(ix, iz);
        self.hardness[i] = h.clamp(0.0, 1.0);
    }
    /// Set moisture at cell (ix, iz).
    pub fn set_moisture(&mut self, ix: usize, iz: usize, m: f64) {
        let i = self.index(ix, iz);
        self.moisture[i] = m.clamp(0.0, 1.0);
    }
    /// Compute the effective friction coefficient at cell (ix, iz).
    ///
    /// Harder and drier surfaces have higher friction.
    /// Formula: μ = μ_rock * hardness * (1 - 0.7 * moisture)
    pub fn friction_at(&self, ix: usize, iz: usize) -> f64 {
        let mu_rock = 0.8_f64;
        let i = self.index(ix, iz);
        let h = self.hardness[i];
        let m = self.moisture[i];
        mu_rock * h * (1.0 - 0.7 * m)
    }
    /// Bilinear interpolation of friction at world-space (x, z).
    pub fn friction_interp(&self, x: f64, z: f64, scale_x: f64, scale_z: f64) -> f64 {
        let cx = (x / scale_x).max(0.0);
        let cz = (z / scale_z).max(0.0);
        let ix = (cx as usize).min(self.nx.saturating_sub(2));
        let iz = (cz as usize).min(self.nz.saturating_sub(2));
        let tx = cx - ix as f64;
        let tz = cz - iz as f64;
        let f00 = self.friction_at(ix, iz);
        let f10 = self.friction_at((ix + 1).min(self.nx - 1), iz);
        let f01 = self.friction_at(ix, (iz + 1).min(self.nz - 1));
        let f11 = self.friction_at((ix + 1).min(self.nx - 1), (iz + 1).min(self.nz - 1));
        bilinear_height_interp(f00, f10, f01, f11, tx, tz)
    }
    /// Apply erosion: reduce hardness proportional to impact force (impulse).
    pub fn erode(&mut self, ix: usize, iz: usize, impulse: f64) {
        let i = self.index(ix, iz);
        self.hardness[i] = (self.hardness[i] - impulse * 0.01).max(0.0);
        self.moisture[i] = (self.moisture[i] + impulse * 0.005).min(1.0);
    }
}
