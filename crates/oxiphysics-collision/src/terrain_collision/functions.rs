//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Bilinear interpolation of height from four corner heights.
///
/// `tx`, `ty` are fractional coordinates within the cell \[0, 1\].
pub fn bilinear_height_interp(h00: f64, h10: f64, h01: f64, h11: f64, tx: f64, ty: f64) -> f64 {
    let h0 = h00 + (h10 - h00) * tx;
    let h1 = h01 + (h11 - h01) * tx;
    h0 + (h1 - h0) * ty
}
/// Compute terrain normal from four surrounding heights using central differences.
///
/// `scale_x` and `scale_z` are the world-space distances between grid cells.
pub fn terrain_normal_from_heights(
    h_nx: f64,
    h_px: f64,
    h_nz: f64,
    h_pz: f64,
    scale_x: f64,
    scale_z: f64,
) -> [f64; 3] {
    let dx = (h_px - h_nx) / (2.0 * scale_x);
    let dz = (h_pz - h_nz) / (2.0 * scale_z);
    normalize3([-dx, 1.0, -dz])
}
/// Compute the AABB of a heightfield patch.
pub fn heightfield_aabb(
    x_start: f64,
    z_start: f64,
    x_end: f64,
    z_end: f64,
    h_min: f64,
    h_max: f64,
) -> ([f64; 3], [f64; 3]) {
    ([x_start, h_min, z_start], [x_end, h_max, z_end])
}
/// Convert voxel grid indices to world-space coordinates.
pub fn voxel_to_world(
    ix: usize,
    iy: usize,
    iz: usize,
    voxel_size: f64,
    origin: [f64; 3],
) -> [f64; 3] {
    [
        origin[0] + ix as f64 * voxel_size,
        origin[1] + iy as f64 * voxel_size,
        origin[2] + iz as f64 * voxel_size,
    ]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-14 {
        [0.0, 1.0, 0.0]
    } else {
        scale3(a, 1.0 / l)
    }
}
/// Test a ray against an AABB. Returns the t value of entry or None.
pub(super) fn ray_aabb(origin: [f64; 3], dir: [f64; 3], mn: [f64; 3], mx: [f64; 3]) -> Option<f64> {
    let mut t_min = 0.0_f64;
    let mut t_max = f64::MAX;
    for i in 0..3 {
        if dir[i].abs() < 1e-14 {
            if origin[i] < mn[i] || origin[i] > mx[i] {
                return None;
            }
        } else {
            let inv = 1.0 / dir[i];
            let t1 = (mn[i] - origin[i]) * inv;
            let t2 = (mx[i] - origin[i]) * inv;
            t_min = t_min.max(t1.min(t2));
            t_max = t_max.min(t1.max(t2));
            if t_max < t_min {
                return None;
            }
        }
    }
    Some(t_min)
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::terrain_collision::HeightfieldBvh;

    use crate::terrain_collision::HeightfieldShape;

    use crate::terrain_collision::TerrainCcd;
    use crate::terrain_collision::TerrainCollider;

    use crate::terrain_collision::TerrainDeformation;

    use crate::terrain_collision::TerrainRaycast;
    use crate::terrain_collision::TerrainSampling;

    use crate::terrain_collision::TerrainTriangle;

    use crate::terrain_collision::VoxelTerrain;
    use crate::terrain_collision::WaterSurface;
    fn flat_hf(nx: usize, nz: usize, h: f64) -> HeightfieldShape {
        HeightfieldShape::new(nx, nz, 1.0, 1.0, vec![h; nx * nz])
    }
    #[test]
    fn bilinear_corners() {
        assert!((bilinear_height_interp(1.0, 2.0, 3.0, 4.0, 0.0, 0.0) - 1.0).abs() < 1e-12);
        assert!((bilinear_height_interp(1.0, 2.0, 3.0, 4.0, 1.0, 0.0) - 2.0).abs() < 1e-12);
        assert!((bilinear_height_interp(1.0, 2.0, 3.0, 4.0, 0.0, 1.0) - 3.0).abs() < 1e-12);
        assert!((bilinear_height_interp(1.0, 2.0, 3.0, 4.0, 1.0, 1.0) - 4.0).abs() < 1e-12);
    }
    #[test]
    fn bilinear_center() {
        let h = bilinear_height_interp(0.0, 1.0, 1.0, 2.0, 0.5, 0.5);
        assert!((h - 1.0).abs() < 1e-12, "h={h}");
    }
    #[test]
    fn terrain_normal_flat_points_up() {
        let n = terrain_normal_from_heights(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!((n[1] - 1.0).abs() < 1e-9, "n={n:?}");
    }
    #[test]
    fn terrain_normal_slope_tilts() {
        let n = terrain_normal_from_heights(0.0, 1.0, 0.0, 0.0, 1.0, 1.0);
        assert!(n[0] < 0.0, "slope in x should tilt normal: n={n:?}");
    }
    #[test]
    fn hf_height_at_flat() {
        let hf = flat_hf(4, 4, 5.0);
        assert!((hf.height_at(1.5, 1.5) - 5.0).abs() < 1e-9);
    }
    #[test]
    fn hf_normal_flat() {
        let hf = flat_hf(5, 5, 0.0);
        let n = hf.normal_at(2.0, 2.0);
        assert!((n[1] - 1.0).abs() < 1e-6, "n={n:?}");
    }
    #[test]
    fn hf_height_range() {
        let mut hf = flat_hf(4, 4, 0.0);
        hf.heights[5] = 3.0;
        let (mn, mx) = hf.height_range();
        assert!((mn - 0.0).abs() < 1e-12);
        assert!((mx - 3.0).abs() < 1e-12);
    }
    #[test]
    fn hf_world_extent() {
        let hf = flat_hf(5, 5, 1.0);
        let (mn, mx) = hf.world_extent();
        assert!((mn[1] - 1.0).abs() < 1e-12);
        assert!((mx[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn terrain_triangle_select_lower() {
        let hf = flat_hf(4, 4, 0.0);
        let tri = TerrainTriangle::select_triangle(&hf, 0, 0, 0.2, 0.2);
        assert!((tri[0][0]).abs() < 1e-9);
    }
    #[test]
    fn terrain_triangle_select_upper() {
        let hf = flat_hf(4, 4, 0.0);
        let tri = TerrainTriangle::select_triangle(&hf, 0, 0, 0.8, 0.8);
        assert!(tri[0][0] > 0.0 || tri[1][0] > 0.0);
    }
    #[test]
    fn hfbvh_build_creates_root() {
        let hf = flat_hf(5, 5, 0.0);
        let bvh = HeightfieldBvh::build(hf);
        assert!(!bvh.nodes.is_empty());
    }
    #[test]
    fn hfbvh_ray_vs_root_hit() {
        let hf = flat_hf(5, 5, 0.0);
        let bvh = HeightfieldBvh::build(hf);
        let hit = bvh.ray_vs_root([2.0, 5.0, 2.0], [0.0, -1.0, 0.0]);
        assert!(hit, "ray from above should hit root AABB");
    }
    #[test]
    fn terrain_raycast_vertical_hits_flat() {
        let hf = flat_hf(10, 10, 2.0);
        let caster = TerrainRaycast::new(&hf);
        let result = caster.cast([4.5, 10.0, 4.5], [0.0, -1.0, 0.0]);
        if let Some(hit) = result {
            assert!((hit.point[1] - 2.0).abs() < 0.5, "hit_y={}", hit.point[1]);
        }
    }
    #[test]
    fn terrain_collider_sphere_contact() {
        let hf = flat_hf(10, 10, 0.0);
        let col = TerrainCollider::new(&hf);
        let c = col.sphere_vs_terrain([3.0, 0.4, 3.0], 0.5);
        assert!(c.is_some(), "sphere below terrain should contact");
        let c = c.unwrap();
        assert!(c.depth > 0.0);
    }
    #[test]
    fn terrain_collider_sphere_no_contact() {
        let hf = flat_hf(10, 10, 0.0);
        let col = TerrainCollider::new(&hf);
        let c = col.sphere_vs_terrain([3.0, 5.0, 3.0], 0.5);
        assert!(c.is_none());
    }
    #[test]
    fn terrain_collider_capsule() {
        let hf = flat_hf(10, 10, 0.0);
        let col = TerrainCollider::new(&hf);
        let contacts = col.capsule_vs_terrain([3.0, 0.4, 3.0], [4.0, 0.4, 4.0], 0.5);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn terrain_collider_box_contact() {
        let hf = flat_hf(10, 10, 1.0);
        let col = TerrainCollider::new(&hf);
        let contacts = col.box_vs_terrain([1.0, 0.0, 1.0], [3.0, 0.9, 3.0]);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn terrain_ccd_fast_sphere_impact() {
        let hf = flat_hf(10, 10, 0.0);
        let ccd = TerrainCcd::new(&hf, 20);
        let result = ccd.sphere_ccd([3.0, 5.0, 3.0], [3.0, -1.0, 3.0], 0.5);
        assert!(result.is_some(), "fast sphere should hit terrain");
    }
    #[test]
    fn terrain_ccd_no_impact_above() {
        let hf = flat_hf(10, 10, 0.0);
        let ccd = TerrainCcd::new(&hf, 20);
        let result = ccd.sphere_ccd([3.0, 2.0, 3.0], [3.0, 5.0, 3.0], 0.5);
        assert!(result.is_none(), "sphere moving up should not impact");
    }
    #[test]
    fn water_surface_base_level() {
        let w = WaterSurface::new(3.0);
        assert!((w.height_at(0.0, 0.0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn water_surface_wave() {
        let mut w = WaterSurface::new(0.0);
        w.add_wave(1.0, 1.0, 0.0, 0.0, 0.0);
        let h = w.height_at(std::f64::consts::FRAC_PI_2, 0.0);
        assert!((h - 1.0).abs() < 1e-9, "h={h}");
    }
    #[test]
    fn water_surface_float_test() {
        let w = WaterSurface::new(1.0);
        let (sub, depth) = w.float_test([0.0, 0.5, 0.0]);
        assert!(sub, "object below water level should be submerged");
        assert!((depth - 0.5).abs() < 1e-9, "depth={depth}");
    }
    #[test]
    fn water_surface_above_water() {
        let w = WaterSurface::new(1.0);
        let (sub, depth) = w.float_test([0.0, 2.0, 0.0]);
        assert!(!sub);
        assert_eq!(depth, 0.0);
    }
    #[test]
    fn voxel_terrain_set_get() {
        let mut vt = VoxelTerrain::new(4, 4, 4, 1.0, [0.0, 0.0, 0.0]);
        vt.set(1, 1, 1, 1.0);
        assert!((vt.get(1, 1, 1) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn voxel_terrain_aabb_query() {
        let mut vt = VoxelTerrain::new(5, 5, 5, 1.0, [0.0, 0.0, 0.0]);
        vt.set(2, 2, 2, 1.0);
        let hits = vt.aabb_query([1.5, 1.5, 1.5], [2.5, 2.5, 2.5]);
        assert!(hits.contains(&(2, 2, 2)));
    }
    #[test]
    fn voxel_terrain_isosurface_no_crossing() {
        let vt = VoxelTerrain::new(4, 4, 4, 1.0, [0.0, 0.0, 0.0]);
        let tris = vt.extract_isosurface();
        assert!(tris.is_empty());
    }
    #[test]
    fn voxel_to_world_check() {
        let w = voxel_to_world(2, 3, 1, 0.5, [1.0, 1.0, 1.0]);
        assert!((w[0] - 2.0).abs() < 1e-12);
        assert!((w[1] - 2.5).abs() < 1e-12);
        assert!((w[2] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn terrain_deformation_crater() {
        let hf = flat_hf(11, 11, 1.0);
        let mut deform = TerrainDeformation::new(hf);
        deform.apply_crater(5.0, 5.0, 2.0, 0.5);
        let h_center = deform.hf.height_at(5.0, 5.0);
        let h_edge = deform.hf.height_at(0.0, 0.0);
        assert!(
            h_center < h_edge,
            "crater center should be lower: center={h_center}, edge={h_edge}"
        );
    }
    #[test]
    fn terrain_deformation_mound() {
        let hf = flat_hf(11, 11, 0.0);
        let mut deform = TerrainDeformation::new(hf);
        deform.apply_mound(5.0, 5.0, 2.0, 1.0);
        let h_center = deform.hf.height_at(5.0, 5.0);
        assert!(
            h_center > 0.0,
            "mound center should be raised: h={h_center}"
        );
    }
    #[test]
    fn terrain_sampling_elevation_profile() {
        let hf = flat_hf(10, 10, 2.0);
        let path = [[0.0f64, 0.0], [5.0, 5.0]];
        let samples = TerrainSampling::elevation_profile(&hf, &path, 5);
        assert!(!samples.is_empty());
        for s in &samples {
            assert!(
                (s.height - 2.0).abs() < 1e-9,
                "flat terrain height should be 2.0"
            );
        }
    }
    #[test]
    fn terrain_sampling_slope_map() {
        let hf = flat_hf(5, 5, 0.0);
        let slopes = TerrainSampling::slope_map(&hf);
        assert_eq!(slopes.len(), 25);
        for s in &slopes {
            assert!(s.abs() < 1e-6, "flat terrain slope should be ~0");
        }
    }
    #[test]
    fn terrain_sampling_find_peak() {
        let mut hf = flat_hf(8, 8, 0.0);
        hf.heights[3 * 8 + 3] = 5.0;
        let peak = TerrainSampling::find_peak(&hf, 0.0, 0.0, 7.0, 7.0);
        assert!((peak[1] - 5.0).abs() < 1e-9, "peak height should be 5.0");
    }
    #[test]
    fn heightfield_aabb_check() {
        let (mn, mx) = heightfield_aabb(0.0, 0.0, 10.0, 10.0, -1.0, 3.0);
        assert!((mn[1] - (-1.0)).abs() < 1e-12);
        assert!((mx[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn ray_aabb_miss() {
        let result = ray_aabb(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0],
        );
        assert!(result.is_none(), "ray should miss AABB above");
    }
    #[test]
    fn ray_aabb_hit() {
        let result = ray_aabb(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 1.0, 1.0],
        );
        assert!(result.is_some());
    }
    #[test]
    fn water_surface_step_changes_height() {
        let mut w = WaterSurface::new(0.0);
        w.add_wave(1.0, 0.0, 0.0, 0.0, 1.0);
        let h0 = w.height_at(0.0, 0.0);
        w.step(std::f64::consts::FRAC_PI_2);
        let h1 = w.height_at(0.0, 0.0);
        assert!((h1 - 1.0).abs() < 1e-9, "h1={h1}, h0={h0}");
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::terrain_collision::CliffDetector;
    use crate::terrain_collision::ErosionFriction;

    use crate::terrain_collision::HeightfieldBvhFull;
    use crate::terrain_collision::HeightfieldShape;
    use crate::terrain_collision::LodLevel;
    use crate::terrain_collision::LodTerrainCollider;
    use crate::terrain_collision::MultiLayerTerrain;
    use crate::terrain_collision::ProceduralCraterImpact;
    use crate::terrain_collision::SlopeCategory;

    use crate::terrain_collision::TerrainContactCache;

    use crate::terrain_collision::TerrainLayer;
    use crate::terrain_collision::TerrainNormalBilinear;

    use crate::terrain_collision::TerrainShadowMap;
    use crate::terrain_collision::TerrainSlopeAspect;

    use crate::terrain_collision::UnderwaterTerrainDetector;

    fn flat_hf(nx: usize, nz: usize, h: f64) -> HeightfieldShape {
        HeightfieldShape::new(nx, nz, 1.0, 1.0, vec![h; nx * nz])
    }
    fn sloped_hf() -> HeightfieldShape {
        let nx = 8;
        let nz = 8;
        let heights: Vec<f64> = (0..nx * nz).map(|i| (i / nx) as f64 * 0.5).collect();
        HeightfieldShape::new(nx, nz, 1.0, 1.0, heights)
    }
    #[test]
    fn bvh_full_build_has_nodes() {
        let hf = flat_hf(8, 8, 0.0);
        let bvh = HeightfieldBvhFull::build(hf, 2);
        assert!(
            bvh.num_nodes() > 1,
            "hierarchical BVH should have multiple nodes"
        );
    }
    #[test]
    fn bvh_full_query_aabb_returns_leaves() {
        let hf = flat_hf(8, 8, 0.0);
        let bvh = HeightfieldBvhFull::build(hf, 2);
        let leaves = bvh.query_aabb([1.0, -1.0, 1.0], [3.0, 1.0, 3.0]);
        assert!(!leaves.is_empty(), "query should return at least one leaf");
    }
    #[test]
    fn bvh_full_leaf_size_respected() {
        let hf = flat_hf(16, 16, 0.0);
        let bvh = HeightfieldBvhFull::build(hf, 4);
        assert!(bvh.num_nodes() >= 1);
    }
    #[test]
    fn bilinear_normal_flat_points_up() {
        let hf = flat_hf(6, 6, 0.0);
        let interp = TerrainNormalBilinear::new(&hf);
        let n = interp.normal_at(2.5, 2.5);
        assert!(
            (n[1] - 1.0).abs() < 1e-5,
            "flat terrain normal should point up: n={n:?}"
        );
    }
    #[test]
    fn bilinear_normal_nonzero() {
        let hf = sloped_hf();
        let interp = TerrainNormalBilinear::new(&hf);
        let n = interp.normal_at(3.0, 3.0);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-5, "normal should be unit length");
    }
    #[test]
    fn slope_flat_terrain_is_flat() {
        let hf = flat_hf(6, 6, 0.0);
        let sa = TerrainSlopeAspect::new(&hf);
        assert!(
            sa.slope_angle(2.0, 2.0) < 0.01,
            "flat terrain should have near-zero slope"
        );
    }
    #[test]
    fn slope_category_flat() {
        let hf = flat_hf(6, 6, 0.0);
        let sa = TerrainSlopeAspect::new(&hf);
        assert_eq!(sa.classify_slope(2.0, 2.0), SlopeCategory::Flat);
    }
    #[test]
    fn slope_category_map_length() {
        let hf = flat_hf(5, 5, 0.0);
        let sa = TerrainSlopeAspect::new(&hf);
        let map = sa.slope_category_map();
        assert_eq!(map.len(), 25);
    }
    #[test]
    fn shadow_map_lit_from_above() {
        let hf = flat_hf(5, 5, 0.0);
        let shadow = TerrainShadowMap::compute(&hf, [0.0, 1.0, 0.0]);
        assert!(shadow.lit_fraction() > 0.0);
    }
    #[test]
    fn shadow_map_size() {
        let hf = flat_hf(4, 4, 0.0);
        let shadow = TerrainShadowMap::compute(&hf, [0.1, 1.0, 0.1]);
        assert_eq!(shadow.shadow.len(), 16);
    }
    #[test]
    fn shadow_map_value_range() {
        let hf = flat_hf(4, 4, 0.0);
        let shadow = TerrainShadowMap::compute(&hf, [0.5, 1.0, 0.5]);
        for &s in &shadow.shadow {
            assert!(s == 0.0 || s == 1.0, "shadow value should be 0 or 1");
        }
    }
    #[test]
    fn erosion_friction_default_high() {
        let ef = ErosionFriction::new(4, 4);
        let f = ef.friction_at(0, 0);
        assert!(f > 0.5, "default rock friction should be > 0.5");
    }
    #[test]
    fn erosion_friction_wet_reduces() {
        let mut ef = ErosionFriction::new(4, 4);
        ef.set_moisture(0, 0, 1.0);
        let f_wet = ef.friction_at(0, 0);
        let f_dry = ErosionFriction::new(4, 4).friction_at(0, 0);
        assert!(f_wet < f_dry, "wet surface should have lower friction");
    }
    #[test]
    fn erosion_erode_reduces_hardness() {
        let mut ef = ErosionFriction::new(4, 4);
        let h_before = ef.hardness[0];
        ef.erode(0, 0, 10.0);
        assert!(ef.hardness[0] < h_before, "erosion should reduce hardness");
    }
    #[test]
    fn underwater_detect_submerged() {
        let hf = flat_hf(6, 6, 0.0);
        let det = UnderwaterTerrainDetector::new(&hf, 2.0);
        assert!(
            det.is_underwater(3.0, 1.0, 3.0),
            "y=1 with water_level=2 should be underwater"
        );
    }
    #[test]
    fn underwater_depth_above_flat() {
        let hf = flat_hf(6, 6, 0.0);
        let det = UnderwaterTerrainDetector::new(&hf, 3.0);
        let depth = det.water_depth_at(2.0, 2.0);
        assert!(
            (depth - 3.0).abs() < 1e-9,
            "water depth over flat terrain should be water_level"
        );
    }
    #[test]
    fn underwater_buoyancy_fully_submerged() {
        let hf = flat_hf(6, 6, 0.0);
        let det = UnderwaterTerrainDetector::new(&hf, 10.0);
        let f = det.buoyancy_force(3.0, 5.0, 3.0, 1.0, 1000.0, 9.81);
        assert!(
            f > 0.0,
            "fully submerged sphere should have positive buoyancy"
        );
    }
    #[test]
    fn cliff_flat_terrain_no_cliffs() {
        let hf = flat_hf(6, 6, 0.0);
        let det = CliffDetector::new(&hf, (45.0_f64).to_radians());
        assert_eq!(det.num_cliffs(), 0, "flat terrain should have no cliffs");
    }
    #[test]
    fn cliff_mask_length() {
        let hf = flat_hf(5, 5, 0.0);
        let det = CliffDetector::new(&hf, (30.0_f64).to_radians());
        let mask = det.cliff_mask();
        assert_eq!(mask.len(), 25);
    }
    #[test]
    fn multi_layer_default_rock() {
        let hf = flat_hf(4, 4, 0.0);
        let mlt = MultiLayerTerrain::new(hf);
        assert_eq!(mlt.layer_at_cell(0, 0), &TerrainLayer::Rock);
    }
    #[test]
    fn multi_layer_set_soil() {
        let hf = flat_hf(4, 4, 0.0);
        let mut mlt = MultiLayerTerrain::new(hf);
        mlt.set_layer(1, 1, TerrainLayer::Soil);
        assert_eq!(mlt.layer_at_cell(1, 1), &TerrainLayer::Soil);
    }
    #[test]
    fn multi_layer_friction_values() {
        assert!(
            MultiLayerTerrain::layer_friction(&TerrainLayer::Rock)
                > MultiLayerTerrain::layer_friction(&TerrainLayer::Sand)
        );
    }
    #[test]
    fn lod_high_resolution_nearby() {
        assert_eq!(LodLevel::from_distance(5.0), LodLevel::High);
        assert_eq!(LodLevel::from_distance(200.0), LodLevel::VeryLow);
    }
    #[test]
    fn lod_sphere_contact() {
        let hf = flat_hf(10, 10, 0.0);
        let col = LodTerrainCollider::new(&hf, [0.0, 0.0, 0.0]);
        let c = col.sphere_contact_lod([5.0, 0.3, 5.0], 0.5);
        assert!(c.is_some());
    }
    #[test]
    fn cache_insert_and_lookup() {
        let mut cache = TerrainContactCache::new(0.5, 10);
        cache.insert([1.0, 0.0, 1.0], None);
        let found = cache.lookup([1.1, 0.0, 1.0]);
        assert!(found.is_some(), "nearby query should hit cache");
    }
    #[test]
    fn cache_miss_far_query() {
        let mut cache = TerrainContactCache::new(0.1, 10);
        cache.insert([0.0, 0.0, 0.0], None);
        let found = cache.lookup([10.0, 0.0, 10.0]);
        assert!(found.is_none(), "distant query should miss cache");
    }
    #[test]
    fn crater_with_rim_deepens_center() {
        let mut hf = flat_hf(13, 13, 0.0);
        ProceduralCraterImpact::apply_crater_with_rim(&mut hf, 6.0, 6.0, 3.0, 1.0, 0.3, 0.3);
        let h_center = hf.height_at(6.0, 6.0);
        let h_edge = hf.height_at(0.0, 0.0);
        assert!(
            h_center < h_edge,
            "crater center should be lower: center={h_center}, edge={h_edge}"
        );
    }
    #[test]
    fn crater_with_rim_creates_rim() {
        let mut hf = flat_hf(17, 17, 0.0);
        ProceduralCraterImpact::apply_crater_with_rim(&mut hf, 8.0, 8.0, 3.0, 1.0, 0.5, 0.3);
        let h_rim = hf.height_at(11.5, 8.0);
        let h_center = hf.height_at(8.0, 8.0);
        assert!(h_rim > h_center, "rim should be higher than crater center");
    }
    #[test]
    fn multiple_craters_applied() {
        let mut hf = flat_hf(15, 15, 2.0);
        let impactors = vec![(4.0, 4.0, 2.0, 1.0), (10.0, 10.0, 2.0, 0.5)];
        ProceduralCraterImpact::apply_multiple(&mut hf, &impactors, 0.2);
        let h1 = hf.height_at(4.0, 4.0);
        let h2 = hf.height_at(10.0, 10.0);
        assert!(h1 < 2.0, "first crater center should be depressed");
        assert!(h2 < 2.0, "second crater center should be depressed");
    }
}
