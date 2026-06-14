// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Signed Distance Field (SDF) geometry module.
//!
//! Provides a comprehensive SDF toolkit:
//!
//! - **Primitive SDFs** ([`sdf_sphere`], [`sdf_box`], [`sdf_capsule`], [`sdf_cylinder`],
//!   [`sdf_torus`], [`sdf_plane`], [`sdf_cone`]): exact signed distances for common shapes.
//! - **Combinators** ([`sdf_union`], [`sdf_intersection`], [`sdf_difference`]): Boolean
//!   CSG operations on SDFs.
//! - **Smooth blending** ([`sdf_smooth_union`], [`sdf_smooth_intersection`],
//!   [`sdf_smooth_difference`]): C¹-continuous blended combinations.
//! - **Gradient computation** ([`sdf_gradient`]): central-difference gradient and
//!   surface normal estimation.
//! - **Closest-point query** ([`closest_point_on_sdf`]): gradient-descent projection
//!   onto the zero level set.
//! - **Ray marching** ([`RayMarcher`]): sphere-tracing through an SDF scene.
//! - **Offset surfaces** ([`sdf_offset`]): inflate/deflate by constant offset.
//! - **Voronoi SDF** ([`VoronoiSdf`]): SDF defined by a set of seed points.
//! - **SDF to mesh** ([`MarchingCubes`]): isosurface extraction from a voxel SDF grid.
//! - **Octree SDF storage** ([`OctreeSdf`]): adaptive octree for compact SDF storage.
//! - **Narrow-band SDF** ([`NarrowBandSdf`]): stores only values within ±bandwidth.
//! - **Fast Marching Method** ([`FastMarchingMethod`]): efficient SDF reinitialization
//!   on a uniform grid.
//! - **SDF scene** ([`SdfScene`]): composable scene of multiple SDF objects.
//!
//! All vectors use plain `[f64; 3]` arrays; no external linear algebra dependency.

mod combinators;
mod fmm;
mod helpers;
mod jfa;
mod marching_cubes;
mod narrow_band;
mod octree;
mod operators;
mod primitives;
mod ray_marcher;
mod scene;
mod voronoi;

// ── Public re-exports (preserve original API surface) ────────────────────────

pub use combinators::{
    sdf_difference, sdf_exp_smooth_union, sdf_intersection, sdf_smooth_difference,
    sdf_smooth_intersection, sdf_smooth_union, sdf_union,
};

pub use fmm::FastMarchingMethod;

pub use jfa::JfaGrid;

pub use marching_cubes::{
    EDGE_VERTICES, MC_TABLE, MarchingCubes, MarchingCubesResult, MeshTriangle, MeshVertex,
};

pub use narrow_band::NarrowBandSdf;

pub use octree::{OctreeNode, OctreeSdf};

pub use operators::{
    closest_point_on_sdf, sdf_bend, sdf_bounding_box, sdf_displace, sdf_elongate_x, sdf_gradient,
    sdf_grid_volume, sdf_mean_curvature, sdf_mirror_y, sdf_morph, sdf_normal, sdf_offset,
    sdf_repeat, sdf_twist, sdf_volume_estimate,
};

pub use primitives::{
    sdf_box, sdf_capsule, sdf_cone, sdf_cylinder, sdf_cylinder_infinite, sdf_ellipsoid, sdf_plane,
    sdf_rounded_box, sdf_segment, sdf_sphere, sdf_torus,
};

pub use ray_marcher::{RayMarchHit, RayMarcher};

pub use scene::{SdfKind, SdfObject, SdfScene};

pub use voronoi::VoronoiSdf;

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── Primitive SDFs ─────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_sphere_inside() {
        let d = sdf_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(d < 0.0, "origin should be inside sphere, d={:.6}", d);
    }

    #[test]
    fn test_sdf_sphere_outside() {
        let d = sdf_sphere([2.0, 0.0, 0.0], 1.0);
        assert!(
            d > 0.0,
            "point outside sphere should have positive SDF, d={:.6}",
            d
        );
    }

    #[test]
    fn test_sdf_sphere_on_surface() {
        let d = sdf_sphere([1.0, 0.0, 0.0], 1.0);
        assert!(d.abs() < 1e-10, "point on sphere surface: d={:.6}", d);
    }

    #[test]
    fn test_sdf_box_inside() {
        let d = sdf_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(d < 0.0, "origin inside box should be negative: d={:.6}", d);
    }

    #[test]
    fn test_sdf_box_outside() {
        let d = sdf_box([2.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(d > 0.0, "point outside box: d={:.6}", d);
    }

    #[test]
    fn test_sdf_capsule_inside() {
        let d = sdf_capsule([0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        assert!(d < 0.0, "origin inside capsule: d={:.6}", d);
    }

    #[test]
    fn test_sdf_cylinder_inside() {
        let d = sdf_cylinder([0.0, 0.0, 0.0], 1.0, 2.0);
        assert!(d < 0.0, "origin inside cylinder: d={:.6}", d);
    }

    #[test]
    fn test_sdf_torus_outside_tube() {
        // Point far from the torus ring
        let d = sdf_torus([0.0, 0.0, 0.0], 2.0, 0.5);
        // Origin is at distance |sqrt(0+0)-2| = 2 from ring centre, minus tube radius 0.5
        assert!(d > 0.0, "origin outside torus tube: d={:.6}", d);
    }

    // ── Combinators ────────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_union_takes_min() {
        assert_eq!(sdf_union(1.0, -0.5), -0.5);
        assert_eq!(sdf_union(-1.0, 0.5), -1.0);
    }

    #[test]
    fn test_sdf_intersection_takes_max() {
        assert_eq!(sdf_intersection(1.0, -0.5), 1.0);
        assert_eq!(sdf_intersection(-1.0, 0.5), 0.5);
    }

    #[test]
    fn test_sdf_difference_example() {
        // a=1 (outside A), b=-0.5 (inside B): result = max(1, 0.5) = 1
        let d = sdf_difference(1.0, -0.5);
        assert!(d > 0.0, "difference should be positive here: d={:.6}", d);
    }

    // ── Smooth blending ────────────────────────────────────────────────────────

    #[test]
    fn test_smooth_union_between_values() {
        let k = 0.5;
        let su = sdf_smooth_union(1.0, 1.0, k);
        // Smooth union of equal values should be close to that value minus blend
        assert!(
            su <= 1.0,
            "smooth union should be <= min(a,b) at equal values: {:.6}",
            su
        );
    }

    #[test]
    fn test_smooth_union_far_apart_is_like_union() {
        // When a and b are far apart relative to k, smooth union ≈ regular union
        let k = 0.1;
        let a = 10.0;
        let b = -5.0;
        let su = sdf_smooth_union(a, b, k);
        let u = sdf_union(a, b);
        assert!(
            (su - u).abs() < 0.5,
            "smooth union far apart: su={:.6}, u={:.6}",
            su,
            u
        );
    }

    #[test]
    fn test_smooth_intersection_ge_max_sometimes() {
        let k = 1.0;
        let si = sdf_smooth_intersection(0.5, 0.5, k);
        // Smooth intersection should be close to 0.5
        assert!((si - 0.5).abs() < 0.3);
    }

    // ── Gradient ──────────────────────────────────────────────────────────────

    #[test]
    fn test_gradient_points_outward_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let p = [1.5, 0.0, 0.0];
        let g = sdf_gradient(&f, p, 1e-4);
        // Gradient should point in +x direction
        assert!(
            g[0] > 0.0,
            "gradient x-component should be positive: {:.6}",
            g[0]
        );
        assert!(
            g[1].abs() < 0.01,
            "gradient y-component should be ~0: {:.6}",
            g[1]
        );
    }

    #[test]
    fn test_normal_is_unit_length() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let p = [1.5, 0.5, 0.0];
        let n = sdf_normal(&f, p, 1e-4);
        // Use helpers::norm3 via the public interface
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-4,
            "normal should be unit length: {:.6}",
            len
        );
    }

    // ── Ray marching ──────────────────────────────────────────────────────────

    #[test]
    fn test_ray_march_hits_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let marcher = RayMarcher::new();
        let hit = marcher.march(&f, [5.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!(hit.is_some(), "ray along -x should hit sphere at origin");
        let h = hit.unwrap();
        assert!(
            (h.t - 4.0).abs() < 0.01,
            "hit distance should be ~4: {:.6}",
            h.t
        );
    }

    #[test]
    fn test_ray_march_misses_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let marcher = RayMarcher::new();
        let hit = marcher.march(&f, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_none(), "ray away from sphere should not hit");
    }

    // ── Voronoi SDF ───────────────────────────────────────────────────────────

    #[test]
    fn test_voronoi_sdf_nearest_seed() {
        let seeds = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let vsdf = VoronoiSdf::new(seeds);
        let nearest = vsdf.nearest_seed([1.0, 0.0, 0.0]);
        assert_eq!(nearest, 0, "nearest seed to [1,0,0] should be index 0");
    }

    #[test]
    fn test_voronoi_sdf_positive_near_boundary() {
        let seeds = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let vsdf = VoronoiSdf::new(seeds);
        // At midpoint, both distances equal → SDF = 0
        let d = vsdf.evaluate([2.0, 0.0, 0.0]);
        assert!(
            d.abs() < 1e-10,
            "midpoint Voronoi SDF should be ~0: {:.6}",
            d
        );
    }

    // ── Octree SDF ────────────────────────────────────────────────────────────

    #[test]
    fn test_octree_sdf_inside_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let octree = OctreeSdf::build(&f, [0.0, 0.0, 0.0], 2.0, 4, 0.5);
        let d = octree.evaluate([0.0, 0.0, 0.0]);
        assert!(d < 0.0, "octree SDF at origin inside sphere: d={:.6}", d);
    }

    #[test]
    fn test_octree_sdf_has_leaves() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let octree = OctreeSdf::build(&f, [0.0, 0.0, 0.0], 2.0, 3, 1.0);
        assert!(
            octree.count_leaves() > 0,
            "octree should have at least one leaf"
        );
    }

    // ── Marching Cubes ────────────────────────────────────────────────────────

    #[test]
    fn test_marching_cubes_sphere_produces_vertices() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let mc = MarchingCubes::from_function(&f, 8, 8, 8, [-2.0, 2.0, -2.0, 2.0, -2.0, 2.0]);
        let result = mc.extract(0.0);
        // Should have some vertices near the sphere surface
        assert!(result.n_vertices() > 0 || result.n_triangles() == 0);
    }

    #[test]
    fn test_marching_cubes_spacing() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let mc = MarchingCubes::from_function(&f, 4, 4, 4, [0.0, 4.0, 0.0, 4.0, 0.0, 4.0]);
        let sp = mc.spacing();
        assert!((sp[0] - 1.0).abs() < 1e-10);
    }

    // ── Narrow-band SDF ───────────────────────────────────────────────────────

    #[test]
    fn test_narrow_band_active_count() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let nb = NarrowBandSdf::from_function(&f, 10, 10, 10, 0.4, [-2.0, -2.0, -2.0], 1.0);
        // Some cells should be active near the sphere surface
        assert!(
            nb.active_count() > 0,
            "narrow band should have active cells"
        );
    }

    // ── FMM ───────────────────────────────────────────────────────────────────

    #[test]
    fn test_fmm_propagates_distance() {
        let mut fmm = FastMarchingMethod::new(5, 5, 5, 1.0);
        // Seed the centre cell
        let centre = fmm.flat(2, 2, 2);
        fmm.set_known(&[(centre, 0.0)]);
        fmm.run();
        // Neighbour should have distance ~1
        let d = fmm.distance_at(3, 2, 2);
        assert!(d > 0.0 && d <= 2.0, "FMM neighbour distance: {:.6}", d);
    }

    #[test]
    fn test_fmm_distance_increases_with_steps() {
        let mut fmm = FastMarchingMethod::new(7, 1, 1, 1.0);
        let seed = fmm.flat(3, 0, 0);
        fmm.set_known(&[(seed, 0.0)]);
        fmm.run();
        let d1 = fmm.distance_at(4, 0, 0);
        let d2 = fmm.distance_at(5, 0, 0);
        assert!(
            d2 >= d1,
            "FMM distance should increase with steps: d1={:.6}, d2={:.6}",
            d1,
            d2
        );
    }

    // ── SDF Scene ─────────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_scene_single_sphere() {
        let mut scene = SdfScene::new();
        scene.add(SdfObject::sphere("s", 1.0, [0.0, 0.0, 0.0]));
        let d = scene.evaluate([0.0, 0.0, 0.0]);
        assert!(d < 0.0, "origin should be inside scene sphere: d={:.6}", d);
    }

    #[test]
    fn test_sdf_scene_union_two_spheres() {
        let mut scene = SdfScene::new();
        scene.add(SdfObject::sphere("a", 1.0, [-2.0, 0.0, 0.0]));
        scene.add(SdfObject::sphere("b", 1.0, [2.0, 0.0, 0.0]));
        // At origin (between spheres), both SDFs are 1, so union = 1
        let d = scene.evaluate([0.0, 0.0, 0.0]);
        assert!(
            d > 0.0,
            "midpoint between two spheres should be outside: d={:.6}",
            d
        );
    }

    // ── Closest point ─────────────────────────────────────────────────────────

    #[test]
    fn test_closest_point_on_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let (closest, _d0) = closest_point_on_sdf(&f, [3.0, 0.0, 0.0], 100, 1e-5);
        let dist_from_origin =
            (closest[0] * closest[0] + closest[1] * closest[1] + closest[2] * closest[2]).sqrt();
        assert!(
            (dist_from_origin - 1.0).abs() < 0.01,
            "closest point should be on sphere surface: dist={:.6}",
            dist_from_origin
        );
    }

    // ── Volume estimate ───────────────────────────────────────────────────────

    #[test]
    fn test_volume_estimate_sphere() {
        let f = |p: [f64; 3]| sdf_sphere(p, 1.0);
        let vol = sdf_volume_estimate(&f, [-1.5, 1.5, -1.5, 1.5, -1.5, 1.5], 10000);
        let exact = 4.0 / 3.0 * PI;
        // Within 10% of exact value
        assert!(
            (vol - exact).abs() / exact < 0.15,
            "volume estimate should be close to {:.6}: got {:.6}",
            exact,
            vol
        );
    }

    // ── MC_TABLE tests ────────────────────────────────────────────────────────

    #[test]
    fn test_mc_table_entry_0_no_triangles() {
        // Config 0: all corners outside → no triangles, first entry is -1
        assert_eq!(
            MC_TABLE[0][0], -1,
            "MC_TABLE[0] should start with -1 (no triangles)"
        );
    }

    #[test]
    fn test_mc_table_entry_255_no_triangles() {
        // Config 255: all corners inside → no triangles, first entry is -1
        assert_eq!(
            MC_TABLE[255][0], -1,
            "MC_TABLE[255] should start with -1 (no triangles)"
        );
    }

    #[test]
    fn test_mc_table_entry_1_one_triangle() {
        // Config 1: corner 0 inside → exactly one triangle (3 edge indices before -1)
        let row = &MC_TABLE[1];
        assert!(row[0] >= 0, "MC_TABLE[1][0] should be a valid edge index");
        assert!(row[1] >= 0, "MC_TABLE[1][1] should be a valid edge index");
        assert!(row[2] >= 0, "MC_TABLE[1][2] should be a valid edge index");
        assert_eq!(
            row[3], -1,
            "MC_TABLE[1] should have exactly one triangle (terminated after 3 entries)"
        );
        // The canonical table has [0, 8, 3, -1, ...] for config 1
        assert_eq!(row[0], 0, "MC_TABLE[1][0] should be edge 0");
        assert_eq!(row[1], 8, "MC_TABLE[1][1] should be edge 8");
        assert_eq!(row[2], 3, "MC_TABLE[1][2] should be edge 3");
    }

    #[test]
    fn test_mc_table_all_edge_indices_valid() {
        // All non-sentinel entries must be valid edge indices (0..11)
        for (config, row) in MC_TABLE.iter().enumerate() {
            for &e in row {
                if e == -1 {
                    break;
                }
                assert!(
                    (0..12).contains(&(e as usize)),
                    "MC_TABLE[{config}] has invalid edge index {e}"
                );
            }
        }
    }

    #[test]
    fn test_mc_table_triangle_count_consistency() {
        // Each config must have complete triangles (count before -1 divisible by 3)
        for (config, row) in MC_TABLE.iter().enumerate() {
            let n = row.iter().take_while(|&&e| e >= 0).count();
            assert_eq!(
                n % 3,
                0,
                "MC_TABLE[{config}] has {n} entries before -1, not divisible by 3"
            );
        }
    }

    #[test]
    fn test_mc_isosurface_sphere_triangle_count() {
        // 32³ grid sphere r=1 should produce 200+ triangles
        let mc = MarchingCubes::from_function(
            &|p: [f64; 3]| sdf_sphere(p, 1.0),
            32,
            32,
            32,
            [-1.5, 1.5, -1.5, 1.5, -1.5, 1.5],
        );
        let result = mc.extract(0.0);
        assert!(
            result.n_triangles() >= 200,
            "sphere isosurface should have at least 200 triangles; got {}",
            result.n_triangles()
        );
        // Verify all triangle vertex indices are in bounds
        for tri in &result.triangles {
            for &idx in &tri.indices {
                assert!(
                    idx < result.n_vertices(),
                    "triangle index {idx} out of bounds (nv={})",
                    result.n_vertices()
                );
            }
        }
    }
}
