// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geometric shape types for collision detection and physics simulation.
//!
//! Provides primitives (sphere, box, capsule, cylinder, cone) and complex
//! shapes (convex hull, triangle mesh, compound, height field) with support
//! point queries, volume/inertia computation, and ray casting.
#![warn(missing_docs)]
#![allow(ambiguous_glob_reexports)]
#![allow(dead_code)]

mod error;
pub use error::*;

pub mod shape;
pub use shape::{RayHit, Shape};

pub mod sphere;
pub use sphere::Sphere;

pub mod box_shape;
pub use box_shape::BoxShape;

pub mod capsule;
pub use capsule::Capsule;

pub mod cylinder;
pub use cylinder::Cylinder;

pub mod cone;
pub use cone::Cone;

pub mod convex_hull;
pub use convex_hull::*;

pub mod triangle_mesh;
pub use triangle_mesh::TriangleMesh;

pub mod heightfield;
pub use heightfield::HeightField;

pub mod compound;
pub use compound::Compound;

pub mod torus;
pub use torus::Torus;

pub mod mesh_ops;
pub use mesh_ops::TriMesh;

pub mod quickhull;
pub use quickhull::ConvexHull3DVec;

pub mod swept;

pub mod csg;

pub mod parametric;

pub mod mesh_repair;

pub mod spatial_hash;

pub mod voronoi;

pub mod decimation;

pub mod offset_surface;
pub use offset_surface::*;

pub mod point_cloud;
pub use point_cloud::*;

pub mod mesh_quality;

pub mod remesh;

pub mod boolean_ops;

pub mod geodesic;

pub mod subdivision;

pub mod mesh_param;
pub use mesh_param::{
    LscmParameterization, ParamTriMesh, TutteParameterization, boundary_vertices, build_halfedge,
    fill_holes, is_boundary_vertex, loop_subdivision, midpoint_subdivision,
    remove_duplicate_vertices, texture_distortion, uv_overlap_check, uv_stretch, vertex_neighbors,
};

pub use voronoi::{
    DelaunayTriangle, DelaunayTriangulation, LegacyVoronoiCell, Point2D, VoronoiCell,
    VoronoiDiagram, VoronoiSite, bowyer_watson, circumcircle, delaunay_to_voronoi, in_circumcircle,
    is_delaunay, lloyd_relaxation, nearest_site, power_diagram, voronoi_area,
};

#[cfg(test)]
mod prop_tests {
    use super::*;
    use oxiphysics_core::math::Vec3;
    use proptest::prelude::*;

    fn positive_f64() -> impl Strategy<Value = f64> {
        0.01_f64..50.0_f64
    }

    fn coord_f64() -> impl Strategy<Value = f64> {
        -50.0_f64..50.0_f64
    }

    proptest! {
        #[test]
        fn prop_sphere_bbox_contains_center(radius in positive_f64()) {
            let sphere = Sphere::new(radius);
            let bb = sphere.bounding_box();
            prop_assert!(
                bb.contains_point(&Vec3::zeros()),
                "sphere bbox does not contain center for radius={}", radius
            );
        }

        #[test]
        fn prop_box_support_within_extents(
            hx in positive_f64(),
            hy in positive_f64(),
            hz in positive_f64(),
            dx in coord_f64(),
            dy in coord_f64(),
            dz in coord_f64(),
        ) {
            let half_extents = Vec3::new(hx, hy, hz);
            let b = BoxShape::new(half_extents);
            let dir = Vec3::new(dx, dy, dz);
            let sp = b.support_point(&dir);
            let eps = 1e-9;
            prop_assert!(sp.x.abs() <= hx + eps, "|sp.x|={} > hx={}", sp.x.abs(), hx);
            prop_assert!(sp.y.abs() <= hy + eps, "|sp.y|={} > hy={}", sp.y.abs(), hy);
            prop_assert!(sp.z.abs() <= hz + eps, "|sp.z|={} > hz={}", sp.z.abs(), hz);
        }

        #[test]
        fn prop_sphere_volume_positive(radius in positive_f64()) {
            let sphere = Sphere::new(radius);
            prop_assert!(sphere.volume() > 0.0, "sphere volume not positive for radius={}", radius);
        }

        #[test]
        fn prop_sphere_inertia_symmetric(
            radius in positive_f64(),
            mass in positive_f64(),
        ) {
            let sphere = Sphere::new(radius);
            let it = sphere.inertia_tensor(mass);
            for i in 0..3 {
                for j in 0..3 {
                    let diff = (it[(i, j)] - it[(j, i)]).abs();
                    prop_assert!(diff < 1e-9, "inertia not symmetric at ({},{})={} vs ({},{})={}", i, j, it[(i,j)], j, i, it[(j,i)]);
                }
            }
        }

        #[test]
        fn prop_sphere_raycast_hit_on_surface(
            radius in positive_f64(),
            // Ray from far away along X axis
            offset in 1.01_f64..50.0_f64,
        ) {
            let sphere = Sphere::new(radius);
            let origin = Vec3::new(-(offset + radius + 1.0), 0.0, 0.0);
            let dir = Vec3::new(1.0, 0.0, 0.0);
            if let Some(hit) = sphere.ray_cast(&origin, &dir, 1000.0) {
                let dist_from_center = hit.point.norm();
                prop_assert!(
                    (dist_from_center - radius).abs() < 1e-6,
                    "hit point not on sphere surface: dist={}, radius={}", dist_from_center, radius
                );
            }
        }

        #[test]
        fn prop_capsule_volume_ge_sphere(
            radius in positive_f64(),
            half_height in positive_f64(),
        ) {
            let capsule = Capsule::new(radius, half_height);
            let sphere = Sphere::new(radius);
            prop_assert!(
                capsule.volume() >= sphere.volume() - 1e-9,
                "capsule volume {} < sphere volume {}", capsule.volume(), sphere.volume()
            );
        }
    }
}
pub mod bspline;
pub mod implicit;
pub mod implicit_surfaces;
pub mod level_set;
pub mod mesh_processing;
pub mod origami;
pub mod sphere_packing;
pub mod topology;

pub mod architectural_geometry;

pub mod fractal_geometry;

pub mod computational_geometry;

pub mod implicit_geometry;

pub mod robot_geometry;

pub mod medical_geometry;

pub mod procedural_geometry;

pub mod spline_geometry;

pub mod discrete_geometry;

pub mod terrain_processing;

pub mod cell_complex;
pub mod convex_decomposition;
pub mod voxel_grid;
pub use voxel_grid::VoxelGrid;
pub mod vhacd;
pub use vhacd::{VHacdConfig, VHacdError, VHacdVoxel};
pub mod geodesic_geometry;
pub mod medial_axis;
pub mod mesh_boolean;
pub mod mesh_repair_ext;
pub mod mesh_simplification;
pub mod nurbs_geometry;
pub mod offset_geometry;
pub mod signed_distance_field;
pub mod topology_geometry;
