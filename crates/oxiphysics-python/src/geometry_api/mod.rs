// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Python geometry bindings.
//!
//! Provides Python-friendly geometry primitives, mesh operations, spatial data
//! structures, and CSG operations using plain `f64` arrays — no nalgebra.

pub mod mesh;
pub mod shapes;
pub mod spatial;
pub mod transform;

// Re-export everything needed by lib.rs and other crate modules
pub use mesh::{PyCsg, PyTriangleMesh, compute_aabb};
pub use shapes::{ConvexHullError, PyConvexHull, PyShape, PyShapeWrapper};
pub use spatial::{PyAabb, PyBvh, PyBvhNode, PyHeightfield, PySpatialHash};
pub use transform::{
    PyGeometryTransform, PyMeshQuality, PyPointCloud, convex_hull_from_mesh, mesh_to_obj_string,
    obj_string_to_mesh, py_compute_aabb, py_convex_hull_from_mesh, py_mesh_to_obj_string,
    py_obj_string_to_mesh,
};

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

/// Register all `geometry` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_geometry_module(
    parent: &pyo3::Bound<'_, pyo3::types::PyModule>,
) -> pyo3::PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = pyo3::types::PyModule::new(parent.py(), "geometry")?;
    child.add_class::<PyShapeWrapper>()?;
    child.add_class::<PyConvexHull>()?;
    child.add_class::<PyTriangleMesh>()?;
    child.add_class::<PyCsg>()?;
    child.add_class::<PyHeightfield>()?;
    child.add_class::<PySpatialHash>()?;
    child.add_class::<PyAabb>()?;
    child.add_class::<PyBvhNode>()?;
    child.add_class::<PyBvh>()?;
    child.add_class::<PyPointCloud>()?;
    child.add_class::<PyGeometryTransform>()?;
    child.add_class::<PyMeshQuality>()?;
    child.add_function(pyo3::wrap_pyfunction!(py_compute_aabb, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(py_mesh_to_obj_string, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(py_obj_string_to_mesh, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(py_convex_hull_from_mesh, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::mesh::{PyTriangleMesh, compute_aabb};
    use super::shapes::{ConvexHullError, PyConvexHull, PyShape};
    use super::spatial::{PyAabb, PyBvh, PyHeightfield, PySpatialHash};
    use super::transform::{
        PyGeometryTransform, PyMeshQuality, PyPointCloud, convex_hull_from_mesh,
        mesh_to_obj_string, obj_string_to_mesh,
    };
    use oxiphysics::geometry::convex_hull::{ConvexHull3D, gjk_intersect};
    use std::f64::consts::PI;

    use super::mesh::PyCsg;

    // --- Helper ---

    fn unit_box_mesh() -> PyTriangleMesh {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let idx = vec![
            0, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6, 0, 5, 1, 0, 4, 5, 2, 7, 3, 2, 6, 7, 0, 3, 7, 0, 7,
            4, 1, 5, 6, 1, 6, 2,
        ];
        PyTriangleMesh::from_raw_internal(verts, idx)
    }

    // --- PyShape tests ---

    #[test]
    fn test_sphere_volume() {
        let s = PyShape::Sphere(1.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((s.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_box_volume() {
        let b = PyShape::Box([1.0, 2.0, 3.0]);
        assert!((b.volume() - 48.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_volume() {
        let c = PyShape::Cylinder {
            radius: 1.0,
            half_height: 1.0,
        };
        assert!((c.volume() - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_torus_volume() {
        let t = PyShape::Torus {
            major_r: 3.0,
            minor_r: 1.0,
        };
        let expected = 2.0 * PI * PI * 3.0 * 1.0;
        assert!((t.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_shape_aabb_sphere() {
        let (mn, mx) = PyShape::Sphere(2.0).aabb();
        assert_eq!(mn, [-2.0, -2.0, -2.0]);
        assert_eq!(mx, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_capsule_surface_area() {
        let c = PyShape::Capsule {
            radius: 1.0,
            half_height: 1.0,
        };
        let expected = 4.0 * PI + 2.0 * PI * 2.0;
        assert!((c.surface_area() - expected).abs() < 1e-10);
    }

    // --- PyConvexHull tests ---

    /// Unit right tetrahedron with its right-angle corner at the origin.
    fn unit_tetra() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }

    /// Axis-aligned cube `[lo, hi]^3` as eight corner points.
    fn cube_corners(lo: f64, hi: f64) -> Vec<[f64; 3]> {
        let mut pts = Vec::new();
        for &x in &[lo, hi] {
            for &y in &[lo, hi] {
                for &z in &[lo, hi] {
                    pts.push([x, y, z]);
                }
            }
        }
        pts
    }

    fn translated(pts: &[[f64; 3]], offset: [f64; 3]) -> Vec<[f64; 3]> {
        pts.iter()
            .map(|p| [p[0] + offset[0], p[1] + offset[1], p[2] + offset[2]])
            .collect()
    }

    #[test]
    fn test_convex_hull_from_points_empty_is_an_error() {
        assert_eq!(
            PyConvexHull::from_points_internal(&[]).err(),
            Some(ConvexHullError::TooFewPoints(0))
        );
    }

    #[test]
    fn test_convex_hull_from_points_basic() {
        let hull = PyConvexHull::from_points_internal(&unit_tetra()).expect("tetrahedron hull");
        assert_eq!(hull.vertices.len(), 4);
        assert_eq!(hull.faces.len(), 12);
    }

    #[test]
    fn test_convex_hull_tetrahedron_volume_is_one_sixth_not_its_bounding_box() {
        let hull = PyConvexHull::from_points_internal(&unit_tetra()).expect("tetrahedron hull");
        // The bounding box of this tetrahedron has volume 1; the tetrahedron 1/6.
        assert!(
            (hull.volume() - 1.0 / 6.0).abs() < 1e-12,
            "volume = {}",
            hull.volume()
        );
        let expected_area = 3.0 * 0.5 + 0.5 * 3.0_f64.sqrt();
        assert!(
            (hull.surface_area() - expected_area).abs() < 1e-12,
            "area = {}",
            hull.surface_area()
        );
    }

    #[test]
    fn test_convex_hull_volume_and_area_of_cube() {
        let hull = PyConvexHull::from_points_internal(&cube_corners(-1.0, 1.0)).expect("cube hull");
        assert!(
            (hull.volume() - 8.0).abs() < 1e-12,
            "volume = {}",
            hull.volume()
        );
        assert!(
            (hull.surface_area() - 24.0).abs() < 1e-12,
            "area = {}",
            hull.surface_area()
        );
    }

    #[test]
    fn test_convex_hull_drops_interior_points() {
        let mut pts = cube_corners(0.0, 2.0);
        pts.push([1.0, 1.0, 1.0]);
        pts.push([0.5, 1.5, 1.0]);
        pts.push([1.9, 0.1, 1.0]);
        let hull = PyConvexHull::from_points_internal(&pts).expect("cube hull");
        assert_eq!(
            hull.vertices.len(),
            8,
            "interior points are not hull vertices"
        );
        assert!((hull.volume() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn test_convex_hull_volume_is_translation_invariant() {
        let far = translated(&unit_tetra(), [1.0e6, -2.0e6, 3.0e6]);
        let hull = PyConvexHull::from_points_internal(&far).expect("tetrahedron hull");
        assert!(
            (hull.volume() - 1.0 / 6.0).abs() < 1e-6,
            "volume = {}",
            hull.volume()
        );
    }

    #[test]
    fn test_convex_hull_rejects_degenerate_input() {
        let coplanar = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        ];
        assert_eq!(
            PyConvexHull::from_points_internal(&coplanar).err(),
            Some(ConvexHullError::Degenerate)
        );
        let too_few = &unit_tetra()[..3];
        assert_eq!(
            PyConvexHull::from_points_internal(too_few).err(),
            Some(ConvexHullError::TooFewPoints(3))
        );
        let mut with_nan = unit_tetra();
        with_nan[2][1] = f64::NAN;
        assert_eq!(
            PyConvexHull::from_points_internal(&with_nan).err(),
            Some(ConvexHullError::NonFinitePoint(2))
        );
    }

    #[test]
    fn test_convex_hull_from_points_python_entry_point() {
        let rows: Vec<Vec<f64>> = unit_tetra().iter().map(|p| p.to_vec()).collect();
        let hull = PyConvexHull::from_points(rows).expect("tetrahedron hull");
        assert!((hull.volume() - 1.0 / 6.0).abs() < 1e-12);
        // A row with fewer than three coordinates is an error, not a silently dropped point.
        let mut short = unit_tetra().iter().map(|p| p.to_vec()).collect::<Vec<_>>();
        short[1] = vec![1.0, 0.0];
        assert!(PyConvexHull::from_points(short).is_err());
        // Fewer than four points cannot enclose a volume.
        assert!(PyConvexHull::from_points(vec![vec![0.0; 3]; 3]).is_err());
    }

    #[test]
    fn test_convex_hull_gjk_overlap() {
        // The second tetrahedron is the first one shifted by (0.25, 0.25, 0.25),
        // so they share the region around (0.3, 0.3, 0.3).
        let h1 = PyConvexHull::from_points_internal(&unit_tetra()).expect("hull 1");
        let h2 = PyConvexHull::from_points_internal(&translated(&unit_tetra(), [0.25; 3]))
            .expect("hull 2");
        assert!(h1.gjk_overlap(&h2));
        assert!(h2.gjk_overlap(&h1));
        assert!(h1.gjk_overlap(&h1));
    }

    #[test]
    fn test_convex_hull_gjk_no_overlap() {
        let pts1 = vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
        ];
        let pts2 = vec![
            [10.0, 10.0, 10.0],
            [10.1, 10.0, 10.0],
            [10.0, 10.1, 10.0],
            [10.0, 10.0, 10.1],
        ];
        let h1 = PyConvexHull::from_points_internal(&pts1).expect("hull 1");
        let h2 = PyConvexHull::from_points_internal(&pts2).expect("hull 2");
        assert!(!h1.gjk_overlap(&h2));
    }

    #[test]
    fn test_convex_hull_gjk_is_not_a_bounding_sphere_test() {
        // Tetrahedron 2 lies entirely in the half-space x + y + z >= 1.5 while
        // tetrahedron 1 lies in x + y + z <= 1, so the hulls are disjoint. Their
        // bounding spheres still overlap: the centroids are ~0.87 apart while
        // both bounding radii are ~0.83.
        let pts2 = vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
        ];
        let h1 = PyConvexHull::from_points_internal(&unit_tetra()).expect("hull 1");
        let h2 = PyConvexHull::from_points_internal(&pts2).expect("hull 2");
        assert!(!h1.gjk_overlap(&h2));
        assert!(!h2.gjk_overlap(&h1));
    }

    #[test]
    fn test_convex_hull_gjk_overlap_matches_the_engine_gjk_for_every_separation() {
        // Cubes [-1, 1]^3 and the same cube shifted by d along x overlap exactly
        // when d < 2. The bounding spheres (radius sqrt(3)) are disjoint for
        // d > 2 sqrt(3) ~ 3.46, so the sweep covers both the early rejection and
        // the GJK decision, which must agree with the engine's GJK on the hulls.
        let base = PyConvexHull::from_points_internal(&cube_corners(-1.0, 1.0)).expect("cube");
        let engine = |hull: &PyConvexHull| ConvexHull3D {
            vertices: hull.vertices.clone(),
            faces: hull
                .faces
                .chunks_exact(3)
                .map(|tri| [tri[0], tri[1], tri[2]])
                .collect(),
        };
        for shift in [
            0.0, 0.5, 1.0, 1.5, 1.9, 2.1, 2.5, 3.0, 3.4, 3.5, 4.0, 6.0, 50.0,
        ] {
            let moved = PyConvexHull::from_points_internal(&translated(
                &cube_corners(-1.0, 1.0),
                [shift, 0.0, 0.0],
            ))
            .expect("shifted cube");
            let expected = shift < 2.0;
            assert_eq!(
                gjk_intersect(&engine(&base), &engine(&moved)),
                expected,
                "engine GJK at shift {shift}"
            );
            assert_eq!(base.gjk_overlap(&moved), expected, "shift {shift}");
            assert_eq!(moved.gjk_overlap(&base), expected, "shift {shift}, swapped");
        }
    }

    #[test]
    fn test_convex_hull_gjk_detects_containment() {
        let outer =
            PyConvexHull::from_points_internal(&cube_corners(-2.0, 2.0)).expect("outer cube");
        let inner =
            PyConvexHull::from_points_internal(&cube_corners(-0.5, 0.5)).expect("inner cube");
        assert!(outer.gjk_overlap(&inner));
        assert!(inner.gjk_overlap(&outer));
    }

    #[test]
    fn test_convex_hull_gjk_separated_cubes_along_axis() {
        let a = PyConvexHull::from_points_internal(&cube_corners(0.0, 1.0)).expect("cube a");
        let near = PyConvexHull::from_points_internal(&translated(
            &cube_corners(0.0, 1.0),
            [0.9, 0.0, 0.0],
        ))
        .expect("overlapping cube");
        let far = PyConvexHull::from_points_internal(&translated(
            &cube_corners(0.0, 1.0),
            [1.1, 0.0, 0.0],
        ))
        .expect("separated cube");
        assert!(a.gjk_overlap(&near));
        assert!(!a.gjk_overlap(&far));
    }

    // --- PyTriangleMesh tests ---

    #[test]
    fn test_mesh_compute_normals_not_empty() {
        let mesh = unit_box_mesh();
        assert_eq!(mesh.normals.len(), mesh.vertices.len());
    }

    #[test]
    fn test_mesh_triangle_count() {
        let mesh = unit_box_mesh();
        assert_eq!(mesh.triangle_count(), 12);
    }

    #[test]
    fn test_mesh_surface_area_unit_box() {
        let mesh = unit_box_mesh();
        let sa = mesh.compute_surface_area();
        assert!((sa - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_mesh_volume_unit_box() {
        let mesh = unit_box_mesh();
        let vol = mesh.compute_volume();
        assert!((vol - 1.0).abs() < 0.1, "vol={vol}");
    }

    #[test]
    fn test_mesh_repair_removes_degenerate() {
        let mut mesh = PyTriangleMesh::from_raw_internal(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![0, 0, 0, 0, 1, 2],
        );
        mesh.repair();
        assert_eq!(mesh.triangle_count(), 1);
    }

    #[test]
    fn test_mesh_smooth_runs() {
        let mut mesh = unit_box_mesh();
        let before = mesh.vertices.clone();
        mesh.smooth(2);
        assert_eq!(mesh.vertices.len(), before.len());
    }

    // --- PyCsg tests ---

    /// Closed triangle mesh of the cube `[lo, hi]`, wound counter-clockwise seen
    /// from outside as the boolean operations require (signed volume `> 0`).
    fn outward_cube(lo: [f64; 3], hi: [f64; 3]) -> PyTriangleMesh {
        let verts = vec![
            [lo[0], lo[1], lo[2]],
            [hi[0], lo[1], lo[2]],
            [hi[0], hi[1], lo[2]],
            [lo[0], hi[1], lo[2]],
            [lo[0], lo[1], hi[2]],
            [hi[0], lo[1], hi[2]],
            [hi[0], hi[1], hi[2]],
            [lo[0], hi[1], hi[2]],
        ];
        let idx = vec![
            0, 2, 1, 0, 3, 2, // -z
            4, 5, 6, 4, 6, 7, // +z
            0, 1, 5, 0, 5, 4, // -y
            3, 6, 2, 3, 7, 6, // +y
            0, 4, 7, 0, 7, 3, // -x
            1, 2, 6, 1, 6, 5, // +x
        ];
        PyTriangleMesh::from_raw_internal(verts, idx)
    }

    #[test]
    fn test_csg_operations_return_the_exact_volumes_of_overlapping_cubes() {
        // B is A shifted by 0.5 along x: the overlap is a 0.5 x 1 x 1 slab whose
        // faces cut the cube faces of both operands.
        let a = outward_cube([0.0; 3], [1.0; 3]);
        let b = outward_cube([0.5, 0.0, 0.0], [1.5, 1.0, 1.0]);
        for (name, result, expected) in [
            ("union", PyCsg::union(&a, &b), 1.5),
            ("intersection", PyCsg::intersection(&a, &b), 0.5),
            ("A minus B", PyCsg::subtraction(&a, &b), 0.5),
            ("B minus A", PyCsg::subtraction(&b, &a), 0.5),
        ] {
            let volume = result.compute_volume();
            assert!(
                (volume - expected).abs() < 1e-9,
                "{name}: volume {volume}, expected {expected}"
            );
            assert!(
                result.triangle_count() >= 12,
                "{name} must be a closed surface"
            );
        }
    }

    #[test]
    fn test_csg_with_a_cube_inside_a_larger_cube() {
        let small = outward_cube([0.0; 3], [1.0; 3]);
        let big = outward_cube([-0.5; 3], [1.5; 3]);
        let volume = |mesh: &PyTriangleMesh| mesh.compute_volume();
        assert!((volume(&PyCsg::union(&small, &big)) - 8.0).abs() < 1e-9);
        assert!((volume(&PyCsg::intersection(&small, &big)) - 1.0).abs() < 1e-9);
        // The big cube with a cavity: 8 - 1.
        assert!((volume(&PyCsg::subtraction(&big, &small)) - 7.0).abs() < 1e-9);
        // Nothing of the small cube is left outside the big one.
        let nothing = PyCsg::subtraction(&small, &big);
        assert_eq!(nothing.triangle_count(), 0);
        assert!(volume(&nothing).abs() < 1e-12);
    }

    #[test]
    fn test_csg_with_disjoint_cubes() {
        let a = outward_cube([0.0; 3], [1.0; 3]);
        let far = outward_cube([5.0; 3], [6.0; 3]);
        let volume = |mesh: &PyTriangleMesh| mesh.compute_volume();
        assert!((volume(&PyCsg::union(&a, &far)) - 2.0).abs() < 1e-9);
        assert!((volume(&PyCsg::subtraction(&a, &far)) - 1.0).abs() < 1e-9);
        assert_eq!(PyCsg::intersection(&a, &far).triangle_count(), 0);
    }

    // --- PyHeightfield tests ---

    #[test]
    fn test_heightfield_flat() {
        let hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        assert!((hf.height_at(5.0, 5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_heightfield_set_and_query() {
        let mut hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        hf.set_height(3, 4, 5.0);
        let h = hf.height_at(3.0, 4.0);
        assert!((h - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_heightfield_normal_flat() {
        let hf = PyHeightfield::new(10, 10, 1.0, 1.0, 1.0);
        let n = hf.normal_at(5.0, 5.0);
        assert!((n[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_heightfield_raycast_hits() {
        let mut hf = PyHeightfield::new(20, 20, 1.0, 1.0, 1.0);
        for r in 0..20 {
            for c in 0..20 {
                hf.set_height(c, r, 1.0);
            }
        }
        let origin = vec![10.0, 10.0, 10.0];
        let dir = vec![0.0, -1.0, 0.0];
        let hit = hf.raycast(origin, dir);
        assert!(hit.is_some());
    }

    // --- PySpatialHash tests ---

    #[test]
    fn test_spatial_hash_insert_query() {
        let mut sh = PySpatialHash::new(1.0);
        let id = sh.insert(vec![0.5, 0.5, 0.5]);
        let res = sh.query(vec![0.5, 0.5, 0.5]);
        assert!(res.contains(&id));
    }

    #[test]
    fn test_spatial_hash_sphere_query() {
        let mut sh = PySpatialHash::new(1.0);
        sh.insert(vec![0.0, 0.0, 0.0]);
        sh.insert(vec![5.0, 5.0, 5.0]);
        let res = sh.sphere_query(vec![0.0, 0.0, 0.0], 1.0);
        assert_eq!(res.len(), 1);
    }

    #[test]
    fn test_spatial_hash_remove() {
        let mut sh = PySpatialHash::new(1.0);
        let id = sh.insert(vec![0.0, 0.0, 0.0]);
        sh.remove(id);
        let res = sh.sphere_query(vec![0.0, 0.0, 0.0], 0.5);
        assert!(!res.contains(&id));
    }

    // --- PyBvh tests ---

    #[test]
    fn test_bvh_build_empty() {
        let bvh = PyBvh::build(vec![]);
        assert!(bvh.nodes.is_empty());
    }

    #[test]
    fn test_bvh_query_aabb_hit() {
        let aabbs = vec![
            PyAabb::new_internal([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            PyAabb::new_internal([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let bvh = PyBvh::build(aabbs);
        let q = PyAabb::new_internal([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let hits = bvh.query_aabb(&q);
        assert!(hits.contains(&0));
        assert!(!hits.contains(&1));
    }

    #[test]
    fn test_bvh_raycast() {
        let aabbs = vec![PyAabb::new_internal([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])];
        let bvh = PyBvh::build(aabbs);
        let hits = bvh.raycast(vec![-5.0, 0.5, 0.5], vec![1.0, 0.0, 0.0]);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].1, 0);
    }

    #[test]
    fn test_bvh_nearest_neighbor() {
        let aabbs = vec![
            PyAabb::new_internal([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            PyAabb::new_internal([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]),
        ];
        let bvh = PyBvh::build(aabbs);
        let nn = bvh.nearest_neighbor(vec![0.5, 0.5, 0.5]);
        assert_eq!(nn, Some(0));
    }

    // --- PyPointCloud tests ---

    #[test]
    fn test_point_cloud_compute_normals() {
        let mut pc = PyPointCloud::from_points(vec![
            vec![1.0, 0.0, 0.0],
            vec![-1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, -1.0, 0.0],
        ]);
        pc.compute_normals(3);
        assert_eq!(pc.normals.len(), 4);
    }

    #[test]
    fn test_point_cloud_simplify() {
        let mut pc =
            PyPointCloud::from_points((0..100).map(|i| vec![i as f64 * 0.01, 0.0, 0.0]).collect());
        pc.simplify(0.1);
        assert!(pc.points.len() < 100);
    }

    #[test]
    fn test_point_cloud_poisson_stub() {
        let pc = PyPointCloud::from_points(vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]]);
        let mesh = pc.poisson_reconstruct();
        assert!(mesh.vertices.is_empty());
    }

    // --- PyGeometryTransform tests ---

    #[test]
    fn test_transform_identity() {
        let t = PyGeometryTransform::identity();
        let p = [1.0, 2.0, 3.0];
        let out = t.apply_point(p.to_vec());
        for i in 0..3 {
            assert!((out[i] - p[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_transform_translate() {
        let mut tr = PyGeometryTransform::identity();
        tr.translation = [1.0, 2.0, 3.0];
        let out = tr.apply_point(vec![0.0, 0.0, 0.0]);
        assert!((out[0] - 1.0).abs() < 1e-10);
        assert!((out[1] - 2.0).abs() < 1e-10);
        assert!((out[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_transform_scale() {
        let mut t = PyGeometryTransform::identity();
        t.scale = 2.0;
        let out = t.apply_point(vec![1.0, 1.0, 1.0]);
        for &out_i in &out {
            assert!((out_i - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_transform_rotate_180() {
        use super::shapes::normalize3;
        let mut t = PyGeometryTransform::identity();
        t.axis = normalize3([0.0, 0.0, 1.0]);
        t.angle = PI;
        let out = t.apply_point(vec![1.0, 0.0, 0.0]);
        assert!((out[0] - (-1.0)).abs() < 1e-10, "x={}", out[0]);
        assert!((out[1]).abs() < 1e-10, "y={}", out[1]);
    }

    #[test]
    fn test_transform_apply_to_mesh() {
        let mut mesh = unit_box_mesh();
        let mut t = PyGeometryTransform::identity();
        t.translation = [1.0, 0.0, 0.0];
        t.apply_to_mesh(&mut mesh);
        for v in &mesh.vertices {
            assert!(v[0] >= 1.0 - 1e-10 && v[0] <= 2.0 + 1e-10);
        }
    }

    // --- PyMeshQuality tests ---

    #[test]
    fn test_mesh_quality_compute() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        assert_eq!(q.aspect_ratios.len(), mesh.triangle_count());
        assert_eq!(q.areas.len(), mesh.triangle_count());
    }

    #[test]
    fn test_mesh_quality_areas_positive() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        for &a in &q.areas {
            assert!(a >= 0.0);
        }
    }

    #[test]
    fn test_mesh_quality_worst_elements() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        let worst = q.worst_elements(3);
        assert!(worst.len() <= 3);
    }

    #[test]
    fn test_mesh_quality_mean_aspect_ratio() {
        let mesh = unit_box_mesh();
        let q = PyMeshQuality::compute(&mesh);
        let mar = q.mean_aspect_ratio();
        assert!(mar >= 1.0);
    }

    // --- Helper function tests ---

    #[test]
    fn test_compute_aabb_empty() {
        let (mn, mx) = compute_aabb(&[]);
        assert_eq!(mn, [0.0; 3]);
        assert_eq!(mx, [0.0; 3]);
    }

    #[test]
    fn test_compute_aabb_points() {
        let pts = vec![[-1.0, 2.0, 3.0], [4.0, -5.0, 6.0]];
        let (mn, mx) = compute_aabb(&pts);
        assert_eq!(mn, [-1.0, -5.0, 3.0]);
        assert_eq!(mx, [4.0, 2.0, 6.0]);
    }

    #[test]
    fn test_obj_roundtrip() {
        let mesh = unit_box_mesh();
        let obj = mesh_to_obj_string(&mesh);
        let parsed = obj_string_to_mesh(&obj);
        assert_eq!(parsed.vertices.len(), mesh.vertices.len());
        assert_eq!(parsed.triangle_count(), mesh.triangle_count());
    }

    #[test]
    fn test_convex_hull_from_mesh_fn() {
        let mesh = unit_box_mesh();
        let hull = convex_hull_from_mesh(&mesh).expect("unit box hull");
        assert_eq!(hull.vertices.len(), 8);
        assert!(
            (hull.volume() - 1.0).abs() < 1e-12,
            "volume = {}",
            hull.volume()
        );
    }

    #[test]
    fn test_convex_hull_from_mesh_rejects_flat_mesh() {
        let flat = PyTriangleMesh::from_raw_internal(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            vec![0, 1, 2, 0, 2, 3],
        );
        assert_eq!(
            convex_hull_from_mesh(&flat).err(),
            Some(ConvexHullError::Degenerate)
        );
    }
}
