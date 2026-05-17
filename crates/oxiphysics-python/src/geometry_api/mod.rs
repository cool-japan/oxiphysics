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
pub use shapes::{PyConvexHull, PyShape, PyShapeWrapper};
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
    use super::shapes::{PyConvexHull, PyShape};
    use super::spatial::{PyAabb, PyBvh, PyHeightfield, PySpatialHash};
    use super::transform::{
        PyGeometryTransform, PyMeshQuality, PyPointCloud, convex_hull_from_mesh,
        mesh_to_obj_string, obj_string_to_mesh,
    };
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

    #[test]
    fn test_convex_hull_from_points_empty() {
        let hull = PyConvexHull::from_points_internal(&[]);
        assert!(hull.vertices.is_empty());
    }

    #[test]
    fn test_convex_hull_from_points_basic() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let hull = PyConvexHull::from_points_internal(&pts);
        assert!(!hull.vertices.is_empty());
        assert!(!hull.faces.is_empty());
    }

    #[test]
    fn test_convex_hull_volume_positive() {
        let pts: Vec<[f64; 3]> = vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ];
        let hull = PyConvexHull::from_points_internal(&pts);
        assert!(hull.volume() > 0.0);
    }

    #[test]
    fn test_convex_hull_gjk_overlap() {
        let pts1 = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let pts2 = vec![
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
        ];
        let h1 = PyConvexHull::from_points_internal(&pts1);
        let h2 = PyConvexHull::from_points_internal(&pts2);
        assert!(h1.gjk_overlap(&h2));
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
        let h1 = PyConvexHull::from_points_internal(&pts1);
        let h2 = PyConvexHull::from_points_internal(&pts2);
        assert!(!h1.gjk_overlap(&h2));
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

    #[test]
    fn test_csg_union_vertex_count() {
        let a = unit_box_mesh();
        let b = unit_box_mesh();
        let u = PyCsg::union(&a, &b);
        assert_eq!(u.vertices.len(), a.vertices.len() + b.vertices.len());
    }

    #[test]
    fn test_csg_intersection_subset() {
        let a = unit_box_mesh();
        let mut b_big = unit_box_mesh();
        for v in &mut b_big.vertices {
            v[0] *= 2.0;
            v[1] *= 2.0;
            v[2] *= 2.0;
            v[0] -= 0.5;
            v[1] -= 0.5;
            v[2] -= 0.5;
        }
        let inter = PyCsg::intersection(&a, &b_big);
        assert!(inter.vertices.len() <= a.vertices.len() + 8);
    }

    #[test]
    fn test_csg_subtraction_removes_some() {
        let a = unit_box_mesh();
        let b = unit_box_mesh();
        let sub = PyCsg::subtraction(&a, &b);
        assert!(sub.triangle_count() <= a.triangle_count());
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
        let hull = convex_hull_from_mesh(&mesh);
        assert!(!hull.vertices.is_empty());
    }
}
