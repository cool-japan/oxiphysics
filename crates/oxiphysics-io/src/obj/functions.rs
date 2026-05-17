//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::items_after_test_module)]
#![allow(clippy::field_reassign_with_default)]

use super::types::{
    MeshInstance, MeshTransform, ObjFace, ObjGroup, ObjMaterial, ObjMesh, ObjMeshStats,
};

#[cfg(test)]
mod tests {

    use crate::obj::types::*;
    use oxiphysics_core::math::Vec3;
    fn make_default_face(vis: Vec<usize>) -> ObjFace {
        ObjFace {
            vertex_indices: vis,
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        }
    }
    #[test]
    fn test_obj_write_and_read_roundtrip() {
        let path = "/tmp/oxiphy_test.obj";
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        ObjWriter::write_legacy(path, &verts, &tris, None).unwrap();
        let (read_verts, read_tris) = ObjReader::read(path).unwrap();
        assert_eq!(read_verts.len(), 3);
        assert_eq!(read_tris.len(), 1);
        assert_eq!(read_tris[0], [0, 1, 2]);
        assert!((read_verts[1].x - 1.0).abs() < 1e-10);
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_obj_write_wavefront() {
        let path = "/tmp/oxiphy_test_wavefront.obj";
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        ObjWriter::write_legacy(path, &verts, &tris, None).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.lines().any(|l| l.starts_with("v ")));
        assert!(content.lines().any(|l| l.starts_with("f ")));
        let v_count = content.lines().filter(|l| l.starts_with("v ")).count();
        let f_count = content.lines().filter(|l| l.starts_with("f ")).count();
        assert_eq!(v_count, 3);
        assert_eq!(f_count, 1);
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_obj_reader_handles_vertex_face_parsing() {
        let path = "/tmp/oxiphy_test_parse.obj";
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let norms = vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let tris = vec![[0, 1, 2]];
        ObjWriter::write_legacy(path, &verts, &tris, Some(&norms)).unwrap();
        let (read_verts, read_tris) = ObjReader::read(path).unwrap();
        assert_eq!(read_verts.len(), 3);
        assert_eq!(read_tris.len(), 1);
        assert_eq!(read_tris[0], [0, 1, 2]);
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_obj_mesh_write_read_vertex_roundtrip() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        let s = ObjWriter::write(&mesh);
        let parsed = ObjReader::from_str(&s).unwrap();
        assert_eq!(parsed.vertices.len(), 3);
        assert!((parsed.vertices[0][0] - 1.0).abs() < 1e-10);
        assert!((parsed.vertices[2][2] - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_obj_mesh_face_format() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        let s = ObjWriter::write(&mesh);
        assert!(s.contains("f 1 2 3"), "face line not found in: {s}");
    }
    #[test]
    fn test_obj_mesh_normal_export() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.normals = vec![[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: Some(vec![0, 1, 2]),
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let s = ObjWriter::write(&mesh);
        assert!(s.contains("vn"), "normals not exported: {s}");
        assert!(s.contains("//"), "face should use v//vn format: {s}");
    }
    #[test]
    fn test_multi_object_groups() {
        let g1 = ObjGroup {
            name: "body".into(),
            face_start: 0,
            face_count: 4,
        };
        let g2 = ObjGroup {
            name: "wheel".into(),
            face_start: 4,
            face_count: 2,
        };
        assert_eq!(g1.name, "body");
        assert_eq!(g2.face_start, 4);
        assert_eq!(g1.face_count + g2.face_count, 6);
    }
    #[test]
    fn test_triangle_soup_count() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        mesh.faces.push(make_default_face(vec![1, 3, 2]));
        let soup = mesh.to_triangle_soup();
        assert_eq!(soup.len(), 2);
    }
    #[test]
    fn test_mtl_writer_output() {
        let mat = ObjMaterial {
            name: "Red".into(),
            kd: [1.0, 0.0, 0.0],
            ks: [0.5, 0.5, 0.5],
            ns: 32.0,
            ka: [0.1, 0.0, 0.0],
            dissolve: 1.0,
            map_kd: None,
        };
        let s = MtlWriter::write(&[mat]);
        assert!(s.contains("newmtl Red"), "material name missing: {s}");
        assert!(s.contains("Kd 1"), "diffuse colour missing: {s}");
        assert!(s.contains("Ns 32"), "shininess missing: {s}");
    }
    #[test]
    fn test_triangle_soup_quad_triangulation() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        mesh.faces.push(make_default_face(vec![0, 1, 2, 3]));
        let soup = mesh.to_triangle_soup();
        assert_eq!(
            soup.len(),
            2,
            "quad should triangulate to 2 triangles, got {}",
            soup.len()
        );
    }
    #[test]
    fn test_group_parsing() {
        let data = "\
v 0 0 0
v 1 0 0
v 0 1 0
v 1 1 0
g group1
f 1 2 3
g group2
f 2 4 3
";
        let mesh = ObjReader::from_str(data).unwrap();
        assert_eq!(mesh.groups.len(), 2);
        assert_eq!(mesh.groups[0].name, "group1");
        assert_eq!(mesh.groups[0].face_count, 1);
        assert_eq!(mesh.groups[1].name, "group2");
        assert_eq!(mesh.groups[1].face_count, 1);
    }
    #[test]
    fn test_smoothing_group_parsing() {
        let data = "\
v 0 0 0
v 1 0 0
v 0 1 0
v 1 1 0
s 1
f 1 2 3
s 2
f 2 4 3
";
        let mesh = ObjReader::from_str(data).unwrap();
        assert_eq!(mesh.faces[0].smoothing_group, 1);
        assert_eq!(mesh.faces[1].smoothing_group, 2);
    }
    #[test]
    fn test_smoothing_group_off() {
        let data = "\
v 0 0 0
v 1 0 0
v 0 1 0
s off
f 1 2 3
";
        let mesh = ObjReader::from_str(data).unwrap();
        assert_eq!(mesh.faces[0].smoothing_group, 0);
    }
    #[test]
    fn test_material_parsing() {
        let data = "\
v 0 0 0
v 1 0 0
v 0 1 0
v 1 1 0
usemtl Red
f 1 2 3
usemtl Blue
f 2 4 3
";
        let mesh = ObjReader::from_str(data).unwrap();
        assert_eq!(mesh.faces[0].material.as_deref(), Some("Red"));
        assert_eq!(mesh.faces[1].material.as_deref(), Some("Blue"));
    }
    #[test]
    fn test_faces_in_group() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        mesh.faces.push(make_default_face(vec![1, 3, 2]));
        mesh.faces.push(make_default_face(vec![0, 3, 2]));
        mesh.groups.push(ObjGroup {
            name: "A".into(),
            face_start: 0,
            face_count: 2,
        });
        mesh.groups.push(ObjGroup {
            name: "B".into(),
            face_start: 2,
            face_count: 1,
        });
        assert_eq!(mesh.faces_in_group("A").len(), 2);
        assert_eq!(mesh.faces_in_group("B").len(), 1);
        assert_eq!(mesh.faces_in_group("C").len(), 0);
    }
    #[test]
    fn test_faces_in_smoothing_group() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 4];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 1,
            material: None,
        });
        mesh.faces.push(ObjFace {
            vertex_indices: vec![1, 3, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 2,
            material: None,
        });
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 3, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 1,
            material: None,
        });
        assert_eq!(mesh.faces_in_smoothing_group(1).len(), 2);
        assert_eq!(mesh.faces_in_smoothing_group(2).len(), 1);
    }
    #[test]
    fn test_faces_with_material() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 4];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: Some("Red".into()),
        });
        mesh.faces.push(ObjFace {
            vertex_indices: vec![1, 3, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: Some("Blue".into()),
        });
        assert_eq!(mesh.faces_with_material("Red").len(), 1);
        assert_eq!(mesh.faces_with_material("Blue").len(), 1);
        assert_eq!(mesh.faces_with_material("Green").len(), 0);
    }
    #[test]
    fn test_triangle_count() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 5];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        mesh.faces.push(make_default_face(vec![0, 1, 2, 3]));
        mesh.faces.push(make_default_face(vec![0, 1, 2, 3, 4]));
        assert_eq!(mesh.triangle_count(), 6);
    }
    #[test]
    fn test_face_normal() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        let n = mesh.face_normal(0).unwrap();
        assert!((n[0] - 0.0).abs() < 1e-10);
        assert!((n[1] - 0.0).abs() < 1e-10);
        assert!((n[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bounding_box() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0], [0.0, 0.0, 0.0]];
        let (min, max) = mesh.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-10);
        assert!((min[1] - (-2.0)).abs() < 1e-10);
        assert!((min[2] - (-3.0)).abs() < 1e-10);
        assert!((max[0] - 1.0).abs() < 1e-10);
        assert!((max[1] - 2.0).abs() < 1e-10);
        assert!((max[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_bounding_box_empty() {
        let mesh = ObjMesh::default();
        assert!(mesh.bounding_box().is_none());
    }
    #[test]
    fn test_texture_coordinate_roundtrip() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.uvs = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: Some(vec![0, 1, 2]),
            smoothing_group: 0,
            material: None,
        });
        let s = ObjWriter::write(&mesh);
        let parsed = ObjReader::from_str(&s).unwrap();
        assert_eq!(parsed.uvs.len(), 3);
        assert!((parsed.uvs[2][1] - 1.0).abs() < 1e-10);
        assert!(parsed.faces[0].uv_indices.is_some());
    }
    #[test]
    fn test_curve_struct() {
        let curve = ObjCurve {
            name: "test_curve".into(),
            degree: 3,
            control_points: vec![0, 1, 2, 3],
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        };
        assert_eq!(curve.degree, 3);
        assert_eq!(curve.control_points.len(), 4);
        assert_eq!(curve.knots.len(), 8);
    }
    #[test]
    fn test_surface_struct() {
        let surface = ObjSurface {
            name: "test_surface".into(),
            degree_u: 2,
            degree_v: 2,
            control_points: vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            n_u: 3,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        };
        assert_eq!(surface.degree_u, 2);
        assert_eq!(surface.n_u, 3);
        assert_eq!(surface.control_points.len(), 9);
    }
    #[test]
    fn test_material_basic() {
        let mat = ObjMaterial::basic("test_mat", [0.5, 0.5, 0.5]);
        assert_eq!(mat.name, "test_mat");
        assert!((mat.kd[0] - 0.5).abs() < 1e-10);
        assert!((mat.dissolve - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_write_with_groups_roundtrip() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        mesh.faces.push(make_default_face(vec![0, 1, 2]));
        mesh.faces.push(make_default_face(vec![1, 3, 2]));
        mesh.groups.push(ObjGroup {
            name: "grp1".into(),
            face_start: 0,
            face_count: 1,
        });
        mesh.groups.push(ObjGroup {
            name: "grp2".into(),
            face_start: 1,
            face_count: 1,
        });
        let s = ObjWriter::write_with_groups(&mesh, true);
        assert!(s.contains("g grp1"));
        assert!(s.contains("g grp2"));
        let parsed = ObjReader::from_str(&s).unwrap();
        assert_eq!(parsed.groups.len(), 2);
    }
    #[test]
    fn test_write_with_uvs() {
        let path = "/tmp/oxiphy_test_uvs.obj";
        let verts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let uvs = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let tris = vec![[0, 1, 2]];
        ObjWriter::write_with_uvs(path, &verts, &uvs, &tris).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.lines().any(|l| l.starts_with("vt ")));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_mtl_writer_with_texture() {
        let mat = ObjMaterial {
            name: "Textured".into(),
            kd: [1.0, 1.0, 1.0],
            ks: [0.0; 3],
            ns: 1.0,
            ka: [0.1, 0.1, 0.1],
            dissolve: 0.8,
            map_kd: Some("diffuse.png".into()),
        };
        let s = MtlWriter::write(&[mat]);
        assert!(s.contains("map_Kd diffuse.png"));
        assert!(s.contains("d 0.8"));
    }
    #[test]
    fn test_write_smoothing_groups() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 1,
            material: None,
        });
        mesh.faces.push(ObjFace {
            vertex_indices: vec![1, 3, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 2,
            material: None,
        });
        let s = ObjWriter::write(&mesh);
        assert!(s.contains("s 1"), "smoothing group 1 missing: {s}");
        assert!(s.contains("s 2"), "smoothing group 2 missing: {s}");
    }
    #[test]
    fn test_write_material_headers() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 4];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: Some("Mat1".into()),
        });
        mesh.faces.push(ObjFace {
            vertex_indices: vec![1, 3, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: Some("Mat2".into()),
        });
        let s = ObjWriter::write(&mesh);
        assert!(s.contains("usemtl Mat1"));
        assert!(s.contains("usemtl Mat2"));
    }

    // ── G4: QEM decimation tests ──────────────────────────────────────────────

    /// Build a UV-sphere-like triangulated mesh with `n_lat` latitude bands and
    /// `n_lon` longitude bands.  Returns an [`ObjMesh`].
    fn make_sphere_mesh(n_lat: usize, n_lon: usize) -> ObjMesh {
        use std::f64::consts::PI;
        let mut vertices: Vec<[f64; 3]> = Vec::new();
        // Stack of latitude rings (poles included via degenerate rings).
        for lat in 0..=n_lat {
            let theta = PI * lat as f64 / n_lat as f64;
            for lon in 0..n_lon {
                let phi = 2.0 * PI * lon as f64 / n_lon as f64;
                vertices.push([
                    theta.sin() * phi.cos(),
                    theta.cos(),
                    theta.sin() * phi.sin(),
                ]);
            }
        }
        let mut faces: Vec<ObjFace> = Vec::new();
        let idx = |lat: usize, lon: usize| lat * n_lon + (lon % n_lon);
        for lat in 0..n_lat {
            for lon in 0..n_lon {
                let a = idx(lat, lon);
                let b = idx(lat + 1, lon);
                let c = idx(lat + 1, lon + 1);
                let d = idx(lat, lon + 1);
                faces.push(ObjFace {
                    vertex_indices: vec![a, b, c],
                    normal_indices: None,
                    uv_indices: None,
                    smoothing_group: 0,
                    material: None,
                });
                faces.push(ObjFace {
                    vertex_indices: vec![a, c, d],
                    normal_indices: None,
                    uv_indices: None,
                    smoothing_group: 0,
                    material: None,
                });
            }
        }
        ObjMesh {
            vertices,
            normals: Vec::new(),
            uvs: Vec::new(),
            faces,
            groups: Vec::new(),
        }
    }

    #[test]
    fn test_decimate_below_target_returns_clone() {
        let mesh = make_sphere_mesh(4, 8); // 64 faces
        let result = ObjLod::decimate(&mesh, 200);
        assert_eq!(result.faces.len(), mesh.faces.len());
    }

    #[test]
    fn test_decimate_qem_reduces_face_count() {
        let mesh = make_sphere_mesh(8, 16); // 256 faces
        let target = 64;
        let result = ObjLod::decimate(&mesh, target);
        assert!(
            result.faces.len() <= target,
            "expected ≤{target} faces, got {}",
            result.faces.len()
        );
        assert!(
            !result.faces.is_empty(),
            "decimated mesh must have at least one face"
        );
    }

    #[test]
    fn test_decimate_vertex_count_decreases() {
        let mesh = make_sphere_mesh(8, 16); // 256 faces, (8+1)*16=144 vertices
        let original_verts = mesh.vertices.len();
        let result = ObjLod::decimate(&mesh, 32);
        assert!(
            result.vertices.len() < original_verts,
            "QEM should reduce vertex count ({} vs {})",
            result.vertices.len(),
            original_verts
        );
    }
}
/// Instantiate `mesh` with a given transform, returning a new transformed `ObjMesh`.
#[allow(dead_code)]
pub fn instantiate_mesh(mesh: &ObjMesh, instance: &MeshInstance) -> ObjMesh {
    let mut out = mesh.clone();
    for v in &mut out.vertices {
        let p = instance.transform.apply(*v);
        *v = p;
    }
    for n in &mut out.normals {
        let rot_only = MeshTransform {
            translation: [0.0; 3],
            scale: 1.0,
            ..instance.transform.clone()
        };
        *n = rot_only.apply(*n);
    }
    out
}
/// Weld (deduplicate) vertices that are closer than `tolerance`.
///
/// Returns the de-duplicated mesh with remapped face indices.
#[allow(dead_code)]
pub fn weld_vertices(mesh: &ObjMesh, tolerance: f64) -> ObjMesh {
    let tol2 = tolerance * tolerance;
    let mut new_verts: Vec<[f64; 3]> = Vec::new();
    let mut remap: Vec<usize> = Vec::with_capacity(mesh.vertices.len());
    for &v in &mesh.vertices {
        let found = new_verts.iter().position(|&u| {
            let dx = u[0] - v[0];
            let dy = u[1] - v[1];
            let dz = u[2] - v[2];
            dx * dx + dy * dy + dz * dz <= tol2
        });
        if let Some(idx) = found {
            remap.push(idx);
        } else {
            remap.push(new_verts.len());
            new_verts.push(v);
        }
    }
    let mut out = ObjMesh {
        vertices: new_verts,
        normals: mesh.normals.clone(),
        uvs: mesh.uvs.clone(),
        groups: mesh.groups.clone(),
        ..Default::default()
    };
    for face in &mesh.faces {
        let new_vis: Vec<usize> = face.vertex_indices.iter().map(|&i| remap[i]).collect();
        out.faces.push(ObjFace {
            vertex_indices: new_vis,
            normal_indices: face.normal_indices.clone(),
            uv_indices: face.uv_indices.clone(),
            smoothing_group: face.smoothing_group,
            material: face.material.clone(),
        });
    }
    out
}
/// Merge two `ObjMesh` instances together.
///
/// Vertices and faces from `b` are appended to `a` with adjusted indices.
#[allow(dead_code)]
pub fn merge_obj_meshes(a: &ObjMesh, b: &ObjMesh) -> ObjMesh {
    let mut out = a.clone();
    let v_offset = a.vertices.len();
    let n_offset = a.normals.len();
    let uv_offset = a.uvs.len();
    out.vertices.extend_from_slice(&b.vertices);
    out.normals.extend_from_slice(&b.normals);
    out.uvs.extend_from_slice(&b.uvs);
    for face in &b.faces {
        let new_vis: Vec<usize> = face.vertex_indices.iter().map(|&i| i + v_offset).collect();
        let new_ns = face
            .normal_indices
            .as_ref()
            .map(|ns| ns.iter().map(|&i| i + n_offset).collect::<Vec<_>>());
        let new_uvs = face
            .uv_indices
            .as_ref()
            .map(|uvs| uvs.iter().map(|&i| i + uv_offset).collect::<Vec<_>>());
        out.faces.push(ObjFace {
            vertex_indices: new_vis,
            normal_indices: new_ns,
            uv_indices: new_uvs,
            smoothing_group: face.smoothing_group,
            material: face.material.clone(),
        });
    }
    let face_offset = a.faces.len();
    for g in &b.groups {
        out.groups.push(ObjGroup {
            name: g.name.clone(),
            face_start: g.face_start + face_offset,
            face_count: g.face_count,
        });
    }
    out
}
/// Recompute per-vertex normals using area-weighted face normals.
///
/// The result is stored in `mesh.normals` and all face `normal_indices`
/// are set to match the corresponding `vertex_indices`.
#[allow(dead_code)]
pub fn recompute_normals(mesh: &mut ObjMesh) {
    let n = mesh.vertices.len();
    let mut accum = vec![[0.0_f64; 3]; n];
    let mut weights = vec![0.0_f64; n];
    for face in &mesh.faces {
        let vis = &face.vertex_indices;
        if vis.len() < 3 {
            continue;
        }
        for i in 1..(vis.len() - 1) {
            let v0 = mesh.vertices[vis[0]];
            let v1 = mesh.vertices[vis[i]];
            let v2 = mesh.vertices[vis[i + 1]];
            let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let nx = e1[1] * e2[2] - e1[2] * e2[1];
            let ny = e1[2] * e2[0] - e1[0] * e2[2];
            let nz = e1[0] * e2[1] - e1[1] * e2[0];
            let area = (nx * nx + ny * ny + nz * nz).sqrt() * 0.5;
            for &vi in &[vis[0], vis[i], vis[i + 1]] {
                accum[vi][0] += nx;
                accum[vi][1] += ny;
                accum[vi][2] += nz;
                weights[vi] += area;
            }
        }
    }
    mesh.normals = accum
        .iter()
        .zip(weights.iter())
        .map(|(n, &w)| {
            if w < 1e-30 {
                [0.0, 0.0, 1.0]
            } else {
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if len < 1e-30 {
                    [0.0, 0.0, 1.0]
                } else {
                    [n[0] / len, n[1] / len, n[2] / len]
                }
            }
        })
        .collect();
    for face in &mut mesh.faces {
        let n_idx: Vec<usize> = face.vertex_indices.clone();
        face.normal_indices = Some(n_idx);
    }
}
/// Parse a Wavefront MTL file content string into a list of materials.
#[allow(dead_code)]
pub fn parse_mtl(data: &str) -> Vec<ObjMaterial> {
    let mut materials: Vec<ObjMaterial> = Vec::new();
    let mut current: Option<ObjMaterial> = None;
    for raw in data.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = line.splitn(2, ' ').collect();
        if tokens.is_empty() {
            continue;
        }
        match tokens[0] {
            "newmtl" => {
                if let Some(mat) = current.take() {
                    materials.push(mat);
                }
                let name = tokens
                    .get(1)
                    .map(|s| s.trim().to_string())
                    .unwrap_or_default();
                current = Some(ObjMaterial {
                    name,
                    kd: [0.8, 0.8, 0.8],
                    ks: [0.0; 3],
                    ns: 1.0,
                    ka: [0.0; 3],
                    dissolve: 1.0,
                    map_kd: None,
                });
            }
            "Kd" | "kd" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    let p: Vec<&str> = rest.split_whitespace().collect();
                    if p.len() >= 3 {
                        mat.kd[0] = p[0].parse().unwrap_or(0.8);
                        mat.kd[1] = p[1].parse().unwrap_or(0.8);
                        mat.kd[2] = p[2].parse().unwrap_or(0.8);
                    }
                }
            }
            "Ks" | "ks" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    let p: Vec<&str> = rest.split_whitespace().collect();
                    if p.len() >= 3 {
                        mat.ks[0] = p[0].parse().unwrap_or(0.0);
                        mat.ks[1] = p[1].parse().unwrap_or(0.0);
                        mat.ks[2] = p[2].parse().unwrap_or(0.0);
                    }
                }
            }
            "Ka" | "ka" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    let p: Vec<&str> = rest.split_whitespace().collect();
                    if p.len() >= 3 {
                        mat.ka[0] = p[0].parse().unwrap_or(0.0);
                        mat.ka[1] = p[1].parse().unwrap_or(0.0);
                        mat.ka[2] = p[2].parse().unwrap_or(0.0);
                    }
                }
            }
            "Ns" | "ns" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    mat.ns = rest.trim().parse().unwrap_or(1.0);
                }
            }
            "d" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    mat.dissolve = rest.trim().parse().unwrap_or(1.0);
                }
            }
            "Tr" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    let tr: f64 = rest.trim().parse().unwrap_or(0.0);
                    mat.dissolve = 1.0 - tr;
                }
            }
            "map_Kd" | "map_kd" => {
                if let Some(ref mut mat) = current
                    && let Some(rest) = tokens.get(1)
                {
                    mat.map_kd = Some(rest.trim().to_string());
                }
            }
            _ => {}
        }
    }
    if let Some(mat) = current {
        materials.push(mat);
    }
    materials
}
/// Compute mesh statistics for an `ObjMesh`.
#[allow(dead_code)]
pub fn compute_mesh_stats(mesh: &ObjMesh) -> ObjMeshStats {
    let mut mat_names: Vec<&str> = Vec::new();
    let mut faces_with_normals = 0;
    let mut faces_with_uvs = 0;
    let mut surface_area = 0.0_f64;
    for face in &mesh.faces {
        if face.normal_indices.is_some() {
            faces_with_normals += 1;
        }
        if face.uv_indices.is_some() {
            faces_with_uvs += 1;
        }
        if let Some(ref m) = face.material
            && !mat_names.contains(&m.as_str())
        {
            mat_names.push(m.as_str());
        }
        let vis = &face.vertex_indices;
        for i in 1..(vis.len().saturating_sub(1)) {
            if vis[0] < mesh.vertices.len()
                && vis[i] < mesh.vertices.len()
                && vis[i + 1] < mesh.vertices.len()
            {
                let v0 = mesh.vertices[vis[0]];
                let v1 = mesh.vertices[vis[i]];
                let v2 = mesh.vertices[vis[i + 1]];
                let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
                let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
                let cx = e1[1] * e2[2] - e1[2] * e2[1];
                let cy = e1[2] * e2[0] - e1[0] * e2[2];
                let cz = e1[0] * e2[1] - e1[1] * e2[0];
                surface_area += (cx * cx + cy * cy + cz * cz).sqrt() * 0.5;
            }
        }
    }
    ObjMeshStats {
        vertex_count: mesh.vertices.len(),
        face_count: mesh.faces.len(),
        triangle_count: mesh.triangle_count(),
        material_count: mat_names.len(),
        group_count: mesh.groups.len(),
        faces_with_normals,
        faces_with_uvs,
        surface_area,
        bbox: mesh.bounding_box(),
    }
}
