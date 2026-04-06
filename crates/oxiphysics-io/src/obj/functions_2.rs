//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::field_reassign_with_default)]

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use crate::obj::types::*;
#[cfg(test)]
mod tests_expanded {
    use super::*;
    #[test]
    fn test_vertex_color_lerp() {
        let c0 = ObjVertexColor::rgb(0.0, 0.0, 0.0);
        let c1 = ObjVertexColor::rgb(1.0, 1.0, 1.0);
        let mid = c0.lerp(c1, 0.5);
        assert!((mid.r - 0.5).abs() < 1e-10);
        assert!((mid.g - 0.5).abs() < 1e-10);
        assert!((mid.b - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_vertex_color_to_array() {
        let c = ObjVertexColor::rgba(0.1, 0.2, 0.3, 0.4);
        let arr = c.to_array();
        assert!((arr[0] - 0.1).abs() < 1e-10);
        assert!((arr[3] - 0.4).abs() < 1e-10);
    }
    #[test]
    fn test_transform_identity() {
        let t = MeshTransform::identity();
        let p = [1.0_f64, 2.0, 3.0];
        let out = t.apply(p);
        assert!((out[0] - 1.0).abs() < 1e-10);
        assert!((out[1] - 2.0).abs() < 1e-10);
        assert!((out[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_translation() {
        let t = MeshTransform::from_translation(1.0, 2.0, 3.0);
        let p = [0.0_f64, 0.0, 0.0];
        let out = t.apply(p);
        assert!((out[0] - 1.0).abs() < 1e-10);
        assert!((out[1] - 2.0).abs() < 1e-10);
        assert!((out[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_scale() {
        let t = MeshTransform {
            translation: [0.0; 3],
            scale: 2.0,
            axis: [0.0, 0.0, 1.0],
            angle: 0.0,
        };
        let p = [1.0_f64, 1.0, 1.0];
        let out = t.apply(p);
        assert!((out[0] - 2.0).abs() < 1e-10);
        assert!((out[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_transform_rotation_z_90() {
        let t = MeshTransform {
            translation: [0.0; 3],
            scale: 1.0,
            axis: [0.0, 0.0, 1.0],
            angle: std::f64::consts::FRAC_PI_2,
        };
        let p = [1.0_f64, 0.0, 0.0];
        let out = t.apply(p);
        assert!((out[0] - 0.0).abs() < 1e-10, "x={}", out[0]);
        assert!((out[1] - 1.0).abs() < 1e-10, "y={}", out[1]);
        assert!(out[2].abs() < 1e-10);
    }
    #[test]
    fn test_weld_vertices_exact() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let welded = weld_vertices(&mesh, 1e-9);
        assert_eq!(welded.vertices.len(), 2, "should collapse duplicate vertex");
    }
    #[test]
    fn test_weld_vertices_no_duplicates() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let welded = weld_vertices(&mesh, 1e-9);
        assert_eq!(welded.vertices.len(), 3);
    }
    #[test]
    fn test_weld_vertices_tolerance() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [0.0001, 0.0, 0.0], [1.0, 0.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let welded = weld_vertices(&mesh, 0.001);
        assert_eq!(welded.vertices.len(), 2);
    }
    #[test]
    fn test_merge_obj_meshes_vertex_count() {
        let mut a = ObjMesh::default();
        a.vertices = vec![[0.0; 3]; 3];
        a.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let mut b = ObjMesh::default();
        b.vertices = vec![[1.0; 3]; 4];
        b.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let merged = merge_obj_meshes(&a, &b);
        assert_eq!(merged.vertices.len(), 7);
        assert_eq!(merged.faces.len(), 2);
        assert_eq!(merged.faces[1].vertex_indices[0], 3);
    }
    #[test]
    fn test_merge_obj_meshes_empty_a() {
        let a = ObjMesh::default();
        let mut b = ObjMesh::default();
        b.vertices = vec![[1.0, 2.0, 3.0]];
        let merged = merge_obj_meshes(&a, &b);
        assert_eq!(merged.vertices.len(), 1);
    }
    #[test]
    fn test_recompute_normals_flat_triangle() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        recompute_normals(&mut mesh);
        assert_eq!(mesh.normals.len(), 3);
        for n in &mesh.normals {
            assert!(n[2] > 0.9, "normal should face +z: {:?}", n);
        }
    }
    #[test]
    fn test_recompute_normals_sets_face_indices() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        recompute_normals(&mut mesh);
        assert!(mesh.faces[0].normal_indices.is_some());
        assert_eq!(mesh.faces[0].normal_indices.as_ref().unwrap(), &[0, 1, 2]);
    }
    #[test]
    fn test_parse_mtl_basic() {
        let data = "newmtl Mat1\nKd 1.0 0.0 0.0\nKs 0.5 0.5 0.5\nNs 32.0\nd 0.8\n";
        let mats = parse_mtl(data);
        assert_eq!(mats.len(), 1);
        assert_eq!(mats[0].name, "Mat1");
        assert!((mats[0].kd[0] - 1.0).abs() < 1e-10);
        assert!((mats[0].ns - 32.0).abs() < 1e-10);
        assert!((mats[0].dissolve - 0.8).abs() < 1e-10);
    }
    #[test]
    fn test_parse_mtl_multiple_materials() {
        let data = "newmtl A\nKd 1 0 0\nnewmtl B\nKd 0 1 0\n";
        let mats = parse_mtl(data);
        assert_eq!(mats.len(), 2);
        assert_eq!(mats[0].name, "A");
        assert_eq!(mats[1].name, "B");
        assert!((mats[1].kd[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_parse_mtl_map_kd() {
        let data = "newmtl Mat\nmap_Kd diffuse.png\n";
        let mats = parse_mtl(data);
        assert_eq!(mats[0].map_kd.as_deref(), Some("diffuse.png"));
    }
    #[test]
    fn test_parse_mtl_transparency() {
        let data = "newmtl Glass\nTr 0.5\n";
        let mats = parse_mtl(data);
        assert!((mats[0].dissolve - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_parse_mtl_empty() {
        let mats = parse_mtl("");
        assert!(mats.is_empty());
    }
    #[test]
    fn test_parse_mtl_comments() {
        let data = "# This is a comment\nnewmtl MyMat\n# another comment\nKd 0.5 0.5 0.5\n";
        let mats = parse_mtl(data);
        assert_eq!(mats.len(), 1);
        assert_eq!(mats[0].name, "MyMat");
    }
    #[test]
    fn test_lod_select() {
        let mut lod = ObjLod::new();
        let m0 = ObjMesh::default();
        let m1 = ObjMesh::default();
        let m2 = ObjMesh::default();
        lod.push(m0, 10.0);
        lod.push(m1, 50.0);
        lod.push(m2, 200.0);
        assert_eq!(lod.num_levels(), 3);
        assert!(lod.select(5.0).is_some());
        assert!(lod.select(30.0).is_some());
        assert!(lod.select(500.0).is_some());
    }
    #[test]
    fn test_lod_decimate_keeps_target() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 6];
        for i in 0..6 {
            mesh.faces.push(ObjFace {
                vertex_indices: vec![i % 3, (i + 1) % 3, (i + 2) % 3],
                normal_indices: None,
                uv_indices: None,
                smoothing_group: 0,
                material: None,
            });
        }
        let dec = ObjLod::decimate(&mesh, 3);
        assert!(dec.faces.len() <= 3 + 1, "decimated to roughly 3 faces");
    }
    #[test]
    fn test_lod_decimate_no_change_if_under_target() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 3];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let dec = ObjLod::decimate(&mesh, 100);
        assert_eq!(
            dec.faces.len(),
            1,
            "no decimation when already under target"
        );
    }
    #[test]
    fn test_scene_add_mesh() {
        let mut scene = ObjScene::new();
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 3];
        let idx = scene.add_mesh(mesh);
        assert_eq!(idx, 0);
        assert_eq!(scene.meshes.len(), 1);
    }
    #[test]
    fn test_scene_total_vertices() {
        let mut scene = ObjScene::new();
        let mut m1 = ObjMesh::default();
        m1.vertices = vec![[0.0; 3]; 4];
        let mut m2 = ObjMesh::default();
        m2.vertices = vec![[1.0; 3]; 6];
        scene.add_mesh(m1);
        scene.add_mesh(m2);
        assert_eq!(scene.total_vertices(), 10);
    }
    #[test]
    fn test_scene_flatten() {
        let mut scene = ObjScene::new();
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let mi = scene.add_mesh(mesh);
        let t = MeshTransform::from_translation(10.0, 0.0, 0.0);
        scene.add_node("node1", Some(mi), t);
        let flat = scene.flatten();
        assert_eq!(flat.vertices.len(), 3);
        assert!((flat.vertices[0][0] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_scene_add_child() {
        let mut scene = ObjScene::new();
        let ni0 = scene.add_node("parent", None, MeshTransform::identity());
        let ni1 = scene.add_node("child", None, MeshTransform::identity());
        scene.add_child(ni0, ni1);
        assert_eq!(scene.nodes[ni0].children, vec![ni1]);
    }
    #[test]
    fn test_compute_mesh_stats_basic() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let stats = compute_mesh_stats(&mesh);
        assert_eq!(stats.vertex_count, 3);
        assert_eq!(stats.face_count, 1);
        assert_eq!(stats.triangle_count, 1);
        assert!(stats.surface_area > 0.0);
    }
    #[test]
    fn test_compute_mesh_stats_materials() {
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
        let stats = compute_mesh_stats(&mesh);
        assert_eq!(stats.material_count, 2);
    }
    #[test]
    fn test_compute_mesh_stats_bbox() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]];
        let stats = compute_mesh_stats(&mesh);
        let (min, max) = stats.bbox.unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-10);
        assert!((max[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vcmesh_parse_basic() {
        let data = "v 0 0 0 1.0 0.0 0.0\nv 1 0 0 0.0 1.0 0.0\nv 0 1 0 0.0 0.0 1.0\nf 1 2 3\n";
        let vcm = ObjVertexColorMesh::from_str(data).unwrap();
        assert_eq!(vcm.mesh.vertices.len(), 3);
        assert_eq!(vcm.colors.len(), 3);
        assert!((vcm.colors[0].r - 1.0).abs() < 1e-10);
        assert!((vcm.colors[1].g - 1.0).abs() < 1e-10);
        assert!((vcm.colors[2].b - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vcmesh_roundtrip() {
        let data = "v 0.0 0.0 0.0 0.5 0.5 0.5\nv 1.0 0.0 0.0 0.5 0.5 0.5\nv 0.0 1.0 0.0 0.5 0.5 0.5\nf 1 2 3\n";
        let vcm = ObjVertexColorMesh::from_str(data).unwrap();
        let s = vcm.to_obj_str();
        let reparsed = ObjVertexColorMesh::from_str(&s).unwrap();
        assert_eq!(reparsed.mesh.vertices.len(), 3);
        assert_eq!(reparsed.mesh.faces.len(), 1);
    }
    #[test]
    fn test_vcmesh_no_color_defaults_white() {
        let data = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let vcm = ObjVertexColorMesh::from_str(data).unwrap();
        assert_eq!(vcm.colors.len(), 3);
        assert!((vcm.colors[0].r - 1.0).abs() < 1e-10);
        assert!((vcm.colors[0].g - 1.0).abs() < 1e-10);
        assert!((vcm.colors[0].b - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_instantiate_mesh_translation() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0, 0.0, 0.0]];
        let inst = MeshInstance {
            name: "inst1".into(),
            transform: MeshTransform::from_translation(5.0, 6.0, 7.0),
        };
        let out = instantiate_mesh(&mesh, &inst);
        assert!((out.vertices[0][0] - 5.0).abs() < 1e-10);
        assert!((out.vertices[0][1] - 6.0).abs() < 1e-10);
        assert!((out.vertices[0][2] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_instantiate_mesh_preserves_face_count() {
        let mut mesh = ObjMesh::default();
        mesh.vertices = vec![[0.0; 3]; 3];
        mesh.faces.push(ObjFace {
            vertex_indices: vec![0, 1, 2],
            normal_indices: None,
            uv_indices: None,
            smoothing_group: 0,
            material: None,
        });
        let inst = MeshInstance {
            name: "x".into(),
            transform: MeshTransform::identity(),
        };
        let out = instantiate_mesh(&mesh, &inst);
        assert_eq!(out.faces.len(), 1);
    }
}
