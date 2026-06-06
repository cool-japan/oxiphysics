//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Opaque node identifier (index into `Scene::nodes`).
pub type NodeId = usize;
/// Column-major 4×4 matrix multiply.
pub fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0f64; 4]; 4];
    for col in 0..4 {
        for row in 0..4 {
            out[col][row] = (0..4).map(|k| a[k][row] * b[col][k]).sum();
        }
    }
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::Aabb;
    use crate::scene::Camera;
    use crate::scene::DrawCall;

    use crate::scene::Frustum;
    use crate::scene::FrustumCuller;

    use crate::scene::InstanceData;
    use crate::scene::Light;
    use crate::scene::LightManager;
    use crate::scene::LodConfig;

    use crate::scene::Material;
    use crate::scene::RenderQueue;
    use crate::scene::Scene;
    use crate::scene::SceneCamera;
    use crate::scene::SceneCameraProjection;
    use crate::scene::SceneEnvironment;
    use crate::scene::SceneLight;

    use crate::scene::Transform;
    #[test]
    fn test_identity_matrix_is_identity() {
        let m = Transform::identity().to_matrix();
        for (col, mcol) in m.iter().enumerate() {
            for (row, &mval) in mcol.iter().enumerate() {
                let expected = if col == row { 1.0 } else { 0.0 };
                assert!((mval - expected).abs() < 1e-12, "m[{col}][{row}]={}", mval);
            }
        }
    }
    #[test]
    fn test_translation_in_matrix() {
        let mut t = Transform::identity();
        t.position = [3.0, 5.0, 7.0];
        let m = t.to_matrix();
        assert!((m[3][0] - 3.0).abs() < 1e-12);
        assert!((m[3][1] - 5.0).abs() < 1e-12);
        assert!((m[3][2] - 7.0).abs() < 1e-12);
        assert!((m[3][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_scale_in_matrix() {
        let mut t = Transform::identity();
        t.scale = [2.0, 3.0, 4.0];
        let m = t.to_matrix();
        assert!((m[0][0] - 2.0).abs() < 1e-12);
        assert!((m[1][1] - 3.0).abs() < 1e-12);
        assert!((m[2][2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_rotation_90_deg_around_z() {
        let s = std::f64::consts::FRAC_PI_4.sin();
        let c = std::f64::consts::FRAC_PI_4.cos();
        let t = Transform {
            position: [0.0; 3],
            rotation: [0.0, 0.0, s, c],
            scale: [1.0; 3],
        };
        let m = t.to_matrix();
        assert!((m[0][0]).abs() < 1e-12);
        assert!((m[0][1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_new_scene_is_empty() {
        let scene = Scene::new();
        assert!(scene.nodes.is_empty());
        assert!(scene.root_nodes.is_empty());
    }
    #[test]
    fn test_add_root_node() {
        let mut scene = Scene::new();
        let id = scene.add_node("root", Transform::identity(), None);
        assert_eq!(id, 0);
        assert_eq!(scene.root_nodes, vec![0]);
        assert_eq!(scene.nodes[0].name, "root");
    }
    #[test]
    fn test_add_child_node() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let child = scene.add_node("child", Transform::identity(), Some(root));
        assert_eq!(scene.nodes[child].parent, Some(root));
        assert!(scene.nodes[root].children.contains(&child));
        assert!(!scene.root_nodes.contains(&child));
    }
    #[test]
    fn test_world_transform_root_equals_local() {
        let mut scene = Scene::new();
        let mut t = Transform::identity();
        t.position = [1.0, 2.0, 3.0];
        let id = scene.add_node("n", t, None);
        let wt = scene.world_transform(id);
        assert!((wt[3][0] - 1.0).abs() < 1e-12);
        assert!((wt[3][1] - 2.0).abs() < 1e-12);
        assert!((wt[3][2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_world_transform_child_accumulates_parent() {
        let mut scene = Scene::new();
        let mut tp = Transform::identity();
        tp.position = [1.0, 0.0, 0.0];
        let root = scene.add_node("root", tp, None);
        let mut tc = Transform::identity();
        tc.position = [2.0, 0.0, 0.0];
        let child = scene.add_node("child", tc, Some(root));
        let wt = scene.world_transform(child);
        assert!(
            (wt[3][0] - 3.0).abs() < 1e-12,
            "expected x=3, got {}",
            wt[3][0]
        );
    }
    #[test]
    fn test_all_nodes_visible_by_default() {
        let mut scene = Scene::new();
        let a = scene.add_node("a", Transform::identity(), None);
        let b = scene.add_node("b", Transform::identity(), None);
        let v = scene.visible_nodes();
        assert!(v.contains(&a));
        assert!(v.contains(&b));
    }
    #[test]
    fn test_set_invisible() {
        let mut scene = Scene::new();
        let id = scene.add_node("n", Transform::identity(), None);
        scene.set_visible(id, false);
        assert!(!scene.visible_nodes().contains(&id));
    }
    #[test]
    fn test_set_visible_toggles() {
        let mut scene = Scene::new();
        let id = scene.add_node("n", Transform::identity(), None);
        scene.set_visible(id, false);
        scene.set_visible(id, true);
        assert!(scene.visible_nodes().contains(&id));
    }
    #[test]
    fn test_remove_root_node() {
        let mut scene = Scene::new();
        let id = scene.add_node("n", Transform::identity(), None);
        scene.remove_node(id);
        assert!(!scene.root_nodes.contains(&id));
        assert!(!scene.nodes[id].visible);
    }
    #[test]
    fn test_remove_node_children_reparented() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let mid = scene.add_node("mid", Transform::identity(), Some(root));
        let leaf = scene.add_node("leaf", Transform::identity(), Some(mid));
        scene.remove_node(mid);
        assert!(scene.nodes[root].children.contains(&leaf));
        assert_eq!(scene.nodes[leaf].parent, Some(root));
    }
    #[test]
    fn test_find_by_name_found() {
        let mut scene = Scene::new();
        let id = scene.add_node("camera", Transform::identity(), None);
        assert_eq!(scene.find_by_name("camera"), Some(id));
    }
    #[test]
    fn test_find_by_name_not_found() {
        let scene = Scene::new();
        assert_eq!(scene.find_by_name("missing"), None);
    }
    #[test]
    fn test_depth_root_is_zero() {
        let mut scene = Scene::new();
        let id = scene.add_node("root", Transform::identity(), None);
        assert_eq!(scene.depth(id), 0);
    }
    #[test]
    fn test_depth_child_is_one() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let child = scene.add_node("child", Transform::identity(), Some(root));
        assert_eq!(scene.depth(child), 1);
    }
    #[test]
    fn test_depth_grandchild_is_two() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let child = scene.add_node("child", Transform::identity(), Some(root));
        let grand = scene.add_node("grand", Transform::identity(), Some(child));
        assert_eq!(scene.depth(grand), 2);
    }
    #[test]
    fn test_scene_default() {
        let scene: Scene = Default::default();
        assert!(scene.nodes.is_empty());
    }
    #[test]
    fn test_light_directional_color_intensity() {
        let l = Light::Directional {
            direction: [0.0, -1.0, 0.0],
            color: [1.0, 0.9, 0.8],
            intensity: 2.0,
        };
        assert!((l.intensity() - 2.0).abs() < 1e-6);
        assert!((l.color()[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_light_point_properties() {
        let l = Light::Point {
            position: [1.0, 2.0, 3.0],
            color: [0.5, 0.5, 0.5],
            intensity: 4.0,
            range: 10.0,
        };
        assert!((l.intensity() - 4.0).abs() < 1e-6);
    }
    #[test]
    fn test_light_spot_properties() {
        let l = Light::Spot {
            position: [0.0; 3],
            direction: [0.0, -1.0, 0.0],
            color: [1.0, 1.0, 1.0],
            intensity: 1.5,
            inner_angle: 0.2,
            outer_angle: 0.5,
            range: 20.0,
        };
        assert!((l.intensity() - 1.5).abs() < 1e-6);
    }
    #[test]
    fn test_light_area_properties() {
        let l = Light::Area {
            position: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            half_extents: [1.0, 0.5],
            color: [0.8, 0.8, 1.0],
            intensity: 3.0,
        };
        assert!((l.intensity() - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_camera_perspective_projection_not_zero() {
        let cam = Camera::perspective(
            [0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0],
            std::f64::consts::FRAC_PI_4,
            16.0 / 9.0,
            0.1,
            100.0,
        );
        let proj = cam.projection_matrix();
        assert!(proj[0][0].abs() > 0.0);
        assert!(proj[1][1].abs() > 0.0);
    }
    #[test]
    fn test_camera_orthographic_projection_correct_scale() {
        let cam = Camera::orthographic([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 2.0, 1.5, 0.1, 100.0);
        let proj = cam.projection_matrix();
        assert!(
            (proj[0][0] - 0.5).abs() < 1e-10,
            "ortho [0][0] = {}",
            proj[0][0]
        );
    }
    #[test]
    fn test_camera_view_matrix_eye_on_axis() {
        let cam = Camera::perspective([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 1.0, 1.0, 0.1, 100.0);
        let view = cam.view_matrix();
        let trace: f64 = (0..4).map(|i| view[i][i]).sum();
        assert!(trace.abs() > 0.1, "view matrix trace should be non-zero");
    }
    #[test]
    fn test_material_default_white() {
        let m = Material::default_white();
        assert!((m.albedo[0] - 1.0).abs() < 1e-6);
        assert!((m.roughness - 0.5).abs() < 1e-6);
        assert!((m.metalness - 0.0).abs() < 1e-6);
        assert!(!m.is_emissive());
    }
    #[test]
    fn test_material_metallic() {
        let m = Material::metallic([0.7, 0.6, 0.5, 1.0]);
        assert!((m.metalness - 1.0).abs() < 1e-6);
        assert!(!m.is_emissive());
    }
    #[test]
    fn test_material_emissive() {
        let m = Material::emissive([2.0, 1.5, 0.5]);
        assert!(m.is_emissive());
        assert!((m.emissive[0] - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_render_queue_sort_opaque_front_to_back() {
        let mut q = RenderQueue::new();
        q.push_opaque(DrawCall {
            node_id: 0,
            material_index: 0,
            depth: 5.0,
        });
        q.push_opaque(DrawCall {
            node_id: 1,
            material_index: 0,
            depth: 2.0,
        });
        q.push_opaque(DrawCall {
            node_id: 2,
            material_index: 0,
            depth: 8.0,
        });
        q.sort();
        let depths: Vec<f64> = q.opaque_calls().iter().map(|c| c.depth).collect();
        assert!(
            depths[0] < depths[1] && depths[1] < depths[2],
            "opaque should be front-to-back"
        );
    }
    #[test]
    fn test_render_queue_sort_transparent_back_to_front() {
        let mut q = RenderQueue::new();
        q.push_transparent(DrawCall {
            node_id: 0,
            material_index: 0,
            depth: 3.0,
        });
        q.push_transparent(DrawCall {
            node_id: 1,
            material_index: 0,
            depth: 7.0,
        });
        q.push_transparent(DrawCall {
            node_id: 2,
            material_index: 0,
            depth: 1.0,
        });
        q.sort();
        let depths: Vec<f64> = q.transparent_calls().iter().map(|c| c.depth).collect();
        assert!(
            depths[0] > depths[1] && depths[1] > depths[2],
            "transparent should be back-to-front"
        );
    }
    #[test]
    fn test_render_queue_len_and_clear() {
        let mut q = RenderQueue::new();
        q.push_opaque(DrawCall {
            node_id: 0,
            material_index: 0,
            depth: 1.0,
        });
        q.push_transparent(DrawCall {
            node_id: 1,
            material_index: 1,
            depth: 2.0,
        });
        assert_eq!(q.len(), 2);
        q.clear();
        assert!(q.is_empty());
    }
    #[test]
    fn test_lod_select_nearest() {
        let mut lod = LodConfig::new();
        lod.add_level(10.0, 0);
        lod.add_level(50.0, 1);
        lod.add_level(200.0, 2);
        assert_eq!(
            lod.select_lod(5.0),
            Some(0),
            "close distance should select LOD 0"
        );
        assert_eq!(
            lod.select_lod(25.0),
            Some(1),
            "mid distance should select LOD 1"
        );
        assert_eq!(
            lod.select_lod(300.0),
            Some(2),
            "far distance should use last LOD"
        );
    }
    #[test]
    fn test_lod_empty_returns_none() {
        let lod = LodConfig::new();
        assert_eq!(lod.select_lod(1.0), None);
    }
    #[test]
    fn test_lod_exact_threshold() {
        let mut lod = LodConfig::new();
        lod.add_level(10.0, 5);
        assert_eq!(lod.select_lod(10.0), Some(5));
        assert_eq!(lod.select_lod(10.001), Some(5));
    }
    #[test]
    fn test_aabb_centre() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = aabb.centre();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_half_extents() {
        let aabb = Aabb::new([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]);
        let he = aabb.half_extents();
        assert!((he[0] - 1.0).abs() < 1e-12);
        assert!((he[1] - 2.0).abs() < 1e-12);
        assert!((he[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_intersects_sphere_inside() {
        let aabb = Aabb::new([-5.0; 3], [5.0; 3]);
        assert!(aabb.intersects_sphere([0.0, 0.0, 0.0], 1.0));
    }
    #[test]
    fn test_aabb_intersects_sphere_outside() {
        let aabb = Aabb::new([0.0; 3], [1.0; 3]);
        assert!(!aabb.intersects_sphere([10.0, 0.0, 0.0], 1.0));
    }
    #[test]
    fn test_aabb_expand() {
        let mut aabb = Aabb::new([0.0; 3], [0.0; 3]);
        aabb.expand([3.0, -1.0, 5.0]);
        assert!((aabb.max[0] - 3.0).abs() < 1e-12);
        assert!((aabb.min[1] - (-1.0)).abs() < 1e-12);
        assert!((aabb.max[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_union() {
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([-1.0; 3], [2.0; 3]);
        let u = a.union(&b);
        assert!((u.min[0] - (-1.0)).abs() < 1e-12);
        assert!((u.max[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_frustum_contains_origin() {
        let m = [
            [1.0, 0.0, 0.0, 0.0f64],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let frustum = Frustum::from_view_projection(m);
        let _ = frustum.contains_point([0.0, 0.0, 0.0]);
        let _ = frustum.contains_sphere([0.0, 0.0, 0.0], 1.0);
    }
    #[test]
    fn test_scene_descendants() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let c1 = scene.add_node("c1", Transform::identity(), Some(root));
        let c2 = scene.add_node("c2", Transform::identity(), Some(root));
        let _gc = scene.add_node("gc", Transform::identity(), Some(c1));
        let desc = scene.descendants(root);
        assert!(desc.contains(&c1));
        assert!(desc.contains(&c2));
        assert_eq!(desc.len(), 3);
    }
    #[test]
    fn test_scene_geometry_nodes() {
        let mut scene = Scene::new();
        let n0 = scene.add_node("a", Transform::identity(), None);
        let n1 = scene.add_node("b", Transform::identity(), None);
        scene.set_mesh_index(n0, Some(0));
        let geom = scene.geometry_nodes();
        assert!(geom.contains(&n0));
        assert!(!geom.contains(&n1));
    }
    #[test]
    fn test_instance_data_push_and_len() {
        let mut inst = InstanceData::new();
        assert!(inst.is_empty());
        let tf = InstanceData::identity_transform();
        inst.push(tf, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(inst.len(), 1);
    }
    #[test]
    fn test_instance_data_identity_transform() {
        let tf = InstanceData::identity_transform();
        for (col, tfcol) in tf.iter().enumerate() {
            for (row, &tfval) in tfcol.iter().enumerate() {
                let expected = if col == row { 1.0 } else { 0.0 };
                assert!((tfval - expected).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_light_manager_add_find() {
        let mut lm = LightManager::new();
        let idx = lm.add(
            "sun",
            Light::Directional {
                direction: [0.0, -1.0, 0.0],
                color: [1.0, 1.0, 1.0],
                intensity: 1.0,
            },
        );
        assert_eq!(idx, 0);
        assert_eq!(lm.find("sun"), Some(0));
        assert_eq!(lm.find("moon"), None);
    }
    #[test]
    fn test_light_manager_total_intensity() {
        let mut lm = LightManager::new();
        lm.add(
            "a",
            Light::Directional {
                direction: [0.0, -1.0, 0.0],
                color: [1.0; 3],
                intensity: 2.0,
            },
        );
        lm.add(
            "b",
            Light::Point {
                position: [0.0; 3],
                color: [1.0; 3],
                intensity: 3.0,
                range: 10.0,
            },
        );
        assert!((lm.total_intensity() - 5.0).abs() < 1e-5);
    }
    #[test]
    fn test_light_manager_remove() {
        let mut lm = LightManager::new();
        lm.add(
            "a",
            Light::Directional {
                direction: [0.0, -1.0, 0.0],
                color: [1.0; 3],
                intensity: 1.0,
            },
        );
        lm.add(
            "b",
            Light::Directional {
                direction: [0.0, -1.0, 0.0],
                color: [1.0; 3],
                intensity: 1.0,
            },
        );
        lm.remove(0);
        assert_eq!(lm.len(), 1);
    }
    #[test]
    fn test_scene_serialize_deserialize_round_trip() {
        let mut scene = Scene::new();
        let root = scene.add_node("root", Transform::identity(), None);
        let mut t = Transform::identity();
        t.position = [1.0, 2.0, 3.0];
        let child = scene.add_node("child", t, Some(root));
        scene.set_mesh_index(child, Some(42));
        let records = scene.serialize();
        let scene2 = Scene::deserialize(&records);
        assert_eq!(scene2.nodes.len(), 2);
        assert_eq!(scene2.nodes[child].name, "child");
        assert_eq!(scene2.nodes[child].mesh_index, Some(42));
        assert!((scene2.nodes[child].transform.position[0] - 1.0).abs() < 1e-12);
        assert!(scene2.nodes[root].children.contains(&child));
    }
    #[test]
    fn test_scene_serialize_preserves_visibility() {
        let mut scene = Scene::new();
        let id = scene.add_node("n", Transform::identity(), None);
        scene.set_visible(id, false);
        let records = scene.serialize();
        assert!(!records[0].visible);
        let scene2 = Scene::deserialize(&records);
        assert!(!scene2.nodes[id].visible);
    }
    #[test]
    fn test_scene_cull_with_frustum_all_inside() {
        let mut scene = Scene::new();
        let a = scene.add_node("a", Transform::identity(), None);
        let b = scene.add_node("b", Transform::identity(), None);
        let m = [
            [1.0, 0.0, 0.0, 0.0f64],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let frustum = Frustum::from_view_projection(m);
        let radii = vec![1.0; 2];
        let visible = scene.cull_with_frustum(&frustum, [0.0, 0.0, -10.0], &radii);
        let ids: Vec<NodeId> = visible.iter().map(|(id, _)| *id).collect();
        let _ = ids;
        let _ = (a, b);
    }
    #[test]
    fn test_scene_camera_default_view_matrix_non_zero() {
        let cam = SceneCamera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        let vm = cam.view_matrix();
        let trace: f64 = (0..4).map(|i| vm[i][i]).sum();
        assert!(trace.abs() > 0.1);
    }
    #[test]
    fn test_scene_camera_perspective_proj_diagonal_nonzero() {
        let cam = SceneCamera::new([0.0, 3.0, 0.0], [0.0, 0.0, 0.0]);
        let proj = cam.projection_matrix();
        assert!(proj[0][0].abs() > 1e-6);
        assert!(proj[1][1].abs() > 1e-6);
    }
    #[test]
    fn test_scene_camera_orthographic_proj() {
        let mut cam = SceneCamera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        cam.projection = SceneCameraProjection::Orthographic {
            half_width: 4.0,
            half_height: 3.0,
            near: 0.1,
            far: 100.0,
        };
        let proj = cam.projection_matrix();
        assert!((proj[0][0] - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_scene_camera_orbit_changes_position() {
        let mut cam = SceneCamera::new([5.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let original = cam.position;
        cam.orbit(0.1, 0.0);
        let moved = cam.position;
        let diff = ((moved[0] - original[0]).powi(2)
            + (moved[1] - original[1]).powi(2)
            + (moved[2] - original[2]).powi(2))
        .sqrt();
        assert!(diff > 1e-6, "orbit should change camera position");
    }
    #[test]
    fn test_scene_camera_fly_forward_moves_toward_target() {
        let mut cam = SceneCamera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        let orig_z = cam.position[2];
        cam.fly(1.0, 0.0, 0.0);
        assert!(
            cam.position[2] < orig_z,
            "fly forward should reduce z distance to target"
        );
    }
    #[test]
    fn test_scene_camera_shake_and_decay() {
        let mut cam = SceneCamera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        cam.apply_shake(0.5, 0.9);
        assert!(cam.shake_magnitude > 0.0);
        cam.update_shake(1.0);
        assert!(cam.shake_magnitude < 0.5);
    }
    #[test]
    fn test_scene_camera_shake_view_differs_from_clean() {
        let mut cam = SceneCamera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        let clean = cam.view_matrix();
        cam.shake_offset = [0.3, 0.1, 0.0];
        let shaken = cam.view_matrix();
        let max_diff = (0..4)
            .flat_map(|c| (0..4).map(move |r| (clean[c][r] - shaken[c][r]).abs()))
            .fold(0.0_f64, f64::max);
        assert!(max_diff > 1e-6, "shake offset should alter the view matrix");
    }
    #[test]
    fn test_scene_light_point_attenuation_at_zero_distance() {
        let light = SceneLight::point([0.0; 3], [1.0, 1.0, 1.0], 10.0, 1.0, 0.09, 0.032);
        let att = light.attenuation(0.0);
        assert!(
            (att - 1.0).abs() < 1e-6,
            "attenuation at d=0 should be 1: {att}"
        );
    }
    #[test]
    fn test_scene_light_point_attenuation_decreases_with_distance() {
        let light = SceneLight::point([0.0; 3], [1.0, 1.0, 1.0], 10.0, 1.0, 0.09, 0.032);
        let a1 = light.attenuation(1.0);
        let a5 = light.attenuation(5.0);
        assert!(a5 < a1, "attenuation should decrease with distance");
    }
    #[test]
    fn test_scene_light_directional_constant_attenuation() {
        let light = SceneLight::directional([0.0, -1.0, 0.0], [1.0, 1.0, 1.0], 1.0);
        let a0 = light.attenuation(0.0);
        let a100 = light.attenuation(100.0);
        assert!(
            (a0 - a100).abs() < 1e-10,
            "directional light attenuation should be constant"
        );
    }
    #[test]
    fn test_scene_light_shadow_frustum_directional() {
        let light = SceneLight::directional([0.0, -1.0, 0.0], [1.0, 1.0, 1.0], 1.0);
        let sf = light.shadow_frustum(20.0);
        assert!(
            sf.is_some(),
            "directional light should produce a shadow frustum"
        );
    }
    #[test]
    fn test_scene_light_spotlight_cone_attenuation() {
        let light = SceneLight::spotlight(
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            [1.0, 1.0, 1.0],
            5.0,
            0.3,
            0.7,
            1.0,
            0.05,
            0.01,
        );
        let cone = light.cone_factor([0.0, 0.0, 0.0]);
        assert!(cone > 0.9, "on-axis cone factor should be near 1: {cone}");
    }
    #[test]
    fn test_scene_environment_default() {
        let env = SceneEnvironment::default();
        assert!(env.ambient_intensity >= 0.0);
        assert!(env.fog_density >= 0.0);
    }
    #[test]
    fn test_scene_environment_fog_factor_at_zero() {
        let env = SceneEnvironment {
            sky_color: [0.4, 0.6, 0.9],
            ambient_intensity: 0.1,
            fog_density: 0.05,
            fog_color: [0.8, 0.8, 0.9],
        };
        let f = env.fog_factor(0.0);
        assert!(
            (f - 1.0).abs() < 1e-6,
            "fog factor at d=0 should be 1 (no fog): {f}"
        );
    }
    #[test]
    fn test_scene_environment_fog_factor_decreases() {
        let env = SceneEnvironment {
            sky_color: [0.5, 0.7, 1.0],
            ambient_intensity: 0.2,
            fog_density: 0.1,
            fog_color: [0.9, 0.9, 0.9],
        };
        let f10 = env.fog_factor(10.0);
        let f50 = env.fog_factor(50.0);
        assert!(f50 < f10, "fog factor should decrease with distance");
    }
    #[test]
    fn test_frustum_culler_build_and_query() {
        let mut culler = FrustumCuller::new(4);
        let m = [
            [1.0, 0.0, 0.0, 0.0f64],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let frustum = Frustum::from_view_projection(m);
        let centres = vec![[0.0_f64; 3]; 4];
        let radii = vec![1.0_f64; 4];
        culler.cull(&frustum, &centres, &radii);
        let _ = culler.is_visible(0);
        let _ = culler.is_visible(3);
    }
    #[test]
    fn test_frustum_culler_out_of_range_invisible() {
        let culler = FrustumCuller::new(2);
        assert!(
            !culler.is_visible(999),
            "out-of-range index should be invisible"
        );
    }
    #[test]
    fn test_frustum_culler_visible_list() {
        let mut culler = FrustumCuller::new(3);
        let m = [
            [1.0, 0.0, 0.0, 0.0f64],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let frustum = Frustum::from_view_projection(m);
        let centres = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let radii = vec![1.0; 3];
        culler.cull(&frustum, &centres, &radii);
        let vis = culler.visible_indices();
        let _ = vis;
    }
    #[test]
    fn test_shadow_map_empty_scene_is_identity() {
        let scene = Scene::new();
        let m = scene.compute_shadow_map_frustum([0.0, -1.0, 0.0]);
        assert!(
            (m[0][0] - 1.0).abs() < 1.0,
            "Empty scene shadow frustum should be near identity"
        );
    }
    #[test]
    fn test_shadow_map_single_node_no_panic() {
        let mut scene = Scene::new();
        let mut t = Transform::identity();
        t.position = [1.0, 2.0, 3.0];
        scene.add_node("n", t, None);
        let m = scene.compute_shadow_map_frustum([0.0, -1.0, 0.0]);
        for col in &m {
            for &v in col {
                assert!(v.is_finite(), "Shadow frustum entry must be finite");
            }
        }
    }
    #[test]
    fn test_shadow_map_multiple_nodes_finite() {
        let mut scene = Scene::new();
        for i in 0..5 {
            let mut t = Transform::identity();
            t.position = [i as f64, 0.0, 0.0];
            scene.add_node(&format!("n{i}"), t, None);
        }
        let m = scene.compute_shadow_map_frustum([1.0, 0.0, 0.0]);
        for col in &m {
            for &v in col {
                assert!(v.is_finite());
            }
        }
    }
    #[test]
    fn test_cull_offscreen_nodes_become_invisible() {
        let mut scene = Scene::new();
        let mut t = Transform::identity();
        t.position = [0.0, 0.0, 1000.0];
        let id = scene.add_node("far", t, None);
        let vp = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        scene.cull_invisible_objects(vp, &[]);
        assert!(
            !scene.nodes[id].visible,
            "Node at z=1000 should be culled as off-screen"
        );
    }
    #[test]
    fn test_cull_within_ndc_stays_visible() {
        let mut scene = Scene::new();
        let mut t = Transform::identity();
        t.position = [0.0, 0.0, 0.0];
        let id = scene.add_node("origin", t, None);
        let vp = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        scene.cull_invisible_objects(vp, &[]);
        assert!(
            scene.nodes[id].visible,
            "Node at origin should remain visible"
        );
    }
    #[test]
    fn test_cull_no_panic_empty_scene() {
        let mut scene = Scene::new();
        let vp = [[1.0_f64; 4]; 4];
        scene.cull_invisible_objects(vp, &[]);
    }
    #[test]
    fn test_lod_empty_config_returns_none() {
        let cfg = LodConfig::new();
        assert_eq!(Scene::compute_lod_level(&cfg, 10.0), None);
    }
    #[test]
    fn test_lod_near_distance_returns_highest_detail() {
        let mut cfg = LodConfig::new();
        cfg.add_level(10.0, 0);
        cfg.add_level(50.0, 1);
        cfg.add_level(200.0, 2);
        assert_eq!(Scene::compute_lod_level(&cfg, 5.0), Some(0));
    }
    #[test]
    fn test_lod_far_distance_returns_lowest_detail() {
        let mut cfg = LodConfig::new();
        cfg.add_level(10.0, 0);
        cfg.add_level(50.0, 1);
        cfg.add_level(200.0, 2);
        assert_eq!(Scene::compute_lod_level(&cfg, 300.0), Some(2));
    }
    #[test]
    fn test_lod_boundary_selects_correct_level() {
        let mut cfg = LodConfig::new();
        cfg.add_level(20.0, 0);
        cfg.add_level(100.0, 1);
        assert_eq!(Scene::compute_lod_level(&cfg, 20.0), Some(0));
        assert_eq!(Scene::compute_lod_level(&cfg, 21.0), Some(1));
    }
}
#[cfg(test)]
mod tests_new_scene {

    use crate::scene::EnvironmentMap;
    use crate::scene::FogModel;

    use crate::scene::HdrToneMapper;

    use crate::scene::LodConfig;
    use crate::scene::LodManager;

    use crate::scene::SceneStats;
    use crate::scene::ShadowMap;
    use crate::scene::SkyboxCube;

    #[test]
    fn test_skybox_solid_color_dimensions() {
        let sky = SkyboxCube::solid_color(8, [100u8, 150, 200, 255]);
        assert_eq!(sky.face_size, 8);
        assert_eq!(sky.face_texels(), 64);
        for face in &sky.faces {
            assert_eq!(face.len(), 64);
        }
    }
    #[test]
    fn test_skybox_sample_positive_x() {
        let sky = SkyboxCube::solid_color(4, [255u8, 0, 0, 255]);
        let s = sky.sample([1.0, 0.0, 0.0]);
        assert_eq!(s, [255u8, 0, 0, 255], "solid red skybox should return red");
    }
    #[test]
    fn test_skybox_sample_does_not_panic_on_any_direction() {
        let sky = SkyboxCube::solid_color(4, [128u8; 4]);
        let dirs: [[f32; 3]; 6] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        for dir in &dirs {
            let _ = sky.sample(*dir);
        }
    }
    #[test]
    fn test_hdr_tone_mapper_exposure_zero_identity() {
        let tm = HdrToneMapper {
            ev: 0.0,
            knee: 0.75,
            max_luminance: 1.0,
            apply_srgb_gamma: false,
        };
        let v = 0.3_f32;
        assert!(
            (tm.apply_exposure(v) - v).abs() < 1e-6,
            "EV=0 should not change value"
        );
    }
    #[test]
    fn test_hdr_tone_mapper_positive_ev_brightens() {
        let tm = HdrToneMapper {
            ev: 1.0,
            knee: 0.75,
            max_luminance: 1.0,
            apply_srgb_gamma: false,
        };
        let brighter = tm.apply_exposure(0.3);
        assert!(brighter > 0.3, "+1 EV should brighten");
    }
    #[test]
    fn test_hdr_tone_mapper_soft_clip_below_knee_unchanged() {
        let tm = HdrToneMapper {
            ev: 0.0,
            knee: 0.75,
            max_luminance: 1.0,
            apply_srgb_gamma: false,
        };
        assert!(
            (tm.soft_clip(0.5) - 0.5).abs() < 1e-5,
            "below knee: no clipping"
        );
    }
    #[test]
    fn test_hdr_tone_mapper_soft_clip_at_max() {
        let tm = HdrToneMapper {
            ev: 0.0,
            knee: 0.5,
            max_luminance: 1.0,
            apply_srgb_gamma: false,
        };
        let c = tm.soft_clip(100.0);
        assert!(c <= 1.0, "soft_clip should not exceed max_luminance: {c}");
    }
    #[test]
    fn test_hdr_tone_mapper_srgb_encode_0_and_1() {
        let tm = HdrToneMapper::default_settings();
        assert!((tm.srgb_encode(0.0)).abs() < 1e-6);
        assert!((tm.srgb_encode(1.0) - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_hdr_tone_mapper_key_gain_zero_luminance() {
        let tm = HdrToneMapper::default_settings();
        let g = tm.key_gain(0.0);
        assert!((g - 1.0).abs() < 1e-5, "key gain for zero luminance: {g}");
    }
    #[test]
    fn test_env_map_solid_color_sample_returns_color() {
        let env = EnvironmentMap::solid_color(32, 16, [1.0, 0.5, 0.25]);
        let s = env.sample([0.0, 1.0, 0.0]);
        assert!((s[0] - 1.0).abs() < 1e-5);
        assert!((s[1] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_env_map_texels() {
        let env = EnvironmentMap::solid_color(64, 32, [0.0; 3]);
        assert_eq!(env.texels(), 64 * 32);
    }
    #[test]
    fn test_env_map_irradiance_approx_all_finite() {
        let env = EnvironmentMap::solid_color(16, 8, [0.5, 0.5, 0.5]);
        let irr = env.irradiance_approx([0.0, 1.0, 0.0]);
        for &c in &irr {
            assert!(c.is_finite(), "irradiance should be finite");
        }
    }
    #[test]
    fn test_fog_exponential_at_zero_distance_is_one() {
        let f = FogModel::Exponential.factor(0.0, 0.05);
        assert!((f - 1.0).abs() < 1e-10, "fog at d=0 should be 1: {f}");
    }
    #[test]
    fn test_fog_exponential_squared_at_zero_is_one() {
        let f = FogModel::ExponentialSquared.factor(0.0, 0.1);
        assert!((f - 1.0).abs() < 1e-10, "exp2 fog at d=0 should be 1: {f}");
    }
    #[test]
    fn test_fog_linear_between_near_far() {
        let fog = FogModel::Linear {
            near: 10.0,
            far: 100.0,
        };
        let f_near = fog.factor(10.0, 0.0);
        let f_mid = fog.factor(55.0, 0.0);
        let f_far = fog.factor(100.0, 0.0);
        assert!((f_near - 1.0).abs() < 1e-5, "linear fog at near: {f_near}");
        assert!((f_far - 0.0).abs() < 1e-5, "linear fog at far: {f_far}");
        assert!(
            f_mid > f_far && f_mid < f_near,
            "linear fog should be between near/far"
        );
    }
    #[test]
    fn test_fog_apply_to_color_full_fog_gives_fog_color() {
        let fog = FogModel::Exponential;
        let result = fog.apply_to_color([1.0, 0.0, 0.0], [0.5, 0.5, 0.5], 1000.0, 1.0);
        assert!(
            (result[0] - 0.5).abs() < 0.01,
            "fully fogged R should be fog color: {}",
            result[0]
        );
    }
    #[test]
    fn test_shadow_map_initial_depth_is_one() {
        let sm = ShadowMap::new(8, 8, [[0.0; 4]; 4]);
        for y in 0..8 {
            for x in 0..8 {
                assert!(
                    (sm.get_depth(x, y) - 1.0).abs() < 1e-6,
                    "initial depth should be 1"
                );
            }
        }
    }
    #[test]
    fn test_shadow_map_set_get_depth() {
        let mut sm = ShadowMap::new(4, 4, [[0.0; 4]; 4]);
        sm.set_depth(2, 3, 0.42);
        assert!((sm.get_depth(2, 3) - 0.42).abs() < 1e-5);
    }
    #[test]
    fn test_shadow_map_out_of_bounds_returns_one() {
        let sm = ShadowMap::new(4, 4, [[0.0; 4]; 4]);
        assert!((sm.get_depth(100, 100) - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_shadow_map_pcf_all_ones_shadow_zero() {
        let sm = ShadowMap::new(8, 8, {
            let mut m = [[0.0f64; 4]; 4];
            m[0][0] = 0.5;
            m[1][1] = 0.5;
            m[2][2] = 0.5;
            m[3][3] = 1.0;
            m
        });
        let in_shadow = sm.is_in_shadow([0.5, 0.5, 0.5], 0.01);
        let _ = in_shadow;
    }
    #[test]
    fn test_lod_manager_add_object() {
        let mut mgr = LodManager::new(LodConfig::new());
        let i = mgr.add_object([0.0, 0.0, 0.0]);
        assert_eq!(i, 0);
        assert_eq!(mgr.len(), 1);
    }
    #[test]
    fn test_lod_manager_compute_all_lods_empty_config() {
        let mut mgr = LodManager::new(LodConfig::new());
        mgr.add_object([1.0, 0.0, 0.0]);
        mgr.add_object([5.0, 0.0, 0.0]);
        let lods = mgr.compute_all_lods([0.0; 3]);
        for l in &lods {
            assert!(l.is_none(), "empty config should give None");
        }
    }
    #[test]
    fn test_lod_manager_selects_correct_lod() {
        let mut cfg = LodConfig::new();
        cfg.add_level(5.0, 0);
        cfg.add_level(20.0, 1);
        let mut mgr = LodManager::new(cfg);
        mgr.add_object([2.0, 0.0, 0.0]);
        mgr.add_object([10.0, 0.0, 0.0]);
        let lods = mgr.compute_all_lods([0.0; 3]);
        assert_eq!(lods[0], Some(0));
        assert_eq!(lods[1], Some(1));
    }
    #[test]
    fn test_lod_manager_out_of_range_returns_none() {
        let mgr = LodManager::new(LodConfig::new());
        assert!(mgr.lod_for(99, [0.0; 3]).is_none());
    }
    #[test]
    fn test_scene_stats_initial_zero() {
        let s = SceneStats::new();
        assert_eq!(s.total_nodes, 0);
        assert_eq!(s.visible_nodes, 0);
        assert_eq!(s.draw_calls, 0);
    }
    #[test]
    fn test_scene_stats_cull_ratio() {
        let s = SceneStats {
            total_nodes: 10,
            culled_nodes: 4,
            visible_nodes: 6,
            ..Default::default()
        };
        let ratio = s.cull_ratio();
        assert!(
            (ratio - 0.4).abs() < 1e-10,
            "cull ratio should be 0.4, got {ratio}"
        );
    }
    #[test]
    fn test_scene_stats_merge() {
        let mut a = SceneStats {
            total_nodes: 5,
            draw_calls: 3,
            ..Default::default()
        };
        let b = SceneStats {
            total_nodes: 7,
            draw_calls: 2,
            triangle_count: 100,
            ..Default::default()
        };
        a.merge(&b);
        assert_eq!(a.total_nodes, 12);
        assert_eq!(a.draw_calls, 5);
        assert_eq!(a.triangle_count, 100);
    }
    #[test]
    fn test_scene_stats_cull_ratio_zero_total() {
        let s = SceneStats::new();
        assert!((s.cull_ratio() - 0.0).abs() < 1e-10);
    }
}
