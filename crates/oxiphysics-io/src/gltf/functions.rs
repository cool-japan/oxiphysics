//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::type_complexity)]
use super::types::{GltfScene, ValidationIssue};

/// Escape special characters in a JSON string value.
pub fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
/// Decode an array of f32 values from a little-endian byte buffer.
pub fn decode_f32_le(bytes: &[u8]) -> Vec<f32> {
    let mut out = Vec::with_capacity(bytes.len() / 4);
    let mut i = 0;
    while i + 3 < bytes.len() {
        let v = f32::from_le_bytes([bytes[i], bytes[i + 1], bytes[i + 2], bytes[i + 3]]);
        out.push(v);
        i += 4;
    }
    out
}
/// Decode an array of u16 values from a little-endian byte buffer.
pub fn decode_u16_le(bytes: &[u8]) -> Vec<u16> {
    let mut out = Vec::with_capacity(bytes.len() / 2);
    let mut i = 0;
    while i + 1 < bytes.len() {
        let v = u16::from_le_bytes([bytes[i], bytes[i + 1]]);
        out.push(v);
        i += 2;
    }
    out
}
/// Decode an array of u32 values from a little-endian byte buffer.
pub fn decode_u32_le(bytes: &[u8]) -> Vec<u32> {
    let mut out = Vec::with_capacity(bytes.len() / 4);
    let mut i = 0;
    while i + 3 < bytes.len() {
        let v = u32::from_le_bytes([bytes[i], bytes[i + 1], bytes[i + 2], bytes[i + 3]]);
        out.push(v);
        i += 4;
    }
    out
}
/// Decode VEC3 data from a buffer of f32 values.
pub fn decode_vec3(floats: &[f32]) -> Vec<[f32; 3]> {
    let mut out = Vec::with_capacity(floats.len() / 3);
    let mut i = 0;
    while i + 2 < floats.len() {
        out.push([floats[i], floats[i + 1], floats[i + 2]]);
        i += 3;
    }
    out
}
/// Decode VEC4 data from a buffer of f32 values.
pub fn decode_vec4(floats: &[f32]) -> Vec<[f32; 4]> {
    let mut out = Vec::with_capacity(floats.len() / 4);
    let mut i = 0;
    while i + 3 < floats.len() {
        out.push([floats[i], floats[i + 1], floats[i + 2], floats[i + 3]]);
        i += 4;
    }
    out
}
/// Write a triangle mesh as a Wavefront OBJ string.
///
/// Produces `v` lines for positions, `vn` lines for normals, and `f` lines
/// using the `v//vn` format (1-indexed). Indices are interpreted as flat
/// per-vertex data (index i -> positions\[i\] and normals\[i\]).
pub fn write_obj(positions: &[[f32; 3]], normals: &[[f32; 3]], indices: &[u32]) -> String {
    let mut out = String::new();
    out.push_str("# OxiPhysics OBJ export\n");
    for p in positions {
        out.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
    }
    for n in normals {
        out.push_str(&format!("vn {} {} {}\n", n[0], n[1], n[2]));
    }
    let mut i = 0;
    while i + 2 < indices.len() {
        let a = indices[i] + 1;
        let b = indices[i + 1] + 1;
        let c = indices[i + 2] + 1;
        out.push_str(&format!("f {a}//{a} {b}//{b} {c}//{c}\n"));
        i += 3;
    }
    out
}
/// Parse a Wavefront OBJ string.
///
/// Returns `(positions, normals, indices)`. Only `v`, `vn`, and `f` records
/// are processed. Face records must use the `v//vn` or `v/vt/vn` format;
/// bare `v` faces use vertex index as normal index.
pub fn parse_obj(content: &str) -> Result<(Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>), String> {
    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut normals: Vec<[f32; 3]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("vn ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let x = parts[1].parse::<f32>().map_err(|e| e.to_string())?;
                let y = parts[2].parse::<f32>().map_err(|e| e.to_string())?;
                let z = parts[3].parse::<f32>().map_err(|e| e.to_string())?;
                normals.push([x, y, z]);
            }
        } else if line.starts_with("v ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                let x = parts[1].parse::<f32>().map_err(|e| e.to_string())?;
                let y = parts[2].parse::<f32>().map_err(|e| e.to_string())?;
                let z = parts[3].parse::<f32>().map_err(|e| e.to_string())?;
                positions.push([x, y, z]);
            }
        } else if line.starts_with("f ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                for token in &parts[1..4] {
                    let vi = token
                        .split('/')
                        .next()
                        .ok_or_else(|| format!("bad face token: {token}"))?
                        .parse::<u32>()
                        .map_err(|e| e.to_string())?;
                    indices.push(vi - 1);
                }
            }
        }
    }
    Ok((positions, normals, indices))
}
/// Compute flat-shading normals for a triangle mesh.
pub fn compute_flat_normals(positions: &[[f32; 3]], indices: &[u32]) -> Vec<[f32; 3]> {
    let mut normals = vec![[0.0f32; 3]; positions.len()];
    let mut i = 0;
    while i + 2 < indices.len() {
        let a = indices[i] as usize;
        let b = indices[i + 1] as usize;
        let c = indices[i + 2] as usize;
        if a < positions.len() && b < positions.len() && c < positions.len() {
            let edge1 = [
                positions[b][0] - positions[a][0],
                positions[b][1] - positions[a][1],
                positions[b][2] - positions[a][2],
            ];
            let edge2 = [
                positions[c][0] - positions[a][0],
                positions[c][1] - positions[a][1],
                positions[c][2] - positions[a][2],
            ];
            let normal = [
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0],
            ];
            let len =
                (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
            let n = if len > 1e-10 {
                [normal[0] / len, normal[1] / len, normal[2] / len]
            } else {
                [0.0, 0.0, 1.0]
            };
            normals[a] = n;
            normals[b] = n;
            normals[c] = n;
        }
        i += 3;
    }
    normals
}
/// Validate a glTF scene for structural correctness.
///
/// Returns `Ok(())` if no errors are found, or `Err(String)` describing the
/// first error encountered.
pub fn validate_scene(scene: &GltfScene) -> Result<(), String> {
    let issues = validate_scene_detailed(scene);
    if let Some(err) = issues.iter().find(|i| i.is_error) {
        Err(err.message.clone())
    } else {
        Ok(())
    }
}
/// Validate a glTF scene and return all issues found (errors and warnings).
pub fn validate_scene_detailed(scene: &GltfScene) -> Vec<ValidationIssue> {
    let mut issues = Vec::new();
    for (ni, node) in scene.nodes.iter().enumerate() {
        if let Some(mesh_idx) = node.mesh
            && mesh_idx >= scene.meshes.len()
        {
            issues
                    .push(
                        ValidationIssue::error(
                            format!(
                                "Node {ni} ('{}') references mesh {mesh_idx} which does not exist (scene has {} meshes)",
                                node.name, scene.meshes.len()
                            ),
                        ),
                    );
        }
        for &child_idx in &node.children {
            if child_idx >= scene.nodes.len() {
                issues
                    .push(
                        ValidationIssue::error(
                            format!(
                                "Node {ni} ('{}') has child index {child_idx} which does not exist (scene has {} nodes)",
                                node.name, scene.nodes.len()
                            ),
                        ),
                    );
            }
        }
    }
    for (mi, mesh) in scene.meshes.iter().enumerate() {
        for (pi, prim) in mesh.primitives.iter().enumerate() {
            if prim.positions.is_empty() {
                issues.push(ValidationIssue::warning(format!(
                    "Mesh {mi} ('{}'), primitive {pi}: has no positions (degenerate primitive)",
                    mesh.name
                )));
            }
            if !prim.positions.is_empty() {
                for (ii, &idx) in prim.indices.iter().enumerate() {
                    if idx as usize >= prim.positions.len() {
                        issues.push(ValidationIssue::error(format!(
                            "Mesh {mi}, primitive {pi}, index {ii}: index {idx} >= vertex count {}",
                            prim.positions.len()
                        )));
                        break;
                    }
                }
            }
            if !prim.normals.is_empty() && prim.normals.len() != prim.positions.len() {
                issues.push(ValidationIssue::warning(format!(
                    "Mesh {mi}, primitive {pi}: normal count ({}) != position count ({})",
                    prim.normals.len(),
                    prim.positions.len()
                )));
            }
        }
    }
    for (mi, mat) in scene.materials.iter().enumerate() {
        match mat.alpha_mode.as_str() {
            "OPAQUE" | "MASK" | "BLEND" => {}
            other => {
                issues.push(ValidationIssue::error(format!(
                    "Material {mi} ('{}') has invalid alphaMode '{other}'",
                    mat.name
                )));
            }
        }
        for (ci, &c) in mat.base_color_factor.iter().enumerate() {
            if !(0.0..=1.0).contains(&c) {
                issues.push(ValidationIssue::warning(format!(
                    "Material {mi} ('{}') baseColorFactor[{ci}] = {c} is outside [0,1]",
                    mat.name
                )));
            }
        }
    }
    for (ai, anim) in scene.animations.iter().enumerate() {
        for (ci, ch) in anim.channels.iter().enumerate() {
            let target_node = ch.target.node;
            if target_node >= scene.nodes.len() {
                issues.push(ValidationIssue::error(format!(
                    "Animation {ai} ('{}'), channel {ci}: target node {target_node} does not exist",
                    anim.name
                )));
            }
            if ch.input_times.is_empty() {
                issues.push(ValidationIssue::warning(format!(
                    "Animation {ai} ('{}'), channel {ci}: no keyframes",
                    anim.name
                )));
            }
        }
    }
    issues
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::gltf::CameraType;

    use crate::gltf::GlbWriter;
    use crate::gltf::GltfAccessor;
    use crate::gltf::GltfAnimation;
    use crate::gltf::GltfAnimationChannel;
    use crate::gltf::GltfAnimationTarget;
    use crate::gltf::GltfBufferView;
    use crate::gltf::GltfMaterial;
    use crate::gltf::GltfMesh;
    use crate::gltf::GltfNode;
    use crate::gltf::GltfPrimitive;

    use crate::gltf::Joint;
    use crate::gltf::JointChannel;
    use crate::gltf::LightType;
    use crate::gltf::MorphPrimitive;
    use crate::gltf::MorphTarget;

    use crate::gltf::SceneCamera;
    use crate::gltf::SceneLight;
    use crate::gltf::Skeleton;
    use crate::gltf::SkeletonAnimation;

    #[test]
    fn test_gltf_scene_to_json_contains_asset_and_meshes() {
        let mut scene = GltfScene::new();
        let prim = GltfPrimitive {
            positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normals: vec![[0.0, 0.0, 1.0]; 3],
            texcoords: vec![],
            indices: vec![0, 1, 2],
        };
        scene.add_mesh(GltfMesh {
            name: "TestMesh".to_string(),
            primitives: vec![prim],
        });
        let json = scene.to_json();
        assert!(json.contains("\"asset\""), "JSON must contain 'asset' key");
        assert!(
            json.contains("\"meshes\""),
            "JSON must contain 'meshes' key"
        );
        assert!(json.contains("\"version\""), "JSON must contain version");
        assert!(json.contains("2.0"), "version must be 2.0");
    }
    #[test]
    fn test_add_mesh_returns_correct_index() {
        let mut scene = GltfScene::new();
        let idx0 = scene.add_mesh(GltfMesh {
            name: "Mesh0".to_string(),
            primitives: vec![],
        });
        let idx1 = scene.add_mesh(GltfMesh {
            name: "Mesh1".to_string(),
            primitives: vec![],
        });
        assert_eq!(idx0, 0);
        assert_eq!(idx1, 1);
        assert_eq!(scene.mesh_count(), 2);
    }
    #[test]
    fn test_write_obj_single_triangle() {
        let positions = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let normals = vec![[0.0f32, 0.0, 1.0]; 3];
        let indices = vec![0u32, 1, 2];
        let obj = write_obj(&positions, &normals, &indices);
        let v_lines: Vec<&str> = obj.lines().filter(|l| l.starts_with("v ")).collect();
        let vn_lines: Vec<&str> = obj.lines().filter(|l| l.starts_with("vn ")).collect();
        let f_lines: Vec<&str> = obj.lines().filter(|l| l.starts_with("f ")).collect();
        assert_eq!(v_lines.len(), 3, "expected 3 vertex lines");
        assert_eq!(vn_lines.len(), 3, "expected 3 normal lines");
        assert_eq!(f_lines.len(), 1, "expected 1 face line");
        assert!(f_lines[0].contains("1//1"), "face must use v//vn format");
    }
    #[test]
    fn test_parse_obj_roundtrip() {
        let positions = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let normals = vec![[0.0f32, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]];
        let indices = vec![0u32, 1, 2];
        let obj_str = write_obj(&positions, &normals, &indices);
        let (parsed_pos, parsed_nrm, parsed_idx) = parse_obj(&obj_str).unwrap();
        assert_eq!(parsed_pos.len(), 3);
        assert_eq!(parsed_nrm.len(), 3);
        assert_eq!(parsed_idx, vec![0, 1, 2]);
        assert!((parsed_pos[1][0] - 1.0).abs() < 1e-6);
        assert!((parsed_nrm[0][2] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_gltf_node_default_translation() {
        let node = GltfNode::default();
        assert_eq!(node.translation, [0.0, 0.0, 0.0]);
        assert_eq!(node.scale, [1.0, 1.0, 1.0]);
        assert_eq!(node.rotation, [0.0, 0.0, 0.0, 1.0]);
        assert!(node.mesh.is_none());
        assert!(node.children.is_empty());
    }
    #[test]
    fn test_primitive_bounding_box() {
        let prim = GltfPrimitive {
            positions: vec![[-1.0, -2.0, -3.0], [4.0, 5.0, 6.0], [0.0, 0.0, 0.0]],
            normals: vec![[0.0, 0.0, 1.0]; 3],
            texcoords: vec![],
            indices: vec![0, 1, 2],
        };
        let (min, max) = prim.bounding_box();
        assert!((min[0] - (-1.0)).abs() < 1e-6);
        assert!((min[1] - (-2.0)).abs() < 1e-6);
        assert!((min[2] - (-3.0)).abs() < 1e-6);
        assert!((max[0] - 4.0).abs() < 1e-6);
        assert!((max[1] - 5.0).abs() < 1e-6);
        assert!((max[2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_primitive_triangle_and_vertex_count() {
        let prim = GltfPrimitive {
            positions: vec![[0.0; 3]; 6],
            normals: vec![[0.0, 0.0, 1.0]; 6],
            texcoords: vec![],
            indices: vec![0, 1, 2, 3, 4, 5],
        };
        assert_eq!(prim.triangle_count(), 2);
        assert_eq!(prim.vertex_count(), 6);
    }
    #[test]
    fn test_extract_triangles() {
        let prim = GltfPrimitive {
            positions: vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            normals: vec![[0.0, 0.0, 1.0]; 3],
            texcoords: vec![],
            indices: vec![0, 1, 2],
        };
        let tris = prim.extract_triangles();
        assert_eq!(tris.len(), 1);
        assert!((tris[0][0][0] - 1.0).abs() < 1e-6);
        assert!((tris[0][1][1] - 1.0).abs() < 1e-6);
        assert!((tris[0][2][2] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_mesh_total_counts() {
        let mesh = GltfMesh {
            name: "test".into(),
            primitives: vec![
                GltfPrimitive {
                    positions: vec![[0.0; 3]; 3],
                    normals: vec![[0.0, 0.0, 1.0]; 3],
                    texcoords: vec![],
                    indices: vec![0, 1, 2],
                },
                GltfPrimitive {
                    positions: vec![[0.0; 3]; 4],
                    normals: vec![[0.0, 0.0, 1.0]; 4],
                    texcoords: vec![],
                    indices: vec![0, 1, 2, 1, 2, 3],
                },
            ],
        };
        assert_eq!(mesh.total_vertex_count(), 7);
        assert_eq!(mesh.total_triangle_count(), 3);
    }
    #[test]
    fn test_accessor_components_and_size() {
        let acc = GltfAccessor {
            buffer_view: 0,
            byte_offset: 0,
            component_type: 5126,
            count: 10,
            element_type: "VEC3".to_string(),
        };
        assert_eq!(acc.components_per_element(), 3);
        assert_eq!(acc.component_size(), 4);
        assert_eq!(acc.byte_length(), 10 * 3 * 4);
    }
    #[test]
    fn test_accessor_scalar_u16() {
        let acc = GltfAccessor {
            buffer_view: 0,
            byte_offset: 0,
            component_type: 5123,
            count: 5,
            element_type: "SCALAR".to_string(),
        };
        assert_eq!(acc.components_per_element(), 1);
        assert_eq!(acc.component_size(), 2);
        assert_eq!(acc.byte_length(), 10);
    }
    #[test]
    fn test_accessor_mat4() {
        let acc = GltfAccessor {
            buffer_view: 0,
            byte_offset: 0,
            component_type: 5126,
            count: 1,
            element_type: "MAT4".to_string(),
        };
        assert_eq!(acc.components_per_element(), 16);
        assert_eq!(acc.byte_length(), 64);
    }
    #[test]
    fn test_decode_f32_le() {
        let val: f32 = 3.125;
        let bytes = val.to_le_bytes();
        let decoded = decode_f32_le(&bytes);
        assert_eq!(decoded.len(), 1);
        assert!((decoded[0] - 3.125).abs() < 1e-5);
    }
    #[test]
    fn test_decode_u16_le() {
        let val: u16 = 12345;
        let bytes = val.to_le_bytes();
        let decoded = decode_u16_le(&bytes);
        assert_eq!(decoded, vec![12345]);
    }
    #[test]
    fn test_decode_u32_le() {
        let val: u32 = 999999;
        let bytes = val.to_le_bytes();
        let decoded = decode_u32_le(&bytes);
        assert_eq!(decoded, vec![999999]);
    }
    #[test]
    fn test_decode_vec3() {
        let floats = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let vecs = decode_vec3(&floats);
        assert_eq!(vecs.len(), 2);
        assert_eq!(vecs[0], [1.0, 2.0, 3.0]);
        assert_eq!(vecs[1], [4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_decode_vec4() {
        let floats = vec![1.0, 2.0, 3.0, 4.0];
        let vecs = decode_vec4(&floats);
        assert_eq!(vecs.len(), 1);
        assert_eq!(vecs[0], [1.0, 2.0, 3.0, 4.0]);
    }
    #[test]
    fn test_animation_duration() {
        let anim = GltfAnimation {
            name: "walk".into(),
            channels: vec![
                GltfAnimationChannel {
                    target: GltfAnimationTarget {
                        node: 0,
                        path: "translation".into(),
                    },
                    interpolation: "LINEAR".into(),
                    input_times: vec![0.0, 0.5, 1.0],
                    output_values: vec![0.0; 9],
                },
                GltfAnimationChannel {
                    target: GltfAnimationTarget {
                        node: 1,
                        path: "rotation".into(),
                    },
                    interpolation: "LINEAR".into(),
                    input_times: vec![0.0, 2.0],
                    output_values: vec![0.0; 8],
                },
            ],
        };
        assert!((anim.duration() - 2.0).abs() < 1e-6);
        assert_eq!(anim.total_keyframe_count(), 5);
    }
    #[test]
    fn test_animation_empty() {
        let anim = GltfAnimation {
            name: "empty".into(),
            channels: vec![],
        };
        assert!((anim.duration()).abs() < 1e-6);
        assert_eq!(anim.total_keyframe_count(), 0);
    }
    #[test]
    fn test_material_default() {
        let mat = GltfMaterial::default();
        assert!(mat.is_opaque());
        assert!(!mat.is_emissive());
        assert_eq!(mat.alpha_mode, "OPAQUE");
        assert!((mat.metallic_factor - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_material_emissive() {
        let mat = GltfMaterial {
            emissive_factor: [1.0, 0.0, 0.0],
            ..Default::default()
        };
        assert!(mat.is_emissive());
    }
    #[test]
    fn test_material_blend() {
        let mut mat = GltfMaterial {
            alpha_mode: "BLEND".to_string(),
            ..Default::default()
        };
        mat.base_color_factor[3] = 0.5;
        assert!(!mat.is_opaque());
    }
    #[test]
    fn test_scene_graph_traversal() {
        let mut scene = GltfScene::new();
        scene.add_node(GltfNode {
            name: "root".into(),
            translation: [1.0, 0.0, 0.0],
            children: vec![1],
            ..GltfNode::default()
        });
        scene.add_node(GltfNode {
            name: "child".into(),
            translation: [0.0, 2.0, 0.0],
            children: vec![2],
            ..GltfNode::default()
        });
        scene.add_node(GltfNode {
            name: "grandchild".into(),
            translation: [0.0, 0.0, 3.0],
            ..GltfNode::default()
        });
        let mut visited = Vec::new();
        scene.traverse_depth_first(|idx, depth, acc_trans| {
            visited.push((idx, depth, acc_trans));
        });
        assert_eq!(visited.len(), 3);
        assert_eq!(visited[0].0, 0);
        assert_eq!(visited[0].1, 0);
        assert!((visited[0].2[0] - 1.0).abs() < 1e-12);
        assert_eq!(visited[1].0, 1);
        assert_eq!(visited[1].1, 1);
        assert!((visited[1].2[0] - 1.0).abs() < 1e-12);
        assert!((visited[1].2[1] - 2.0).abs() < 1e-12);
        assert_eq!(visited[2].0, 2);
        assert_eq!(visited[2].1, 2);
        assert!((visited[2].2[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_collect_mesh_primitives() {
        let mut scene = GltfScene::new();
        let mesh_idx = scene.add_mesh(GltfMesh {
            name: "m0".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0; 3]; 3],
                normals: vec![[0.0, 0.0, 1.0]; 3],
                texcoords: vec![],
                indices: vec![0, 1, 2],
            }],
        });
        scene.add_node(GltfNode {
            name: "node0".into(),
            mesh: Some(mesh_idx),
            ..GltfNode::default()
        });
        let prims = scene.collect_mesh_primitives();
        assert_eq!(prims.len(), 1);
        assert_eq!(prims[0].0, "node0");
        assert_eq!(prims[0].1.vertex_count(), 3);
    }
    #[test]
    fn test_nodes_using_mesh() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "m".into(),
            primitives: vec![],
        });
        scene.add_node(GltfNode {
            name: "a".into(),
            mesh: Some(0),
            ..GltfNode::default()
        });
        scene.add_node(GltfNode {
            name: "b".into(),
            mesh: None,
            ..GltfNode::default()
        });
        scene.add_node(GltfNode {
            name: "c".into(),
            mesh: Some(0),
            ..GltfNode::default()
        });
        let users = scene.nodes_using_mesh(0);
        assert_eq!(users, vec![0, 2]);
    }
    #[test]
    fn test_scene_total_counts() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "m".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0; 3]; 4],
                normals: vec![[0.0, 0.0, 1.0]; 4],
                texcoords: vec![],
                indices: vec![0, 1, 2, 1, 2, 3],
            }],
        });
        assert_eq!(scene.total_vertex_count(), 4);
        assert_eq!(scene.total_triangle_count(), 2);
    }
    #[test]
    fn test_compute_flat_normals() {
        let positions = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let indices = vec![0u32, 1, 2];
        let normals = compute_flat_normals(&positions, &indices);
        assert_eq!(normals.len(), 3);
        for n in &normals {
            assert!((n[2] - 1.0).abs() < 1e-5);
        }
    }
    #[test]
    fn test_json_with_materials() {
        let mut scene = GltfScene::new();
        scene.add_node(GltfNode::default());
        scene.add_mesh(GltfMesh {
            name: "m".into(),
            primitives: vec![],
        });
        scene.add_material(GltfMaterial {
            name: "red".into(),
            base_color_factor: [1.0, 0.0, 0.0, 1.0],
            ..GltfMaterial::default()
        });
        let json = scene.to_json();
        assert!(json.contains("\"materials\""));
        assert!(json.contains("\"red\""));
        assert!(json.contains("pbrMetallicRoughness"));
    }
    #[test]
    fn test_json_with_children() {
        let mut scene = GltfScene::new();
        scene.add_node(GltfNode {
            name: "parent".into(),
            children: vec![1],
            ..GltfNode::default()
        });
        scene.add_node(GltfNode {
            name: "child".into(),
            ..GltfNode::default()
        });
        let json = scene.to_json();
        assert!(json.contains("\"children\""));
    }
    #[test]
    fn test_buffer_view() {
        let bv = GltfBufferView {
            buffer: 0,
            byte_offset: 100,
            byte_length: 240,
            byte_stride: Some(12),
        };
        assert_eq!(bv.byte_offset, 100);
        assert_eq!(bv.byte_stride, Some(12));
    }
    #[test]
    fn test_morph_target_blend() {
        let base = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let target = MorphTarget {
            name: "smile".into(),
            position_deltas: vec![[0.0, 0.1, 0.0], [0.0, 0.2, 0.0], [0.0, 0.15, 0.0]],
        };
        let blended = target.apply(&base, 0.5);
        assert_eq!(blended.len(), 3);
        assert!((blended[0][1] - 0.05).abs() < 1e-6, "y should be 0.05");
        assert!((blended[1][0] - 1.0).abs() < 1e-6, "x unchanged");
    }
    #[test]
    fn test_morph_target_zero_weight() {
        let base = vec![[1.0f32, 2.0, 3.0]];
        let target = MorphTarget {
            name: "t".into(),
            position_deltas: vec![[10.0, 10.0, 10.0]],
        };
        let result = target.apply(&base, 0.0);
        assert!((result[0][0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_morph_target_full_weight() {
        let base = vec![[1.0f32, 2.0, 3.0]];
        let target = MorphTarget {
            name: "t".into(),
            position_deltas: vec![[1.0, -1.0, 0.5]],
        };
        let result = target.apply(&base, 1.0);
        assert!((result[0][0] - 2.0).abs() < 1e-6);
        assert!((result[0][1] - 1.0).abs() < 1e-6);
        assert!((result[0][2] - 3.5).abs() < 1e-6);
    }
    #[test]
    fn test_primitive_morph_targets_blend() {
        let mut prim = MorphPrimitive {
            base: GltfPrimitive {
                positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                normals: vec![[0.0, 0.0, 1.0]; 3],
                texcoords: vec![],
                indices: vec![0, 1, 2],
            },
            targets: vec![
                MorphTarget {
                    name: "A".into(),
                    position_deltas: vec![[0.0, 1.0, 0.0]; 3],
                },
                MorphTarget {
                    name: "B".into(),
                    position_deltas: vec![[0.0, 0.0, 2.0]; 3],
                },
            ],
            weights: vec![0.5, 0.0],
        };
        let blended = prim.blend();
        assert!((blended[0][1] - 0.5).abs() < 1e-6);
        prim.set_weight(1, 1.0);
        let blended2 = prim.blend();
        assert!((blended2[0][2] - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_joint_default() {
        let j = Joint::default();
        assert_eq!(j.translation, [0.0, 0.0, 0.0]);
        assert_eq!(j.rotation, [0.0, 0.0, 0.0, 1.0]);
        assert_eq!(j.scale, [1.0, 1.0, 1.0]);
        assert!(j.children.is_empty());
    }
    #[test]
    fn test_skeleton_add_and_count() {
        let mut skel = Skeleton::new();
        skel.add_joint(Joint {
            name: "root".into(),
            ..Joint::default()
        });
        skel.add_joint(Joint {
            name: "hip".into(),
            ..Joint::default()
        });
        assert_eq!(skel.joint_count(), 2);
        assert_eq!(skel.joints[0].name, "root");
    }
    #[test]
    fn test_skeleton_animation_export_json() {
        let mut skel = Skeleton::new();
        skel.add_joint(Joint {
            name: "bone0".into(),
            ..Joint::default()
        });
        let anim = SkeletonAnimation {
            name: "run".into(),
            joint_channels: vec![JointChannel {
                joint_index: 0,
                path: "rotation".into(),
                times: vec![0.0, 0.5, 1.0],
                values: vec![
                    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.707, 0.707, 0.0, 0.0, 0.0, 1.0,
                ],
                interpolation: "LINEAR".into(),
            }],
        };
        let json = skel.export_animation_json(&anim);
        assert!(json.contains("run"), "should contain animation name");
        assert!(json.contains("rotation"), "should contain channel path");
        assert!(json.contains("bone0"), "should contain joint name");
    }
    #[test]
    fn test_skeleton_animation_duration() {
        let anim = SkeletonAnimation {
            name: "walk".into(),
            joint_channels: vec![
                JointChannel {
                    joint_index: 0,
                    path: "translation".into(),
                    times: vec![0.0, 1.5],
                    values: vec![0.0; 6],
                    interpolation: "LINEAR".into(),
                },
                JointChannel {
                    joint_index: 1,
                    path: "rotation".into(),
                    times: vec![0.0, 0.5],
                    values: vec![0.0; 8],
                    interpolation: "STEP".into(),
                },
            ],
        };
        assert!((anim.duration() - 1.5).abs() < 1e-6);
    }
    #[test]
    fn test_glb_magic_and_version() {
        let writer = GlbWriter::new();
        let scene = GltfScene::new();
        let bytes = writer.write(&scene);
        assert!(bytes.len() >= 12, "GLB must be at least 12 bytes");
        assert_eq!(&bytes[0..4], b"glTF", "magic bytes must be 'glTF'");
        let version = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
        assert_eq!(version, 2, "GLB version must be 2");
    }
    #[test]
    fn test_glb_header_length_consistent() {
        let writer = GlbWriter::new();
        let scene = GltfScene::new();
        let bytes = writer.write(&scene);
        let total_len = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
        assert_eq!(
            total_len as usize,
            bytes.len(),
            "header length field must equal actual length"
        );
    }
    #[test]
    fn test_glb_with_mesh_round_trip() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "cube".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                normals: vec![[0.0, 0.0, 1.0]; 3],
                texcoords: vec![],
                indices: vec![0, 1, 2],
            }],
        });
        let writer = GlbWriter::new();
        let bytes = writer.write(&scene);
        assert!(bytes.len() > 28, "GLB must have JSON chunk");
        let json_chunk_type = u32::from_le_bytes([bytes[16], bytes[17], bytes[18], bytes[19]]);
        assert_eq!(json_chunk_type, 0x4E4F534A, "first chunk must be JSON");
    }
    #[test]
    fn test_camera_node_json() {
        let mut scene = GltfScene::new();
        let cam = SceneCamera {
            name: "main_cam".into(),
            camera_type: CameraType::Perspective {
                yfov: 0.785,
                znear: 0.01,
                zfar: Some(1000.0),
                aspect_ratio: Some(1.778),
            },
        };
        scene.add_camera(cam);
        let json = scene.to_json_with_hierarchy();
        assert!(json.contains("cameras"), "JSON should contain 'cameras'");
        assert!(
            json.contains("perspective"),
            "JSON should contain 'perspective'"
        );
    }
    #[test]
    fn test_light_node_json() {
        let mut scene = GltfScene::new();
        let light = SceneLight {
            name: "sun".into(),
            light_type: LightType::Directional,
            color: [1.0, 0.98, 0.9],
            intensity: 100_000.0,
        };
        scene.add_light(light);
        let json = scene.to_json_with_hierarchy();
        assert!(
            json.contains("extensions"),
            "should contain 'extensions' for lights"
        );
        assert!(
            json.contains("KHR_lights_punctual"),
            "should contain KHR extension"
        );
        assert!(json.contains("directional"), "should contain light type");
    }
    #[test]
    fn test_validate_empty_scene_passes() {
        let scene = GltfScene::new();
        let result = validate_scene(&scene);
        assert!(result.is_ok(), "empty scene should pass validation");
    }
    #[test]
    fn test_validate_dangling_node_mesh_reference_fails() {
        let mut scene = GltfScene::new();
        scene.add_node(GltfNode {
            name: "bad".into(),
            mesh: Some(5),
            ..GltfNode::default()
        });
        let result = validate_scene(&scene);
        assert!(
            result.is_err(),
            "dangling mesh reference should fail validation"
        );
    }
    #[test]
    fn test_validate_valid_mesh_ref_passes() {
        let mut scene = GltfScene::new();
        let idx = scene.add_mesh(GltfMesh {
            name: "m".into(),
            primitives: vec![],
        });
        scene.add_node(GltfNode {
            name: "good".into(),
            mesh: Some(idx),
            ..GltfNode::default()
        });
        let result = validate_scene(&scene);
        assert!(result.is_ok(), "valid mesh reference should pass");
    }
    #[test]
    fn test_validate_dangling_child_reference_fails() {
        let mut scene = GltfScene::new();
        scene.add_node(GltfNode {
            name: "parent".into(),
            children: vec![99],
            ..GltfNode::default()
        });
        let result = validate_scene(&scene);
        assert!(result.is_err(), "dangling child reference should fail");
    }
    #[test]
    fn test_validate_degenerate_primitive_warns() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "degen".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![],
                normals: vec![],
                texcoords: vec![],
                indices: vec![],
            }],
        });
        let issues = validate_scene_detailed(&scene);
        assert!(
            !issues.is_empty(),
            "degenerate primitive should produce issues"
        );
    }
    #[test]
    fn test_validate_index_out_of_range_fails() {
        let scene_with_bad_prim = GltfScene {
            nodes: vec![],
            meshes: vec![GltfMesh {
                name: "m".into(),
                primitives: vec![GltfPrimitive {
                    positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                    normals: vec![[0.0, 0.0, 1.0]; 2],
                    texcoords: vec![],
                    indices: vec![0, 1, 5],
                }],
            }],
            accessors: vec![],
            buffer_views: vec![],
            animations: vec![],
            materials: vec![],
            cameras: vec![],
            lights: vec![],
        };
        let issues = validate_scene_detailed(&scene_with_bad_prim);
        assert!(
            !issues.is_empty(),
            "out-of-range index should produce issues"
        );
    }
}
/// Produce a minimal but complete glTF 2.0 JSON document from a [`GltfScene`].
///
/// The output is a valid JSON string that can be saved as `.gltf`.
/// Buffer data is embedded as base64 data URIs when `embed_buffers` is true.
pub fn scene_to_gltf_json(scene: &GltfScene, embed_buffers: bool) -> String {
    let mut buffer_data: Vec<u8> = Vec::new();
    let mut buffer_views_json: Vec<String> = Vec::new();
    let mut accessors_json: Vec<String> = Vec::new();
    let mut mesh_json_items: Vec<String> = Vec::new();
    for mesh in &scene.meshes {
        for prim in &mesh.primitives {
            let pos_bv_start = buffer_data.len();
            for p in &prim.positions {
                for c in p {
                    buffer_data.extend_from_slice(&c.to_le_bytes());
                }
            }
            let pos_bv_len = buffer_data.len() - pos_bv_start;
            let pos_bv_idx = buffer_views_json.len();
            buffer_views_json.push(format!(
                r#"{{ "buffer": 0, "byteOffset": {start}, "byteLength": {len}, "target": 34962 }}"#,
                start = pos_bv_start,
                len = pos_bv_len
            ));
            let pos_acc_idx = accessors_json.len();
            accessors_json
                .push(
                    format!(
                        r#"{{ "bufferView": {bv}, "byteOffset": 0, "componentType": 5126, "type": "VEC3", "count": {count} }}"#,
                        bv = pos_bv_idx, count = prim.positions.len()
                    ),
                );
            let idx_bv_start = buffer_data.len();
            while !buffer_data.len().is_multiple_of(4) {
                buffer_data.push(0u8);
            }
            let idx_bv_start = if buffer_data.len() == idx_bv_start {
                idx_bv_start
            } else {
                buffer_data.len()
            };
            let _ = idx_bv_start;
            let idx_real_start = buffer_data.len();
            for &i in &prim.indices {
                buffer_data.extend_from_slice(&i.to_le_bytes());
            }
            let idx_bv_len = buffer_data.len() - idx_real_start;
            let idx_bv_idx = buffer_views_json.len();
            buffer_views_json.push(format!(
                r#"{{ "buffer": 0, "byteOffset": {start}, "byteLength": {len}, "target": 34963 }}"#,
                start = idx_real_start,
                len = idx_bv_len
            ));
            let idx_acc_idx = accessors_json.len();
            accessors_json
                .push(
                    format!(
                        r#"{{ "bufferView": {bv}, "byteOffset": 0, "componentType": 5125, "type": "SCALAR", "count": {count} }}"#,
                        bv = idx_bv_idx, count = prim.indices.len()
                    ),
                );
            let prim_json = format!(
                r#"{{ "attributes": {{ "POSITION": {pa} }}, "indices": {ia} }}"#,
                pa = pos_acc_idx,
                ia = idx_acc_idx
            );
            mesh_json_items.push(format!(
                r#"{{ "name": "{}", "primitives": [{}] }}"#,
                mesh.name, prim_json
            ));
        }
    }
    let buf_json = if embed_buffers {
        let b64 = base64_encode(&buffer_data);
        format!(
            r#"{{ "uri": "data:application/octet-stream;base64,{}", "byteLength": {} }}"#,
            b64,
            buffer_data.len()
        )
    } else {
        format!(
            r#"{{ "uri": "buffer.bin", "byteLength": {} }}"#,
            buffer_data.len()
        )
    };
    let bvs = buffer_views_json.join(", ");
    let accs = accessors_json.join(", ");
    let meshes = mesh_json_items.join(", ");
    format!(
        r#"{{
  "asset": {{ "version": "2.0", "generator": "OxiPhysics" }},
  "scene": 0,
  "scenes": [{{ "nodes": [] }}],
  "nodes": [],
  "meshes": [{meshes}],
  "accessors": [{accs}],
  "bufferViews": [{bvs}],
  "buffers": [{buf}]
}}"#,
        meshes = meshes,
        accs = accs,
        bvs = bvs,
        buf = buf_json,
    )
}
/// Minimal base64 encoder (no external dep).
pub(super) fn base64_encode(data: &[u8]) -> String {
    pub(super) const CHARS: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    let mut i = 0;
    while i < data.len() {
        let b0 = data[i] as u32;
        let b1 = if i + 1 < data.len() {
            data[i + 1] as u32
        } else {
            0
        };
        let b2 = if i + 2 < data.len() {
            data[i + 2] as u32
        } else {
            0
        };
        let n = (b0 << 16) | (b1 << 8) | b2;
        out.push(CHARS[((n >> 18) & 0x3F) as usize] as char);
        out.push(CHARS[((n >> 12) & 0x3F) as usize] as char);
        if i + 1 < data.len() {
            out.push(CHARS[((n >> 6) & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
        if i + 2 < data.len() {
            out.push(CHARS[(n & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
        i += 3;
    }
    out
}
#[cfg(test)]
mod tests_gltf_additions {
    use super::*;
    use crate::gltf::AccessorType;
    use crate::gltf::AnimationChannelBuilder;
    use crate::gltf::ComponentType;
    use crate::gltf::GlbWriter;
    use crate::gltf::GltfMesh;
    use crate::gltf::GltfPrimitive;
    use crate::gltf::Interpolation;
    use crate::gltf::PbrMaterialBuilder;
    use crate::gltf::TypedAccessor;

    #[test]
    fn test_component_type_byte_size() {
        assert_eq!(ComponentType::UnsignedByte.byte_size(), 1);
        assert_eq!(ComponentType::UnsignedShort.byte_size(), 2);
        assert_eq!(ComponentType::UnsignedInt.byte_size(), 4);
        assert_eq!(ComponentType::Float.byte_size(), 4);
    }
    #[test]
    fn test_accessor_type_num_components() {
        assert_eq!(AccessorType::Scalar.num_components(), 1);
        assert_eq!(AccessorType::Vec3.num_components(), 3);
        assert_eq!(AccessorType::Vec4.num_components(), 4);
        assert_eq!(AccessorType::Mat4.num_components(), 16);
    }
    #[test]
    fn test_accessor_type_as_str() {
        assert_eq!(AccessorType::Vec3.as_str(), "VEC3");
        assert_eq!(AccessorType::Mat4.as_str(), "MAT4");
        assert_eq!(AccessorType::Scalar.as_str(), "SCALAR");
    }
    #[test]
    fn test_component_type_code() {
        assert_eq!(ComponentType::Float.component_type_code(), 5126);
        assert_eq!(ComponentType::UnsignedInt.component_type_code(), 5125);
    }
    #[test]
    fn test_typed_accessor_byte_length() {
        let acc = TypedAccessor::new("pos", 0, 0, ComponentType::Float, AccessorType::Vec3, 10);
        assert_eq!(acc.byte_length(), 120);
    }
    #[test]
    fn test_typed_accessor_to_json_contains_type() {
        let acc = TypedAccessor::new("vel", 1, 0, ComponentType::Float, AccessorType::Vec3, 5);
        let json = acc.to_json();
        assert!(json.contains("VEC3"), "JSON should contain type");
        assert!(json.contains("5126"), "JSON should contain componentType");
    }
    #[test]
    fn test_typed_accessor_decode_f32() {
        let values: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut buffer: Vec<u8> = Vec::new();
        for v in &values {
            buffer.extend_from_slice(&v.to_le_bytes());
        }
        let acc = TypedAccessor::new("v", 0, 0, ComponentType::Float, AccessorType::Vec3, 2);
        let decoded = acc.decode_f32(&buffer).unwrap();
        assert_eq!(decoded.len(), 6);
        assert!((decoded[0] - 1.0).abs() < 1e-6);
        assert!((decoded[5] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_typed_accessor_decode_u32() {
        let indices: Vec<u32> = vec![0, 1, 2, 3, 4, 5];
        let mut buffer: Vec<u8> = Vec::new();
        for i in &indices {
            buffer.extend_from_slice(&i.to_le_bytes());
        }
        let acc = TypedAccessor::new(
            "idx",
            0,
            0,
            ComponentType::UnsignedInt,
            AccessorType::Scalar,
            6,
        );
        let decoded = acc.decode_u32(&buffer).unwrap();
        assert_eq!(decoded, indices);
    }
    #[test]
    fn test_typed_accessor_wrong_type_returns_none() {
        let buffer = vec![0u8; 24];
        let acc = TypedAccessor::new(
            "idx",
            0,
            0,
            ComponentType::UnsignedInt,
            AccessorType::Scalar,
            6,
        );
        assert!(acc.decode_f32(&buffer).is_none());
    }
    #[test]
    fn test_pbr_builder_default() {
        let mat = PbrMaterialBuilder::default();
        assert_eq!(mat.metallic_factor, 0.0);
        assert_eq!(mat.roughness_factor, 0.5);
        assert_eq!(mat.alpha_mode, "OPAQUE");
    }
    #[test]
    fn test_pbr_builder_chain() {
        let mat = PbrMaterialBuilder::new("metal")
            .base_color(0.8, 0.8, 0.8, 1.0)
            .metallic_roughness(0.9, 0.1)
            .emissive(0.0, 0.0, 0.0)
            .double_sided(true)
            .alpha_mode("BLEND");
        assert!((mat.metallic_factor - 0.9).abs() < 1e-6);
        assert!(mat.double_sided);
        assert_eq!(mat.alpha_mode, "BLEND");
    }
    #[test]
    fn test_pbr_to_json_contains_keys() {
        let mat = PbrMaterialBuilder::new("glass")
            .base_color(0.2, 0.4, 0.8, 0.5)
            .alpha_mode("BLEND");
        let json = mat.to_json();
        assert!(
            json.contains("pbrMetallicRoughness"),
            "should contain PBR key"
        );
        assert!(
            json.contains("baseColorFactor"),
            "should contain baseColorFactor"
        );
        assert!(json.contains("BLEND"), "should contain alpha mode");
        assert!(json.contains("glass"), "should contain name");
    }
    #[test]
    fn test_pbr_build_to_gltf_material() {
        let mat = PbrMaterialBuilder::new("copper")
            .metallic_roughness(0.8, 0.3)
            .build();
        assert_eq!(mat.name, "copper");
        assert!((mat.metallic_factor - 0.8).abs() < 1e-6);
        assert!((mat.roughness_factor - 0.3).abs() < 1e-6);
    }
    #[test]
    fn test_interpolation_as_str() {
        assert_eq!(Interpolation::Linear.as_str(), "LINEAR");
        assert_eq!(Interpolation::Step.as_str(), "STEP");
        assert_eq!(Interpolation::CubicSpline.as_str(), "CUBICSPLINE");
    }
    #[test]
    fn test_animation_channel_builder_push_and_duration() {
        let ch = AnimationChannelBuilder::new(0, "translation", Interpolation::Linear)
            .push(0.0, vec![0.0, 0.0, 0.0])
            .push(0.5, vec![0.0, 1.0, 0.0])
            .push(1.0, vec![0.0, 2.0, 0.0]);
        assert_eq!(ch.len(), 3);
        assert!((ch.duration() - 1.0).abs() < 1e-6);
        assert!(!ch.is_empty());
    }
    #[test]
    fn test_animation_channel_times_and_values() {
        let ch = AnimationChannelBuilder::new(1, "rotation", Interpolation::Step)
            .push(0.0, vec![0.0, 0.0, 0.0, 1.0])
            .push(1.0, vec![0.0, 0.707, 0.0, 0.707]);
        let times = ch.times();
        let vals = ch.values_flat();
        assert_eq!(times.len(), 2);
        assert_eq!(vals.len(), 8);
        assert!((times[0] - 0.0).abs() < 1e-6);
        assert!((times[1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_animation_channel_json_fragments() {
        let ch = AnimationChannelBuilder::new(0, "scale", Interpolation::Linear)
            .push(0.0, vec![1.0, 1.0, 1.0])
            .push(1.0, vec![2.0, 2.0, 2.0]);
        let (sampler, channel) = ch.to_json_fragments(0, 10, 11);
        assert!(
            sampler.contains("LINEAR"),
            "sampler should contain interpolation"
        );
        assert!(channel.contains("scale"), "channel should contain path");
        assert!(
            channel.contains("\"node\": 0"),
            "channel should contain node index"
        );
    }
    #[test]
    fn test_animation_channel_empty() {
        let ch = AnimationChannelBuilder::new(0, "translation", Interpolation::Linear);
        assert!(ch.is_empty());
        assert_eq!(ch.len(), 0);
        assert!((ch.duration() - 0.0).abs() < 1e-6);
    }
    #[test]
    fn test_base64_encode_hello() {
        let encoded = base64_encode(b"Hello");
        assert_eq!(encoded, "SGVsbG8=");
    }
    #[test]
    fn test_base64_encode_empty() {
        assert_eq!(base64_encode(b""), "");
    }
    #[test]
    fn test_base64_encode_length_multiple_3() {
        let encoded = base64_encode(b"Man");
        assert_eq!(encoded, "TWFu");
    }
    #[test]
    fn test_scene_to_gltf_json_empty_scene() {
        let scene = GltfScene::new();
        let json = scene_to_gltf_json(&scene, false);
        assert!(
            json.contains("\"version\": \"2.0\""),
            "should contain asset version"
        );
        assert!(json.contains("OxiPhysics"), "should contain generator");
    }
    #[test]
    fn test_scene_to_gltf_json_with_mesh() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "triangle".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                normals: vec![[0.0, 0.0, 1.0]; 3],
                texcoords: vec![],
                indices: vec![0, 1, 2],
            }],
        });
        let json = scene_to_gltf_json(&scene, false);
        assert!(json.contains("triangle"), "should contain mesh name");
        assert!(json.contains("VEC3"), "should contain VEC3 accessor type");
        assert!(json.contains("SCALAR"), "should contain SCALAR for indices");
    }
    #[test]
    fn test_scene_to_gltf_json_embed_buffer() {
        let mut scene = GltfScene::new();
        scene.add_mesh(GltfMesh {
            name: "cube".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0, 0.0, 0.0]],
                normals: vec![[0.0, 1.0, 0.0]],
                texcoords: vec![],
                indices: vec![0],
            }],
        });
        let json = scene_to_gltf_json(&scene, true);
        assert!(
            json.contains("data:application/octet-stream;base64,"),
            "embedded buffer should use data URI"
        );
    }

    // -----------------------------------------------------------------------
    // J1: glTF 2.0 Binary (GLB) buffer-packing tests
    // -----------------------------------------------------------------------

    /// Build a minimal scene with one triangle mesh for GLB tests.
    fn make_triangle_scene() -> GltfScene {
        let mut scene = GltfScene::new();
        let mesh_idx = scene.add_mesh(GltfMesh {
            name: "triangle".into(),
            primitives: vec![GltfPrimitive {
                positions: vec![[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                normals: vec![[0.0_f32, 0.0, 1.0]; 3],
                texcoords: vec![[0.0_f32, 0.0], [1.0, 0.0], [0.0, 1.0]],
                indices: vec![0u32, 1, 2],
            }],
        });
        scene.add_node(crate::gltf::GltfNode {
            name: "node0".into(),
            mesh: Some(mesh_idx),
            ..crate::gltf::GltfNode::default()
        });
        scene
    }

    #[test]
    fn test_write_glb_header_magic() {
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);
        assert!(bytes.len() >= 12, "GLB must have at least a 12-byte header");
        assert_eq!(&bytes[0..4], b"glTF", "magic bytes must be 'glTF'");
        let version = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);
        assert_eq!(version, 2, "GLB version must be 2");
    }

    #[test]
    fn test_write_glb_total_length_consistent() {
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);
        let total_len = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
        assert_eq!(
            total_len as usize,
            bytes.len(),
            "total_length field must match actual byte length"
        );
    }

    #[test]
    fn test_write_glb_chunk_types() {
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);
        // JSON chunk starts at byte 12
        let json_type = u32::from_le_bytes([bytes[16], bytes[17], bytes[18], bytes[19]]);
        assert_eq!(
            json_type, 0x4E4F534A,
            "first chunk must be JSON (0x4E4F534A)"
        );
        // JSON chunk length
        let json_chunk_len =
            u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]) as usize;
        // BIN chunk starts after header(12) + json_chunk_header(8) + json_data
        let bin_start = 12 + 8 + json_chunk_len;
        assert!(
            bin_start + 8 <= bytes.len(),
            "must have room for BIN chunk header"
        );
        let bin_type = u32::from_le_bytes([
            bytes[bin_start + 4],
            bytes[bin_start + 5],
            bytes[bin_start + 6],
            bytes[bin_start + 7],
        ]);
        assert_eq!(
            bin_type, 0x004E4942,
            "second chunk must be BIN\\0 (0x004E4942)"
        );
    }

    #[test]
    fn test_write_glb_vertex_positions_roundtrip() {
        let expected_positions: [[f32; 3]; 3] = [
            [0.0_f32, 0.0, 0.0],
            [1.0_f32, 0.0, 0.0],
            [0.0_f32, 1.0, 0.0],
        ];
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);

        // Extract BIN chunk offset and data.
        let json_chunk_len =
            u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]) as usize;
        let bin_start = 12 + 8 + json_chunk_len;
        let bin_chunk_len = u32::from_le_bytes([
            bytes[bin_start],
            bytes[bin_start + 1],
            bytes[bin_start + 2],
            bytes[bin_start + 3],
        ]) as usize;
        let bin_data = &bytes[bin_start + 8..bin_start + 8 + bin_chunk_len];

        // Position data starts at byteOffset=0 in the BIN chunk.
        let n_verts = expected_positions.len();
        let pos_size = n_verts * 12; // 3 × f32 × 4 bytes
        assert!(
            bin_data.len() >= pos_size,
            "BIN buffer must be large enough for position data"
        );

        for (i, expected) in expected_positions.iter().enumerate() {
            let off = i * 12;
            let x = f32::from_le_bytes([
                bin_data[off],
                bin_data[off + 1],
                bin_data[off + 2],
                bin_data[off + 3],
            ]);
            let y = f32::from_le_bytes([
                bin_data[off + 4],
                bin_data[off + 5],
                bin_data[off + 6],
                bin_data[off + 7],
            ]);
            let z = f32::from_le_bytes([
                bin_data[off + 8],
                bin_data[off + 9],
                bin_data[off + 10],
                bin_data[off + 11],
            ]);
            assert!(
                (x - expected[0]).abs() < f32::EPSILON,
                "position[{i}].x mismatch: {x} vs {}",
                expected[0]
            );
            assert!(
                (y - expected[1]).abs() < f32::EPSILON,
                "position[{i}].y mismatch: {y} vs {}",
                expected[1]
            );
            assert!(
                (z - expected[2]).abs() < f32::EPSILON,
                "position[{i}].z mismatch: {z} vs {}",
                expected[2]
            );
        }
    }

    #[test]
    fn test_write_glb_bufferview_offsets_layout() {
        // Verify that the JSON contains the expected bufferView byteOffset sequence.
        // For a triangle (3 verts, 3 indices):
        //   pos_size  = 3 * 12 = 36 bytes  → byteOffset=0
        //   norm_size = 3 * 12 = 36 bytes  → byteOffset=36
        //   uv_size   = 3 * 8  = 24 bytes  → byteOffset=72
        //   idx starts at 96 (96 is already 4-byte aligned)
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);

        // Extract JSON chunk string.
        let json_chunk_len =
            u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]) as usize;
        let json_bytes = &bytes[20..20 + json_chunk_len];
        // trim trailing spaces (padding)
        let json_str = std::str::from_utf8(json_bytes)
            .expect("valid UTF-8")
            .trim_end();

        // Check expected offsets appear in the JSON.
        assert!(
            json_str.contains("\"byteOffset\": 0"),
            "position bufferView byteOffset must be 0"
        );
        assert!(
            json_str.contains("\"byteOffset\": 36"),
            "normal bufferView byteOffset must be 36"
        );
        assert!(
            json_str.contains("\"byteOffset\": 72"),
            "uv bufferView byteOffset must be 72"
        );
        assert!(
            json_str.contains("\"byteOffset\": 96"),
            "index bufferView byteOffset must be 96"
        );
    }

    #[test]
    fn test_write_glb_json_contains_accessors() {
        let scene = make_triangle_scene();
        let writer = GlbWriter::new();
        let bytes = writer.write_glb(&scene);
        let json_chunk_len =
            u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]) as usize;
        let json_str = std::str::from_utf8(&bytes[20..20 + json_chunk_len])
            .expect("UTF-8")
            .trim_end();
        assert!(
            json_str.contains("\"accessors\""),
            "JSON must contain accessors array"
        );
        assert!(
            json_str.contains("\"bufferViews\""),
            "JSON must contain bufferViews array"
        );
        assert!(
            json_str.contains("\"buffers\""),
            "JSON must contain buffers array"
        );
        assert!(
            json_str.contains("VEC3"),
            "accessors must include VEC3 type"
        );
        assert!(
            json_str.contains("VEC2"),
            "accessors must include VEC2 type for UV"
        );
        assert!(
            json_str.contains("SCALAR"),
            "accessors must include SCALAR type for indices"
        );
        // POSITION accessor must have min/max bounds
        assert!(
            json_str.contains("\"min\""),
            "POSITION accessor must have min bounds"
        );
        assert!(
            json_str.contains("\"max\""),
            "POSITION accessor must have max bounds"
        );
    }
}
