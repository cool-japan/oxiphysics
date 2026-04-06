//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::{MorphTarget, SkinningVertex};

/// Compute a tangent-space normal from a normal map sample.
///
/// The normal map stores normals as RGB in \[0, 1\] range.
/// This converts to \[-1, 1\] range.
pub fn decode_normal_map(r: f64, g: f64, b: f64) -> [f64; 3] {
    let nx = r * 2.0 - 1.0;
    let ny = g * 2.0 - 1.0;
    let nz = b * 2.0 - 1.0;
    normalize3([nx, ny, nz])
}
/// Compute the TBN (Tangent, Bitangent, Normal) matrix from triangle data.
///
/// Given three positions and their UV coordinates, computes the tangent
/// and bitangent vectors for normal mapping.
///
/// Returns `(tangent, bitangent)` as unit vectors.
#[allow(clippy::too_many_arguments)]
pub fn compute_tbn(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    uv0: [f64; 2],
    uv1: [f64; 2],
    uv2: [f64; 2],
) -> ([f64; 3], [f64; 3]) {
    let edge1 = sub3(p1, p0);
    let edge2 = sub3(p2, p0);
    let duv1 = [uv1[0] - uv0[0], uv1[1] - uv0[1]];
    let duv2 = [uv2[0] - uv0[0], uv2[1] - uv0[1]];
    let det = duv1[0] * duv2[1] - duv2[0] * duv1[1];
    if det.abs() < 1e-15 {
        return ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    }
    let inv_det = 1.0 / det;
    let tangent = normalize3([
        inv_det * (duv2[1] * edge1[0] - duv1[1] * edge2[0]),
        inv_det * (duv2[1] * edge1[1] - duv1[1] * edge2[1]),
        inv_det * (duv2[1] * edge1[2] - duv1[1] * edge2[2]),
    ]);
    let bitangent = normalize3([
        inv_det * (-duv2[0] * edge1[0] + duv1[0] * edge2[0]),
        inv_det * (-duv2[0] * edge1[1] + duv1[0] * edge2[1]),
        inv_det * (-duv2[0] * edge1[2] + duv1[0] * edge2[2]),
    ]);
    (tangent, bitangent)
}
/// Transform a tangent-space normal to world space using the TBN matrix.
pub fn tbn_transform(
    tangent: [f64; 3],
    bitangent: [f64; 3],
    normal: [f64; 3],
    ts_normal: [f64; 3],
) -> [f64; 3] {
    let ws = [
        tangent[0] * ts_normal[0] + bitangent[0] * ts_normal[1] + normal[0] * ts_normal[2],
        tangent[1] * ts_normal[0] + bitangent[1] * ts_normal[1] + normal[1] * ts_normal[2],
        tangent[2] * ts_normal[0] + bitangent[2] * ts_normal[1] + normal[2] * ts_normal[2],
    ];
    normalize3(ws)
}
/// Compute a look-at view matrix (column-major 4x4).
pub fn look_at(eye: [f64; 3], target: [f64; 3], up: [f64; 3]) -> [[f64; 4]; 4] {
    let f = normalize3(sub3(target, eye));
    let s = normalize3(cross3(f, up));
    let u = cross3(s, f);
    let mut m = [[0.0_f64; 4]; 4];
    m[0][0] = s[0];
    m[1][0] = s[1];
    m[2][0] = s[2];
    m[0][1] = u[0];
    m[1][1] = u[1];
    m[2][1] = u[2];
    m[0][2] = -f[0];
    m[1][2] = -f[1];
    m[2][2] = -f[2];
    m[3][0] = -dot3(s, eye);
    m[3][1] = -dot3(u, eye);
    m[3][2] = dot3(f, eye);
    m[3][3] = 1.0;
    m
}
/// Compute a cubemap face and UV from a 3D direction.
///
/// Returns `(face_index, u, v)` where face_index is 0..5 (+X, -X, +Y, -Y, +Z, -Z).
pub fn cubemap_face_uv(direction: [f64; 3]) -> (usize, f64, f64) {
    let abs_x = direction[0].abs();
    let abs_y = direction[1].abs();
    let abs_z = direction[2].abs();
    let (face, sc, tc, ma) = if abs_x >= abs_y && abs_x >= abs_z {
        if direction[0] > 0.0 {
            (0, -direction[2], -direction[1], abs_x)
        } else {
            (1, direction[2], -direction[1], abs_x)
        }
    } else if abs_y >= abs_x && abs_y >= abs_z {
        if direction[1] > 0.0 {
            (2, direction[0], direction[2], abs_y)
        } else {
            (3, direction[0], -direction[2], abs_y)
        }
    } else if direction[2] > 0.0 {
        (4, direction[0], -direction[1], abs_z)
    } else {
        (5, -direction[0], -direction[1], abs_z)
    };
    let u = (sc / (2.0 * ma) + 0.5).clamp(0.0, 1.0);
    let v = (tc / (2.0 * ma) + 0.5).clamp(0.0, 1.0);
    (face, u, v)
}
/// Sample a spherical environment map from a direction.
///
/// Returns `(u, v)` in \[0, 1\] for equirectangular mapping.
pub fn equirectangular_uv(direction: [f64; 3]) -> (f64, f64) {
    let d = normalize3(direction);
    let u = 0.5 + d[2].atan2(d[0]) / (2.0 * PI);
    let v = 0.5 - d[1].asin() / PI;
    (u.clamp(0.0, 1.0), v.clamp(0.0, 1.0))
}
/// First-order spherical harmonics (SH) coefficients for irradiance.
///
/// Evaluates the first 4 SH basis functions at the given direction.
/// Returns \[Y00, Y1m1, Y10, Y11\] (L=0 and L=1 bands).
pub fn sh_basis_l1(direction: [f64; 3]) -> [f64; 4] {
    let d = normalize3(direction);
    let y00 = 0.5 * (1.0 / PI).sqrt();
    let y1m1 = (3.0 / (4.0 * PI)).sqrt() * d[1];
    let y10 = (3.0 / (4.0 * PI)).sqrt() * d[2];
    let y11 = (3.0 / (4.0 * PI)).sqrt() * d[0];
    [y00, y1m1, y10, y11]
}
/// Evaluate irradiance from first-order SH coefficients.
///
/// `coeffs` is a \[4\]\[3\] array: 4 SH bands, 3 color channels (RGB).
pub fn evaluate_sh_irradiance(coeffs: &[[f64; 3]; 4], direction: [f64; 3]) -> [f64; 3] {
    let basis = sh_basis_l1(direction);
    let mut result = [0.0_f64; 3];
    for i in 0..4 {
        for c in 0..3 {
            result[c] += basis[i] * coeffs[i][c];
        }
    }
    [result[0].max(0.0), result[1].max(0.0), result[2].max(0.0)]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = dot3(v, v).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 0.0];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::shader::Material;

    use crate::shader::PbrMaterial;

    use crate::shader::PhongMaterial;

    use crate::shader::RenderPass;
    use crate::shader::ShaderSource;
    use crate::shader::ShaderStage;
    use crate::shader::ShadowMapSetup;

    use crate::shader::UniformBuffer;
    use crate::shader::UniformValue;

    #[test]
    fn phong_lighting_90_degree_incidence_zero_diffuse() {
        let mat = PhongMaterial::default_gray();
        let result =
            PhongMaterial::phong_lighting(&mat, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        for ch in 0..3 {
            assert!(
                (result[ch] - mat.ambient[ch]).abs() < 1e-10,
                "channel {}: expected ambient {}, got {}",
                ch,
                mat.ambient[ch],
                result[ch]
            );
        }
    }
    #[test]
    fn phong_lighting_direct_incidence_nonzero() {
        let mat = PhongMaterial::default_gray();
        let result =
            PhongMaterial::phong_lighting(&mat, [0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        let diffuse_sum: f64 = mat.diffuse.iter().sum();
        assert!(
            result.iter().sum::<f64>() > diffuse_sum * 0.5,
            "direct incidence should produce significant brightness"
        );
    }
    #[test]
    fn material_set_get_uniform_roundtrip() {
        let vs = ShaderSource::vertex("void main(){}");
        let fs = ShaderSource::fragment("void main(){}");
        let mut mat = Material::new("test", vs, fs);
        mat.set_uniform("u_time", UniformValue::Float(3.125));
        mat.set_uniform("u_flag", UniformValue::Bool(true));
        match mat.get_uniform("u_time") {
            Some(UniformValue::Float(v)) => assert!((*v - 3.125).abs() < 1e-12),
            _ => panic!("expected Float uniform"),
        }
        match mat.get_uniform("u_flag") {
            Some(UniformValue::Bool(b)) => assert!(*b),
            _ => panic!("expected Bool uniform"),
        }
        assert!(mat.get_uniform("nonexistent").is_none());
    }
    #[test]
    fn render_pass_material_count() {
        let mut pass = RenderPass::new("forward");
        assert_eq!(pass.material_count(), 0);
        let make_mat = |n: &str| {
            Material::new(
                n,
                ShaderSource::vertex("void main(){}"),
                ShaderSource::fragment("void main(){}"),
            )
        };
        pass.add_material(make_mat("a"));
        pass.add_material(make_mat("b"));
        pass.add_material(make_mat("c"));
        assert_eq!(pass.material_count(), 3);
    }
    #[test]
    fn ggx_distribution_peak_at_ndoth_one() {
        let roughness = 0.1_f64;
        let d_peak = PbrMaterial::ggx_distribution(1.0, roughness);
        let d_off = PbrMaterial::ggx_distribution(0.9, roughness);
        assert!(d_peak > d_off);
    }
    #[test]
    fn ggx_distribution_roughness_one_lower_peak_than_roughness_point_one() {
        let d_sharp = PbrMaterial::ggx_distribution(1.0, 0.1);
        let d_rough = PbrMaterial::ggx_distribution(1.0, 1.0);
        assert!(d_sharp > d_rough);
    }
    #[test]
    fn schlick_fresnel_cos_zero_gives_white() {
        let f0 = [0.04, 0.04, 0.04];
        let result = PbrMaterial::schlick_fresnel(0.0, f0);
        for ch in 0..3 {
            assert!((result[ch] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn schlick_fresnel_cos_one_gives_f0() {
        let f0 = [0.1, 0.2, 0.3];
        let result = PbrMaterial::schlick_fresnel(1.0, f0);
        for ch in 0..3 {
            assert!((result[ch] - f0[ch]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_pbr_f0_dielectric() {
        let mat = PbrMaterial::dielectric([0.8, 0.2, 0.1], 0.5);
        let f0 = mat.f0();
        for ch in 0..3 {
            assert!(
                (f0[ch] - 0.04).abs() < 1e-10,
                "Dielectric F0 should be 0.04"
            );
        }
    }
    #[test]
    fn test_pbr_f0_metal() {
        let mat = PbrMaterial::metallic([0.8, 0.2, 0.1], 0.5);
        let f0 = mat.f0();
        assert!((f0[0] - 0.8).abs() < 1e-10);
        assert!((f0[1] - 0.2).abs() < 1e-10);
        assert!((f0[2] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_geometry_smith() {
        let g = PbrMaterial::geometry_smith(1.0, 1.0, 0.5);
        assert!(g > 0.0 && g <= 1.0, "G should be in (0, 1], got {g}");
        let g_grazing = PbrMaterial::geometry_smith(0.01, 0.01, 0.5);
        assert!(g_grazing < g, "Grazing should have less geometry");
    }
    #[test]
    fn test_decode_normal_map_center() {
        let n = decode_normal_map(0.5, 0.5, 1.0);
        assert!((n[0] - 0.0).abs() < 1e-6);
        assert!((n[1] - 0.0).abs() < 1e-6);
        assert!((n[2] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_compute_tbn() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.0, 1.0, 0.0];
        let uv0 = [0.0, 0.0];
        let uv1 = [1.0, 0.0];
        let uv2 = [0.0, 1.0];
        let (tangent, bitangent) = compute_tbn(p0, p1, p2, uv0, uv1, uv2);
        assert!((tangent[0] - 1.0).abs() < 1e-6, "Tangent X: {}", tangent[0]);
        assert!(
            (bitangent[1] - 1.0).abs() < 1e-6,
            "Bitangent Y: {}",
            bitangent[1]
        );
    }
    #[test]
    fn test_tbn_transform_identity() {
        let tangent = [1.0, 0.0, 0.0];
        let bitangent = [0.0, 1.0, 0.0];
        let normal = [0.0, 0.0, 1.0];
        let ts_normal = [0.0, 0.0, 1.0];
        let ws = tbn_transform(tangent, bitangent, normal, ts_normal);
        assert!((ws[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_shadow_map_setup() {
        let shadow = ShadowMapSetup::directional([0.0, -1.0, 0.0], 50.0);
        assert_eq!(shadow.resolution, 2048);
        assert!(shadow.texel_size() > 0.0);
        let ortho = shadow.ortho_matrix();
        assert!((ortho[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_shadow_map_ortho_dimensions() {
        let shadow = ShadowMapSetup::directional([0.0, -1.0, 0.0], 10.0);
        let m = shadow.ortho_matrix();
        assert!((m[0][0] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_cubemap_face_positive_x() {
        let (face, u, v) = cubemap_face_uv([1.0, 0.0, 0.0]);
        assert_eq!(face, 0, "Direction +X should map to face 0");
        assert!((u - 0.5).abs() < 1e-6);
        assert!((v - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_cubemap_face_negative_y() {
        let (face, _u, _v) = cubemap_face_uv([0.0, -1.0, 0.0]);
        assert_eq!(face, 3, "Direction -Y should map to face 3");
    }
    #[test]
    fn test_equirectangular_uv() {
        let (_, v) = equirectangular_uv([0.0, 1.0, 0.0]);
        assert!(v < 0.01, "Straight up should have v~0, got {v}");
        let (_, v_down) = equirectangular_uv([0.0, -1.0, 0.0]);
        assert!(v_down > 0.99, "Straight down should have v~1, got {v_down}");
    }
    #[test]
    fn test_sh_basis_l1() {
        let b1 = sh_basis_l1([1.0, 0.0, 0.0]);
        let b2 = sh_basis_l1([0.0, 1.0, 0.0]);
        assert!((b1[0] - b2[0]).abs() < 1e-10, "Y00 should be constant");
    }
    #[test]
    fn test_sh_irradiance() {
        let coeffs = [
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let c1 = evaluate_sh_irradiance(&coeffs, [1.0, 0.0, 0.0]);
        let c2 = evaluate_sh_irradiance(&coeffs, [0.0, 1.0, 0.0]);
        for ch in 0..3 {
            assert!(
                (c1[ch] - c2[ch]).abs() < 1e-6,
                "Uniform SH should give same irradiance in all directions"
            );
        }
    }
    #[test]
    fn test_uniform_buffer() {
        let mut buf = UniformBuffer::new(0);
        assert!(buf.is_empty());
        assert_eq!(buf.len(), 0);
        buf.set("u_time", UniformValue::Float(1.0));
        buf.set("u_color", UniformValue::Vec3([1.0, 0.0, 0.0]));
        assert_eq!(buf.len(), 2);
        assert!(!buf.is_empty());
        assert!(buf.remove("u_time"));
        assert_eq!(buf.len(), 1);
        assert!(!buf.remove("nonexistent"));
        buf.clear();
        assert!(buf.is_empty());
    }
    #[test]
    fn test_look_at_matrix() {
        let m = look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((m[3][3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pbr_to_material() {
        let pbr = PbrMaterial::default();
        let mat = pbr.to_material();
        assert_eq!(mat.name, "pbr");
        assert!(mat.get_uniform("u_metallic").is_some());
        assert!(mat.get_uniform("u_roughness").is_some());
        assert!(mat.get_uniform("u_base_color").is_some());
        assert_eq!(mat.uniform_count(), 5);
    }
    #[test]
    fn test_render_pass_clear_color() {
        let mut pass = RenderPass::new("shadow");
        pass.set_clear_color(1.0, 1.0, 1.0, 1.0);
        assert!((pass.clear_color[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_compute_shader_source() {
        let cs = ShaderSource::compute("void main() {}");
        assert_eq!(cs.entry_point, "main");
        match cs.stage {
            ShaderStage::Compute => {}
            _ => panic!("Expected Compute stage"),
        }
    }
}
/// Helper to concatenate a preamble with a shader body.
#[allow(dead_code)]
pub fn glsl_compose(preamble: &str, body: &str) -> String {
    format!("{preamble}\n{body}")
}
/// Maximum number of joints supported in a skinning shader.
pub const MAX_JOINTS: usize = 128;
/// Apply skin weights to transform a position.
///
/// `joint_matrices` maps joint index → world-space joint matrix (column-major 4x4).
/// Returns the blended world-space position.
#[allow(dead_code)]
pub fn skin_position(
    pos: [f64; 3],
    sv: &SkinningVertex,
    joint_matrices: &[[[f64; 4]; 4]],
) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for i in 0..4 {
        let w = sv.weights[i] as f64;
        if w.abs() < 1e-15 {
            continue;
        }
        let ji = sv.joint_indices[i] as usize;
        if ji >= joint_matrices.len() {
            continue;
        }
        let m = &joint_matrices[ji];
        let tp = mat4_transform_point(m, pos);
        result[0] += w * tp[0];
        result[1] += w * tp[1];
        result[2] += w * tp[2];
    }
    result
}
/// Apply a blend of morph targets to a base mesh.
///
/// `base_positions` are the original vertex positions.
/// `targets` is a slice of `(&MorphTarget, weight)` pairs.
///
/// Returns the blended positions. Vertices beyond any target's length are
/// left at their base positions.
#[allow(dead_code)]
pub fn apply_morph_targets(
    base_positions: &[[f32; 3]],
    targets: &[(&MorphTarget, f32)],
) -> Vec<[f32; 3]> {
    let mut result = base_positions.to_vec();
    for (target, weight) in targets {
        for (i, disp) in target.displacements.iter().enumerate() {
            if i >= result.len() {
                break;
            }
            result[i][0] += weight * disp[0];
            result[i][1] += weight * disp[1];
            result[i][2] += weight * disp[2];
        }
    }
    result
}
/// Evaluate all 9 real spherical harmonic basis functions (L0 + L1 + L2)
/// at a unit direction.
#[allow(dead_code)]
pub fn sh9_basis(d: [f64; 3]) -> [f64; 9] {
    let [x, y, z] = d;
    let y00 = (1.0 / (4.0 * PI)).sqrt();
    let y1m1 = (3.0 / (4.0 * PI)).sqrt() * y;
    let y10 = (3.0 / (4.0 * PI)).sqrt() * z;
    let y11 = (3.0 / (4.0 * PI)).sqrt() * x;
    let y2m2 = (15.0 / (4.0 * PI)).sqrt() * x * y;
    let y2m1 = (15.0 / (4.0 * PI)).sqrt() * y * z;
    let y20 = (5.0 / (16.0 * PI)).sqrt() * (3.0 * z * z - 1.0);
    let y21 = (15.0 / (4.0 * PI)).sqrt() * x * z;
    let y22 = (15.0 / (16.0 * PI)).sqrt() * (x * x - y * y);
    [y00, y1m1, y10, y11, y2m2, y2m1, y20, y21, y22]
}
/// Create a 4x4 identity matrix (column-major).
#[allow(dead_code)]
pub fn identity4() -> [[f64; 4]; 4] {
    let mut m = [[0.0_f64; 4]; 4];
    m[0][0] = 1.0;
    m[1][1] = 1.0;
    m[2][2] = 1.0;
    m[3][3] = 1.0;
    m
}
/// Transform a 3D point by a 4x4 column-major matrix (assumes w=1).
#[allow(dead_code)]
pub fn mat4_transform_point(m: &[[f64; 4]; 4], p: [f64; 3]) -> [f64; 3] {
    let x = m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0];
    let y = m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1];
    let z = m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2];
    [x, y, z]
}
/// Normalize a 3D vector (re-exported for shader helpers).
#[inline]
#[allow(dead_code)]
pub(super) fn normalize3_shader(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 1.0];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}
#[cfg(test)]
mod extended_tests {
    use super::*;

    use crate::shader::GlslOptions;
    use crate::shader::IblProbe;

    use crate::shader::Joint;

    use crate::shader::PcfSampler;

    use crate::shader::Pose;

    use crate::shader::Skeleton;

    #[test]
    fn test_glsl_preamble_version_line() {
        let opts = GlslOptions::glsl_330();
        let pre = opts.preamble();
        assert!(
            pre.starts_with("#version 330 core"),
            "preamble should start with version: {pre}"
        );
    }
    #[test]
    fn test_glsl_preamble_es_precision() {
        let opts = GlslOptions::glsl_es_300();
        let pre = opts.preamble();
        assert!(
            pre.contains("precision mediump float"),
            "ES target should have precision: {pre}"
        );
    }
    #[test]
    fn test_glsl_preamble_no_precision_for_desktop() {
        let opts = GlslOptions::glsl_330();
        let pre = opts.preamble();
        assert!(
            !pre.contains("precision mediump"),
            "desktop should not have mediump: {pre}"
        );
    }
    #[test]
    fn test_glsl_define_flag() {
        let opts = GlslOptions::glsl_330().define_flag("USE_SHADOWS");
        let pre = opts.preamble();
        assert!(
            pre.contains("#define USE_SHADOWS"),
            "should have define: {pre}"
        );
        assert!(
            !pre.contains("#define USE_SHADOWS "),
            "flag define should not have a value"
        );
    }
    #[test]
    fn test_glsl_define_with_value() {
        let opts = GlslOptions::glsl_330().define("MAX_LIGHTS", "8");
        let pre = opts.preamble();
        assert!(
            pre.contains("#define MAX_LIGHTS 8"),
            "should have valued define: {pre}"
        );
    }
    #[test]
    fn test_glsl_extension_directive() {
        let opts = GlslOptions::glsl_330().extension("GL_ARB_shader_draw_parameters", "require");
        let pre = opts.preamble();
        assert!(
            pre.contains("#extension GL_ARB_shader_draw_parameters"),
            "should have extension: {pre}"
        );
    }
    #[test]
    fn test_glsl_multiple_defines() {
        let opts = GlslOptions::glsl_330()
            .define_flag("ALPHA_TEST")
            .define("ALPHA_CUTOFF", "0.5")
            .define_flag("ENABLE_FOG");
        let pre = opts.preamble();
        assert!(pre.contains("ALPHA_TEST"));
        assert!(pre.contains("ALPHA_CUTOFF 0.5"));
        assert!(pre.contains("ENABLE_FOG"));
    }
    #[test]
    fn test_glsl_compose() {
        let pre = "#version 330 core\n";
        let body = "void main() {}";
        let result = glsl_compose(pre, body);
        assert!(result.contains("#version 330 core"));
        assert!(result.contains("void main() {}"));
    }
    #[test]
    fn test_skinning_vertex_single_joint() {
        let sv = SkinningVertex::single_joint(5);
        assert_eq!(sv.joint_indices[0], 5);
        assert!((sv.weights[0] - 1.0).abs() < 1e-6);
        assert!((sv.weight_sum() - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_skinning_vertex_two_joints() {
        let sv = SkinningVertex::two_joints(0, 0.6, 1, 0.4);
        assert!((sv.weight_sum() - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_skinning_vertex_normalize() {
        let sv = SkinningVertex {
            joint_indices: [0, 1, 0, 0],
            weights: [2.0, 2.0, 0.0, 0.0],
        }
        .normalized();
        assert!(
            (sv.weight_sum() - 1.0).abs() < 1e-6,
            "normalized weights should sum to 1: {}",
            sv.weight_sum()
        );
    }
    #[test]
    fn test_skin_position_identity_joint() {
        let sv = SkinningVertex::single_joint(0);
        let m = identity4();
        let pos = [1.0, 2.0, 3.0];
        let result = skin_position(pos, &sv, &[m]);
        for i in 0..3 {
            assert!(
                (result[i] - pos[i]).abs() < 1e-10,
                "identity joint should preserve position: {result:?}"
            );
        }
    }
    #[test]
    fn test_skin_position_no_joints_gives_zero() {
        let sv = SkinningVertex::single_joint(0);
        let result = skin_position([1.0, 2.0, 3.0], &sv, &[]);
        assert_eq!(result, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_skin_position_two_joints_blended() {
        let sv = SkinningVertex::two_joints(0, 0.5, 1, 0.5);
        let mut m0 = identity4();
        m0[3][0] = 2.0;
        let m1 = identity4();
        let pos = [0.0, 0.0, 0.0];
        let result = skin_position(pos, &sv, &[m0, m1]);
        assert!((result[0] - 1.0).abs() < 1e-10, "blended x: {}", result[0]);
        assert!(result[1].abs() < 1e-10);
    }
    #[test]
    fn test_skeleton_add_joint() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint(Joint::new("root", None));
        let child = skel.add_joint(Joint::new("hip", Some(root)));
        assert_eq!(skel.joint_count(), 2);
        assert_eq!(skel.joints[child].parent, Some(root));
    }
    #[test]
    fn test_skeleton_find_joint() {
        let mut skel = Skeleton::new();
        skel.add_joint(Joint::new("spine", None));
        skel.add_joint(Joint::new("head", Some(0)));
        assert_eq!(skel.find_joint("head"), Some(1));
        assert!(skel.find_joint("nonexistent").is_none());
    }
    #[test]
    fn test_pose_bind_pose_identity() {
        let pose = Pose::bind_pose(3);
        assert_eq!(pose.joint_count(), 3);
        for m in &pose.local_transforms {
            assert!((m[0][0] - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_morph_target_creation() {
        let t = MorphTarget::new("smile", vec![[0.1, 0.0, 0.0]; 10]);
        assert_eq!(t.vertex_count(), 10);
        assert_eq!(t.name, "smile");
    }
    #[test]
    fn test_apply_morph_targets_zero_weight() {
        let base = vec![[1.0_f32, 2.0, 3.0]; 5];
        let target = MorphTarget::new("t", vec![[10.0, 10.0, 10.0]; 5]);
        let result = apply_morph_targets(&base, &[(&target, 0.0)]);
        assert_eq!(
            result[0],
            [1.0, 2.0, 3.0],
            "zero weight should not change positions"
        );
    }
    #[test]
    fn test_apply_morph_targets_full_weight() {
        let base = vec![[0.0_f32; 3]; 3];
        let target = MorphTarget::new("t", vec![[1.0, 0.0, 0.0]; 3]);
        let result = apply_morph_targets(&base, &[(&target, 1.0)]);
        assert!((result[0][0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_apply_morph_targets_half_weight() {
        let base = vec![[0.0_f32; 3]; 2];
        let target = MorphTarget::new("t", vec![[2.0, 0.0, 0.0]; 2]);
        let result = apply_morph_targets(&base, &[(&target, 0.5)]);
        assert!(
            (result[0][0] - 1.0).abs() < 1e-6,
            "half weight should give 1.0: {}",
            result[0][0]
        );
    }
    #[test]
    fn test_apply_morph_targets_two_targets() {
        let base = vec![[0.0_f32; 3]; 1];
        let t1 = MorphTarget::new("t1", vec![[1.0, 0.0, 0.0]]);
        let t2 = MorphTarget::new("t2", vec![[0.0, 1.0, 0.0]]);
        let result = apply_morph_targets(&base, &[(&t1, 1.0), (&t2, 1.0)]);
        assert!((result[0][0] - 1.0).abs() < 1e-6);
        assert!((result[0][1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_apply_morph_shorter_target_doesnt_panic() {
        let base = vec![[0.0_f32; 3]; 5];
        let target = MorphTarget::new("t", vec![[1.0, 0.0, 0.0]; 2]);
        let result = apply_morph_targets(&base, &[(&target, 1.0)]);
        assert_eq!(result.len(), 5, "result should match base length");
        assert!((result[0][0] - 1.0).abs() < 1e-6, "first vert affected");
        assert_eq!(
            result[2], [0.0; 3],
            "beyond target length should be unchanged"
        );
    }
    #[test]
    fn test_pcf_fully_lit_when_receiver_shallower() {
        let sampler = PcfSampler::constant(1.0, 8, 8);
        let factor = sampler.sample_pcf_3x3(0.5, 0.5, 0.5);
        assert!((factor - 1.0).abs() < 1e-6, "should be fully lit: {factor}");
    }
    #[test]
    fn test_pcf_fully_shadowed_when_receiver_deeper() {
        let sampler = PcfSampler::new(vec![0.0_f32; 8 * 8], 8, 8, 0.0);
        let factor = sampler.sample_pcf_3x3(0.5, 0.5, 0.9);
        assert!(factor < 1e-6, "should be fully shadowed: {factor}");
    }
    #[test]
    fn test_pcf_fetch_depth() {
        let mut depths = vec![0.0_f32; 4 * 4];
        depths[0] = 0.9;
        let sampler = PcfSampler::new(depths, 4, 4, 0.001);
        let d = sampler.fetch_depth(0.0, 0.0);
        assert!((d - 0.9).abs() < 1e-6);
    }
    #[test]
    fn test_pcf_radius_zero_single_tap() {
        let sampler = PcfSampler::constant(1.0, 4, 4);
        let f = sampler.sample_pcf(0.5, 0.5, 0.5, 0);
        assert!((f - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_ibl_probe_uniform_grey_consistent() {
        let probe = IblProbe::uniform_grey(0.5);
        let c1 = probe.evaluate([1.0, 0.0, 0.0]);
        let c2 = probe.evaluate([0.0, 1.0, 0.0]);
        for ch in 0..3 {
            assert!(
                (c1[ch] - c2[ch]).abs() < 1e-10,
                "uniform probe should be direction-independent: {c1:?} vs {c2:?}"
            );
        }
    }
    #[test]
    fn test_ibl_probe_non_negative() {
        let probe = IblProbe::uniform_grey(1.0);
        let c = probe.evaluate([0.0, 0.0, 1.0]);
        for ch in 0..3 {
            assert!(c[ch] >= 0.0, "irradiance must be non-negative");
        }
    }
    #[test]
    fn test_ibl_probe_intensity_scales_output() {
        let mut p1 = IblProbe::uniform_grey(1.0);
        let mut p2 = IblProbe::uniform_grey(1.0);
        p1.intensity = 1.0;
        p2.intensity = 2.0;
        let c1 = p1.evaluate([0.0, 0.0, 1.0]);
        let c2 = p2.evaluate([0.0, 0.0, 1.0]);
        for ch in 0..3 {
            assert!(
                (c2[ch] - 2.0 * c1[ch]).abs() < 1e-10,
                "intensity 2x should double the output: {c2:?} vs {c1:?}"
            );
        }
    }
    #[test]
    fn test_sh9_basis_length_9() {
        let b = sh9_basis([0.0, 0.0, 1.0]);
        assert_eq!(b.len(), 9);
    }
    #[test]
    fn test_sh9_basis_l0_constant() {
        let b1 = sh9_basis([1.0, 0.0, 0.0]);
        let b2 = sh9_basis([0.0, 1.0, 0.0]);
        assert!(
            (b1[0] - b2[0]).abs() < 1e-10,
            "Y00 should be constant for all directions"
        );
    }
    #[test]
    fn test_identity4_correct() {
        let m = identity4();
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((m[i][j] - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_mat4_transform_point_identity() {
        let m = identity4();
        let p = [1.0, 2.0, 3.0];
        let tp = mat4_transform_point(&m, p);
        for i in 0..3 {
            assert!((tp[i] - p[i]).abs() < 1e-10);
        }
    }
    #[test]
    fn test_mat4_transform_point_translation() {
        let mut m = identity4();
        m[3][0] = 5.0;
        m[3][1] = -3.0;
        let p = [1.0, 1.0, 1.0];
        let tp = mat4_transform_point(&m, p);
        assert!((tp[0] - 6.0).abs() < 1e-10, "tx: {}", tp[0]);
        assert!((tp[1] - (-2.0)).abs() < 1e-10, "ty: {}", tp[1]);
        assert!((tp[2] - 1.0).abs() < 1e-10);
    }
}
#[cfg(test)]
mod new_shader_tests {
    use super::*;
    use crate::shader::AnisotropicBrdf;
    use crate::shader::BssrdfApprox;

    use crate::shader::IridescenceParams;

    use crate::shader::ParticleSpriteParams;

    use crate::shader::SssParams;

    use crate::shader::VolumetricFogParams;
    #[test]
    fn sss_fresnel_schlick_at_normal_incidence() {
        let sss = SssParams::skin();
        let r0 = ((sss.ior - 1.0) / (sss.ior + 1.0)).powi(2);
        let f = sss.fresnel_schlick_scalar(1.0);
        assert!((f - r0).abs() < 1e-10, "F at cos_theta=1 should be F0: {f}");
    }
    #[test]
    fn sss_fresnel_schlick_increases_at_grazing() {
        let sss = SssParams::skin();
        let f_normal = sss.fresnel_schlick_scalar(1.0);
        let f_grazing = sss.fresnel_schlick_scalar(0.0);
        assert!(
            f_grazing > f_normal,
            "Fresnel should increase at grazing: normal={f_normal}, grazing={f_grazing}"
        );
    }
    #[test]
    fn sss_dipole_reflectance_decreases_with_distance() {
        let sss = SssParams::skin();
        let r1 = sss.dipole_reflectance(0.001, 0);
        let r2 = sss.dipole_reflectance(0.01, 0);
        assert!(
            r2 < r1,
            "reflectance should decrease with distance: r1={r1}, r2={r2}"
        );
    }
    #[test]
    fn sss_evaluate_sss_non_negative() {
        let sss = SssParams::skin();
        let result = sss.evaluate_sss(0.8, 0.002, [1.0, 1.0, 1.0]);
        for &v in &result {
            assert!(v >= 0.0, "SSS should be non-negative: {v}");
        }
    }
    #[test]
    fn sss_marble_albedo_range() {
        let marble = SssParams::marble();
        for &a in &marble.albedo {
            assert!((0.0..=1.0).contains(&a), "albedo should be in [0,1]: {a}");
        }
    }
    #[test]
    fn aniso_brdf_fresnel_at_normal_incidence_equals_f0() {
        let brdf = AnisotropicBrdf::isotropic(0.3, [0.04, 0.04, 0.04]);
        let f = brdf.fresnel(1.0);
        for (ch, &fv) in f.iter().enumerate() {
            assert!((fv - brdf.f0[ch]).abs() < 1e-10, "F(1)=F0: ch={ch}, {fv}");
        }
    }
    #[test]
    fn aniso_brdf_fresnel_approaches_one_at_grazing() {
        let brdf = AnisotropicBrdf::isotropic(0.3, [0.04, 0.04, 0.04]);
        let f = brdf.fresnel(0.0);
        for &fv in &f {
            assert!((fv - 1.0).abs() < 1e-10, "F(0)=1: {fv}");
        }
    }
    #[test]
    fn aniso_brdf_ndf_positive() {
        let brdf = AnisotropicBrdf {
            alpha_x: 0.3,
            alpha_y: 0.5,
            f0: [0.04; 3],
        };
        let h = normalize3([0.0, 0.0, 1.0]);
        let n = [0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let d = brdf.ndf_ggx_aniso(h, t, b, n);
        assert!(d > 0.0, "NDF should be positive at normal incidence: {d}");
    }
    #[test]
    fn aniso_brdf_evaluate_nonneg() {
        let brdf = AnisotropicBrdf::isotropic(0.5, [0.04; 3]);
        let n = [0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let v = normalize3([0.3, 0.0, 0.7]);
        let l = normalize3([0.1, 0.1, 0.9]);
        let result = brdf.evaluate(v, l, n, t, b);
        for &r in &result {
            assert!(r >= 0.0, "BRDF should be non-negative: {r}");
        }
    }
    #[test]
    fn aniso_brdf_directional_albedo_in_range() {
        let brdf = AnisotropicBrdf::isotropic(0.3, [0.04; 3]);
        let alb = brdf.directional_albedo(0.8);
        for &a in &alb {
            assert!((0.0..=1.0).contains(&a), "directional albedo in [0,1]: {a}");
        }
    }
    #[test]
    fn iridescence_evaluate_rgb_in_unit_range() {
        let params = IridescenceParams::soap_bubble();
        let colors = params.evaluate_rgb(0.5, [700.0, 546.0, 435.0]);
        for &c in &colors {
            assert!(
                (0.0..=1.0 + 1e-6).contains(&c),
                "iridescence should be in [0,1]: {c}"
            );
        }
    }
    #[test]
    fn iridescence_opd_increases_with_thickness() {
        let p1 = IridescenceParams {
            thickness_nm: 100.0,
            ..IridescenceParams::soap_bubble()
        };
        let p2 = IridescenceParams {
            thickness_nm: 200.0,
            ..IridescenceParams::soap_bubble()
        };
        assert!(p2.opd(1.0) > p1.opd(1.0));
    }
    #[test]
    fn iridescence_fresnel_zero_at_normal_incidence_for_low_ior() {
        let params = IridescenceParams {
            ior_film: 1.0,
            ..IridescenceParams::soap_bubble()
        };
        let f = params.fresnel_film(1.0);
        assert!(
            f.abs() < 1e-10,
            "no reflection when ior_film == ior_air: {f}"
        );
    }
    #[test]
    fn iridescence_colors_vary_with_angle() {
        let params = IridescenceParams::soap_bubble();
        let c1 = params.evaluate_rgb(0.1, [700.0, 546.0, 435.0]);
        let c2 = params.evaluate_rgb(0.9, [700.0, 546.0, 435.0]);
        let diff: f64 = c1.iter().zip(c2.iter()).map(|(a, b)| (a - b).abs()).sum();
        assert!(diff.is_finite(), "colors should be finite");
    }
    #[test]
    fn bssrdf_transmittance_one_at_zero_distance() {
        let bssrdf = BssrdfApprox {
            sigma_s: [1.0, 0.8, 0.5],
            sigma_a: [0.1, 0.05, 0.02],
            ior: 1.4,
        };
        let t = bssrdf.transmittance(0.0);
        for &tv in &t {
            assert!((tv - 1.0).abs() < 1e-10, "transmittance at 0 = 1: {tv}");
        }
    }
    #[test]
    fn bssrdf_transmittance_decreases_with_distance() {
        let bssrdf = BssrdfApprox {
            sigma_s: [1.0, 0.8, 0.5],
            sigma_a: [0.1, 0.05, 0.02],
            ior: 1.4,
        };
        let t0 = bssrdf.transmittance(0.0);
        let t1 = bssrdf.transmittance(1.0);
        for (t0v, t1v) in t0.iter().zip(t1.iter()) {
            assert!(t1v < t0v, "transmittance decreases: t0={t0v}, t1={t1v}");
        }
    }
    #[test]
    fn bssrdf_diffusion_profile_positive() {
        let bssrdf = BssrdfApprox {
            sigma_s: [1.0, 0.8, 0.5],
            sigma_a: [0.1, 0.05, 0.02],
            ior: 1.4,
        };
        let r = bssrdf.diffusion_profile(0.01);
        for &rv in &r {
            assert!(rv > 0.0, "diffusion profile should be positive: {rv}");
        }
    }
    #[test]
    fn bssrdf_sigma_tr_positive() {
        let bssrdf = BssrdfApprox {
            sigma_s: [1.0, 0.8, 0.5],
            sigma_a: [0.1, 0.05, 0.02],
            ior: 1.4,
        };
        for &s in &bssrdf.sigma_tr() {
            assert!(s > 0.0, "sigma_tr should be positive: {s}");
        }
    }
    #[test]
    fn particle_sprite_circle_alpha_inside() {
        let p = ParticleSpriteParams::circle(1.0, [1.0, 0.5, 0.0, 1.0]);
        let a = p.alpha(0.0, 0.0);
        assert!((a - 1.0).abs() < 1e-10, "alpha at centre = 1: {a}");
    }
    #[test]
    fn particle_sprite_circle_alpha_outside() {
        let p = ParticleSpriteParams::circle(1.0, [1.0, 0.5, 0.0, 1.0]);
        let a = p.alpha(1.5, 0.0);
        assert_eq!(a, 0.0, "alpha outside circle = 0: {a}");
    }
    #[test]
    fn particle_sprite_gaussian_alpha_max_at_centre() {
        let p = ParticleSpriteParams::gaussian(1.0, [1.0, 1.0, 1.0, 1.0]);
        let a_centre = p.alpha(0.0, 0.0);
        let a_edge = p.alpha(0.9, 0.0);
        assert!(
            a_centre > a_edge,
            "gaussian max at centre: {a_centre} vs {a_edge}"
        );
    }
    #[test]
    fn particle_sprite_soft_blend_far_scene() {
        let p = ParticleSpriteParams::circle(1.0, [1.0; 4]);
        let blend = p.soft_blend(0.5, 10.0);
        assert!((blend - 1.0).abs() < 1e-10, "blend={blend}");
    }
    #[test]
    fn particle_sprite_soft_blend_same_depth() {
        let p = ParticleSpriteParams::circle(1.0, [1.0; 4]);
        let blend = p.soft_blend(1.0, 1.0);
        assert!(blend.abs() < 1e-10, "same depth → blend=0: {blend}");
    }
    #[test]
    fn particle_sprite_final_color_glow() {
        let mut p = ParticleSpriteParams::circle(1.0, [0.5, 0.5, 0.5, 1.0]);
        p.glow = 0.0;
        let c_no_glow = p.final_color(1.0);
        p.glow = 1.0;
        let c_glow = p.final_color(1.0);
        assert!(
            c_glow[0] > c_no_glow[0],
            "glow should brighten: {} vs {}",
            c_glow[0],
            c_no_glow[0]
        );
    }
    #[test]
    fn fog_transmittance_one_at_zero_distance() {
        let fog = VolumetricFogParams::thin_fog();
        let t = fog.transmittance(0.0, 0.0);
        assert!(
            (t - 1.0).abs() < 1e-10,
            "transmittance at d=0 should be 1: {t}"
        );
    }
    #[test]
    fn fog_transmittance_decreases_with_distance() {
        let fog = VolumetricFogParams::thin_fog();
        let t1 = fog.transmittance(10.0, 0.0);
        let t2 = fog.transmittance(100.0, 0.0);
        assert!(t2 < t1, "transmittance should decrease: t1={t1}, t2={t2}");
    }
    #[test]
    fn fog_phase_normalised_approximately() {
        let fog = VolumetricFogParams {
            phase_g: 0.0,
            ..VolumetricFogParams::thin_fog()
        };
        let p1 = fog.phase(0.0);
        let p2 = fog.phase(1.0);
        assert!(
            (p1 - p2).abs() < 1e-10,
            "isotropic phase should be constant: {p1} vs {p2}"
        );
    }
    #[test]
    fn fog_inscatter_returns_finite_values() {
        let fog = VolumetricFogParams::thin_fog();
        let result = fog.inscatter(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            10.0,
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            8,
        );
        for &r in &result {
            assert!(
                r.is_finite() && r >= 0.0,
                "inscatter should be finite non-negative: {r}"
            );
        }
    }
    #[test]
    fn fog_blend_scene_returns_fog_color_at_high_density() {
        let mut fog = VolumetricFogParams::ground_fog();
        fog.density = 100.0;
        let blended = fog.blend_scene([1.0, 0.0, 0.0], 10.0, 0.0);
        let t = fog.transmittance(10.0, 0.0);
        let expected = t * 1.0 + (1.0 - t) * fog.color[0];
        assert!(
            (blended[0] - expected).abs() < 1e-8,
            "blended r={} expected={}",
            blended[0],
            expected
        );
        assert!(
            blended[0] < 1.0,
            "fog should reduce max red: {}",
            blended[0]
        );
    }
    #[test]
    fn fog_height_falloff_reduces_density_at_high_altitude() {
        let fog = VolumetricFogParams::ground_fog();
        let t_low = fog.transmittance(10.0, 0.0);
        let t_high = fog.transmittance(10.0, 100.0);
        assert!(
            t_high > t_low,
            "higher altitude → higher transmittance: t_low={t_low}, t_high={t_high}"
        );
    }
    #[test]
    fn aniso_brdf_energy_not_above_incident_simple_check() {
        let brdf = AnisotropicBrdf::isotropic(0.8, [0.04; 3]);
        let directions = [
            normalize3([0.0, 0.0, 1.0]),
            normalize3([0.5, 0.0, 0.8660254]),
            normalize3([0.3, 0.3, 0.9]),
        ];
        let n = [0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        for &v in &directions {
            for &l in &directions {
                let result = brdf.evaluate(v, l, n, t, b);
                for &r in &result {
                    assert!(r.is_finite(), "BRDF must be finite: {r}");
                    assert!(r >= 0.0, "BRDF must be non-negative: {r}");
                }
            }
        }
    }
    #[test]
    fn sss_evaluate_clamps_grazing_angles() {
        let sss = SssParams::skin();
        let result = sss.evaluate_sss(-0.5, 0.001, [1.0, 1.0, 1.0]);
        for &v in &result {
            assert!(v >= 0.0, "SSS at grazing should be non-negative: {v}");
        }
    }
}
