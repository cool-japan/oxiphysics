//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::type_complexity)]
#![allow(clippy::needless_range_loop)]
use crate::primitives::{Color, RenderMesh, Vertex};
use oxiphysics_core::math::Vec3;

use super::types::MeshData;

/// Convert a `Vec3` (f64) to an `[f32; 3]` array.
pub(super) fn v3_to_f32(v: Vec3) -> [f32; 3] {
    [v.x as f32, v.y as f32, v.z as f32]
}
/// Generate a UV sphere mesh.
///
/// `subdivisions` controls the number of latitude/longitude segments.
pub fn sphere_mesh(center: Vec3, radius: f64, subdivisions: usize, color: Color) -> RenderMesh {
    let stacks = subdivisions.max(3);
    let slices = subdivisions.max(3) * 2;
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=stacks {
        let phi = std::f64::consts::PI * (i as f64) / (stacks as f64);
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        for j in 0..=slices {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / (slices as f64);
            let sin_theta = theta.sin();
            let cos_theta = theta.cos();
            let nx = sin_phi * cos_theta;
            let ny = cos_phi;
            let nz = sin_phi * sin_theta;
            let pos = center + Vec3::new(nx * radius, ny * radius, nz * radius);
            vertices.push(Vertex {
                position: v3_to_f32(pos),
                normal: [nx as f32, ny as f32, nz as f32],
                color,
            });
        }
    }
    for i in 0..stacks {
        for j in 0..slices {
            let row0 = i * (slices + 1);
            let row1 = (i + 1) * (slices + 1);
            indices.push((row0 + j) as u32);
            indices.push((row1 + j) as u32);
            indices.push((row1 + j + 1) as u32);
            indices.push((row0 + j) as u32);
            indices.push((row1 + j + 1) as u32);
            indices.push((row0 + j + 1) as u32);
        }
    }
    RenderMesh { vertices, indices }
}
/// Generate a box mesh with 24 vertices (unique normals per face) and 36 indices.
pub fn box_mesh(center: Vec3, half_extents: Vec3, color: Color) -> RenderMesh {
    let hx = half_extents.x;
    let hy = half_extents.y;
    let hz = half_extents.z;
    let face_data: [([f64; 3], [[f64; 3]; 4]); 6] = [
        (
            [1.0, 0.0, 0.0],
            [[hx, -hy, -hz], [hx, hy, -hz], [hx, hy, hz], [hx, -hy, hz]],
        ),
        (
            [-1.0, 0.0, 0.0],
            [
                [-hx, -hy, hz],
                [-hx, hy, hz],
                [-hx, hy, -hz],
                [-hx, -hy, -hz],
            ],
        ),
        (
            [0.0, 1.0, 0.0],
            [[-hx, hy, -hz], [hx, hy, -hz], [hx, hy, hz], [-hx, hy, hz]],
        ),
        (
            [0.0, -1.0, 0.0],
            [
                [-hx, -hy, hz],
                [hx, -hy, hz],
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
            ],
        ),
        (
            [0.0, 0.0, 1.0],
            [[-hx, -hy, hz], [hx, -hy, hz], [hx, hy, hz], [-hx, hy, hz]],
        ),
        (
            [0.0, 0.0, -1.0],
            [
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
                [-hx, hy, -hz],
                [hx, hy, -hz],
            ],
        ),
    ];
    let mut vertices = Vec::with_capacity(24);
    let mut indices = Vec::with_capacity(36);
    for (normal, positions) in &face_data {
        let base = vertices.len() as u32;
        for pos in positions {
            let p = center + Vec3::new(pos[0], pos[1], pos[2]);
            vertices.push(Vertex {
                position: v3_to_f32(p),
                normal: [normal[0] as f32, normal[1] as f32, normal[2] as f32],
                color,
            });
        }
        indices.push(base);
        indices.push(base + 1);
        indices.push(base + 2);
        indices.push(base);
        indices.push(base + 2);
        indices.push(base + 3);
    }
    RenderMesh { vertices, indices }
}
/// Generate a plane mesh (quad) centered at `center` with the given `normal`.
///
/// The quad has side length `2 * half_size`.
pub fn plane_mesh(center: Vec3, normal: Vec3, half_size: f64, color: Color) -> RenderMesh {
    let n = normal.normalize();
    let up = if n.dot(&Vec3::new(0.0, 1.0, 0.0)).abs() < 0.99 {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let tangent = n.cross(&up).normalize();
    let bitangent = n.cross(&tangent);
    let corners = [
        center + (-tangent - bitangent) * half_size,
        center + (tangent - bitangent) * half_size,
        center + (tangent + bitangent) * half_size,
        center + (-tangent + bitangent) * half_size,
    ];
    let nf = v3_to_f32(n);
    let vertices: Vec<Vertex> = corners
        .iter()
        .map(|&p| Vertex {
            position: v3_to_f32(p),
            normal: nf,
            color,
        })
        .collect();
    let indices = vec![0, 1, 2, 0, 2, 3];
    RenderMesh { vertices, indices }
}
/// Normalize a 3-element f64 vector in place. Returns the magnitude.
#[allow(dead_code)]
pub(super) fn normalize3_f64(v: &mut [f64; 3]) -> f64 {
    let m = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if m > 1e-30 {
        v[0] /= m;
        v[1] /= m;
        v[2] /= m;
    }
    m
}
/// Generate a UV sphere mesh.
///
/// Returns `(vertices, triangles, normals)` where normals are per-vertex.
#[allow(dead_code)]
pub fn generate_uv_sphere(
    center: [f64; 3],
    radius: f64,
    rings: usize,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let r = rings.max(2);
    let s = sectors.max(3);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    for i in 0..=r {
        let phi = std::f64::consts::PI * (i as f64) / (r as f64);
        let sp = phi.sin();
        let cp = phi.cos();
        for j in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / (s as f64);
            let nx = sp * theta.cos();
            let ny = cp;
            let nz = sp * theta.sin();
            verts.push([
                center[0] + radius * nx,
                center[1] + radius * ny,
                center[2] + radius * nz,
            ]);
            normals.push([nx, ny, nz]);
        }
    }
    for i in 0..r {
        for j in 0..s {
            let a = i * (s + 1) + j;
            let b = a + s + 1;
            tris.push([a, b, b + 1]);
            tris.push([a, b + 1, a + 1]);
        }
    }
    (verts, tris, normals)
}
/// Generate a cylinder mesh.
#[allow(dead_code)]
pub fn generate_cylinder(
    center: [f64; 3],
    radius: f64,
    height: f64,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let s = sectors.max(3);
    let half_h = height / 2.0;
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    for &y_off in &[-half_h, half_h] {
        for i in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
            let nx = theta.cos();
            let nz = theta.sin();
            verts.push([
                center[0] + radius * nx,
                center[1] + y_off,
                center[2] + radius * nz,
            ]);
            normals.push([nx, 0.0, nz]);
        }
    }
    for i in 0..s {
        let a = i;
        let b = i + s + 1;
        tris.push([a, b, b + 1]);
        tris.push([a, b + 1, a + 1]);
    }
    let bot_center_idx = verts.len();
    verts.push([center[0], center[1] - half_h, center[2]]);
    normals.push([0.0, -1.0, 0.0]);
    let top_center_idx = verts.len();
    verts.push([center[0], center[1] + half_h, center[2]]);
    normals.push([0.0, 1.0, 0.0]);
    for i in 0..s {
        tris.push([bot_center_idx, i + 1, i]);
        tris.push([top_center_idx, s + 1 + i, s + 1 + i + 1]);
    }
    (verts, tris, normals)
}
/// Generate a cone mesh.
#[allow(dead_code)]
pub fn generate_cone(
    center: [f64; 3],
    radius: f64,
    height: f64,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let s = sectors.max(3);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
        let nx = theta.cos();
        let nz = theta.sin();
        verts.push([center[0] + radius * nx, center[1], center[2] + radius * nz]);
        let slant = (radius * radius + height * height).sqrt();
        normals.push([nx * height / slant, radius / slant, nz * height / slant]);
    }
    let apex_idx = verts.len();
    verts.push([center[0], center[1] + height, center[2]]);
    normals.push([0.0, 1.0, 0.0]);
    for i in 0..s {
        tris.push([i, i + 1, apex_idx]);
    }
    let base_center_idx = verts.len();
    verts.push([center[0], center[1], center[2]]);
    normals.push([0.0, -1.0, 0.0]);
    for i in 0..s {
        tris.push([base_center_idx, i + 1, i]);
    }
    (verts, tris, normals)
}
/// Generate a torus mesh.
#[allow(dead_code)]
pub fn generate_torus(
    center: [f64; 3],
    major_r: f64,
    minor_r: f64,
    major_segs: usize,
    minor_segs: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let ms = major_segs.max(3);
    let ns = minor_segs.max(3);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    for i in 0..=ms {
        let u = 2.0 * std::f64::consts::PI * (i as f64) / (ms as f64);
        let cu = u.cos();
        let su = u.sin();
        for j in 0..=ns {
            let v = 2.0 * std::f64::consts::PI * (j as f64) / (ns as f64);
            let cv = v.cos();
            let sv = v.sin();
            let x = (major_r + minor_r * cv) * cu;
            let y = minor_r * sv;
            let z = (major_r + minor_r * cv) * su;
            verts.push([center[0] + x, center[1] + y, center[2] + z]);
            normals.push([cv * cu, sv, cv * su]);
        }
    }
    for i in 0..ms {
        for j in 0..ns {
            let a = i * (ns + 1) + j;
            let b = a + ns + 1;
            tris.push([a, b, b + 1]);
            tris.push([a, b + 1, a + 1]);
        }
    }
    (verts, tris, normals)
}
/// Generate a subdivided plane mesh.
#[allow(dead_code)]
pub fn generate_plane(
    center: [f64; 3],
    normal: [f64; 3],
    size: f64,
    subdivisions: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let n_div = subdivisions.max(1);
    let mut n = normal;
    normalize3_f64(&mut n);
    let up = if n[1].abs() < 0.99 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut t = [
        n[1] * up[2] - n[2] * up[1],
        n[2] * up[0] - n[0] * up[2],
        n[0] * up[1] - n[1] * up[0],
    ];
    normalize3_f64(&mut t);
    let b = [
        n[1] * t[2] - n[2] * t[1],
        n[2] * t[0] - n[0] * t[2],
        n[0] * t[1] - n[1] * t[0],
    ];
    let half = size / 2.0;
    let step = size / n_div as f64;
    let mut verts = Vec::new();
    let mut normals_out = Vec::new();
    let mut tris = Vec::new();
    for iy in 0..=n_div {
        for ix in 0..=n_div {
            let u = -half + ix as f64 * step;
            let v = -half + iy as f64 * step;
            verts.push([
                center[0] + u * t[0] + v * b[0],
                center[1] + u * t[1] + v * b[1],
                center[2] + u * t[2] + v * b[2],
            ]);
            normals_out.push(n);
        }
    }
    let w = n_div + 1;
    for iy in 0..n_div {
        for ix in 0..n_div {
            let a = iy * w + ix;
            tris.push([a, a + w, a + w + 1]);
            tris.push([a, a + w + 1, a + 1]);
        }
    }
    (verts, tris, normals_out)
}
/// Generate a capsule mesh (hemisphere + cylinder + hemisphere).
#[allow(dead_code)]
pub fn generate_capsule(
    p0: [f64; 3],
    p1: [f64; 3],
    radius: f64,
    rings: usize,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let r = rings.max(2);
    let s = sectors.max(3);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    for i in 0..=r {
        let phi = std::f64::consts::PI * 0.5 + std::f64::consts::PI * 0.5 * (i as f64) / (r as f64);
        let sp = phi.sin();
        let cp = phi.cos();
        for j in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / (s as f64);
            let nx = sp * theta.cos();
            let ny = cp;
            let nz = sp * theta.sin();
            verts.push([
                p0[0] + radius * nx,
                p0[1] + radius * ny,
                p0[2] + radius * nz,
            ]);
            normals.push([nx, ny, nz]);
        }
    }
    let bottom_count = verts.len();
    for i in 0..=r {
        let phi = std::f64::consts::PI * 0.5 * (i as f64) / (r as f64);
        let sp = phi.sin();
        let cp = phi.cos();
        for j in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * (j as f64) / (s as f64);
            let nx = sp * theta.cos();
            let ny = cp;
            let nz = sp * theta.sin();
            verts.push([
                p1[0] + radius * nx,
                p1[1] + radius * ny,
                p1[2] + radius * nz,
            ]);
            normals.push([nx, ny, nz]);
        }
    }
    for i in 0..r {
        for j in 0..s {
            let a = i * (s + 1) + j;
            let b = a + s + 1;
            tris.push([a, b, b + 1]);
            tris.push([a, b + 1, a + 1]);
        }
    }
    for i in 0..r {
        for j in 0..s {
            let a = bottom_count + i * (s + 1) + j;
            let b = a + s + 1;
            tris.push([a, b, b + 1]);
            tris.push([a, b + 1, a + 1]);
        }
    }
    (verts, tris, normals)
}
/// Generate an arrow mesh from `start` to `end` for force/velocity visualization.
///
/// The arrow consists of a cylindrical shaft and a cone tip.
pub fn arrow_mesh(start: Vec3, end: Vec3, shaft_radius: f64, color: Color) -> RenderMesh {
    let dir = end - start;
    let length = dir.norm();
    if length < 1e-12 {
        return RenderMesh::new();
    }
    let axis = dir / length;
    let up = if axis.dot(&Vec3::new(0.0, 1.0, 0.0)).abs() < 0.99 {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let tangent = axis.cross(&up).normalize();
    let bitangent = axis.cross(&tangent);
    let segments = 8;
    let cone_length = length * 0.2;
    let cone_radius = shaft_radius * 2.0;
    let shaft_length = length - cone_length;
    let shaft_end = start + axis * shaft_length;
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for i in 0..segments {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (segments as f64);
        let cos_a = a0.cos();
        let sin_a = a0.sin();
        let offset = tangent * (shaft_radius * cos_a) + bitangent * (shaft_radius * sin_a);
        let n = (tangent * cos_a + bitangent * sin_a).normalize();
        let nf = v3_to_f32(n);
        vertices.push(Vertex {
            position: v3_to_f32(start + offset),
            normal: nf,
            color,
        });
        vertices.push(Vertex {
            position: v3_to_f32(shaft_end + offset),
            normal: nf,
            color,
        });
    }
    for i in 0..segments {
        let i0 = (i * 2) as u32;
        let i1 = (i * 2 + 1) as u32;
        let i2 = (((i + 1) % segments) * 2) as u32;
        let i3 = (((i + 1) % segments) * 2 + 1) as u32;
        indices.extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
    }
    let cone_base_start = vertices.len() as u32;
    let tip = end;
    let tip_normal = v3_to_f32(axis);
    let cone_center_idx = vertices.len() as u32;
    vertices.push(Vertex {
        position: v3_to_f32(shaft_end),
        normal: v3_to_f32(-axis),
        color,
    });
    for i in 0..segments {
        let a0 = 2.0 * std::f64::consts::PI * (i as f64) / (segments as f64);
        let cos_a = a0.cos();
        let sin_a = a0.sin();
        let offset = tangent * (cone_radius * cos_a) + bitangent * (cone_radius * sin_a);
        let n = (tangent * cos_a + bitangent * sin_a).normalize();
        vertices.push(Vertex {
            position: v3_to_f32(shaft_end + offset),
            normal: v3_to_f32(n * 0.894 + axis * 0.447),
            color,
        });
    }
    let tip_idx = vertices.len() as u32;
    vertices.push(Vertex {
        position: v3_to_f32(tip),
        normal: tip_normal,
        color,
    });
    for i in 0..segments {
        let base_i = cone_base_start + 1 + (i as u32);
        let base_next = cone_base_start + 1 + (((i + 1) % segments) as u32);
        indices.extend_from_slice(&[base_i, base_next, tip_idx]);
    }
    for i in 0..segments {
        let base_i = cone_base_start + 1 + (i as u32);
        let base_next = cone_base_start + 1 + (((i + 1) % segments) as u32);
        indices.extend_from_slice(&[cone_center_idx, base_next, base_i]);
    }
    RenderMesh { vertices, indices }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn assert_unit_normals(normals: &[[f64; 3]]) {
        for (i, n) in normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (mag - 1.0).abs() < 0.02,
                "normal {i} magnitude = {mag}, expected ~1.0"
            );
        }
    }
    #[test]
    fn test_uv_sphere_vertex_count() {
        let rings = 4;
        let sectors = 6;
        let (verts, _tris, normals) = generate_uv_sphere([0.0; 3], 1.0, rings, sectors);
        let expected_verts = (rings + 1) * (sectors + 1);
        assert_eq!(verts.len(), expected_verts);
        assert_eq!(normals.len(), expected_verts);
    }
    #[test]
    fn test_uv_sphere_triangle_count() {
        let rings = 4;
        let sectors = 6;
        let (_verts, tris, _normals) = generate_uv_sphere([0.0; 3], 1.0, rings, sectors);
        let expected_tris = rings * sectors * 2;
        assert_eq!(tris.len(), expected_tris);
    }
    #[test]
    fn test_uv_sphere_normals_unit() {
        let (_verts, _tris, normals) = generate_uv_sphere([0.0; 3], 1.0, 4, 6);
        assert_unit_normals(&normals);
    }
    #[test]
    fn test_cylinder_vertex_count() {
        let sectors = 8;
        let (verts, _tris, _normals) = generate_cylinder([0.0; 3], 1.0, 2.0, sectors);
        let expected = 2 * (sectors + 1) + 2;
        assert_eq!(verts.len(), expected);
    }
    #[test]
    fn test_cone_has_apex() {
        let (verts, _tris, _normals) = generate_cone([0.0; 3], 1.0, 3.0, 6);
        let apex = verts.iter().find(|v| (v[1] - 3.0).abs() < 1e-10);
        assert!(apex.is_some(), "cone should have apex vertex at height");
    }
    #[test]
    fn test_torus_vertex_count() {
        let ms = 8;
        let ns = 6;
        let (verts, tris, normals) = generate_torus([0.0; 3], 2.0, 0.5, ms, ns);
        let expected_verts = (ms + 1) * (ns + 1);
        assert_eq!(verts.len(), expected_verts);
        assert_eq!(normals.len(), expected_verts);
        let expected_tris = ms * ns * 2;
        assert_eq!(tris.len(), expected_tris);
    }
    #[test]
    fn test_plane_vertex_count() {
        let sub = 3;
        let (verts, tris, normals) = generate_plane([0.0; 3], [0.0, 1.0, 0.0], 2.0, sub);
        let expected_verts = (sub + 1) * (sub + 1);
        assert_eq!(verts.len(), expected_verts);
        assert_eq!(normals.len(), expected_verts);
        let expected_tris = sub * sub * 2;
        assert_eq!(tris.len(), expected_tris);
    }
    #[test]
    fn test_capsule_non_empty() {
        let (verts, tris, normals) = generate_capsule([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5, 3, 6);
        assert!(!verts.is_empty());
        assert!(!tris.is_empty());
        assert_eq!(verts.len(), normals.len());
    }
    #[test]
    fn test_torus_normals_unit() {
        let (_verts, _tris, normals) = generate_torus([0.0; 3], 2.0, 0.5, 8, 6);
        assert_unit_normals(&normals);
    }
}
/// Generate a flat disc mesh in the XZ plane centered at `center`.
///
/// `inner_radius` may be 0.0 for a solid disc, or positive for an annulus.
/// Returns `(vertices, triangles, normals)` where all normals point along +Y.
#[allow(dead_code)]
pub fn generate_disc(
    center: [f64; 3],
    inner_radius: f64,
    outer_radius: f64,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let s = sectors.max(3);
    let r_inner = inner_radius.max(0.0);
    let r_outer = outer_radius.max(r_inner + 1e-9);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    if r_inner < 1e-12 {
        let center_idx = verts.len();
        verts.push(center);
        normals.push([0.0, 1.0, 0.0]);
        let ring_base = verts.len();
        for i in 0..s {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
            verts.push([
                center[0] + r_outer * theta.cos(),
                center[1],
                center[2] + r_outer * theta.sin(),
            ]);
            normals.push([0.0, 1.0, 0.0]);
        }
        for i in 0..s {
            tris.push([center_idx, ring_base + i, ring_base + (i + 1) % s]);
        }
    } else {
        let inner_base = verts.len();
        for i in 0..s {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
            verts.push([
                center[0] + r_inner * theta.cos(),
                center[1],
                center[2] + r_inner * theta.sin(),
            ]);
            normals.push([0.0, 1.0, 0.0]);
        }
        let outer_base = verts.len();
        for i in 0..s {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
            verts.push([
                center[0] + r_outer * theta.cos(),
                center[1],
                center[2] + r_outer * theta.sin(),
            ]);
            normals.push([0.0, 1.0, 0.0]);
        }
        for i in 0..s {
            let i0 = inner_base + i;
            let i1 = inner_base + (i + 1) % s;
            let o0 = outer_base + i;
            let o1 = outer_base + (i + 1) % s;
            tris.push([i0, o0, o1]);
            tris.push([i0, o1, i1]);
        }
    }
    (verts, tris, normals)
}
/// Generate a terrain (heightfield) mesh from a height function `h(x, z) -> y`.
///
/// The terrain spans `[x_min..x_max] × [z_min..z_max]` with `nx` × `nz` cells.
/// Normals are computed from central differences.
/// Returns `(vertices, triangles, normals)`.
#[allow(dead_code)]
pub fn generate_terrain<F>(
    x_min: f64,
    x_max: f64,
    z_min: f64,
    z_max: f64,
    nx: usize,
    nz: usize,
    height_fn: F,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>)
where
    F: Fn(f64, f64) -> f64,
{
    let nx = nx.max(1);
    let nz = nz.max(1);
    let dx = (x_max - x_min) / nx as f64;
    let dz = (z_max - z_min) / nz as f64;
    let rows = nz + 1;
    let cols = nx + 1;
    let mut verts: Vec<[f64; 3]> = Vec::with_capacity(rows * cols);
    let mut normals: Vec<[f64; 3]> = Vec::with_capacity(rows * cols);
    let mut tris: Vec<[usize; 3]> = Vec::new();
    for iz in 0..rows {
        for ix in 0..cols {
            let x = x_min + ix as f64 * dx;
            let z = z_min + iz as f64 * dz;
            let y = height_fn(x, z);
            verts.push([x, y, z]);
            normals.push([0.0, 1.0, 0.0]);
        }
    }
    for iz in 0..nz {
        for ix in 0..nx {
            let a = iz * cols + ix;
            let b = a + 1;
            let c = a + cols;
            let d = c + 1;
            tris.push([a, c, d]);
            tris.push([a, d, b]);
        }
    }
    for iz in 0..rows {
        for ix in 0..cols {
            let idx = iz * cols + ix;
            let xp = if ix + 1 < cols {
                verts[(iz) * cols + ix + 1][1]
            } else {
                verts[idx][1]
            };
            let xm = if ix > 0 {
                verts[(iz) * cols + ix - 1][1]
            } else {
                verts[idx][1]
            };
            let zp = if iz + 1 < rows {
                verts[(iz + 1) * cols + ix][1]
            } else {
                verts[idx][1]
            };
            let zm = if iz > 0 {
                verts[(iz - 1) * cols + ix][1]
            } else {
                verts[idx][1]
            };
            let nx_comp = -(xp - xm) / (2.0 * dx.max(1e-12));
            let nz_comp = -(zp - zm) / (2.0 * dz.max(1e-12));
            let ny_comp = 1.0_f64;
            let mag = (nx_comp * nx_comp + ny_comp * ny_comp + nz_comp * nz_comp)
                .sqrt()
                .max(1e-30);
            normals[idx] = [nx_comp / mag, ny_comp / mag, nz_comp / mag];
        }
    }
    (verts, tris, normals)
}
/// Generate an arrow glyph mesh from `origin` pointing in `direction` for
/// vector field visualization.
///
/// The arrow is an octagonal prism shaft plus a cone tip.
/// Returns `(vertices, triangles, normals)`.
#[allow(dead_code)]
pub fn generate_arrow_glyph(
    origin: [f64; 3],
    direction: [f64; 3],
    shaft_radius: f64,
    tip_radius: f64,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let s = sectors.max(3);
    let len =
        (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
            .sqrt();
    if len < 1e-12 {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let axis = [direction[0] / len, direction[1] / len, direction[2] / len];
    let up = if axis[1].abs() < 0.99 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut tangent = [
        axis[1] * up[2] - axis[2] * up[1],
        axis[2] * up[0] - axis[0] * up[2],
        axis[0] * up[1] - axis[1] * up[0],
    ];
    normalize3_f64(&mut tangent);
    let bitangent = [
        axis[1] * tangent[2] - axis[2] * tangent[1],
        axis[2] * tangent[0] - axis[0] * tangent[2],
        axis[0] * tangent[1] - axis[1] * tangent[0],
    ];
    let shaft_len = len * 0.75;
    let tip_len = len - shaft_len;
    let shaft_end = [
        origin[0] + axis[0] * shaft_len,
        origin[1] + axis[1] * shaft_len,
        origin[2] + axis[2] * shaft_len,
    ];
    let tip = [
        origin[0] + direction[0],
        origin[1] + direction[1],
        origin[2] + direction[2],
    ];
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    let shaft_bot_base = verts.len();
    for i in 0..s {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        let nx = cos_t * tangent[0] + sin_t * bitangent[0];
        let ny = cos_t * tangent[1] + sin_t * bitangent[1];
        let nz = cos_t * tangent[2] + sin_t * bitangent[2];
        verts.push([
            origin[0] + shaft_radius * (cos_t * tangent[0] + sin_t * bitangent[0]),
            origin[1] + shaft_radius * (cos_t * tangent[1] + sin_t * bitangent[1]),
            origin[2] + shaft_radius * (cos_t * tangent[2] + sin_t * bitangent[2]),
        ]);
        normals.push([nx, ny, nz]);
    }
    let shaft_top_base = verts.len();
    for i in 0..s {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        let nx = cos_t * tangent[0] + sin_t * bitangent[0];
        let ny = cos_t * tangent[1] + sin_t * bitangent[1];
        let nz = cos_t * tangent[2] + sin_t * bitangent[2];
        verts.push([
            shaft_end[0] + shaft_radius * (cos_t * tangent[0] + sin_t * bitangent[0]),
            shaft_end[1] + shaft_radius * (cos_t * tangent[1] + sin_t * bitangent[1]),
            shaft_end[2] + shaft_radius * (cos_t * tangent[2] + sin_t * bitangent[2]),
        ]);
        normals.push([nx, ny, nz]);
    }
    for i in 0..s {
        let b0 = shaft_bot_base + i;
        let b1 = shaft_bot_base + (i + 1) % s;
        let t0 = shaft_top_base + i;
        let t1 = shaft_top_base + (i + 1) % s;
        tris.push([b0, t0, t1]);
        tris.push([b0, t1, b1]);
    }
    let cone_base_idx = verts.len();
    let slant = (tip_radius * tip_radius + tip_len * tip_len)
        .sqrt()
        .max(1e-12);
    for i in 0..s {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        let nt = [
            cos_t * tangent[0] + sin_t * bitangent[0],
            cos_t * tangent[1] + sin_t * bitangent[1],
            cos_t * tangent[2] + sin_t * bitangent[2],
        ];
        let nx = nt[0] * tip_len / slant + axis[0] * tip_radius / slant;
        let ny = nt[1] * tip_len / slant + axis[1] * tip_radius / slant;
        let nz = nt[2] * tip_len / slant + axis[2] * tip_radius / slant;
        verts.push([
            shaft_end[0] + tip_radius * (cos_t * tangent[0] + sin_t * bitangent[0]),
            shaft_end[1] + tip_radius * (cos_t * tangent[1] + sin_t * bitangent[1]),
            shaft_end[2] + tip_radius * (cos_t * tangent[2] + sin_t * bitangent[2]),
        ]);
        normals.push([nx, ny, nz]);
    }
    let apex_idx = verts.len();
    verts.push(tip);
    normals.push(axis);
    for i in 0..s {
        let b0 = cone_base_idx + i;
        let b1 = cone_base_idx + (i + 1) % s;
        tris.push([b0, b1, apex_idx]);
    }
    let bot_center_idx = verts.len();
    verts.push(origin);
    normals.push([-axis[0], -axis[1], -axis[2]]);
    for i in 0..s {
        tris.push([
            bot_center_idx,
            shaft_bot_base + (i + 1) % s,
            shaft_bot_base + i,
        ]);
    }
    (verts, tris, normals)
}
/// Generate an icosphere by subdividing an icosahedron `subdivisions` times and
/// projecting all vertices onto a sphere of `radius`.
///
/// This is an alias to `generate_geodesic_sphere` using a more descriptive name.
/// Returns `(vertices, triangles, normals)`.
#[allow(dead_code)]
#[inline]
pub fn generate_icosphere(
    center: [f64; 3],
    radius: f64,
    subdivisions: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    generate_geodesic_sphere(center, radius, subdivisions)
}
/// Generate a geodesic sphere by subdividing an icosahedron `subdivisions` times.
///
/// Returns `(vertices, triangles, normals)`.
#[allow(dead_code)]
pub fn generate_geodesic_sphere(
    center: [f64; 3],
    radius: f64,
    subdivisions: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let raw: [[f64; 3]; 12] = [
        [-1.0, phi, 0.0],
        [1.0, phi, 0.0],
        [-1.0, -phi, 0.0],
        [1.0, -phi, 0.0],
        [0.0, -1.0, phi],
        [0.0, 1.0, phi],
        [0.0, -1.0, -phi],
        [0.0, 1.0, -phi],
        [phi, 0.0, -1.0],
        [phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
    ];
    let mut verts: Vec<[f64; 3]> = raw
        .iter()
        .map(|v| {
            let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            [v[0] / mag, v[1] / mag, v[2] / mag]
        })
        .collect();
    let mut tris: Vec<[usize; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..subdivisions {
        let mut new_tris = Vec::new();
        let mut midpoint_cache: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        let get_midpoint = |a: usize,
                            b: usize,
                            verts: &mut Vec<[f64; 3]>,
                            cache: &mut std::collections::HashMap<(usize, usize), usize>|
         -> usize {
            let key = (a.min(b), a.max(b));
            if let Some(&idx) = cache.get(&key) {
                return idx;
            }
            let va = verts[a];
            let vb = verts[b];
            let mut mid = [
                (va[0] + vb[0]) * 0.5,
                (va[1] + vb[1]) * 0.5,
                (va[2] + vb[2]) * 0.5,
            ];
            let mag = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2])
                .sqrt()
                .max(1e-30);
            mid[0] /= mag;
            mid[1] /= mag;
            mid[2] /= mag;
            let idx = verts.len();
            verts.push(mid);
            cache.insert(key, idx);
            idx
        };
        for &[a, b, c] in &tris {
            let ab = get_midpoint(a, b, &mut verts, &mut midpoint_cache);
            let bc = get_midpoint(b, c, &mut verts, &mut midpoint_cache);
            let ca = get_midpoint(c, a, &mut verts, &mut midpoint_cache);
            new_tris.push([a, ab, ca]);
            new_tris.push([b, bc, ab]);
            new_tris.push([c, ca, bc]);
            new_tris.push([ab, bc, ca]);
        }
        tris = new_tris;
    }
    let normals: Vec<[f64; 3]> = verts.to_vec();
    let out_verts: Vec<[f64; 3]> = verts
        .iter()
        .map(|v| {
            [
                center[0] + v[0] * radius,
                center[1] + v[1] * radius,
                center[2] + v[2] * radius,
            ]
        })
        .collect();
    (out_verts, tris, normals)
}
/// Generate a frustum (truncated cone) mesh.
///
/// `r_bottom` and `r_top` are the bottom and top radii, `height` is the height.
#[allow(dead_code)]
pub fn generate_frustum(
    center: [f64; 3],
    r_bottom: f64,
    r_top: f64,
    height: f64,
    sectors: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let s = sectors.max(3);
    let half_h = height / 2.0;
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    let mut tris = Vec::new();
    let slant = ((r_bottom - r_top).powi(2) + height.powi(2))
        .sqrt()
        .max(1e-30);
    let ny_comp = (r_bottom - r_top) / slant;
    let r_comp = height / slant;
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (s as f64);
        let ct = theta.cos();
        let st = theta.sin();
        verts.push([
            center[0] + r_bottom * ct,
            center[1] - half_h,
            center[2] + r_bottom * st,
        ]);
        normals.push([r_comp * ct, ny_comp, r_comp * st]);
        verts.push([
            center[0] + r_top * ct,
            center[1] + half_h,
            center[2] + r_top * st,
        ]);
        normals.push([r_comp * ct, ny_comp, r_comp * st]);
    }
    for i in 0..s {
        let b0 = i * 2;
        let b1 = i * 2 + 1;
        let t0 = (i + 1) * 2;
        let t1 = (i + 1) * 2 + 1;
        tris.push([b0, t0, b1]);
        tris.push([b1, t0, t1]);
    }
    let bot_center = verts.len();
    verts.push([center[0], center[1] - half_h, center[2]]);
    normals.push([0.0, -1.0, 0.0]);
    for i in 0..s {
        let a = i * 2;
        let b = ((i + 1) % s) * 2;
        tris.push([bot_center, b, a]);
    }
    let top_center = verts.len();
    verts.push([center[0], center[1] + half_h, center[2]]);
    normals.push([0.0, 1.0, 0.0]);
    for i in 0..s {
        let a = i * 2 + 1;
        let b = ((i + 1) % s) * 2 + 1;
        tris.push([top_center, a, b]);
    }
    (verts, tris, normals)
}
/// Bicubic Bernstein basis function B_{i,3}(t).
pub(super) fn bernstein3(i: usize, t: f64) -> f64 {
    let binom = [1.0, 3.0, 3.0, 1.0];
    binom[i] * t.powi(i as i32) * (1.0 - t).powi((3 - i) as i32)
}
/// Generate a parametric bicubic Bézier patch surface.
///
/// `control_pts[4][4]` are 16 control points in 3D.
/// `u_div` and `v_div` are tessellation subdivisions.
#[allow(dead_code)]
pub fn generate_bicubic_patch(
    control_pts: &[[f64; 3]; 16],
    u_div: usize,
    v_div: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<[f64; 3]>) {
    let u_n = u_div.max(1);
    let v_n = v_div.max(1);
    let mut verts = Vec::new();
    let mut normals = Vec::new();
    for iv in 0..=v_n {
        let v = iv as f64 / v_n as f64;
        for iu in 0..=u_n {
            let u = iu as f64 / u_n as f64;
            let mut p = [0.0_f64; 3];
            for j in 0..4 {
                let bv = bernstein3(j, v);
                for i in 0..4 {
                    let bu = bernstein3(i, u);
                    let cp = control_pts[j * 4 + i];
                    for k in 0..3 {
                        p[k] += bu * bv * cp[k];
                    }
                }
            }
            verts.push(p);
            normals.push([0.0, 1.0, 0.0]);
        }
    }
    let row = u_n + 1;
    let mut tris = Vec::new();
    for iv in 0..v_n {
        for iu in 0..u_n {
            let a = iv * row + iu;
            let b = a + row;
            tris.push([a, b, b + 1]);
            tris.push([a, b + 1, a + 1]);
        }
    }
    (verts, tris, normals)
}
/// Generate a 3-axis gizmo mesh (RGB arrows for X, Y, Z).
///
/// Returns three meshes for X (red), Y (green), Z (blue) axes.
#[allow(dead_code)]
pub fn generate_axes_gizmo(
    origin: [f64; 3],
    length: f64,
    shaft_radius: f64,
    _sectors: usize,
) -> [(Vec<[f64; 3]>, Vec<[usize; 3]>); 3] {
    use crate::primitives::Color;
    let colors = [Color::red(), Color::green(), Color::blue()];
    let directions: [[f64; 3]; 3] = [[length, 0.0, 0.0], [0.0, length, 0.0], [0.0, 0.0, length]];
    std::array::from_fn(|axis| {
        let start = Vec3::new(origin[0], origin[1], origin[2]);
        let dir = directions[axis];
        let end = Vec3::new(origin[0] + dir[0], origin[1] + dir[1], origin[2] + dir[2]);
        let mesh = arrow_mesh(start, end, shaft_radius, colors[axis]);
        let verts: Vec<[f64; 3]> = mesh
            .vertices
            .iter()
            .map(|v| {
                [
                    v.position[0] as f64,
                    v.position[1] as f64,
                    v.position[2] as f64,
                ]
            })
            .collect();
        let tris: Vec<[usize; 3]> = mesh
            .indices
            .chunks(3)
            .map(|c| [c[0] as usize, c[1] as usize, c[2] as usize])
            .collect();
        (verts, tris)
    })
}
#[cfg(test)]
mod expanded_mesh_tests {
    use super::*;
    #[test]
    fn test_geodesic_sphere_vertex_count_grows_with_subdivisions() {
        let (v0, t0, _) = generate_geodesic_sphere([0.0; 3], 1.0, 0);
        let (v1, t1, _) = generate_geodesic_sphere([0.0; 3], 1.0, 1);
        assert!(
            v1.len() > v0.len(),
            "subdivision 1 should have more verts than 0"
        );
        assert!(
            t1.len() > t0.len(),
            "subdivision 1 should have more tris than 0"
        );
        assert_eq!(t0.len(), 20);
    }
    #[test]
    fn test_geodesic_sphere_normals_unit() {
        let (_, _, normals) = generate_geodesic_sphere([0.0; 3], 1.0, 1);
        for (i, n) in normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((mag - 1.0).abs() < 1e-10, "normal {i} mag = {mag}");
        }
    }
    #[test]
    fn test_geodesic_sphere_verts_on_surface() {
        let r = 2.5;
        let (verts, _, _) = generate_geodesic_sphere([1.0, 2.0, 3.0], r, 1);
        for (i, v) in verts.iter().enumerate() {
            let d = ((v[0] - 1.0).powi(2) + (v[1] - 2.0).powi(2) + (v[2] - 3.0).powi(2)).sqrt();
            assert!((d - r).abs() < 1e-10, "vertex {i} not on sphere: dist={d}");
        }
    }
    #[test]
    fn test_frustum_vertex_count() {
        let sectors = 6;
        let (verts, _tris, _normals) = generate_frustum([0.0; 3], 1.0, 0.5, 2.0, sectors);
        let expected = 2 * (sectors + 1) + 2;
        assert_eq!(
            verts.len(),
            expected,
            "frustum vertex count = {}",
            verts.len()
        );
    }
    #[test]
    fn test_frustum_cone_when_top_radius_zero() {
        let (_verts, tris, _normals) = generate_frustum([0.0; 3], 1.0, 0.0, 2.0, 6);
        assert!(!tris.is_empty(), "frustum (cone) should have triangles");
    }
    #[test]
    fn test_bicubic_patch_vertex_count() {
        let cp: [[f64; 3]; 16] = std::array::from_fn(|i| {
            let u = (i % 4) as f64;
            let v = (i / 4) as f64;
            [u, v, 0.0]
        });
        let u_div = 4;
        let v_div = 4;
        let (verts, tris, _) = generate_bicubic_patch(&cp, u_div, v_div);
        assert_eq!(verts.len(), (u_div + 1) * (v_div + 1));
        assert_eq!(tris.len(), u_div * v_div * 2);
    }
    #[test]
    fn test_bicubic_patch_flat_plane() {
        let cp: [[f64; 3]; 16] = std::array::from_fn(|i| {
            let u = (i % 4) as f64 / 3.0;
            let v = (i / 4) as f64 / 3.0;
            [u, v, 0.0]
        });
        let (verts, _, _) = generate_bicubic_patch(&cp, 4, 4);
        for (i, v) in verts.iter().enumerate() {
            assert!(v[2].abs() < 1e-12, "flat patch z[{i}] = {}", v[2]);
        }
    }
    #[test]
    fn test_axes_gizmo_three_axes() {
        let axes = generate_axes_gizmo([0.0; 3], 1.0, 0.05, 6);
        for (i, (verts, tris)) in axes.iter().enumerate() {
            assert!(!verts.is_empty(), "axis {i} should have vertices");
            assert!(!tris.is_empty(), "axis {i} should have triangles");
        }
    }
    #[test]
    fn test_generate_capsule_vertex_count() {
        let rings = 3;
        let sectors = 6;
        let (verts, tris, normals) =
            generate_capsule([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5, rings, sectors);
        assert!(!verts.is_empty());
        assert!(!tris.is_empty());
        assert_eq!(verts.len(), normals.len());
    }
    #[test]
    fn test_generate_cylinder_triangle_count() {
        let s = 8;
        let (_verts, tris, _normals) = generate_cylinder([0.0; 3], 1.0, 2.0, s);
        assert_eq!(tris.len(), 4 * s);
    }
}
#[cfg(test)]
mod new_mesh_tests {
    use super::*;
    fn assert_unit_normals_f64(normals: &[[f64; 3]], tol: f64) {
        for (i, n) in normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (mag - 1.0).abs() < tol,
                "normal {i} magnitude = {mag:.6}, expected ~1.0"
            );
        }
    }
    #[test]
    fn test_disc_solid_vertex_count() {
        let sectors = 8;
        let (verts, tris, normals) = generate_disc([0.0; 3], 0.0, 1.0, sectors);
        assert_eq!(verts.len(), sectors + 1);
        assert_eq!(normals.len(), sectors + 1);
        assert_eq!(tris.len(), sectors);
    }
    #[test]
    fn test_disc_solid_normals_y_up() {
        let (_, _, normals) = generate_disc([0.0; 3], 0.0, 2.0, 12);
        for (i, n) in normals.iter().enumerate() {
            assert!(
                (n[1] - 1.0).abs() < 1e-10,
                "disc normal[{i}] should be [0,1,0], got {:?}",
                n
            );
        }
    }
    #[test]
    fn test_disc_solid_verts_on_plane() {
        let center = [1.0, 3.0, -2.0];
        let r = 2.5;
        let (verts, _, _) = generate_disc(center, 0.0, r, 16);
        for (i, v) in verts.iter().enumerate() {
            assert!(
                (v[1] - center[1]).abs() < 1e-10,
                "vert {i} y={} != center y={}",
                v[1],
                center[1]
            );
            if i > 0 {
                let dx = v[0] - center[0];
                let dz = v[2] - center[2];
                let d = (dx * dx + dz * dz).sqrt();
                assert!((d - r).abs() < 1e-10, "vert {i} radius={d} != {r}");
            }
        }
    }
    #[test]
    fn test_annulus_vertex_count() {
        let sectors = 10;
        let (verts, tris, normals) = generate_disc([0.0; 3], 0.5, 1.5, sectors);
        assert_eq!(verts.len(), 2 * sectors);
        assert_eq!(normals.len(), 2 * sectors);
        assert_eq!(tris.len(), 2 * sectors);
    }
    #[test]
    fn test_annulus_inner_verts_at_inner_radius() {
        let r_inner = 1.0;
        let r_outer = 3.0;
        let sectors = 8;
        let center = [0.0; 3];
        let (verts, _, _) = generate_disc(center, r_inner, r_outer, sectors);
        for i in 0..sectors {
            let v = verts[i];
            let dx = v[0] - center[0];
            let dz = v[2] - center[2];
            let d = (dx * dx + dz * dz).sqrt();
            assert!(
                (d - r_inner).abs() < 1e-10,
                "inner vert {i} radius={d} != {r_inner}"
            );
        }
    }
    #[test]
    fn test_annulus_outer_verts_at_outer_radius() {
        let r_inner = 1.0;
        let r_outer = 3.0;
        let sectors = 8;
        let center = [0.0; 3];
        let (verts, _, _) = generate_disc(center, r_inner, r_outer, sectors);
        for i in sectors..2 * sectors {
            let v = verts[i];
            let dx = v[0] - center[0];
            let dz = v[2] - center[2];
            let d = (dx * dx + dz * dz).sqrt();
            assert!(
                (d - r_outer).abs() < 1e-10,
                "outer vert {i} radius={d} != {r_outer}"
            );
        }
    }
    #[test]
    fn test_annulus_normals_y_up() {
        let (_, _, normals) = generate_disc([0.0; 3], 1.0, 2.0, 8);
        assert_unit_normals_f64(&normals, 1e-9);
    }
    #[test]
    fn test_terrain_flat_vertex_count() {
        let nx = 4;
        let nz = 5;
        let (verts, tris, normals) = generate_terrain(0.0, 4.0, 0.0, 5.0, nx, nz, |_x, _z| 0.0);
        assert_eq!(verts.len(), (nx + 1) * (nz + 1));
        assert_eq!(normals.len(), (nx + 1) * (nz + 1));
        assert_eq!(tris.len(), nx * nz * 2);
    }
    #[test]
    fn test_terrain_flat_normals_y_up() {
        let (_, _, normals) = generate_terrain(-1.0, 1.0, -1.0, 1.0, 4, 4, |_x, _z| 0.0);
        assert_unit_normals_f64(&normals, 1e-9);
        for (i, n) in normals.iter().enumerate() {
            assert!(
                (n[1] - 1.0).abs() < 1e-9,
                "flat terrain normal[{i}] y={}",
                n[1]
            );
        }
    }
    #[test]
    fn test_terrain_height_values_applied() {
        let (verts, _, _) = generate_terrain(0.0, 1.0, 0.0, 1.0, 4, 4, |x, z| x + z);
        for v in &verts {
            let expected_y = v[0] + v[2];
            assert!(
                (v[1] - expected_y).abs() < 1e-9,
                "height mismatch: y={} expected={}",
                v[1],
                expected_y
            );
        }
    }
    #[test]
    fn test_terrain_sinusoidal_normals_unit() {
        let (_, _, normals) = generate_terrain(0.0, 4.0, 0.0, 4.0, 8, 8, |x, z| {
            (x * std::f64::consts::PI * 0.5).sin() * (z * std::f64::consts::PI * 0.5).cos()
        });
        assert_unit_normals_f64(&normals, 0.01);
    }
    #[test]
    fn test_terrain_triangle_indices_in_bounds() {
        let nx = 3;
        let nz = 3;
        let (verts, tris, _) = generate_terrain(0.0, 3.0, 0.0, 3.0, nx, nz, |_, _| 0.0);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex index {vi} out of bounds");
            }
        }
    }
    #[test]
    fn test_arrow_glyph_non_empty() {
        let (verts, tris, normals) = generate_arrow_glyph([0.0; 3], [0.0, 1.0, 0.0], 0.05, 0.15, 8);
        assert!(!verts.is_empty());
        assert!(!tris.is_empty());
        assert_eq!(verts.len(), normals.len());
    }
    #[test]
    fn test_arrow_glyph_zero_direction_returns_empty() {
        let (verts, tris, normals) = generate_arrow_glyph([0.0; 3], [0.0, 0.0, 0.0], 0.05, 0.15, 8);
        assert!(verts.is_empty());
        assert!(tris.is_empty());
        assert!(normals.is_empty());
    }
    #[test]
    fn test_arrow_glyph_apex_near_tip() {
        let origin = [0.0; 3];
        let dir = [3.0, 0.0, 0.0];
        let (verts, _, _) = generate_arrow_glyph(origin, dir, 0.1, 0.2, 8);
        let apex = verts
            .iter()
            .max_by(|a, b| a[0].partial_cmp(&b[0]).unwrap())
            .unwrap();
        assert!((apex[0] - 3.0).abs() < 1e-9, "apex x={}", apex[0]);
    }
    #[test]
    fn test_arrow_glyph_triangle_indices_in_bounds() {
        let (verts, tris, _) = generate_arrow_glyph([1.0, 2.0, 3.0], [0.0, 0.0, 2.0], 0.1, 0.3, 6);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} out of bounds");
            }
        }
    }
    #[test]
    fn test_icosphere_is_geodesic_sphere() {
        let (v0, t0, _) = generate_icosphere([0.0; 3], 1.0, 0);
        assert_eq!(t0.len(), 20, "icosahedron should have 20 faces");
        assert_eq!(v0.len(), 12, "icosahedron should have 12 vertices");
    }
    #[test]
    fn test_icosphere_subdivision_grows() {
        let (v1, t1, _) = generate_icosphere([0.0; 3], 1.0, 1);
        let (v2, t2, _) = generate_icosphere([0.0; 3], 1.0, 2);
        assert!(v2.len() > v1.len());
        assert!(t2.len() > t1.len());
    }
    #[test]
    fn test_icosphere_verts_on_sphere() {
        let r = 3.0;
        let (verts, _, _) = generate_icosphere([1.0, 2.0, 3.0], r, 1);
        for (i, v) in verts.iter().enumerate() {
            let d = ((v[0] - 1.0).powi(2) + (v[1] - 2.0).powi(2) + (v[2] - 3.0).powi(2)).sqrt();
            assert!((d - r).abs() < 1e-9, "vertex {i} dist={d} != {r}");
        }
    }
    #[test]
    fn test_icosphere_normals_unit() {
        let (_, _, normals) = generate_icosphere([0.0; 3], 2.0, 2);
        assert_unit_normals_f64(&normals, 1e-10);
    }
    #[test]
    fn test_uv_sphere_verts_on_surface() {
        let r = 1.5;
        let c = [1.0, 2.0, 3.0];
        let (verts, _, _) = generate_uv_sphere(c, r, 6, 8);
        for (i, v) in verts.iter().enumerate() {
            let d = ((v[0] - c[0]).powi(2) + (v[1] - c[1]).powi(2) + (v[2] - c[2]).powi(2)).sqrt();
            assert!((d - r).abs() < 1e-9, "vertex {i} dist={d} != {r}");
        }
    }
    #[test]
    fn test_cylinder_top_bottom_caps_exist() {
        let s = 6;
        let (verts, _, _) = generate_cylinder([0.0; 3], 1.0, 2.0, s);
        let has_top = verts.iter().any(|v| (v[1] - 1.0).abs() < 1e-9);
        let has_bot = verts.iter().any(|v| (v[1] + 1.0).abs() < 1e-9);
        assert!(has_top, "cylinder should have top cap center at y=1");
        assert!(has_bot, "cylinder should have bottom cap center at y=-1");
    }
    #[test]
    fn test_cone_base_cap_center_at_origin() {
        let (verts, _, _) = generate_cone([0.0; 3], 1.0, 2.0, 8);
        let has_base = verts
            .iter()
            .any(|v| v[0].abs() < 1e-10 && v[2].abs() < 1e-10 && v[1].abs() < 1e-10);
        assert!(has_base, "cone base center should be at (0,0,0)");
    }
    #[test]
    fn test_torus_normals_unit_high_seg() {
        let (_, _, normals) = generate_torus([0.0; 3], 3.0, 0.5, 16, 12);
        assert_unit_normals_f64(&normals, 0.01);
    }
    #[test]
    fn test_frustum_side_normals_outward() {
        let (verts, _, normals) = generate_frustum([0.0; 3], 1.0, 0.5, 2.0, 6);
        for (i, (v, n)) in verts.iter().zip(normals.iter()).enumerate() {
            if v[0].abs() < 1e-9 && v[2].abs() < 1e-9 {
                continue;
            }
            let radial_dot = v[0] * n[0] + v[2] * n[2];
            assert!(
                radial_dot >= -1e-6,
                "side normal[{i}] should have outward radial component, dot={radial_dot}"
            );
        }
    }
}
/// Generate a UV sphere as a `MeshData` with UV coordinates.
///
/// `stacks` and `slices` control latitude/longitude subdivisions.
#[allow(dead_code)]
pub fn generate_sphere(center: [f64; 3], radius: f64, stacks: usize, slices: usize) -> MeshData {
    let st = stacks.max(2);
    let sl = slices.max(3);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=st {
        let phi = std::f64::consts::PI * i as f64 / st as f64;
        let sp = phi.sin();
        let cp = phi.cos();
        let v_uv = i as f64 / st as f64;
        for j in 0..=sl {
            let theta = 2.0 * std::f64::consts::PI * j as f64 / sl as f64;
            let nx = sp * theta.cos();
            let ny = cp;
            let nz = sp * theta.sin();
            positions.push([
                center[0] + radius * nx,
                center[1] + radius * ny,
                center[2] + radius * nz,
            ]);
            normals.push([nx, ny, nz]);
            uvs.push([j as f64 / sl as f64, v_uv]);
        }
    }
    for i in 0..st {
        for j in 0..sl {
            let a = i * (sl + 1) + j;
            let b = a + sl + 1;
            indices.push([a, b, b + 1]);
            indices.push([a, b + 1, a + 1]);
        }
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate a box mesh as a `MeshData` with per-face normals and UV coordinates.
#[allow(dead_code)]
pub fn generate_box(center: [f64; 3], half_extents: [f64; 3]) -> MeshData {
    let [cx, cy, cz] = center;
    let [hx, hy, hz] = half_extents;
    let face_data: [([f64; 3], [[f64; 3]; 4], [[f64; 2]; 4]); 6] = [
        (
            [1.0, 0.0, 0.0],
            [[hx, -hy, -hz], [hx, hy, -hz], [hx, hy, hz], [hx, -hy, hz]],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
        (
            [-1.0, 0.0, 0.0],
            [
                [-hx, -hy, hz],
                [-hx, hy, hz],
                [-hx, hy, -hz],
                [-hx, -hy, -hz],
            ],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
        (
            [0.0, 1.0, 0.0],
            [[-hx, hy, -hz], [hx, hy, -hz], [hx, hy, hz], [-hx, hy, hz]],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
        (
            [0.0, -1.0, 0.0],
            [
                [-hx, -hy, hz],
                [hx, -hy, hz],
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
            ],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
        (
            [0.0, 0.0, 1.0],
            [[-hx, -hy, hz], [hx, -hy, hz], [hx, hy, hz], [-hx, hy, hz]],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
        (
            [0.0, 0.0, -1.0],
            [
                [hx, -hy, -hz],
                [-hx, -hy, -hz],
                [-hx, hy, -hz],
                [hx, hy, -hz],
            ],
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ),
    ];
    let mut positions = Vec::with_capacity(24);
    let mut normals = Vec::with_capacity(24);
    let mut uvs = Vec::with_capacity(24);
    let mut indices = Vec::with_capacity(12);
    for (normal, corners, face_uvs) in &face_data {
        let base = positions.len();
        for (pos_off, uv) in corners.iter().zip(face_uvs.iter()) {
            positions.push([cx + pos_off[0], cy + pos_off[1], cz + pos_off[2]]);
            normals.push(*normal);
            uvs.push(*uv);
        }
        indices.push([base, base + 1, base + 2]);
        indices.push([base, base + 2, base + 3]);
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate a cylinder as a `MeshData` with caps and UV coordinates.
#[allow(dead_code)]
pub fn generate_cylinder_mesh(
    center: [f64; 3],
    radius: f64,
    height: f64,
    sectors: usize,
) -> MeshData {
    let s = sectors.max(3);
    let [cx, cy, cz] = center;
    let half_h = height * 0.5;
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for (ring, y_off) in [(-half_h, 0.0_f64), (half_h, 1.0_f64)] {
        for i in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
            let nx = theta.cos();
            let nz = theta.sin();
            positions.push([cx + radius * nx, cy + ring, cz + radius * nz]);
            normals.push([nx, 0.0, nz]);
            uvs.push([i as f64 / s as f64, y_off]);
        }
    }
    for i in 0..s {
        let bot = i;
        let top = i + s + 1;
        indices.push([bot, top, top + 1]);
        indices.push([bot, top + 1, bot + 1]);
    }
    let bot_center = positions.len();
    positions.push([cx, cy - half_h, cz]);
    normals.push([0.0, -1.0, 0.0]);
    uvs.push([0.5, 0.5]);
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
        positions.push([
            cx + radius * theta.cos(),
            cy - half_h,
            cz + radius * theta.sin(),
        ]);
        normals.push([0.0, -1.0, 0.0]);
        uvs.push([0.5 + 0.5 * theta.cos(), 0.5 + 0.5 * theta.sin()]);
    }
    for i in 0..s {
        indices.push([bot_center, bot_center + i + 2, bot_center + i + 1]);
    }
    let top_center = positions.len();
    positions.push([cx, cy + half_h, cz]);
    normals.push([0.0, 1.0, 0.0]);
    uvs.push([0.5, 0.5]);
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
        positions.push([
            cx + radius * theta.cos(),
            cy + half_h,
            cz + radius * theta.sin(),
        ]);
        normals.push([0.0, 1.0, 0.0]);
        uvs.push([0.5 + 0.5 * theta.cos(), 0.5 + 0.5 * theta.sin()]);
    }
    for i in 0..s {
        indices.push([top_center, top_center + i + 1, top_center + i + 2]);
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate a cone as a `MeshData` with a flat base cap.
#[allow(dead_code)]
pub fn generate_cone_mesh(center: [f64; 3], radius: f64, height: f64, sectors: usize) -> MeshData {
    let s = sectors.max(3);
    let [cx, cy, cz] = center;
    let slant = (radius * radius + height * height).sqrt().max(1e-12);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
        let nx = theta.cos() * height / slant;
        let nz = theta.sin() * height / slant;
        let ny = radius / slant;
        positions.push([cx + radius * theta.cos(), cy, cz + radius * theta.sin()]);
        normals.push([nx, ny, nz]);
        uvs.push([i as f64 / s as f64, 0.0]);
    }
    let apex_idx = positions.len();
    positions.push([cx, cy + height, cz]);
    normals.push([0.0, 1.0, 0.0]);
    uvs.push([0.5, 1.0]);
    for i in 0..s {
        indices.push([i, i + 1, apex_idx]);
    }
    let base_center = positions.len();
    positions.push([cx, cy, cz]);
    normals.push([0.0, -1.0, 0.0]);
    uvs.push([0.5, 0.5]);
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
        positions.push([cx + radius * theta.cos(), cy, cz + radius * theta.sin()]);
        normals.push([0.0, -1.0, 0.0]);
        uvs.push([0.5 + 0.5 * theta.cos(), 0.5 + 0.5 * theta.sin()]);
    }
    for i in 0..s {
        indices.push([base_center, base_center + i + 2, base_center + i + 1]);
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate a torus as a `MeshData` with UV coordinates.
#[allow(dead_code)]
pub fn generate_torus_mesh(
    center: [f64; 3],
    major_r: f64,
    minor_r: f64,
    major_segs: usize,
    minor_segs: usize,
) -> MeshData {
    let ms = major_segs.max(3);
    let ns = minor_segs.max(3);
    let [cx, cy, cz] = center;
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=ms {
        let u = 2.0 * std::f64::consts::PI * i as f64 / ms as f64;
        let cu = u.cos();
        let su = u.sin();
        for j in 0..=ns {
            let v = 2.0 * std::f64::consts::PI * j as f64 / ns as f64;
            let cv = v.cos();
            let sv = v.sin();
            let x = (major_r + minor_r * cv) * cu;
            let y = minor_r * sv;
            let z = (major_r + minor_r * cv) * su;
            positions.push([cx + x, cy + y, cz + z]);
            normals.push([cv * cu, sv, cv * su]);
            uvs.push([i as f64 / ms as f64, j as f64 / ns as f64]);
        }
    }
    for i in 0..ms {
        for j in 0..ns {
            let a = i * (ns + 1) + j;
            let b = a + ns + 1;
            indices.push([a, b, b + 1]);
            indices.push([a, b + 1, a + 1]);
        }
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate a subdivided plane as a `MeshData` with UV coordinates.
#[allow(dead_code)]
pub fn generate_plane_mesh(
    center: [f64; 3],
    normal: [f64; 3],
    size: f64,
    subdivisions: usize,
) -> MeshData {
    let n_div = subdivisions.max(1);
    let mut n = normal;
    normalize3_f64(&mut n);
    let up = if n[1].abs() < 0.99 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut t = [
        n[1] * up[2] - n[2] * up[1],
        n[2] * up[0] - n[0] * up[2],
        n[0] * up[1] - n[1] * up[0],
    ];
    normalize3_f64(&mut t);
    let b = [
        n[1] * t[2] - n[2] * t[1],
        n[2] * t[0] - n[0] * t[2],
        n[0] * t[1] - n[1] * t[0],
    ];
    let half = size * 0.5;
    let step = size / n_div as f64;
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for iy in 0..=n_div {
        for ix in 0..=n_div {
            let u_off = -half + ix as f64 * step;
            let v_off = -half + iy as f64 * step;
            positions.push([
                center[0] + u_off * t[0] + v_off * b[0],
                center[1] + u_off * t[1] + v_off * b[1],
                center[2] + u_off * t[2] + v_off * b[2],
            ]);
            normals.push(n);
            uvs.push([ix as f64 / n_div as f64, iy as f64 / n_div as f64]);
        }
    }
    let w = n_div + 1;
    for iy in 0..n_div {
        for ix in 0..n_div {
            let a = iy * w + ix;
            indices.push([a, a + w, a + w + 1]);
            indices.push([a, a + w + 1, a + 1]);
        }
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
/// Generate an arrow for vector field visualization as a `MeshData`.
///
/// The arrow consists of a cylindrical shaft (75% of length) and a cone tip.
#[allow(dead_code)]
pub fn generate_arrow(
    origin: [f64; 3],
    direction: [f64; 3],
    shaft_radius: f64,
    tip_radius: f64,
    sectors: usize,
) -> MeshData {
    let s = sectors.max(3);
    let len =
        (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
            .sqrt();
    if len < 1e-12 {
        return MeshData::empty();
    }
    let axis = [direction[0] / len, direction[1] / len, direction[2] / len];
    let up = if axis[1].abs() < 0.99 {
        [0.0_f64, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut tangent = [
        axis[1] * up[2] - axis[2] * up[1],
        axis[2] * up[0] - axis[0] * up[2],
        axis[0] * up[1] - axis[1] * up[0],
    ];
    normalize3_f64(&mut tangent);
    let bitangent = [
        axis[1] * tangent[2] - axis[2] * tangent[1],
        axis[2] * tangent[0] - axis[0] * tangent[2],
        axis[0] * tangent[1] - axis[1] * tangent[0],
    ];
    let shaft_len = len * 0.75;
    let tip_len = len - shaft_len;
    let shaft_end = [
        origin[0] + axis[0] * shaft_len,
        origin[1] + axis[1] * shaft_len,
        origin[2] + axis[2] * shaft_len,
    ];
    let tip = [
        origin[0] + direction[0],
        origin[1] + direction[1],
        origin[2] + direction[2],
    ];
    let slant = (tip_radius * tip_radius + tip_len * tip_len)
        .sqrt()
        .max(1e-12);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut uvs = Vec::new();
    let mut indices = Vec::new();
    for (ring_center, v_uv) in [(&origin, 0.0_f64), (&shaft_end, 1.0_f64)] {
        for i in 0..=s {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
            let ct = theta.cos();
            let st = theta.sin();
            let nx = ct * tangent[0] + st * bitangent[0];
            let ny = ct * tangent[1] + st * bitangent[1];
            let nz = ct * tangent[2] + st * bitangent[2];
            positions.push([
                ring_center[0] + shaft_radius * (ct * tangent[0] + st * bitangent[0]),
                ring_center[1] + shaft_radius * (ct * tangent[1] + st * bitangent[1]),
                ring_center[2] + shaft_radius * (ct * tangent[2] + st * bitangent[2]),
            ]);
            normals.push([nx, ny, nz]);
            uvs.push([i as f64 / s as f64, v_uv]);
        }
    }
    for i in 0..s {
        let b0 = i;
        let b1 = i + 1;
        let t0 = i + s + 1;
        let t1 = i + s + 2;
        indices.push([b0, t0, t1]);
        indices.push([b0, t1, b1]);
    }
    let cone_base = positions.len();
    for i in 0..=s {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / s as f64;
        let ct = theta.cos();
        let st = theta.sin();
        let nt_x = ct * tangent[0] + st * bitangent[0];
        let nt_y = ct * tangent[1] + st * bitangent[1];
        let nt_z = ct * tangent[2] + st * bitangent[2];
        positions.push([
            shaft_end[0] + tip_radius * (ct * tangent[0] + st * bitangent[0]),
            shaft_end[1] + tip_radius * (ct * tangent[1] + st * bitangent[1]),
            shaft_end[2] + tip_radius * (ct * tangent[2] + st * bitangent[2]),
        ]);
        normals.push([
            nt_x * tip_len / slant + axis[0] * tip_radius / slant,
            nt_y * tip_len / slant + axis[1] * tip_radius / slant,
            nt_z * tip_len / slant + axis[2] * tip_radius / slant,
        ]);
        uvs.push([i as f64 / s as f64, 0.0]);
    }
    let apex_idx = positions.len();
    positions.push(tip);
    normals.push(axis);
    uvs.push([0.5, 1.0]);
    for i in 0..s {
        indices.push([cone_base + i, cone_base + i + 1, apex_idx]);
    }
    MeshData {
        positions,
        normals,
        uvs,
        indices,
    }
}
