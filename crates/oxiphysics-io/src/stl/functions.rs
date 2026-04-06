//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_strip)]
#![allow(clippy::items_after_test_module)]
use std::collections::HashMap;

use super::types::{
    EdgeKey, StlColor, StlMesh, StlStatistics, StlTriangle, StlValidation, StlValidationReport,
};

/// Compute the outward normal of a triangle using the right-hand rule.
pub(super) fn compute_normal(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> [f32; 3] {
    let a = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let b = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    normalize3_f32(cross3_f32(a, b))
}
/// Cross product of two 3-vectors.
pub(super) fn cross3_f32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Normalise a 3-vector; returns the zero vector if the magnitude is tiny.
pub(super) fn normalize3_f32(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-30 {
        [0.0; 3]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
/// Area of a single triangle.
pub(super) fn triangle_area(tri: &StlTriangle) -> f32 {
    let a = [
        tri.v1[0] - tri.v0[0],
        tri.v1[1] - tri.v0[1],
        tri.v1[2] - tri.v0[2],
    ];
    let b = [
        tri.v2[0] - tri.v0[0],
        tri.v2[1] - tri.v0[1],
        tri.v2[2] - tri.v0[2],
    ];
    let c = cross3_f32(a, b);
    0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}
/// Read 3 little-endian f32 values starting at `offset` in `data`.
pub(super) fn read_f32x3(data: &[u8], offset: usize) -> Result<[f32; 3], String> {
    if offset + 12 > data.len() {
        return Err(format!(
            "Buffer too short to read 3 f32 values at offset {}",
            offset
        ));
    }
    let x = f32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ]);
    let y = f32::from_le_bytes([
        data[offset + 4],
        data[offset + 5],
        data[offset + 6],
        data[offset + 7],
    ]);
    let z = f32::from_le_bytes([
        data[offset + 8],
        data[offset + 9],
        data[offset + 10],
        data[offset + 11],
    ]);
    Ok([x, y, z])
}
/// Parse three floats that follow a prefix on the same line.
pub(super) fn parse_vec3_from_line(line: &str, prefix: &str) -> Result<[f32; 3], String> {
    let rest = line
        .strip_prefix(prefix)
        .ok_or_else(|| format!("Expected prefix '{}' in line: {}", prefix, line))?
        .trim();
    let parts: Vec<&str> = rest.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(format!(
            "Expected 3 floats after '{}', got: {}",
            prefix, rest
        ));
    }
    Ok([
        parts[0].parse::<f32>().map_err(|e| e.to_string())?,
        parts[1].parse::<f32>().map_err(|e| e.to_string())?,
        parts[2].parse::<f32>().map_err(|e| e.to_string())?,
    ])
}
/// Parse a "vertex x y z" line.
pub(super) fn parse_vertex_line(line: &str) -> Result<[f32; 3], String> {
    parse_vec3_from_line(line, "vertex")
}
/// Validate an STL mesh and return a validation report.
#[allow(dead_code)]
pub fn validate_stl(mesh: &StlMesh) -> StlValidation {
    let mut degenerate_count = 0;
    let mut nan_inf_count = 0;
    let mut normals_valid = true;
    for tri in &mesh.triangles {
        let has_bad = [&tri.v0, &tri.v1, &tri.v2, &tri.normal]
            .iter()
            .any(|v| v.iter().any(|f| f.is_nan() || f.is_infinite()));
        if has_bad {
            nan_inf_count += 1;
            continue;
        }
        let area = triangle_area(tri);
        if area < 1e-10 {
            degenerate_count += 1;
        }
        let n = &tri.normal;
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if (len - 1.0).abs() > 0.01 && len > 1e-10 {
            normals_valid = false;
        }
    }
    let mut edge_counts: HashMap<EdgeKey, usize> = HashMap::new();
    for tri in &mesh.triangles {
        for &(a, b) in &[(tri.v0, tri.v1), (tri.v1, tri.v2), (tri.v2, tri.v0)] {
            let key = EdgeKey::new(a, b);
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }
    let non_manifold_edges = edge_counts.values().filter(|&&c| c != 2).count();
    let is_watertight = mesh.triangles.is_empty() || non_manifold_edges == 0;
    StlValidation {
        degenerate_count,
        non_manifold_edges,
        is_watertight,
        nan_inf_count,
        normals_valid,
    }
}
/// Recalculate all triangle normals from vertex positions.
#[allow(dead_code)]
pub fn recalculate_normals(mesh: &mut StlMesh) {
    for tri in &mut mesh.triangles {
        tri.normal = compute_normal(tri.v0, tri.v1, tri.v2);
    }
}
/// Flip all triangle normals (reverses winding order).
#[allow(dead_code)]
pub fn flip_normals(mesh: &mut StlMesh) {
    for tri in &mut mesh.triangles {
        std::mem::swap(&mut tri.v1, &mut tri.v2);
        tri.normal = [-tri.normal[0], -tri.normal[1], -tri.normal[2]];
    }
}
/// Remove degenerate (zero-area) triangles from the mesh.
#[allow(dead_code)]
pub fn remove_degenerate_triangles(mesh: &mut StlMesh) {
    mesh.triangles.retain(|tri| triangle_area(tri) >= 1e-10);
}
/// Remove duplicate triangles (exact vertex match).
#[allow(dead_code)]
pub fn remove_duplicate_triangles(mesh: &mut StlMesh) {
    let mut seen: std::collections::HashSet<[u32; 9]> = std::collections::HashSet::new();
    mesh.triangles.retain(|tri| {
        let key = [
            tri.v0[0].to_bits(),
            tri.v0[1].to_bits(),
            tri.v0[2].to_bits(),
            tri.v1[0].to_bits(),
            tri.v1[1].to_bits(),
            tri.v1[2].to_bits(),
            tri.v2[0].to_bits(),
            tri.v2[1].to_bits(),
            tri.v2[2].to_bits(),
        ];
        seen.insert(key)
    });
}
/// Merge two STL meshes into one.
#[allow(dead_code)]
pub fn merge_meshes(a: &StlMesh, b: &StlMesh) -> StlMesh {
    let mut result = StlMesh::new(&format!("{}_merged_{}", a.name, b.name));
    result
        .triangles
        .reserve(a.triangles.len() + b.triangles.len());
    for tri in &a.triangles {
        result.triangles.push(StlTriangle {
            normal: tri.normal,
            v0: tri.v0,
            v1: tri.v1,
            v2: tri.v2,
        });
    }
    for tri in &b.triangles {
        result.triangles.push(StlTriangle {
            normal: tri.normal,
            v0: tri.v0,
            v1: tri.v1,
            v2: tri.v2,
        });
    }
    result
}
/// Compute statistics for an STL mesh.
#[allow(dead_code)]
pub fn compute_statistics(mesh: &StlMesh) -> StlStatistics {
    let n = mesh.triangles.len();
    if n == 0 {
        return StlStatistics {
            triangle_count: 0,
            surface_area: 0.0,
            bb_min: [0.0; 3],
            bb_max: [0.0; 3],
            bb_size: [0.0; 3],
            avg_triangle_area: 0.0,
            min_triangle_area: 0.0,
            max_triangle_area: 0.0,
            avg_edge_length: 0.0,
            approx_unique_vertices: 0,
        };
    }
    let (bb_min, bb_max) = mesh.bounding_box();
    let bb_size = [
        bb_max[0] - bb_min[0],
        bb_max[1] - bb_min[1],
        bb_max[2] - bb_min[2],
    ];
    let mut total_area = 0.0_f32;
    let mut min_area = f32::MAX;
    let mut max_area = 0.0_f32;
    let mut total_edge_len = 0.0_f32;
    let mut vertex_set: std::collections::HashSet<[u32; 3]> = std::collections::HashSet::new();
    for tri in &mesh.triangles {
        let area = triangle_area(tri);
        total_area += area;
        if area < min_area {
            min_area = area;
        }
        if area > max_area {
            max_area = area;
        }
        for &(a, b) in &[(tri.v0, tri.v1), (tri.v1, tri.v2), (tri.v2, tri.v0)] {
            let dx = b[0] - a[0];
            let dy = b[1] - a[1];
            let dz = b[2] - a[2];
            total_edge_len += (dx * dx + dy * dy + dz * dz).sqrt();
        }
        for v in [&tri.v0, &tri.v1, &tri.v2] {
            vertex_set.insert([v[0].to_bits(), v[1].to_bits(), v[2].to_bits()]);
        }
    }
    StlStatistics {
        triangle_count: n,
        surface_area: total_area,
        bb_min,
        bb_max,
        bb_size,
        avg_triangle_area: total_area / n as f32,
        min_triangle_area: min_area,
        max_triangle_area: max_area,
        avg_edge_length: total_edge_len / (n as f32 * 3.0),
        approx_unique_vertices: vertex_set.len(),
    }
}
/// Translate all vertices by the given offset.
#[allow(dead_code)]
pub fn translate(mesh: &mut StlMesh, offset: [f32; 3]) {
    for tri in &mut mesh.triangles {
        for v in [&mut tri.v0, &mut tri.v1, &mut tri.v2] {
            v[0] += offset[0];
            v[1] += offset[1];
            v[2] += offset[2];
        }
    }
}
/// Scale all vertices uniformly from the origin.
#[allow(dead_code)]
pub fn scale_uniform(mesh: &mut StlMesh, factor: f32) {
    for tri in &mut mesh.triangles {
        for v in [&mut tri.v0, &mut tri.v1, &mut tri.v2] {
            v[0] *= factor;
            v[1] *= factor;
            v[2] *= factor;
        }
    }
}
/// Scale each axis independently.
#[allow(dead_code)]
pub fn scale_nonuniform(mesh: &mut StlMesh, factors: [f32; 3]) {
    for tri in &mut mesh.triangles {
        for v in [&mut tri.v0, &mut tri.v1, &mut tri.v2] {
            v[0] *= factors[0];
            v[1] *= factors[1];
            v[2] *= factors[2];
        }
    }
    recalculate_normals(mesh);
}
/// Center the mesh at the origin (translate so bounding box center is at 0,0,0).
#[allow(dead_code)]
pub fn center_at_origin(mesh: &mut StlMesh) {
    let (lo, hi) = mesh.bounding_box();
    let offset = [
        -(lo[0] + hi[0]) / 2.0,
        -(lo[1] + hi[1]) / 2.0,
        -(lo[2] + hi[2]) / 2.0,
    ];
    translate(mesh, offset);
}
/// Rotate mesh about the Z axis by `angle_deg` degrees.
#[allow(dead_code)]
pub fn rotate_z(mesh: &mut StlMesh, angle_deg: f32) {
    let theta = angle_deg.to_radians();
    let c = theta.cos();
    let s = theta.sin();
    for tri in &mut mesh.triangles {
        for v in [&mut tri.v0, &mut tri.v1, &mut tri.v2] {
            let x = v[0];
            let y = v[1];
            v[0] = c * x - s * y;
            v[1] = s * x + c * y;
        }
    }
    recalculate_normals(mesh);
}
#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> ([f32; 3], [f32; 3], [f32; 3]) {
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    }
    fn make_cube() -> StlMesh {
        let mut m = StlMesh::new("cube");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        m.add_triangle([1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        m.add_triangle([0.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]);
        m.add_triangle([1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]);
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        m.add_triangle([1.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        m.add_triangle([0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]);
        m.add_triangle([1.0, 1.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0]);
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        m.add_triangle([0.0, 1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
        m.add_triangle([1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]);
        m.add_triangle([1.0, 1.0, 1.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]);
        m
    }
    #[test]
    fn test_new_empty() {
        let m = StlMesh::new("test");
        assert_eq!(m.name, "test");
        assert!(m.triangles.is_empty());
    }
    #[test]
    fn test_add_triangle_count() {
        let mut m = StlMesh::new("t");
        let (v0, v1, v2) = unit_triangle();
        m.add_triangle(v0, v1, v2);
        assert_eq!(m.triangles.len(), 1);
    }
    #[test]
    fn test_normal_is_unit() {
        let (v0, v1, v2) = unit_triangle();
        let n = compute_normal(v0, v1, v2);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-5, "len={len}");
    }
    #[test]
    fn test_normal_z_up() {
        let n = compute_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((n[2] - 1.0).abs() < 1e-5, "expected +z normal, got {n:?}");
    }
    #[test]
    fn test_cross3() {
        let x = [1.0_f32, 0.0, 0.0];
        let y = [0.0_f32, 1.0, 0.0];
        let z = cross3_f32(x, y);
        assert!((z[2] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_normalize_unit() {
        let v = [3.0_f32, 0.0, 0.0];
        let n = normalize3_f32(v);
        assert!((n[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_normalize_zero_returns_zero() {
        let n = normalize3_f32([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_binary_length() {
        let mut m = StlMesh::new("x");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let bytes = m.to_binary_bytes();
        assert_eq!(bytes.len(), 84 + 50);
    }
    #[test]
    fn test_binary_roundtrip_name() {
        let mut m = StlMesh::new("mysolid");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let bytes = m.to_binary_bytes();
        let m2 = StlMesh::from_binary_bytes(&bytes).unwrap();
        assert_eq!(m2.name, "mysolid");
    }
    #[test]
    fn test_binary_roundtrip_triangle_count() {
        let cube = make_cube();
        let bytes = cube.to_binary_bytes();
        let m2 = StlMesh::from_binary_bytes(&bytes).unwrap();
        assert_eq!(m2.triangles.len(), 12);
    }
    #[test]
    fn test_binary_roundtrip_vertices() {
        let mut m = StlMesh::new("s");
        m.add_triangle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]);
        let bytes = m.to_binary_bytes();
        let m2 = StlMesh::from_binary_bytes(&bytes).unwrap();
        let tri = &m2.triangles[0];
        assert_eq!(tri.v0, [1.0, 2.0, 3.0]);
        assert_eq!(tri.v1, [4.0, 5.0, 6.0]);
        assert_eq!(tri.v2, [7.0, 8.0, 9.0]);
    }
    #[test]
    fn test_binary_too_short_error() {
        let result = StlMesh::from_binary_bytes(&[0u8; 10]);
        assert!(result.is_err());
    }
    #[test]
    fn test_binary_count_mismatch_error() {
        let mut bytes = vec![0u8; 84];
        bytes[80..84].copy_from_slice(&5u32.to_le_bytes());
        let result = StlMesh::from_binary_bytes(&bytes);
        assert!(result.is_err());
    }
    #[test]
    fn test_ascii_contains_solid_name() {
        let m = StlMesh::new("sphere");
        let s = m.to_ascii();
        assert!(s.contains("solid sphere"));
        assert!(s.contains("endsolid sphere"));
    }
    #[test]
    fn test_ascii_roundtrip_single_triangle() {
        let mut m = StlMesh::new("tri");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let s = m.to_ascii();
        let m2 = StlMesh::from_ascii(&s).unwrap();
        assert_eq!(m2.triangles.len(), 1);
    }
    #[test]
    fn test_ascii_roundtrip_vertices() {
        let mut m = StlMesh::new("t");
        m.add_triangle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]);
        let s = m.to_ascii();
        let m2 = StlMesh::from_ascii(&s).unwrap();
        let t = &m2.triangles[0];
        assert!((t.v0[0] - 1.0).abs() < 1e-4);
        assert!((t.v2[2] - 9.0).abs() < 1e-4);
    }
    #[test]
    fn test_ascii_roundtrip_name() {
        let m = StlMesh::new("named_solid");
        let s = m.to_ascii();
        let m2 = StlMesh::from_ascii(&s).unwrap();
        assert_eq!(m2.name, "named_solid");
    }
    #[test]
    fn test_ascii_bad_header_error() {
        let result = StlMesh::from_ascii("not_solid\n");
        assert!(result.is_err());
    }
    #[test]
    fn test_surface_area_unit_triangle() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let area = m.surface_area();
        assert!((area - 0.5).abs() < 1e-5, "area={area}");
    }
    #[test]
    fn test_surface_area_cube() {
        let cube = make_cube();
        let area = cube.surface_area();
        assert!((area - 6.0).abs() < 1e-4, "cube area={area}");
    }
    #[test]
    fn test_surface_area_empty() {
        let m = StlMesh::new("empty");
        assert_eq!(m.surface_area(), 0.0);
    }
    #[test]
    fn test_bounding_box_single_triangle() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let (lo, hi) = m.bounding_box();
        assert_eq!(lo, [0.0, 0.0, 0.0]);
        assert_eq!(hi, [1.0, 1.0, 0.0]);
    }
    #[test]
    fn test_bounding_box_empty() {
        let m = StlMesh::new("empty");
        let (lo, hi) = m.bounding_box();
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [0.0; 3]);
    }
    #[test]
    fn test_bounding_box_cube() {
        let cube = make_cube();
        let (lo, hi) = cube.bounding_box();
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [1.0; 3]);
    }
    #[test]
    fn test_watertight_cube() {
        let cube = make_cube();
        assert!(cube.is_watertight(), "unit cube should be watertight");
    }
    #[test]
    fn test_not_watertight_open() {
        let mut m = StlMesh::new("open");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        assert!(!m.is_watertight());
    }
    #[test]
    fn test_watertight_empty() {
        let m = StlMesh::new("empty");
        assert!(m.is_watertight());
    }
    #[test]
    fn test_validate_cube() {
        let cube = make_cube();
        let v = validate_stl(&cube);
        assert!(v.is_watertight);
        assert_eq!(v.degenerate_count, 0);
        assert_eq!(v.nan_inf_count, 0);
        assert!(v.normals_valid);
        assert_eq!(v.non_manifold_edges, 0);
    }
    #[test]
    fn test_validate_degenerate() {
        let mut m = StlMesh::new("degen");
        m.add_triangle([0., 0., 0.], [0., 0., 0.], [0., 0., 0.]);
        let v = validate_stl(&m);
        assert_eq!(v.degenerate_count, 1);
    }
    #[test]
    fn test_validate_empty() {
        let m = StlMesh::new("empty");
        let v = validate_stl(&m);
        assert!(v.is_watertight);
        assert_eq!(v.degenerate_count, 0);
    }
    #[test]
    fn test_validate_open_mesh() {
        let mut m = StlMesh::new("open");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let v = validate_stl(&m);
        assert!(!v.is_watertight);
        assert!(v.non_manifold_edges > 0);
    }
    #[test]
    fn test_recalculate_normals() {
        let mut m = StlMesh::new("t");
        m.triangles.push(StlTriangle {
            normal: [0.0, 0.0, 0.0],
            v0: [0., 0., 0.],
            v1: [1., 0., 0.],
            v2: [0., 1., 0.],
        });
        recalculate_normals(&mut m);
        let n = m.triangles[0].normal;
        assert!((n[2] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_flip_normals() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let orig_n = m.triangles[0].normal;
        flip_normals(&mut m);
        let flipped_n = m.triangles[0].normal;
        assert!((flipped_n[2] + orig_n[2]).abs() < 1e-5);
    }
    #[test]
    fn test_remove_degenerate_triangles() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        m.add_triangle([0., 0., 0.], [0., 0., 0.], [0., 0., 0.]);
        assert_eq!(m.triangles.len(), 2);
        remove_degenerate_triangles(&mut m);
        assert_eq!(m.triangles.len(), 1);
    }
    #[test]
    fn test_remove_duplicate_triangles() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        assert_eq!(m.triangles.len(), 2);
        remove_duplicate_triangles(&mut m);
        assert_eq!(m.triangles.len(), 1);
    }
    #[test]
    fn test_remove_duplicate_keeps_unique() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        m.add_triangle([1., 1., 1.], [2., 1., 1.], [1., 2., 1.]);
        remove_duplicate_triangles(&mut m);
        assert_eq!(m.triangles.len(), 2);
    }
    #[test]
    fn test_merge_meshes() {
        let mut a = StlMesh::new("a");
        a.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let mut b = StlMesh::new("b");
        b.add_triangle([1., 1., 1.], [2., 1., 1.], [1., 2., 1.]);
        let merged = merge_meshes(&a, &b);
        assert_eq!(merged.triangles.len(), 2);
        assert!(merged.name.contains("merged"));
    }
    #[test]
    fn test_merge_empty() {
        let a = StlMesh::new("a");
        let b = StlMesh::new("b");
        let merged = merge_meshes(&a, &b);
        assert!(merged.triangles.is_empty());
    }
    #[test]
    fn test_statistics_cube() {
        let cube = make_cube();
        let stats = compute_statistics(&cube);
        assert_eq!(stats.triangle_count, 12);
        assert!((stats.surface_area - 6.0).abs() < 0.01);
        assert_eq!(stats.bb_min, [0.0; 3]);
        assert_eq!(stats.bb_max, [1.0; 3]);
        assert!((stats.bb_size[0] - 1.0).abs() < 1e-5);
        assert!(stats.avg_triangle_area > 0.0);
        assert!(stats.avg_edge_length > 0.0);
        assert!(stats.approx_unique_vertices <= 8);
    }
    #[test]
    fn test_statistics_empty() {
        let m = StlMesh::new("empty");
        let stats = compute_statistics(&m);
        assert_eq!(stats.triangle_count, 0);
        assert_eq!(stats.surface_area, 0.0);
        assert_eq!(stats.approx_unique_vertices, 0);
    }
    #[test]
    fn test_statistics_single_triangle() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let stats = compute_statistics(&m);
        assert_eq!(stats.triangle_count, 1);
        assert!((stats.surface_area - 0.5).abs() < 1e-5);
        assert_eq!(stats.approx_unique_vertices, 3);
    }
    #[test]
    fn test_translate() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        translate(&mut m, [10.0, 20.0, 30.0]);
        assert!((m.triangles[0].v0[0] - 10.0).abs() < 1e-5);
        assert!((m.triangles[0].v0[1] - 20.0).abs() < 1e-5);
        assert!((m.triangles[0].v0[2] - 30.0).abs() < 1e-5);
    }
    #[test]
    fn test_scale_uniform() {
        let mut m = StlMesh::new("t");
        m.add_triangle([1., 0., 0.], [2., 0., 0.], [1., 1., 0.]);
        scale_uniform(&mut m, 2.0);
        assert!((m.triangles[0].v0[0] - 2.0).abs() < 1e-5);
        assert!((m.triangles[0].v1[0] - 4.0).abs() < 1e-5);
    }
    #[test]
    fn test_scale_nonuniform() {
        let mut m = StlMesh::new("t");
        m.add_triangle([1., 1., 1.], [2., 1., 1.], [1., 2., 1.]);
        scale_nonuniform(&mut m, [2.0, 3.0, 4.0]);
        assert!((m.triangles[0].v0[0] - 2.0).abs() < 1e-5);
        assert!((m.triangles[0].v0[1] - 3.0).abs() < 1e-5);
        assert!((m.triangles[0].v0[2] - 4.0).abs() < 1e-5);
    }
    #[test]
    fn test_center_at_origin() {
        let mut m = StlMesh::new("t");
        m.add_triangle([2., 2., 2.], [4., 2., 2.], [2., 4., 2.]);
        center_at_origin(&mut m);
        let (lo, hi) = m.bounding_box();
        let cx = (lo[0] + hi[0]) / 2.0;
        let cy = (lo[1] + hi[1]) / 2.0;
        assert!(cx.abs() < 1e-4, "cx={cx}");
        assert!(cy.abs() < 1e-4, "cy={cy}");
    }
    #[test]
    fn test_rotate_z() {
        let mut m = StlMesh::new("t");
        m.add_triangle([1., 0., 0.], [2., 0., 0.], [1., 1., 0.]);
        rotate_z(&mut m, 90.0);
        assert!((m.triangles[0].v0[0]).abs() < 0.01, "x should be ~0");
        assert!((m.triangles[0].v0[1] - 1.0).abs() < 0.01, "y should be ~1");
    }
    #[test]
    fn test_translate_preserves_area() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let area_before = m.surface_area();
        translate(&mut m, [100.0, 200.0, 300.0]);
        let area_after = m.surface_area();
        assert!((area_before - area_after).abs() < 1e-5);
    }
    #[test]
    fn test_scale_uniform_area() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0., 0., 0.], [1., 0., 0.], [0., 1., 0.]);
        let area_before = m.surface_area();
        scale_uniform(&mut m, 3.0);
        let area_after = m.surface_area();
        assert!((area_after - area_before * 9.0).abs() < 0.01);
    }
}
/// Validate an STL mesh for watertightness and manifoldness.
///
/// An edge is represented as a sorted pair of vertex indices (by quantised
/// coordinates).  In a valid closed manifold every edge appears exactly twice.
pub fn validate_mesh(mesh: &StlMesh) -> StlValidationReport {
    let quantise = |v: [f32; 3]| -> [i64; 3] {
        [
            (v[0] * 1_000_000.0) as i64,
            (v[1] * 1_000_000.0) as i64,
            (v[2] * 1_000_000.0) as i64,
        ]
    };
    let mut edge_count: HashMap<([i64; 3], [i64; 3]), usize> = HashMap::new();
    let mut degenerate_count = 0usize;
    for tri in &mesh.triangles {
        let e1 = [
            tri.v1[0] - tri.v0[0],
            tri.v1[1] - tri.v0[1],
            tri.v1[2] - tri.v0[2],
        ];
        let e2 = [
            tri.v2[0] - tri.v0[0],
            tri.v2[1] - tri.v0[1],
            tri.v2[2] - tri.v0[2],
        ];
        let cross = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let area2 = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        if area2 < 1e-12 {
            degenerate_count += 1;
        }
        let va = quantise(tri.v0);
        let vb = quantise(tri.v1);
        let vc = quantise(tri.v2);
        for (a, b) in [(va, vb), (vb, vc), (vc, va)] {
            let key = if a <= b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let boundary_edge_count = edge_count.values().filter(|&&c| c != 2).count();
    let is_watertight = boundary_edge_count == 0;
    let is_manifold = edge_count.values().all(|&c| c == 2);
    StlValidationReport {
        boundary_edge_count,
        is_watertight,
        is_manifold,
        degenerate_count,
    }
}
/// Serialize `mesh` to binary STL bytes with per-triangle color in the
/// attribute field (VisCAM format).  `colors` must be the same length as
/// `mesh.triangles`; if shorter, remaining triangles get attribute = 0.
pub fn to_binary_bytes_with_color(mesh: &StlMesh, colors: &[StlColor]) -> Vec<u8> {
    let n = mesh.triangles.len();
    let mut buf = Vec::with_capacity(84 + n * 50);
    let mut header = [0u8; 80];
    let nb = mesh.name.as_bytes();
    let cl = nb.len().min(80);
    header[..cl].copy_from_slice(&nb[..cl]);
    buf.extend_from_slice(&header);
    buf.extend_from_slice(&(n as u32).to_le_bytes());
    for (idx, tri) in mesh.triangles.iter().enumerate() {
        for &f in &tri.normal {
            buf.extend_from_slice(&f.to_le_bytes());
        }
        for &f in &tri.v0 {
            buf.extend_from_slice(&f.to_le_bytes());
        }
        for &f in &tri.v1 {
            buf.extend_from_slice(&f.to_le_bytes());
        }
        for &f in &tri.v2 {
            buf.extend_from_slice(&f.to_le_bytes());
        }
        let attr: u16 = colors.get(idx).map(|c| c.encode()).unwrap_or(0);
        buf.extend_from_slice(&attr.to_le_bytes());
    }
    buf
}
/// Merge two STL meshes into a new mesh with the given name.
///
/// All triangles from both meshes are copied; no deduplication is performed.
pub fn merge_meshes_named(a: &StlMesh, b: &StlMesh, name: &str) -> StlMesh {
    let mut out = StlMesh::new(name);
    for tri in &a.triangles {
        out.triangles.push(StlTriangle {
            normal: tri.normal,
            v0: tri.v0,
            v1: tri.v1,
            v2: tri.v2,
        });
    }
    for tri in &b.triangles {
        out.triangles.push(StlTriangle {
            normal: tri.normal,
            v0: tri.v0,
            v1: tri.v1,
            v2: tri.v2,
        });
    }
    out
}
/// Recompute all triangle normals from vertex positions using the right-hand rule.
pub fn fix_normals(mesh: &mut StlMesh) {
    for tri in mesh.triangles.iter_mut() {
        tri.normal = compute_normal(tri.v0, tri.v1, tri.v2);
    }
}
/// Flip the winding order of all triangles in `mesh`.
///
/// This reverses the orientation of every triangle (useful when all normals
/// are pointing inward and need to be flipped outward).
pub fn flip_winding(mesh: &mut StlMesh) {
    for tri in mesh.triangles.iter_mut() {
        std::mem::swap(&mut tri.v1, &mut tri.v2);
        tri.normal = compute_normal(tri.v0, tri.v1, tri.v2);
    }
}
/// Attempt to close small holes in `mesh` by inserting degenerate "cap"
/// triangles over boundary edges.
///
/// **Important**: this is a heuristic repair that works only for meshes with
/// a small number of boundary edges.  The cap triangles use the midpoint of
/// each boundary edge as the third vertex.
pub fn repair_close_holes(mesh: &mut StlMesh) {
    let quantise = |v: [f32; 3]| -> [i64; 3] {
        [
            (v[0] * 1_000_000.0) as i64,
            (v[1] * 1_000_000.0) as i64,
            (v[2] * 1_000_000.0) as i64,
        ]
    };
    let dequantise = |k: [i64; 3]| -> [f32; 3] {
        [
            k[0] as f32 / 1_000_000.0,
            k[1] as f32 / 1_000_000.0,
            k[2] as f32 / 1_000_000.0,
        ]
    };
    let mut edge_count: HashMap<([i64; 3], [i64; 3]), usize> = HashMap::new();
    for tri in &mesh.triangles {
        let va = quantise(tri.v0);
        let vb = quantise(tri.v1);
        let vc = quantise(tri.v2);
        for (a, b) in [(va, vb), (vb, vc), (vc, va)] {
            let key = if a <= b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let mut new_tris: Vec<StlTriangle> = Vec::new();
    for ((ka, kb), &cnt) in &edge_count {
        if cnt != 1 {
            continue;
        }
        let va = dequantise(*ka);
        let vb = dequantise(*kb);
        let mid = [
            (va[0] + vb[0]) / 2.0,
            (va[1] + vb[1]) / 2.0,
            (va[2] + vb[2]) / 2.0,
        ];
        let normal = compute_normal(va, vb, mid);
        new_tris.push(StlTriangle {
            normal,
            v0: va,
            v1: vb,
            v2: mid,
        });
    }
    mesh.triangles.extend(new_tris);
}
/// Compute the axis-aligned bounding box of a standalone triangle set.
///
/// Returns `([min_x, min_y, min_z], [max_x, max_y, max_z])` or `None` if
/// `triangles` is empty.
pub fn triangle_bounding_box(triangles: &[StlTriangle]) -> Option<([f32; 3], [f32; 3])> {
    if triangles.is_empty() {
        return None;
    }
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for tri in triangles {
        for v in [tri.v0, tri.v1, tri.v2] {
            for d in 0..3 {
                if v[d] < lo[d] {
                    lo[d] = v[d];
                }
                if v[d] > hi[d] {
                    hi[d] = v[d];
                }
            }
        }
    }
    Some((lo, hi))
}
/// Convert an STL mesh to a Wavefront OBJ string.
///
/// Since STL has no shared-vertex semantics, each triangle is emitted as
/// three independent vertices.  Groups are not emitted; a single `g default`
/// line is written.
pub fn stl_to_obj(mesh: &StlMesh) -> String {
    let mut out = String::new();
    out.push_str("# Converted from STL by OxiPhysics\n");
    out.push_str(&format!("o {}\n", mesh.name));
    out.push_str("g default\n");
    for tri in &mesh.triangles {
        for v in [tri.v0, tri.v1, tri.v2] {
            out.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
        }
    }
    for i in 0..mesh.triangles.len() {
        let base = i * 3 + 1;
        out.push_str(&format!("f {} {} {}\n", base, base + 1, base + 2));
    }
    out
}
/// Parse a minimal Wavefront OBJ string back into an STL mesh.
///
/// Only `v` and `f` records are handled.  The mesh name is taken from the
/// `o` record if present.
pub fn obj_to_stl(obj: &str) -> StlMesh {
    let mut vertices: Vec<[f32; 3]> = Vec::new();
    let mut name = "imported".to_string();
    let mut mesh = StlMesh::new(&name);
    for line in obj.lines() {
        let line = line.trim();
        if line.starts_with("o ") {
            name = line[2..].trim().to_string();
            mesh.name = name.clone();
        } else if line.starts_with("v ") && !line.starts_with("vt") && !line.starts_with("vn") {
            let parts: Vec<f32> = line[2..]
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if parts.len() >= 3 {
                vertices.push([parts[0], parts[1], parts[2]]);
            }
        } else if line.starts_with("f ") {
            let idx: Vec<usize> = line[2..]
                .split_whitespace()
                .filter_map(|tok| tok.split('/').next()?.parse::<usize>().ok())
                .collect();
            if idx.len() >= 3 {
                let v0 = vertices.get(idx[0] - 1).copied().unwrap_or([0.0; 3]);
                let v1 = vertices.get(idx[1] - 1).copied().unwrap_or([0.0; 3]);
                let v2 = vertices.get(idx[2] - 1).copied().unwrap_or([0.0; 3]);
                mesh.add_triangle(v0, v1, v2);
            }
        }
    }
    mesh
}
