//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::items_after_test_module)]

use super::types::{
    ElementTypeStats, ExtendedQualityReport, FullQualityReport, MeshComparison, MeshQualityReport,
    QualityCheckResult, QualitySuggestion, QualityThresholds, TriangleMesh,
};

#[inline]
pub(super) fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn norm(a: [f32; 3]) -> f32 {
    dot(a, a).sqrt()
}
/// Compute the area of a triangle given its three vertex positions.
#[allow(dead_code)]
pub fn triangle_area(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let ab = sub(v1, v0);
    let ac = sub(v2, v0);
    norm(cross(ab, ac)) * 0.5
}
/// Compute the aspect ratio of a triangle: longest edge / shortest altitude.
/// Returns 1.0 for an equilateral triangle.
#[allow(dead_code)]
pub fn triangle_aspect_ratio(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let e0 = norm(sub(v1, v0));
    let e1 = norm(sub(v2, v1));
    let e2 = norm(sub(v0, v2));
    let area = triangle_area(v0, v1, v2);
    if area < 1e-15 {
        return f32::INFINITY;
    }
    let longest = e0.max(e1).max(e2);
    let shortest_altitude = 2.0 * area / longest;
    longest / shortest_altitude
}
/// Compute the skewness of a triangle based on the maximum angle.
///
/// Skewness = (theta_max - theta_equi) / (180 - theta_equi)
/// where theta_equi = 60 degrees. Returns 0 for equilateral, 1 for degenerate.
#[allow(dead_code)]
pub fn triangle_skewness(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let theta_max = triangle_max_angle_deg(v0, v1, v2);
    let theta_equi = 60.0_f32;
    (theta_max - theta_equi) / (180.0 - theta_equi)
}
/// Compute the minimum interior angle of a triangle in degrees.
#[allow(dead_code)]
pub fn triangle_min_angle_deg(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let angles = triangle_angles_deg(v0, v1, v2);
    angles[0].min(angles[1]).min(angles[2])
}
/// Compute the maximum interior angle of a triangle in degrees.
#[allow(dead_code)]
pub fn triangle_max_angle_deg(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let angles = triangle_angles_deg(v0, v1, v2);
    angles[0].max(angles[1]).max(angles[2])
}
/// Compute the inscribed circle radius: r = area / s (s = semi-perimeter).
#[allow(dead_code)]
pub fn triangle_inscribed_radius(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let a = norm(sub(v2, v1));
    let b = norm(sub(v0, v2));
    let c = norm(sub(v1, v0));
    let s = (a + b + c) * 0.5;
    if s < 1e-15 {
        return 0.0;
    }
    triangle_area(v0, v1, v2) / s
}
/// Compute the circumscribed circle radius: R = a*b*c / (4*area).
#[allow(dead_code)]
pub fn triangle_circumscribed_radius(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
    let a = norm(sub(v2, v1));
    let b = norm(sub(v0, v2));
    let c = norm(sub(v1, v0));
    let area = triangle_area(v0, v1, v2);
    if area < 1e-15 {
        return f32::INFINITY;
    }
    (a * b * c) / (4.0 * area)
}
pub(super) fn triangle_angles_deg(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> [f32; 3] {
    let a = norm(sub(v2, v1));
    let b = norm(sub(v0, v2));
    let c = norm(sub(v1, v0));
    let cos_alpha = if b * c > 1e-15 {
        ((b * b + c * c - a * a) / (2.0 * b * c)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let cos_beta = if a * c > 1e-15 {
        ((a * a + c * c - b * b) / (2.0 * a * c)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let cos_gamma = if a * b > 1e-15 {
        ((a * a + b * b - c * c) / (2.0 * a * b)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let rad_to_deg = 180.0 / std::f32::consts::PI;
    [
        cos_alpha.acos() * rad_to_deg,
        cos_beta.acos() * rad_to_deg,
        cos_gamma.acos() * rad_to_deg,
    ]
}
/// Compute a quality report for the given triangle mesh.
#[allow(dead_code)]
pub fn compute_quality_report(mesh: &TriangleMesh) -> MeshQualityReport {
    if mesh.triangles.is_empty() {
        return MeshQualityReport {
            min_area: 0.0,
            max_area: 0.0,
            mean_area: 0.0,
            min_aspect_ratio: 0.0,
            max_aspect_ratio: 0.0,
            mean_skewness: 0.0,
            n_degenerate: 0,
        };
    }
    let n = mesh.triangles.len();
    let mut min_area = f32::INFINITY;
    let mut max_area = f32::NEG_INFINITY;
    let mut sum_area = 0.0_f32;
    let mut min_ar = f32::INFINITY;
    let mut max_ar = f32::NEG_INFINITY;
    let mut sum_skew = 0.0_f32;
    let mut n_degenerate = 0_u32;
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let area = triangle_area(v0, v1, v2);
        if area < 1e-10 {
            n_degenerate += 1;
        }
        min_area = min_area.min(area);
        max_area = max_area.max(area);
        sum_area += area;
        let ar = triangle_aspect_ratio(v0, v1, v2);
        min_ar = min_ar.min(ar);
        max_ar = max_ar.max(ar);
        sum_skew += triangle_skewness(v0, v1, v2);
    }
    MeshQualityReport {
        min_area,
        max_area,
        mean_area: sum_area / n as f32,
        min_aspect_ratio: min_ar,
        max_aspect_ratio: max_ar,
        mean_skewness: sum_skew / n as f32,
        n_degenerate,
    }
}
/// Compute the axis-aligned bounding box of a mesh.
/// Returns (min_corner, max_corner).
#[allow(dead_code)]
pub fn mesh_bounding_box(mesh: &TriangleMesh) -> ([f32; 3], [f32; 3]) {
    if mesh.vertices.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = [f32::INFINITY; 3];
    let mut mx = [f32::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for i in 0..3 {
            if v[i] < mn[i] {
                mn[i] = v[i];
            }
            if v[i] > mx[i] {
                mx[i] = v[i];
            }
        }
    }
    (mn, mx)
}
/// Compute the total surface area of the mesh.
#[allow(dead_code)]
pub fn mesh_surface_area(mesh: &TriangleMesh) -> f32 {
    mesh.triangles
        .iter()
        .map(|tri| {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            triangle_area(v0, v1, v2)
        })
        .sum()
}
/// Compute the average edge length of the mesh.
#[allow(dead_code)]
pub fn mesh_average_edge_length(mesh: &TriangleMesh) -> f32 {
    if mesh.triangles.is_empty() {
        return 0.0;
    }
    let total: f32 = mesh
        .triangles
        .iter()
        .map(|tri| {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            norm(sub(v1, v0)) + norm(sub(v2, v1)) + norm(sub(v0, v2))
        })
        .sum();
    total / (mesh.triangles.len() as f32 * 3.0)
}
/// Compute the signed volume of a closed mesh using the divergence theorem.
#[allow(dead_code)]
pub fn mesh_volume_signed(mesh: &TriangleMesh) -> f32 {
    mesh.triangles
        .iter()
        .map(|tri| {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            let c = cross(v1, v2);
            dot(v0, c) / 6.0
        })
        .sum()
}
/// Compute the number of adjacent triangles per vertex (vertex valence).
#[allow(dead_code)]
pub fn vertex_valence(mesh: &TriangleMesh) -> Vec<u32> {
    let mut valence = vec![0_u32; mesh.vertices.len()];
    for tri in &mesh.triangles {
        valence[tri[0] as usize] += 1;
        valence[tri[1] as usize] += 1;
        valence[tri[2] as usize] += 1;
    }
    valence
}
/// Remove duplicate vertices within the given tolerance, remapping triangle indices.
#[allow(dead_code)]
pub fn remove_duplicate_vertices(mesh: &TriangleMesh, tolerance: f32) -> TriangleMesh {
    let mut new_verts: Vec<[f32; 3]> = Vec::new();
    let mut remap: Vec<u32> = Vec::with_capacity(mesh.vertices.len());
    for v in &mesh.vertices {
        let found = new_verts.iter().enumerate().find(|(_, u)| {
            let dx = v[0] - u[0];
            let dy = v[1] - u[1];
            let dz = v[2] - u[2];
            (dx * dx + dy * dy + dz * dz).sqrt() <= tolerance
        });
        if let Some((idx, _)) = found {
            remap.push(idx as u32);
        } else {
            remap.push(new_verts.len() as u32);
            new_verts.push(*v);
        }
    }
    let new_tris: Vec<[u32; 3]> = mesh
        .triangles
        .iter()
        .map(|tri| {
            [
                remap[tri[0] as usize],
                remap[tri[1] as usize],
                remap[tri[2] as usize],
            ]
        })
        .collect();
    TriangleMesh::from_raw(new_verts, new_tris)
}
/// Remove triangles with area < 1e-10 (degenerate triangles).
#[allow(dead_code)]
pub fn remove_degenerate_triangles(mesh: &TriangleMesh) -> TriangleMesh {
    let tris: Vec<[u32; 3]> = mesh
        .triangles
        .iter()
        .filter(|tri| {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            triangle_area(v0, v1, v2) >= 1e-10
        })
        .copied()
        .collect();
    TriangleMesh::from_raw(mesh.vertices.clone(), tris)
}
/// Flip the normals of a mesh by reversing the winding order of all triangles.
#[allow(dead_code)]
pub fn flip_normals(mesh: &TriangleMesh) -> TriangleMesh {
    let tris: Vec<[u32; 3]> = mesh
        .triangles
        .iter()
        .map(|tri| [tri[0], tri[2], tri[1]])
        .collect();
    TriangleMesh::from_raw(mesh.vertices.clone(), tris)
}
/// Check a mesh against quality thresholds.
#[allow(dead_code)]
pub fn check_quality(mesh: &TriangleMesh, thresholds: &QualityThresholds) -> QualityCheckResult {
    let mut n_bad_ar = 0_u32;
    let mut n_bad_skew = 0_u32;
    let mut n_small_angle = 0_u32;
    let mut n_large_angle = 0_u32;
    let mut n_below_area = 0_u32;
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let area = triangle_area(v0, v1, v2);
        if area < thresholds.min_area {
            n_below_area += 1;
        }
        let ar = triangle_aspect_ratio(v0, v1, v2);
        if ar > thresholds.max_aspect_ratio {
            n_bad_ar += 1;
        }
        let skew = triangle_skewness(v0, v1, v2);
        if skew > thresholds.max_skewness {
            n_bad_skew += 1;
        }
        let min_a = triangle_min_angle_deg(v0, v1, v2);
        if min_a < thresholds.min_angle {
            n_small_angle += 1;
        }
        let max_a = triangle_max_angle_deg(v0, v1, v2);
        if max_a > thresholds.max_angle {
            n_large_angle += 1;
        }
    }
    let passes = n_bad_ar == 0
        && n_bad_skew == 0
        && n_small_angle == 0
        && n_large_angle == 0
        && n_below_area == 0;
    QualityCheckResult {
        n_bad_aspect_ratio: n_bad_ar,
        n_bad_skewness: n_bad_skew,
        n_small_angle,
        n_large_angle,
        n_below_min_area: n_below_area,
        passes,
    }
}
/// Generate quality improvement suggestions for a mesh.
#[allow(dead_code)]
pub fn suggest_improvements(
    mesh: &TriangleMesh,
    thresholds: &QualityThresholds,
) -> Vec<QualitySuggestion> {
    let mut suggestions = Vec::new();
    let check = check_quality(mesh, thresholds);
    if check.n_below_min_area > 0 {
        suggestions.push(QualitySuggestion::RemoveDegenerates {
            count: check.n_below_min_area,
        });
    }
    if check.n_bad_aspect_ratio > 0 {
        let max_ar = mesh
            .triangles
            .iter()
            .map(|tri| {
                let v0 = mesh.vertices[tri[0] as usize];
                let v1 = mesh.vertices[tri[1] as usize];
                let v2 = mesh.vertices[tri[2] as usize];
                triangle_aspect_ratio(v0, v1, v2)
            })
            .fold(0.0_f32, |a, b| if b.is_finite() { a.max(b) } else { a });
        suggestions.push(QualitySuggestion::RefineHighAspectRatio {
            count: check.n_bad_aspect_ratio,
            max_ar,
        });
    }
    if check.n_large_angle > 0 {
        suggestions.push(QualitySuggestion::SplitLargeAngleTriangles {
            count: check.n_large_angle,
        });
    }
    let n_verts = mesh.vertices.len();
    let mut n_dupes = 0_usize;
    for i in 0..n_verts {
        for j in (i + 1)..n_verts {
            let d = sub(mesh.vertices[i], mesh.vertices[j]);
            if norm(d) < 1e-6 {
                n_dupes += 1;
                break;
            }
        }
    }
    if n_dupes > 0 {
        suggestions.push(QualitySuggestion::RemoveDuplicateVertices {
            estimated_duplicates: n_dupes,
        });
    }
    if suggestions.is_empty() {
        suggestions.push(QualitySuggestion::MeshIsGoodQuality);
    }
    suggestions
}
/// Compute element type statistics for a mesh.
#[allow(dead_code)]
pub fn element_type_stats(mesh: &TriangleMesh) -> ElementTypeStats {
    let mut stats = ElementTypeStats {
        total: mesh.triangles.len() as u32,
        ..Default::default()
    };
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let area = triangle_area(v0, v1, v2);
        if area < 1e-10 {
            stats.n_degenerate += 1;
            continue;
        }
        let angles = triangle_angles_deg(v0, v1, v2);
        let min_a = angles[0].min(angles[1]).min(angles[2]);
        let max_a = angles[0].max(angles[1]).max(angles[2]);
        if min_a >= 50.0 && max_a <= 70.0 {
            stats.n_equilateral += 1;
        } else if (max_a - 90.0).abs() < 5.0 {
            stats.n_right += 1;
        } else if max_a > 95.0 {
            stats.n_obtuse += 1;
        } else {
            stats.n_acute += 1;
        }
    }
    stats
}
/// Compare two meshes and return the differences.
#[allow(dead_code)]
pub fn compare_meshes(mesh_a: &TriangleMesh, mesh_b: &TriangleMesh) -> MeshComparison {
    let rep_a = compute_quality_report(mesh_a);
    let rep_b = compute_quality_report(mesh_b);
    MeshComparison {
        vertex_count_diff: mesh_b.vertex_count() as i64 - mesh_a.vertex_count() as i64,
        triangle_count_diff: mesh_b.triangle_count() as i64 - mesh_a.triangle_count() as i64,
        surface_area_diff: mesh_surface_area(mesh_b) - mesh_surface_area(mesh_a),
        avg_edge_length_diff: mesh_average_edge_length(mesh_b) - mesh_average_edge_length(mesh_a),
        mean_skewness_diff: rep_b.mean_skewness - rep_a.mean_skewness,
        fewer_degenerates: rep_b.n_degenerate < rep_a.n_degenerate,
    }
}
/// Compute an extended quality report with histograms.
#[allow(dead_code)]
pub fn compute_extended_quality_report(mesh: &TriangleMesh) -> ExtendedQualityReport {
    let basic = compute_quality_report(mesh);
    let n = mesh.triangles.len();
    if n == 0 {
        return ExtendedQualityReport {
            basic,
            aspect_ratio_histogram: [0; 10],
            min_angle_histogram: [0; 9],
            median_area: 0.0,
            area_std_dev: 0.0,
        };
    }
    let mut areas = Vec::with_capacity(n);
    let mut ar_hist = [0_u32; 10];
    let mut angle_hist = [0_u32; 9];
    let max_ar = basic.max_aspect_ratio.min(100.0);
    let ar_bin_size = (max_ar - 1.0).max(0.01) / 10.0;
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let area = triangle_area(v0, v1, v2);
        areas.push(area);
        let ar = triangle_aspect_ratio(v0, v1, v2);
        if ar.is_finite() {
            let bin = ((ar - 1.0) / ar_bin_size).floor() as usize;
            let bin = bin.min(9);
            ar_hist[bin] += 1;
        }
        let min_a = triangle_min_angle_deg(v0, v1, v2);
        let abin = (min_a / 10.0).floor() as usize;
        let abin = abin.min(8);
        angle_hist[abin] += 1;
    }
    areas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_area = areas[n / 2];
    let mean_area = basic.mean_area;
    let variance: f32 = areas.iter().map(|a| (a - mean_area).powi(2)).sum::<f32>() / n as f32;
    let area_std_dev = variance.sqrt();
    ExtendedQualityReport {
        basic,
        aspect_ratio_histogram: ar_hist,
        min_angle_histogram: angle_hist,
        median_area,
        area_std_dev,
    }
}
/// Compute the Euler characteristic: V - E + F.
/// For a closed surface, this should be 2.
#[allow(dead_code)]
pub fn euler_characteristic(mesh: &TriangleMesh) -> i64 {
    let v = mesh.vertex_count() as i64;
    let f = mesh.triangle_count() as i64;
    let e = (f * 3) / 2;
    v - e + f
}
/// Compute the ratio of inradius to circumradius for each triangle.
/// For an equilateral triangle, this ratio is 0.5.
/// Closer to 0 means lower quality.
#[allow(dead_code)]
pub fn radius_ratio_quality(mesh: &TriangleMesh) -> Vec<f32> {
    mesh.triangles
        .iter()
        .map(|tri| {
            let v0 = mesh.vertices[tri[0] as usize];
            let v1 = mesh.vertices[tri[1] as usize];
            let v2 = mesh.vertices[tri[2] as usize];
            let r_in = triangle_inscribed_radius(v0, v1, v2);
            let r_circum = triangle_circumscribed_radius(v0, v1, v2);
            if r_circum > 1e-15 && r_circum.is_finite() {
                r_in / r_circum
            } else {
                0.0
            }
        })
        .collect()
}
/// Compute mean radius ratio for the mesh (higher = better quality).
#[allow(dead_code)]
pub fn mean_radius_ratio(mesh: &TriangleMesh) -> f32 {
    let ratios = radius_ratio_quality(mesh);
    if ratios.is_empty() {
        return 0.0;
    }
    ratios.iter().sum::<f32>() / ratios.len() as f32
}
#[cfg(test)]
mod tests {
    use super::*;
    fn equilateral() -> ([f32; 3], [f32; 3], [f32; 3]) {
        let v0 = [0.0_f32, 0.0, 0.0];
        let v1 = [2.0, 0.0, 0.0];
        let v2 = [1.0, 3.0_f32.sqrt(), 0.0];
        (v0, v1, v2)
    }
    fn right_triangle() -> ([f32; 3], [f32; 3], [f32; 3]) {
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    }
    fn make_good_mesh() -> TriangleMesh {
        let (v0, v1, v2) = equilateral();
        TriangleMesh::from_raw(vec![v0, v1, v2], vec![[0, 1, 2]])
    }
    fn make_mixed_mesh() -> TriangleMesh {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 3.0_f32.sqrt(), 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [3, 4, 5], [6, 7, 8]];
        TriangleMesh::from_raw(verts, tris)
    }
    #[test]
    fn test_new_mesh_empty() {
        let m = TriangleMesh::new();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
    }
    #[test]
    fn test_from_raw() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        assert_eq!(m.vertex_count(), 3);
        assert_eq!(m.triangle_count(), 1);
    }
    #[test]
    fn test_default_mesh() {
        let m = TriangleMesh::default();
        assert_eq!(m.vertex_count(), 0);
    }
    #[test]
    fn test_triangle_area_right() {
        let (v0, v1, v2) = right_triangle();
        let area = triangle_area(v0, v1, v2);
        assert!((area - 0.5).abs() < 1e-5, "area={area}");
    }
    #[test]
    fn test_triangle_area_equilateral() {
        let (v0, v1, v2) = equilateral();
        let expected = 3.0_f32.sqrt();
        let area = triangle_area(v0, v1, v2);
        assert!(
            (area - expected).abs() < 1e-4,
            "area={area} expected={expected}"
        );
    }
    #[test]
    fn test_triangle_area_degenerate() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [2.0, 0.0, 0.0];
        let area = triangle_area(v0, v1, v2);
        assert!(area < 1e-6, "area should be 0 for collinear points");
    }
    #[test]
    fn test_aspect_ratio_equilateral() {
        let (v0, v1, v2) = equilateral();
        let ar = triangle_aspect_ratio(v0, v1, v2);
        assert!(ar >= 1.0, "ar={ar}");
        assert!(ar < 1.2, "equilateral ar should be near 1.1547, got {ar}");
    }
    #[test]
    fn test_aspect_ratio_degenerate() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [2.0, 0.0, 0.0];
        let ar = triangle_aspect_ratio(v0, v1, v2);
        assert!(ar.is_infinite(), "degenerate should give infinity");
    }
    #[test]
    fn test_skewness_equilateral() {
        let (v0, v1, v2) = equilateral();
        let s = triangle_skewness(v0, v1, v2);
        assert!(s.abs() < 0.01, "equilateral skewness should be ~0, got {s}");
    }
    #[test]
    fn test_skewness_right_triangle() {
        let (v0, v1, v2) = right_triangle();
        let s = triangle_skewness(v0, v1, v2);
        assert!(
            (s - 0.25).abs() < 0.01,
            "right-triangle skewness should be ~0.25, got {s}"
        );
    }
    #[test]
    fn test_angles_equilateral() {
        let (v0, v1, v2) = equilateral();
        let min_a = triangle_min_angle_deg(v0, v1, v2);
        let max_a = triangle_max_angle_deg(v0, v1, v2);
        assert!((min_a - 60.0).abs() < 0.5, "min_angle={min_a}");
        assert!((max_a - 60.0).abs() < 0.5, "max_angle={max_a}");
    }
    #[test]
    fn test_angles_right_triangle() {
        let (v0, v1, v2) = right_triangle();
        let max_a = triangle_max_angle_deg(v0, v1, v2);
        assert!((max_a - 90.0).abs() < 0.5, "max_angle={max_a}");
    }
    #[test]
    fn test_inscribed_radius_equilateral() {
        let (v0, v1, v2) = equilateral();
        let r = triangle_inscribed_radius(v0, v1, v2);
        let expected = 1.0 / 3.0_f32.sqrt();
        assert!((r - expected).abs() < 0.01, "r={r} expected~{expected}");
    }
    #[test]
    fn test_circumscribed_radius_equilateral() {
        let (v0, v1, v2) = equilateral();
        let r = triangle_circumscribed_radius(v0, v1, v2);
        let expected = 2.0 / 3.0_f32.sqrt();
        assert!((r - expected).abs() < 0.01, "R={r} expected~{expected}");
    }
    #[test]
    fn test_inradius_right_triangle() {
        let (v0, v1, v2) = right_triangle();
        let r = triangle_inscribed_radius(v0, v1, v2);
        let s = (1.0 + 1.0 + 2.0_f32.sqrt()) * 0.5;
        let expected = 0.5 / s;
        assert!((r - expected).abs() < 0.001, "r={r} expected={expected}");
    }
    #[test]
    fn test_quality_report_empty() {
        let m = TriangleMesh::new();
        let rep = compute_quality_report(&m);
        assert_eq!(rep.n_degenerate, 0);
    }
    #[test]
    fn test_quality_report_single_triangle() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let rep = compute_quality_report(&m);
        assert!((rep.min_area - 0.5).abs() < 1e-4);
        assert!((rep.max_area - 0.5).abs() < 1e-4);
        assert_eq!(rep.n_degenerate, 0);
    }
    #[test]
    fn test_quality_report_degenerate() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let rep = compute_quality_report(&m);
        assert_eq!(rep.n_degenerate, 1);
    }
    #[test]
    fn test_bounding_box() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, 0.5, 1.0]];
        let m = TriangleMesh::from_raw(verts, vec![]);
        let (mn, mx) = mesh_bounding_box(&m);
        assert!((mn[0] - (-1.0)).abs() < 1e-6);
        assert!((mx[0] - 1.0).abs() < 1e-6);
        assert!((mx[1] - 2.0).abs() < 1e-6);
        assert!((mx[2] - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_surface_area() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        assert!((mesh_surface_area(&m) - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_average_edge_length() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let avg = mesh_average_edge_length(&m);
        let expected = (1.0 + 1.0 + 2.0_f32.sqrt()) / 3.0;
        assert!((avg - expected).abs() < 1e-5, "avg={avg}");
    }
    #[test]
    fn test_vertex_valence() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let val = vertex_valence(&m);
        assert_eq!(val, vec![1, 1, 1]);
    }
    #[test]
    fn test_signed_volume_tetrahedron() {
        let verts = vec![
            [0.0_f32, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let m = TriangleMesh::from_raw(verts, tris);
        let vol = mesh_volume_signed(&m).abs();
        let expected = 1.0 / 6.0;
        assert!((vol - expected).abs() < 0.01, "vol={vol}");
    }
    #[test]
    fn test_remove_duplicate_vertices() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [3, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let m2 = remove_duplicate_vertices(&m, 1e-5);
        assert_eq!(m2.vertex_count(), 3, "should have 3 unique vertices");
        assert_eq!(m2.triangle_count(), 2);
        assert_eq!(m2.triangles[0][0], m2.triangles[1][0]);
    }
    #[test]
    fn test_remove_degenerate_triangles() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3]];
        let m = TriangleMesh::from_raw(verts, tris);
        let m2 = remove_degenerate_triangles(&m);
        assert_eq!(m2.triangle_count(), 1);
    }
    #[test]
    fn test_flip_normals() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let m = TriangleMesh::from_raw(verts, tris);
        let m2 = flip_normals(&m);
        assert_eq!(m2.triangles[0], [0, 2, 1]);
    }
    #[test]
    fn test_quality_thresholds_default() {
        let t = QualityThresholds::default();
        assert!(t.max_aspect_ratio > 1.0);
        assert!(t.max_skewness > 0.0);
        assert!(t.min_angle > 0.0);
        assert!(t.max_angle < 180.0);
    }
    #[test]
    fn test_check_quality_good_mesh() {
        let m = make_good_mesh();
        let t = QualityThresholds::default();
        let r = check_quality(&m, &t);
        assert!(
            r.passes,
            "equilateral triangle should pass default thresholds"
        );
    }
    #[test]
    fn test_check_quality_degenerate() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let m = TriangleMesh::from_raw(verts, vec![[0, 1, 2]]);
        let t = QualityThresholds::default();
        let r = check_quality(&m, &t);
        assert!(!r.passes);
        assert!(r.n_below_min_area > 0);
    }
    #[test]
    fn test_check_quality_strict_thresholds() {
        let (v0, v1, v2) = right_triangle();
        let m = TriangleMesh::from_raw(vec![v0, v1, v2], vec![[0, 1, 2]]);
        let t = QualityThresholds {
            max_angle: 85.0,
            ..QualityThresholds::default()
        };
        let r = check_quality(&m, &t);
        assert!(r.n_large_angle > 0);
        assert!(!r.passes);
    }
    #[test]
    fn test_suggest_good_mesh() {
        let m = make_good_mesh();
        let t = QualityThresholds::default();
        let suggestions = suggest_improvements(&m, &t);
        assert!(
            suggestions
                .iter()
                .any(|s| matches!(s, QualitySuggestion::MeshIsGoodQuality))
        );
    }
    #[test]
    fn test_suggest_degenerate_mesh() {
        let m = make_mixed_mesh();
        let t = QualityThresholds::default();
        let suggestions = suggest_improvements(&m, &t);
        assert!(
            suggestions
                .iter()
                .any(|s| matches!(s, QualitySuggestion::RemoveDegenerates { .. }))
        );
    }
    #[test]
    fn test_element_type_stats_equilateral() {
        let m = make_good_mesh();
        let stats = element_type_stats(&m);
        assert_eq!(stats.total, 1);
        assert_eq!(stats.n_equilateral, 1);
        assert_eq!(stats.n_degenerate, 0);
    }
    #[test]
    fn test_element_type_stats_mixed() {
        let m = make_mixed_mesh();
        let stats = element_type_stats(&m);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.n_degenerate, 1);
        assert_eq!(stats.n_equilateral, 1);
        assert_eq!(stats.n_right, 1);
    }
    #[test]
    fn test_element_type_stats_empty() {
        let m = TriangleMesh::new();
        let stats = element_type_stats(&m);
        assert_eq!(stats.total, 0);
    }
    #[test]
    fn test_compare_same_mesh() {
        let m = make_good_mesh();
        let cmp = compare_meshes(&m, &m);
        assert_eq!(cmp.vertex_count_diff, 0);
        assert_eq!(cmp.triangle_count_diff, 0);
        assert!(cmp.surface_area_diff.abs() < 1e-6);
    }
    #[test]
    fn test_compare_different_meshes() {
        let m1 = make_good_mesh();
        let m2 = make_mixed_mesh();
        let cmp = compare_meshes(&m1, &m2);
        assert!(cmp.vertex_count_diff > 0);
        assert!(cmp.triangle_count_diff > 0);
    }
    #[test]
    fn test_extended_report_empty() {
        let m = TriangleMesh::new();
        let rep = compute_extended_quality_report(&m);
        assert_eq!(rep.median_area, 0.0);
    }
    #[test]
    fn test_extended_report_single() {
        let m = make_good_mesh();
        let rep = compute_extended_quality_report(&m);
        assert!(rep.median_area > 0.0);
        assert!(rep.area_std_dev >= 0.0);
    }
    #[test]
    fn test_extended_report_histograms() {
        let m = make_mixed_mesh();
        let rep = compute_extended_quality_report(&m);
        let total_ar: u32 = rep.aspect_ratio_histogram.iter().sum();
        assert!(total_ar <= 3);
        let total_angle: u32 = rep.min_angle_histogram.iter().sum();
        assert_eq!(total_angle, 3);
    }
    #[test]
    fn test_euler_characteristic_single() {
        let m = make_good_mesh();
        let chi = euler_characteristic(&m);
        assert!(chi > 0);
    }
    #[test]
    fn test_radius_ratio_equilateral() {
        let m = make_good_mesh();
        let ratios = radius_ratio_quality(&m);
        assert_eq!(ratios.len(), 1);
        assert!((ratios[0] - 0.5).abs() < 0.01, "ratio={}", ratios[0]);
    }
    #[test]
    fn test_mean_radius_ratio_good_mesh() {
        let m = make_good_mesh();
        let mean = mean_radius_ratio(&m);
        assert!((mean - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_mean_radius_ratio_empty() {
        let m = TriangleMesh::new();
        assert_eq!(mean_radius_ratio(&m), 0.0);
    }
    #[test]
    fn test_radius_ratio_right_triangle() {
        let (v0, v1, v2) = right_triangle();
        let m = TriangleMesh::from_raw(vec![v0, v1, v2], vec![[0, 1, 2]]);
        let ratios = radius_ratio_quality(&m);
        assert!(ratios[0] < 0.5);
        assert!(ratios[0] > 0.0);
    }
}
/// Compute a full quality report for a mesh.
#[allow(dead_code)]
pub fn compute_full_quality_report(mesh: &TriangleMesh) -> FullQualityReport {
    let basic = compute_quality_report(mesh);
    let stats = element_type_stats(mesh);
    let n = mesh.triangles.len();
    if n == 0 {
        return FullQualityReport {
            basic,
            global_min_angle_deg: 0.0,
            global_max_angle_deg: 0.0,
            mean_angle_deg: 60.0,
            min_jacobian: 0.0,
            max_jacobian: 0.0,
            mean_jacobian: 0.0,
            element_type_counts: Vec::new(),
            n_poorly_shaped: 0,
            total_surface_area: 0.0,
        };
    }
    let mut g_min_angle = f32::INFINITY;
    let mut g_max_angle = f32::NEG_INFINITY;
    let mut sum_angle = 0.0_f32;
    let mut min_jac = f32::INFINITY;
    let mut max_jac = f32::NEG_INFINITY;
    let mut sum_jac = 0.0_f32;
    let mut n_poorly = 0_u32;
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let min_a = triangle_min_angle_deg(v0, v1, v2);
        let max_a = triangle_max_angle_deg(v0, v1, v2);
        let mean_a = (min_a + max_a) / 2.0;
        g_min_angle = g_min_angle.min(min_a);
        g_max_angle = g_max_angle.max(max_a);
        sum_angle += mean_a;
        let area = triangle_area(v0, v1, v2);
        let jac = 2.0 * area;
        min_jac = min_jac.min(jac);
        max_jac = max_jac.max(jac);
        sum_jac += jac;
        let ar = triangle_aspect_ratio(v0, v1, v2);
        if ar > 5.0 || min_a < 15.0 {
            n_poorly += 1;
        }
    }
    let element_type_counts = vec![
        ("equilateral".to_string(), stats.n_equilateral),
        ("right".to_string(), stats.n_right),
        ("obtuse".to_string(), stats.n_obtuse),
        ("acute".to_string(), stats.n_acute),
        ("degenerate".to_string(), stats.n_degenerate),
    ];
    FullQualityReport {
        global_min_angle_deg: if n > 0 { g_min_angle } else { 0.0 },
        global_max_angle_deg: if n > 0 { g_max_angle } else { 0.0 },
        mean_angle_deg: if n > 0 { sum_angle / n as f32 } else { 60.0 },
        min_jacobian: if n > 0 { min_jac } else { 0.0 },
        max_jacobian: if n > 0 { max_jac } else { 0.0 },
        mean_jacobian: if n > 0 { sum_jac / n as f32 } else { 0.0 },
        element_type_counts,
        n_poorly_shaped: n_poorly,
        total_surface_area: mesh_surface_area(mesh),
        basic,
    }
}
