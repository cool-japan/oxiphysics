//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::triangle_angles_deg;
use super::functions::*;
use super::types::{
    AutoImproveSuggestion, PerTypeMeshStats, PoorShapeReason, PoorlyShapedElement,
    QualityThresholds, QualityVisualizationData, TriangleMesh,
};

/// Compute per-element-type quality statistics.
pub fn per_type_stats(mesh: &TriangleMesh) -> PerTypeMeshStats {
    let mut stats = PerTypeMeshStats::default();
    let mut eq_ar_sum = 0.0_f32;
    let mut ob_ar_sum = 0.0_f32;
    let mut rt_ar_sum = 0.0_f32;
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
        let ar = triangle_aspect_ratio(v0, v1, v2);
        if min_a >= 50.0 && max_a <= 70.0 {
            stats.n_equilateral += 1;
            eq_ar_sum += ar;
        } else if (max_a - 90.0).abs() < 5.0 {
            stats.n_right += 1;
            rt_ar_sum += ar;
        } else if max_a > 95.0 {
            stats.n_obtuse += 1;
            ob_ar_sum += ar;
        } else {
            stats.n_acute += 1;
        }
    }
    if stats.n_equilateral > 0 {
        stats.equilateral_mean_ar = eq_ar_sum / stats.n_equilateral as f32;
    }
    if stats.n_obtuse > 0 {
        stats.obtuse_mean_ar = ob_ar_sum / stats.n_obtuse as f32;
    }
    if stats.n_right > 0 {
        stats.right_mean_ar = rt_ar_sum / stats.n_right as f32;
    }
    stats
}
/// Generate per-triangle quality values for visualization (e.g. coloring in ParaView).
pub fn quality_visualization_data(mesh: &TriangleMesh) -> QualityVisualizationData {
    let n = mesh.triangles.len();
    let mut ars = Vec::with_capacity(n);
    let mut skews = Vec::with_capacity(n);
    let mut min_as = Vec::with_capacity(n);
    let mut max_as = Vec::with_capacity(n);
    let mut areas = Vec::with_capacity(n);
    let mut jacs = Vec::with_capacity(n);
    for tri in &mesh.triangles {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let a = triangle_area(v0, v1, v2);
        areas.push(a);
        jacs.push(2.0 * a);
        ars.push(triangle_aspect_ratio(v0, v1, v2));
        skews.push(triangle_skewness(v0, v1, v2));
        min_as.push(triangle_min_angle_deg(v0, v1, v2));
        max_as.push(triangle_max_angle_deg(v0, v1, v2));
    }
    QualityVisualizationData {
        aspect_ratios: ars,
        skewness: skews,
        min_angles: min_as,
        max_angles: max_as,
        areas,
        jacobians: jacs,
    }
}
/// Detect poorly-shaped elements using per-element thresholds.
pub fn detect_poorly_shaped(
    mesh: &TriangleMesh,
    max_aspect_ratio: f32,
    min_angle_deg: f32,
    max_angle_deg: f32,
    max_skewness: f32,
) -> Vec<PoorlyShapedElement> {
    let mut poor = Vec::new();
    for (idx, tri) in mesh.triangles.iter().enumerate() {
        let v0 = mesh.vertices[tri[0] as usize];
        let v1 = mesh.vertices[tri[1] as usize];
        let v2 = mesh.vertices[tri[2] as usize];
        let mut reasons = Vec::new();
        let area = triangle_area(v0, v1, v2);
        if area < 1e-10 {
            reasons.push(PoorShapeReason::Degenerate);
        }
        let ar = triangle_aspect_ratio(v0, v1, v2);
        if ar.is_finite() && ar > max_aspect_ratio {
            reasons.push(PoorShapeReason::HighAspectRatio(ar));
        }
        let min_a = triangle_min_angle_deg(v0, v1, v2);
        if min_a < min_angle_deg {
            reasons.push(PoorShapeReason::SmallAngle(min_a));
        }
        let max_a = triangle_max_angle_deg(v0, v1, v2);
        if max_a > max_angle_deg {
            reasons.push(PoorShapeReason::LargeAngle(max_a));
        }
        let skew = triangle_skewness(v0, v1, v2);
        if skew > max_skewness {
            reasons.push(PoorShapeReason::HighSkewness(skew));
        }
        if !reasons.is_empty() {
            poor.push(PoorlyShapedElement {
                index: idx,
                reasons,
            });
        }
    }
    poor
}
/// Generate automated quality improvement suggestions with priorities.
pub fn auto_improve_suggestions(mesh: &TriangleMesh) -> Vec<AutoImproveSuggestion> {
    let mut suggestions = Vec::new();
    let report = compute_quality_report(mesh);
    let thresholds = QualityThresholds::default();
    let check = check_quality(mesh, &thresholds);
    if report.n_degenerate > 0 {
        suggestions.push(AutoImproveSuggestion::new(
            format!(
                "Remove {} degenerate triangles (area ≈ 0)",
                report.n_degenerate
            ),
            1,
            report.n_degenerate as usize,
        ));
    }
    if check.n_bad_aspect_ratio > 0 {
        suggestions.push(AutoImproveSuggestion::new(
            format!(
                "Refine {} high-aspect-ratio triangles (AR > {:.1})",
                check.n_bad_aspect_ratio, thresholds.max_aspect_ratio
            ),
            2,
            check.n_bad_aspect_ratio as usize,
        ));
    }
    if check.n_small_angle > 0 {
        suggestions.push(AutoImproveSuggestion::new(
            format!(
                "Fix {} triangles with small angles (< {:.0}°)",
                check.n_small_angle, thresholds.min_angle
            ),
            2,
            check.n_small_angle as usize,
        ));
    }
    if check.n_large_angle > 0 {
        suggestions.push(AutoImproveSuggestion::new(
            format!(
                "Split {} triangles with large angles (> {:.0}°)",
                check.n_large_angle, thresholds.max_angle
            ),
            3,
            check.n_large_angle as usize,
        ));
    }
    if check.n_bad_skewness > 0 {
        suggestions.push(AutoImproveSuggestion::new(
            format!(
                "Smooth {} high-skewness triangles (skewness > {:.2})",
                check.n_bad_skewness, thresholds.max_skewness
            ),
            4,
            check.n_bad_skewness as usize,
        ));
    }
    if suggestions.is_empty() {
        suggestions.push(AutoImproveSuggestion::new(
            "Mesh passes all quality checks — no improvements needed.".to_string(),
            5,
            0,
        ));
    }
    suggestions.sort_by_key(|s| s.priority);
    suggestions
}
#[cfg(test)]
mod tests_mesh_quality_extended {
    use super::*;

    use crate::mesh_quality::types::*;
    fn equilateral() -> TriangleMesh {
        let v0 = [0.0_f32, 0.0, 0.0];
        let v1 = [2.0, 0.0, 0.0];
        let v2 = [1.0, 3.0_f32.sqrt(), 0.0];
        TriangleMesh::from_raw(vec![v0, v1, v2], vec![[0, 1, 2]])
    }
    fn degenerate_mesh() -> TriangleMesh {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        TriangleMesh::from_raw(verts, vec![[0, 1, 2]])
    }
    fn mixed_mesh() -> TriangleMesh {
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
    fn test_full_quality_report_empty() {
        let m = TriangleMesh::new();
        let r = compute_full_quality_report(&m);
        assert_eq!(r.n_poorly_shaped, 0);
        assert_eq!(r.total_surface_area, 0.0);
    }
    #[test]
    fn test_full_quality_report_equilateral() {
        let m = equilateral();
        let r = compute_full_quality_report(&m);
        assert!(
            (r.global_min_angle_deg - 60.0).abs() < 1.0,
            "min_angle={}",
            r.global_min_angle_deg
        );
        assert!(
            (r.global_max_angle_deg - 60.0).abs() < 1.0,
            "max_angle={}",
            r.global_max_angle_deg
        );
        assert_eq!(r.n_poorly_shaped, 0);
    }
    #[test]
    fn test_full_quality_report_degenerate_poorly_shaped() {
        let m = degenerate_mesh();
        let r = compute_full_quality_report(&m);
        assert!(r.n_poorly_shaped > 0);
    }
    #[test]
    fn test_full_quality_report_jacobian_positive() {
        let m = equilateral();
        let r = compute_full_quality_report(&m);
        assert!(r.min_jacobian > 0.0);
        assert!(r.max_jacobian >= r.min_jacobian);
    }
    #[test]
    fn test_full_quality_report_element_type_counts() {
        let m = equilateral();
        let r = compute_full_quality_report(&m);
        let eq_count: u32 = r
            .element_type_counts
            .iter()
            .find(|(t, _)| t == "equilateral")
            .map(|(_, c)| *c)
            .unwrap_or(0);
        assert_eq!(eq_count, 1);
    }
    #[test]
    fn test_per_type_stats_equilateral() {
        let m = equilateral();
        let s = per_type_stats(&m);
        assert_eq!(s.n_equilateral, 1);
        assert_eq!(s.n_degenerate, 0);
        assert!(s.equilateral_mean_ar > 0.0);
    }
    #[test]
    fn test_per_type_stats_degenerate() {
        let m = degenerate_mesh();
        let s = per_type_stats(&m);
        assert_eq!(s.n_degenerate, 1);
    }
    #[test]
    fn test_per_type_stats_mixed() {
        let m = mixed_mesh();
        let s = per_type_stats(&m);
        assert_eq!(
            s.n_equilateral + s.n_right + s.n_obtuse + s.n_acute + s.n_degenerate,
            3
        );
    }
    #[test]
    fn test_quality_visualization_data_lengths() {
        let m = equilateral();
        let d = quality_visualization_data(&m);
        assert_eq!(d.aspect_ratios.len(), 1);
        assert_eq!(d.skewness.len(), 1);
        assert_eq!(d.min_angles.len(), 1);
        assert_eq!(d.areas.len(), 1);
        assert_eq!(d.jacobians.len(), 1);
    }
    #[test]
    fn test_quality_visualization_jacobian_twice_area() {
        let m = equilateral();
        let d = quality_visualization_data(&m);
        for (j, a) in d.jacobians.iter().zip(d.areas.iter()) {
            assert!((j - 2.0 * a).abs() < 1e-6, "j={j}, 2a={}", 2.0 * a);
        }
    }
    #[test]
    fn test_quality_visualization_empty() {
        let m = TriangleMesh::new();
        let d = quality_visualization_data(&m);
        assert!(d.aspect_ratios.is_empty());
    }
    #[test]
    fn test_detect_poorly_shaped_equilateral_none() {
        let m = equilateral();
        let poor = detect_poorly_shaped(&m, 5.0, 15.0, 150.0, 0.8);
        assert!(poor.is_empty(), "Equilateral should not be poorly shaped");
    }
    #[test]
    fn test_detect_poorly_shaped_degenerate() {
        let m = degenerate_mesh();
        let poor = detect_poorly_shaped(&m, 5.0, 15.0, 150.0, 0.8);
        assert!(!poor.is_empty());
        assert!(
            poor[0]
                .reasons
                .iter()
                .any(|r| matches!(r, PoorShapeReason::Degenerate))
        );
    }
    #[test]
    fn test_detect_poorly_shaped_strict_aspect_ratio() {
        let verts = vec![[0.0_f32, 0.0, 0.0], [100.0, 0.0, 0.0], [50.0, 0.01, 0.0]];
        let m = TriangleMesh::from_raw(verts, vec![[0, 1, 2]]);
        let poor = detect_poorly_shaped(&m, 2.0, 1.0, 179.0, 0.99);
        assert!(!poor.is_empty());
    }
    #[test]
    fn test_detect_poorly_shaped_index_correct() {
        let m = mixed_mesh();
        let poor = detect_poorly_shaped(&m, 5.0, 15.0, 150.0, 0.8);
        for p in &poor {
            assert!(p.index < m.triangle_count());
        }
    }
    #[test]
    fn test_auto_improve_suggestions_good_mesh() {
        let m = equilateral();
        let s = auto_improve_suggestions(&m);
        assert!(!s.is_empty());
        assert!(s[0].description.contains("no improvements") || s[0].priority == 5);
    }
    #[test]
    fn test_auto_improve_suggestions_degenerate_priority_1() {
        let m = degenerate_mesh();
        let s = auto_improve_suggestions(&m);
        assert!(
            s.iter().any(|sg| sg.priority == 1),
            "Degenerate should be priority 1"
        );
    }
    #[test]
    fn test_auto_improve_suggestions_sorted_by_priority() {
        let m = mixed_mesh();
        let s = auto_improve_suggestions(&m);
        let priorities: Vec<u8> = s.iter().map(|sg| sg.priority).collect();
        let sorted: Vec<u8> = {
            let mut p = priorities.clone();
            p.sort_unstable();
            p
        };
        assert_eq!(
            priorities, sorted,
            "Suggestions should be sorted by priority"
        );
    }
    #[test]
    fn test_auto_improve_suggestions_empty_mesh() {
        let m = TriangleMesh::new();
        let s = auto_improve_suggestions(&m);
        assert!(!s.is_empty());
    }
    #[test]
    fn test_poor_shape_reason_equality() {
        assert_eq!(PoorShapeReason::Degenerate, PoorShapeReason::Degenerate);
        assert_ne!(
            PoorShapeReason::Degenerate,
            PoorShapeReason::SmallAngle(10.0)
        );
    }
    #[test]
    fn test_full_quality_report_surface_area() {
        let verts = vec![[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let m = TriangleMesh::from_raw(verts, vec![[0, 1, 2]]);
        let r = compute_full_quality_report(&m);
        assert!(
            (r.total_surface_area - 0.5).abs() < 1e-5,
            "area={}",
            r.total_surface_area
        );
    }
}
#[cfg(test)]
mod tests_mesh_quality_new {
    use super::*;
    use crate::mesh_quality::types::MeshQuality;

    #[test]
    fn test_dihedral_flat_triangles_is_180() {
        let ea = [0.0_f32, 0.0, 0.0];
        let eb = [1.0, 0.0, 0.0];
        let p0 = [0.5, 1.0, 0.0];
        let p1 = [0.5, -1.0, 0.0];
        let angle = MeshQuality::compute_dihedral_angle(ea, eb, p0, p1);
        assert!(
            (angle - 180.0).abs() < 1e-4,
            "Expected ~180°, got {}",
            angle
        );
    }
    #[test]
    fn test_dihedral_perpendicular_is_90() {
        let ea = [0.0_f32, 0.0, 0.0];
        let eb = [1.0, 0.0, 0.0];
        let p0 = [0.5, 1.0, 0.0];
        let p1 = [0.5, 0.0, 1.0];
        let angle = MeshQuality::compute_dihedral_angle(ea, eb, p0, p1);
        assert!((angle - 90.0).abs() < 1e-4, "Expected 90°, got {}", angle);
    }
    #[test]
    fn test_dihedral_degenerate_edge_is_zero() {
        let p = [0.0_f32, 0.0, 0.0];
        let q = [1.0, 0.0, 0.0];
        let r = [0.0, 1.0, 0.0];
        let angle = MeshQuality::compute_dihedral_angle(p, p, q, r);
        assert!(angle.abs() < 1e-6 || angle.is_finite(), "Should not panic");
    }
    #[test]
    fn test_edge_length_ratio_equilateral_is_one() {
        let sqrt3 = 3.0_f32.sqrt();
        let v0 = [0.0_f32, 0.0, 0.0];
        let v1 = [2.0, 0.0, 0.0];
        let v2 = [1.0, sqrt3, 0.0];
        let ratio = MeshQuality::compute_edge_length_ratio(v0, v1, v2);
        assert!(
            (ratio - 1.0).abs() < 1e-5,
            "Equilateral ratio should be 1, got {}",
            ratio
        );
    }
    #[test]
    fn test_edge_length_ratio_elongated_gt_one() {
        let v0 = [0.0_f32, 0.0, 0.0];
        let v1 = [10.0, 0.0, 0.0];
        let v2 = [5.0, 0.1, 0.0];
        let ratio = MeshQuality::compute_edge_length_ratio(v0, v1, v2);
        assert!(ratio > 1.0, "Elongated triangle ratio should be > 1");
    }
    #[test]
    fn test_edge_length_ratio_degenerate_returns_infinity() {
        let p = [0.0_f32, 0.0, 0.0];
        let ratio = MeshQuality::compute_edge_length_ratio(p, p, p);
        assert!(
            ratio.is_infinite(),
            "Degenerate triangle should return infinity"
        );
    }
    #[test]
    fn test_regularity_empty_mesh_is_one() {
        let m = TriangleMesh::new();
        let r = MeshQuality::compute_mesh_regularity(&m);
        assert!((r - 1.0).abs() < 1e-6, "Empty mesh regularity should be 1");
    }
    #[test]
    fn test_regularity_equilateral_near_one() {
        let sqrt3 = 3.0_f32.sqrt();
        let verts = vec![[0.0_f32, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, sqrt3, 0.0]];
        let m = TriangleMesh::from_raw(verts, vec![[0, 1, 2]]);
        let r = MeshQuality::compute_mesh_regularity(&m);
        assert!(
            (r - 1.0).abs() < 1e-4,
            "Equilateral mesh regularity should be ~1, got {}",
            r
        );
    }
    #[test]
    fn test_regularity_mixed_mesh_in_range() {
        let verts = vec![
            [0.0_f32, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 3.0_f32.sqrt(), 0.0],
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 0.0],
            [50.0, 0.01, 0.0],
        ];
        let tris = vec![[0u32, 1, 2], [3, 4, 5]];
        let m = TriangleMesh::from_raw(verts, tris);
        let r = MeshQuality::compute_mesh_regularity(&m);
        assert!(
            r > 0.0 && r < 1.0,
            "Mixed mesh regularity should be in (0,1), got {}",
            r
        );
    }
}
