//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// A suggestion for improving mesh quality.
#[derive(Debug, Clone)]
pub enum QualitySuggestion {
    /// Remove degenerate triangles.
    RemoveDegenerates {
        /// Number of degenerate triangles.
        count: u32,
    },
    /// Refine triangles with high aspect ratio.
    RefineHighAspectRatio {
        /// Number of high aspect ratio triangles.
        count: u32,
        /// Maximum aspect ratio found.
        max_ar: f32,
    },
    /// Split triangles with large angles.
    SplitLargeAngleTriangles {
        /// Number of large-angle triangles.
        count: u32,
    },
    /// Remove duplicate vertices.
    RemoveDuplicateVertices {
        /// Estimated number of duplicate vertices.
        estimated_duplicates: usize,
    },
    /// Mesh is good quality.
    MeshIsGoodQuality,
}
/// Extended quality report with histogram data.
#[derive(Debug, Clone)]
pub struct ExtendedQualityReport {
    /// Basic quality report.
    pub basic: MeshQualityReport,
    /// Histogram of aspect ratios in 10 bins from 1.0 to max.
    pub aspect_ratio_histogram: [u32; 10],
    /// Histogram of minimum angles in 9 bins (0-10, 10-20, ..., 80-90 degrees).
    pub min_angle_histogram: [u32; 9],
    /// Median area.
    pub median_area: f32,
    /// Standard deviation of areas.
    pub area_std_dev: f32,
}
/// Per-triangle quality values for visualization.
#[derive(Debug, Clone)]
pub struct QualityVisualizationData {
    /// Aspect ratio per triangle.
    pub aspect_ratios: Vec<f32>,
    /// Skewness per triangle.
    pub skewness: Vec<f32>,
    /// Minimum angle per triangle.
    pub min_angles: Vec<f32>,
    /// Maximum angle per triangle.
    pub max_angles: Vec<f32>,
    /// Area per triangle.
    pub areas: Vec<f32>,
    /// Jacobian (2*area) per triangle.
    pub jacobians: Vec<f32>,
}
/// Result of comparing two meshes.
#[derive(Debug, Clone)]
pub struct MeshComparison {
    /// Difference in vertex count (mesh_b - mesh_a).
    pub vertex_count_diff: i64,
    /// Difference in triangle count (mesh_b - mesh_a).
    pub triangle_count_diff: i64,
    /// Difference in total surface area.
    pub surface_area_diff: f32,
    /// Difference in average edge length.
    pub avg_edge_length_diff: f32,
    /// Difference in mean skewness.
    pub mean_skewness_diff: f32,
    /// Whether mesh_b has fewer degenerate triangles.
    pub fewer_degenerates: bool,
}
/// A poorly-shaped element with its index and the reason it is flagged.
#[derive(Debug, Clone)]
pub struct PoorlyShapedElement {
    /// Triangle index.
    pub index: usize,
    /// Reason(s) for being flagged.
    pub reasons: Vec<PoorShapeReason>,
}
/// Reason a triangle is considered poorly shaped.
#[derive(Debug, Clone, PartialEq)]
pub enum PoorShapeReason {
    /// High aspect ratio (> threshold).
    HighAspectRatio(f32),
    /// Small minimum angle (< threshold degrees).
    SmallAngle(f32),
    /// Large maximum angle (> threshold degrees).
    LargeAngle(f32),
    /// Very small area (degenerate).
    Degenerate,
    /// High skewness (> threshold).
    HighSkewness(f32),
}
/// Additional triangle mesh quality metrics.
pub struct MeshQuality;
impl MeshQuality {
    /// Compute the dihedral angle (in degrees) between two adjacent triangles
    /// sharing an edge defined by vertices `ea` and `eb`.
    ///
    /// The first triangle is `(ea, eb, p0)` and the second is `(ea, eb, p1)`.
    /// Returns the angle between their face normals, in the range `[0°, 180°]`.
    pub fn compute_dihedral_angle(ea: [f32; 3], eb: [f32; 3], p0: [f32; 3], p1: [f32; 3]) -> f32 {
        let n0 = cross(sub(eb, ea), sub(p0, ea));
        let n1 = cross(sub(eb, ea), sub(p1, ea));
        let len0 = norm(n0);
        let len1 = norm(n1);
        if len0 < 1e-15 || len1 < 1e-15 {
            return 0.0;
        }
        let cos_theta = (dot(n0, n1) / (len0 * len1)).clamp(-1.0, 1.0);
        cos_theta.acos().to_degrees()
    }
    /// Compute the edge-length ratio (max edge / min edge) for a triangle.
    ///
    /// A value of 1 indicates an equilateral triangle.  Returns `f32::INFINITY`
    /// if the minimum edge length is effectively zero (degenerate triangle).
    pub fn compute_edge_length_ratio(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> f32 {
        let e0 = norm(sub(v1, v0));
        let e1 = norm(sub(v2, v1));
        let e2 = norm(sub(v0, v2));
        let max_e = e0.max(e1).max(e2);
        let min_e = e0.min(e1).min(e2);
        if min_e < 1e-15 {
            f32::INFINITY
        } else {
            max_e / min_e
        }
    }
    /// Compute a global mesh regularity score in `[0, 1]`.
    ///
    /// Regularity is defined as:
    /// ```text
    /// R = 1 / (1 + sigma)
    /// ```
    /// where `sigma` is the standard deviation of the per-triangle aspect ratios.
    /// A perfectly regular mesh (all equilateral triangles) has `sigma = 0` → `R = 1`.
    /// A highly irregular mesh has large `sigma` → `R → 0`.
    ///
    /// Returns `1.0` for an empty mesh.
    pub fn compute_mesh_regularity(mesh: &TriangleMesh) -> f32 {
        let n = mesh.triangle_count();
        if n == 0 {
            return 1.0;
        }
        let ratios: Vec<f32> = mesh
            .triangles
            .iter()
            .map(|t| {
                let v0 = mesh.vertices[t[0] as usize];
                let v1 = mesh.vertices[t[1] as usize];
                let v2 = mesh.vertices[t[2] as usize];
                triangle_aspect_ratio(v0, v1, v2)
            })
            .collect();
        let finite: Vec<f32> = ratios.iter().copied().filter(|v| v.is_finite()).collect();
        if finite.is_empty() {
            return 0.0;
        }
        let mean = finite.iter().sum::<f32>() / finite.len() as f32;
        let variance =
            finite.iter().map(|&r| (r - mean) * (r - mean)).sum::<f32>() / finite.len() as f32;
        let sigma = variance.sqrt();
        1.0 / (1.0 + sigma)
    }
}
/// A triangle mesh consisting of vertices and triangle connectivity.
pub struct TriangleMesh {
    /// 3D vertex positions.
    pub vertices: Vec<[f32; 3]>,
    /// Triangle face indices (each entry is \[v0, v1, v2\]).
    pub triangles: Vec<[u32; 3]>,
}
impl TriangleMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }
    /// Create a mesh from raw vertex and triangle data.
    pub fn from_raw(verts: Vec<[f32; 3]>, tris: Vec<[u32; 3]>) -> Self {
        Self {
            vertices: verts,
            triangles: tris,
        }
    }
    /// Return the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    /// Return the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }
}
/// Full mesh quality report including angles, Jacobians, and element counts.
#[derive(Debug, Clone)]
pub struct FullQualityReport {
    /// Basic quality report.
    pub basic: MeshQualityReport,
    /// Minimum interior angle across all triangles (degrees).
    pub global_min_angle_deg: f32,
    /// Maximum interior angle across all triangles (degrees).
    pub global_max_angle_deg: f32,
    /// Mean interior angle across all triangles (degrees).
    pub mean_angle_deg: f32,
    /// Minimum Jacobian (2D: 2*area, proxy for triangle quality).
    pub min_jacobian: f32,
    /// Maximum Jacobian.
    pub max_jacobian: f32,
    /// Mean Jacobian.
    pub mean_jacobian: f32,
    /// Count of triangles per element type string.
    pub element_type_counts: Vec<(String, u32)>,
    /// Poorly-shaped element count (aspect ratio > 5 OR min angle < 15 deg).
    pub n_poorly_shaped: u32,
    /// Total surface area.
    pub total_surface_area: f32,
}
/// Quality thresholds for evaluating mesh quality.
#[derive(Debug, Clone)]
pub struct QualityThresholds {
    /// Maximum acceptable aspect ratio.
    pub max_aspect_ratio: f32,
    /// Maximum acceptable skewness (0..1).
    pub max_skewness: f32,
    /// Minimum acceptable angle (degrees).
    pub min_angle: f32,
    /// Maximum acceptable angle (degrees).
    pub max_angle: f32,
    /// Minimum acceptable area.
    pub min_area: f32,
}
/// Statistics about element shapes in the mesh.
#[derive(Debug, Clone, Default)]
pub struct ElementTypeStats {
    /// Number of near-equilateral triangles (all angles within 50-70 degrees).
    pub n_equilateral: u32,
    /// Number of right triangles (one angle within 85-95 degrees).
    pub n_right: u32,
    /// Number of obtuse triangles (one angle > 95 degrees).
    pub n_obtuse: u32,
    /// Number of acute triangles (all angles < 85 degrees, non-equilateral).
    pub n_acute: u32,
    /// Number of degenerate triangles (area < threshold).
    pub n_degenerate: u32,
    /// Total number of triangles.
    pub total: u32,
}
/// Quality statistics broken down per element category.
#[derive(Debug, Clone, Default)]
pub struct PerTypeMeshStats {
    /// Mean aspect ratio of equilateral triangles.
    pub equilateral_mean_ar: f32,
    /// Mean aspect ratio of obtuse triangles.
    pub obtuse_mean_ar: f32,
    /// Mean aspect ratio of right triangles.
    pub right_mean_ar: f32,
    /// Count of each category.
    pub n_equilateral: u32,
    /// Count of obtuse triangles.
    pub n_obtuse: u32,
    /// Count of right triangles.
    pub n_right: u32,
    /// Count of acute triangles.
    pub n_acute: u32,
    /// Count of degenerate triangles.
    pub n_degenerate: u32,
}
/// Result of checking a mesh against quality thresholds.
#[derive(Debug, Clone)]
pub struct QualityCheckResult {
    /// Number of triangles exceeding max aspect ratio.
    pub n_bad_aspect_ratio: u32,
    /// Number of triangles exceeding max skewness.
    pub n_bad_skewness: u32,
    /// Number of triangles with angles below min_angle.
    pub n_small_angle: u32,
    /// Number of triangles with angles above max_angle.
    pub n_large_angle: u32,
    /// Number of triangles below min area.
    pub n_below_min_area: u32,
    /// Whether the mesh passes all checks.
    pub passes: bool,
}
/// Summary of mesh quality metrics computed over all triangles.
#[derive(Debug, Clone)]
pub struct MeshQualityReport {
    /// Minimum triangle area.
    pub min_area: f32,
    /// Maximum triangle area.
    pub max_area: f32,
    /// Mean triangle area.
    pub mean_area: f32,
    /// Minimum aspect ratio.
    pub min_aspect_ratio: f32,
    /// Maximum aspect ratio.
    pub max_aspect_ratio: f32,
    /// Mean skewness (0 = perfect, 1 = degenerate).
    pub mean_skewness: f32,
    /// Number of degenerate triangles (area < 1e-10).
    pub n_degenerate: u32,
}
/// Automated quality improvement suggestion with priority.
#[derive(Debug, Clone)]
pub struct AutoImproveSuggestion {
    /// Human-readable description.
    pub description: String,
    /// Priority (1 = highest, 5 = lowest).
    pub priority: u8,
    /// Estimated number of affected elements.
    pub affected: usize,
}
impl AutoImproveSuggestion {
    /// Create a new suggestion.
    pub fn new(description: impl Into<String>, priority: u8, affected: usize) -> Self {
        Self {
            description: description.into(),
            priority,
            affected,
        }
    }
}
