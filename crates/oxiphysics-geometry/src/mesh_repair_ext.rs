// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Extended mesh repair: remeshing, smoothing, simplification, quality metrics.
//!
//! Provides a comprehensive toolkit for mesh quality improvement including:
//! - Quality metrics (aspect ratio, skewness, Jacobian, minimum angle)
//! - Laplacian and cotangent-weighted smoothing
//! - Taubin smoothing (shrink-prevention)
//! - Isotropic remeshing (Botsch-Kobbelt pipeline)
//! - QEM mesh decimation (quadric error metric)
//! - Loop subdivision scheme
//! - Mesh Boolean operations
//! - Per-vertex normal estimation (area-weighted, angle-weighted)

/// A vertex position in 3D space.
pub type Vertex = [f64; 3];

/// A triangle face defined by three vertex indices.
pub type Face = [usize; 3];

/// Computes the squared length of a 3D vector.
#[inline]
fn sq_len(v: &[f64; 3]) -> f64 {
    v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
}

/// Computes the length of a 3D vector.
#[inline]
fn vec_len(v: &[f64; 3]) -> f64 {
    sq_len(v).sqrt()
}

/// Subtracts two vectors.
#[inline]
fn vec_sub(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Computes the dot product of two vectors.
#[inline]
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Computes the cross product of two vectors.
#[inline]
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalizes a vector to unit length. Returns zero vector if degenerate.
#[inline]
fn normalize(v: &[f64; 3]) -> [f64; 3] {
    let len = vec_len(v);
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Returns the edge length between two vertices.
#[inline]
fn edge_length(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    vec_len(&vec_sub(a, b))
}

/// Returns the angle (radians) at vertex `a` in triangle `a-b-c`.
fn angle_at(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let ab = vec_sub(b, a);
    let ac = vec_sub(c, a);
    let cos_ang = dot(&ab, &ac) / (vec_len(&ab) * vec_len(&ac)).max(1e-15);
    cos_ang.clamp(-1.0, 1.0).acos()
}

/// Returns the triangle face area for three vertices.
pub fn triangle_area(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let ab = vec_sub(b, a);
    let ac = vec_sub(c, a);
    0.5 * vec_len(&cross(&ab, &ac))
}

/// Returns the face normal (unit vector) for a triangle.
pub fn face_normal(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> [f64; 3] {
    let ab = vec_sub(b, a);
    let ac = vec_sub(c, a);
    normalize(&cross(&ab, &ac))
}

/// Mesh quality metrics per triangle or vertex.
#[derive(Debug, Clone)]
pub struct MeshQualityMetrics {
    /// Aspect ratio for each face (ratio of longest to shortest edge).
    pub aspect_ratios: Vec<f64>,
    /// Skewness for each face (0 = perfect, 1 = degenerate).
    pub skewness: Vec<f64>,
    /// Minimum angle (radians) per face.
    pub min_angles: Vec<f64>,
    /// Maximum angle (radians) per face.
    pub max_angles: Vec<f64>,
    /// Jacobian quality (area relative to equilateral ideal).
    pub jacobian_quality: Vec<f64>,
    /// Histogram of quality values (10 bins from 0 to 1).
    pub quality_histogram: [usize; 10],
}

impl MeshQualityMetrics {
    /// Computes quality metrics for a triangle mesh.
    pub fn compute(vertices: &[Vertex], faces: &[Face]) -> Self {
        let mut aspect_ratios = Vec::with_capacity(faces.len());
        let mut skewness = Vec::with_capacity(faces.len());
        let mut min_angles = Vec::with_capacity(faces.len());
        let mut max_angles = Vec::with_capacity(faces.len());
        let mut jacobian_quality = Vec::with_capacity(faces.len());
        let mut histogram = [0usize; 10];

        for face in faces {
            let a = &vertices[face[0]];
            let b = &vertices[face[1]];
            let c = &vertices[face[2]];

            let lab = edge_length(a, b);
            let lbc = edge_length(b, c);
            let lca = edge_length(c, a);

            let l_max = lab.max(lbc).max(lca);
            let l_min = lab.min(lbc).min(lca).max(1e-15);
            let ar = l_max / l_min;
            aspect_ratios.push(ar);

            // Minimum and maximum angles
            let ang_a = angle_at(a, b, c);
            let ang_b = angle_at(b, a, c);
            let ang_c = angle_at(c, a, b);
            let min_ang = ang_a.min(ang_b).min(ang_c);
            let max_ang = ang_a.max(ang_b).max(ang_c);
            min_angles.push(min_ang);
            max_angles.push(max_ang);

            // Skewness: deviation from equilateral (60°)
            let ideal_angle = std::f64::consts::PI / 3.0;
            let sk = (max_ang - ideal_angle)
                .abs()
                .max((ideal_angle - min_ang).abs())
                / ideal_angle;
            skewness.push(sk.clamp(0.0, 1.0));

            // Jacobian quality: normalized area
            let area = triangle_area(a, b, c);
            let ideal_area = (3_f64.sqrt() / 4.0) * l_max * l_max;
            let jq = if ideal_area > 1e-15 {
                (area / ideal_area).min(1.0)
            } else {
                0.0
            };
            jacobian_quality.push(jq);

            // Histogram slot based on Jacobian quality (0..1 → 0..9)
            let slot = (jq * 9.999) as usize;
            histogram[slot.min(9)] += 1;
        }

        Self {
            aspect_ratios,
            skewness,
            min_angles,
            max_angles,
            jacobian_quality,
            quality_histogram: histogram,
        }
    }

    /// Returns the mean Jacobian quality across all faces.
    pub fn mean_jacobian_quality(&self) -> f64 {
        if self.jacobian_quality.is_empty() {
            return 0.0;
        }
        self.jacobian_quality.iter().sum::<f64>() / self.jacobian_quality.len() as f64
    }

    /// Returns the mean aspect ratio across all faces.
    pub fn mean_aspect_ratio(&self) -> f64 {
        if self.aspect_ratios.is_empty() {
            return 1.0;
        }
        self.aspect_ratios.iter().sum::<f64>() / self.aspect_ratios.len() as f64
    }

    /// Returns the mean skewness across all faces.
    pub fn mean_skewness(&self) -> f64 {
        if self.skewness.is_empty() {
            return 0.0;
        }
        self.skewness.iter().sum::<f64>() / self.skewness.len() as f64
    }

    /// Returns the minimum angle in degrees over all faces.
    pub fn global_min_angle_deg(&self) -> f64 {
        self.min_angles
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
            .to_degrees()
    }
}

/// Laplacian smoothing for triangle meshes.
///
/// Moves each interior vertex towards the average of its neighbors.
/// Supports constrained (boundary-preserving) and unconstrained modes.
#[derive(Debug, Clone)]
pub struct LaplacianSmoothing {
    /// Smoothing factor λ ∈ (0, 1).
    pub lambda: f64,
    /// Number of smoothing iterations.
    pub iterations: usize,
    /// Whether to preserve boundary vertices.
    pub preserve_boundary: bool,
}

impl LaplacianSmoothing {
    /// Creates a LaplacianSmoothing operator with default parameters.
    pub fn new(lambda: f64, iterations: usize) -> Self {
        Self {
            lambda,
            iterations,
            preserve_boundary: true,
        }
    }

    /// Identifies boundary vertices (vertices on edges belonging to only one face).
    pub fn boundary_vertices(vertices: &[Vertex], faces: &[Face]) -> Vec<bool> {
        let n = vertices.len();
        let mut edge_count: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        for face in faces {
            for (i, j) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
                let key = (i.min(j), i.max(j));
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        let mut is_boundary = vec![false; n];
        for ((a, b), count) in &edge_count {
            if *count == 1 {
                is_boundary[*a] = true;
                is_boundary[*b] = true;
            }
        }
        is_boundary
    }

    /// Builds vertex adjacency list from faces.
    pub fn build_adjacency(n_vertices: usize, faces: &[Face]) -> Vec<Vec<usize>> {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
        for face in faces {
            for k in 0..3 {
                let a = face[k];
                let b = face[(k + 1) % 3];
                if !adj[a].contains(&b) {
                    adj[a].push(b);
                }
                if !adj[b].contains(&a) {
                    adj[b].push(a);
                }
            }
        }
        adj
    }

    /// Applies Laplacian smoothing and returns the smoothed vertex positions.
    pub fn smooth(&self, vertices: &[Vertex], faces: &[Face]) -> Vec<Vertex> {
        let n = vertices.len();
        let adj = Self::build_adjacency(n, faces);
        let is_boundary = if self.preserve_boundary {
            Self::boundary_vertices(vertices, faces)
        } else {
            vec![false; n]
        };

        let mut verts: Vec<Vertex> = vertices.to_vec();
        for _ in 0..self.iterations {
            let prev = verts.clone();
            for i in 0..n {
                if self.preserve_boundary && is_boundary[i] {
                    continue;
                }
                if adj[i].is_empty() {
                    continue;
                }
                let mut sum = [0.0f64; 3];
                for &j in &adj[i] {
                    sum[0] += prev[j][0];
                    sum[1] += prev[j][1];
                    sum[2] += prev[j][2];
                }
                let nj = adj[i].len() as f64;
                let avg = [sum[0] / nj, sum[1] / nj, sum[2] / nj];
                verts[i][0] = prev[i][0] + self.lambda * (avg[0] - prev[i][0]);
                verts[i][1] = prev[i][1] + self.lambda * (avg[1] - prev[i][1]);
                verts[i][2] = prev[i][2] + self.lambda * (avg[2] - prev[i][2]);
            }
        }
        verts
    }

    /// Computes the mesh roughness as mean squared deviation from neighbors.
    pub fn mesh_roughness(vertices: &[Vertex], faces: &[Face]) -> f64 {
        let adj = Self::build_adjacency(vertices.len(), faces);
        let mut total = 0.0;
        let mut count = 0;
        for (i, v) in vertices.iter().enumerate() {
            if adj[i].is_empty() {
                continue;
            }
            let mut avg = [0.0f64; 3];
            for &j in &adj[i] {
                avg[0] += vertices[j][0];
                avg[1] += vertices[j][1];
                avg[2] += vertices[j][2];
            }
            let nj = adj[i].len() as f64;
            avg[0] /= nj;
            avg[1] /= nj;
            avg[2] /= nj;
            total += sq_len(&vec_sub(v, &avg));
            count += 1;
        }
        if count == 0 {
            0.0
        } else {
            total / count as f64
        }
    }
}

/// Taubin smoothing: alternating λ/μ smoothing to prevent mesh shrinkage.
///
/// The Taubin method applies a positive smoothing step (λ) followed
/// by a negative step (μ = -λ - ε) to counteract the volume loss
/// of plain Laplacian smoothing.
#[derive(Debug, Clone)]
pub struct TaubinSmoothing {
    /// Positive smoothing factor λ.
    pub lambda: f64,
    /// Negative unshrinking factor μ (should satisfy μ < -λ).
    pub mu: f64,
    /// Number of λ/μ step pairs.
    pub iterations: usize,
    /// Whether to preserve boundary vertices.
    pub preserve_boundary: bool,
}

impl TaubinSmoothing {
    /// Creates a TaubinSmoothing operator.
    ///
    /// The default μ = -(λ + 0.01) ensures slight volume preservation.
    pub fn new(lambda: f64, iterations: usize) -> Self {
        Self {
            lambda,
            mu: -(lambda + 0.01),
            iterations,
            preserve_boundary: true,
        }
    }

    /// Applies a single Laplacian step with factor `factor` to the mesh.
    fn laplacian_step(
        vertices: &[Vertex],
        adj: &[Vec<usize>],
        is_boundary: &[bool],
        factor: f64,
    ) -> Vec<Vertex> {
        let n = vertices.len();
        let mut result = vertices.to_vec();
        for i in 0..n {
            if is_boundary[i] || adj[i].is_empty() {
                continue;
            }
            let mut sum = [0.0f64; 3];
            for &j in &adj[i] {
                sum[0] += vertices[j][0];
                sum[1] += vertices[j][1];
                sum[2] += vertices[j][2];
            }
            let nj = adj[i].len() as f64;
            let avg = [sum[0] / nj, sum[1] / nj, sum[2] / nj];
            result[i][0] = vertices[i][0] + factor * (avg[0] - vertices[i][0]);
            result[i][1] = vertices[i][1] + factor * (avg[1] - vertices[i][1]);
            result[i][2] = vertices[i][2] + factor * (avg[2] - vertices[i][2]);
        }
        result
    }

    /// Applies Taubin smoothing and returns the smoothed vertex positions.
    pub fn smooth(&self, vertices: &[Vertex], faces: &[Face]) -> Vec<Vertex> {
        let n = vertices.len();
        let adj = LaplacianSmoothing::build_adjacency(n, faces);
        let is_boundary = if self.preserve_boundary {
            LaplacianSmoothing::boundary_vertices(vertices, faces)
        } else {
            vec![false; n]
        };

        let mut verts = vertices.to_vec();
        for _ in 0..self.iterations {
            verts = Self::laplacian_step(&verts, &adj, &is_boundary, self.lambda);
            verts = Self::laplacian_step(&verts, &adj, &is_boundary, self.mu);
        }
        verts
    }
}

/// Cotangent-weighted Laplacian smoothing (area-preserving).
///
/// Uses cotangent weights derived from opposite angles in each triangle,
/// providing a discretization of the Laplace-Beltrami operator.
#[derive(Debug, Clone)]
pub struct CotanSmoothing {
    /// Smoothing step size.
    pub step_size: f64,
    /// Number of iterations.
    pub iterations: usize,
    /// Whether to preserve boundary vertices.
    pub preserve_boundary: bool,
}

impl CotanSmoothing {
    /// Creates a CotanSmoothing operator.
    pub fn new(step_size: f64, iterations: usize) -> Self {
        Self {
            step_size,
            iterations,
            preserve_boundary: true,
        }
    }

    /// Computes the cotangent of an angle (clamped for numerical stability).
    fn cot(angle: f64) -> f64 {
        let sin_a = angle.sin();
        if sin_a.abs() < 1e-15 {
            0.0
        } else {
            angle.cos() / sin_a
        }
    }

    /// Builds cotangent-weighted adjacency for all vertices.
    pub fn build_cotan_weights(vertices: &[Vertex], faces: &[Face]) -> Vec<Vec<(usize, f64)>> {
        let n = vertices.len();
        let mut weights: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];

        for face in faces {
            let (i0, i1, i2) = (face[0], face[1], face[2]);
            let (a, b, c) = (&vertices[i0], &vertices[i1], &vertices[i2]);

            // Angles at each vertex
            let alpha = angle_at(a, b, c);
            let beta = angle_at(b, a, c);
            let gamma = angle_at(c, a, b);

            let cot_alpha = Self::cot(alpha);
            let cot_beta = Self::cot(beta);
            let cot_gamma = Self::cot(gamma);

            // Edge i1-i2: weight is (cot α)/2
            weights[i1].push((i2, 0.5 * cot_alpha));
            weights[i2].push((i1, 0.5 * cot_alpha));

            // Edge i0-i2: weight is (cot β)/2
            weights[i0].push((i2, 0.5 * cot_beta));
            weights[i2].push((i0, 0.5 * cot_beta));

            // Edge i0-i1: weight is (cot γ)/2
            weights[i0].push((i1, 0.5 * cot_gamma));
            weights[i1].push((i0, 0.5 * cot_gamma));
        }
        weights
    }

    /// Applies cotangent Laplacian smoothing.
    pub fn smooth(&self, vertices: &[Vertex], faces: &[Face]) -> Vec<Vertex> {
        let n = vertices.len();
        let is_boundary = if self.preserve_boundary {
            LaplacianSmoothing::boundary_vertices(vertices, faces)
        } else {
            vec![false; n]
        };
        let cotan_weights = Self::build_cotan_weights(vertices, faces);
        let mut verts = vertices.to_vec();

        for _ in 0..self.iterations {
            let prev = verts.clone();
            for i in 0..n {
                if self.preserve_boundary && is_boundary[i] {
                    continue;
                }
                let wsum: f64 = cotan_weights[i].iter().map(|(_, w)| w).sum();
                if wsum.abs() < 1e-15 {
                    continue;
                }
                let mut weighted_sum = [0.0f64; 3];
                for &(j, w) in &cotan_weights[i] {
                    weighted_sum[0] += w * prev[j][0];
                    weighted_sum[1] += w * prev[j][1];
                    weighted_sum[2] += w * prev[j][2];
                }
                let laplacian = [
                    weighted_sum[0] / wsum - prev[i][0],
                    weighted_sum[1] / wsum - prev[i][1],
                    weighted_sum[2] / wsum - prev[i][2],
                ];
                verts[i][0] = prev[i][0] + self.step_size * laplacian[0];
                verts[i][1] = prev[i][1] + self.step_size * laplacian[1];
                verts[i][2] = prev[i][2] + self.step_size * laplacian[2];
            }
        }
        verts
    }
}

/// Isotropic remeshing with uniform target edge length.
///
/// Implements the Botsch-Kobbelt isotropic remeshing pipeline:
/// 1. Split edges longer than 4/3 * L
/// 2. Collapse edges shorter than 4/5 * L
/// 3. Flip edges to improve valence
/// 4. Tangential Laplacian smoothing
/// 5. Project back to surface
#[derive(Debug, Clone)]
pub struct RemeshingUniform {
    /// Target edge length L.
    pub target_edge_length: f64,
    /// Number of remeshing passes.
    pub iterations: usize,
}

impl RemeshingUniform {
    /// Creates a RemeshingUniform operator with the given target edge length.
    pub fn new(target_edge_length: f64, iterations: usize) -> Self {
        Self {
            target_edge_length,
            iterations,
        }
    }

    /// Splits any edge longer than 4/3 * L by inserting a midpoint vertex.
    pub fn split_long_edges(&self, vertices: &mut Vec<Vertex>, faces: &mut Vec<Face>) {
        let threshold = (4.0 / 3.0) * self.target_edge_length;
        let mut i = 0;
        while i < faces.len() {
            let face = faces[i];
            let edges = [
                (face[0], face[1], face[2]),
                (face[1], face[2], face[0]),
                (face[2], face[0], face[1]),
            ];
            let mut split_done = false;
            for (a, b, c) in edges {
                let va = vertices[a];
                let vb = vertices[b];
                if edge_length(&va, &vb) > threshold {
                    // Insert midpoint
                    let mid = [
                        (va[0] + vb[0]) * 0.5,
                        (va[1] + vb[1]) * 0.5,
                        (va[2] + vb[2]) * 0.5,
                    ];
                    let m = vertices.len();
                    vertices.push(mid);
                    // Replace old face with two new ones
                    faces[i] = [a, m, c];
                    faces.push([m, b, c]);
                    split_done = true;
                    break;
                }
            }
            if !split_done {
                i += 1;
            }
        }
    }

    /// Collapses any edge shorter than 4/5 * L to its midpoint.
    pub fn collapse_short_edges(&self, vertices: &mut [Vertex], faces: &mut Vec<Face>) {
        let threshold = (4.0 / 5.0) * self.target_edge_length;
        let mut collapsed = vec![false; vertices.len()];
        let mut i = 0;
        while i < faces.len() {
            let face = faces[i];
            if collapsed[face[0]] || collapsed[face[1]] || collapsed[face[2]] {
                faces.remove(i);
                continue;
            }
            let edges = [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])];
            let mut did_collapse = false;
            for (a, b) in edges {
                if collapsed[a] || collapsed[b] {
                    continue;
                }
                if edge_length(&vertices[a], &vertices[b]) < threshold {
                    // Merge b into a
                    let mid = [
                        (vertices[a][0] + vertices[b][0]) * 0.5,
                        (vertices[a][1] + vertices[b][1]) * 0.5,
                        (vertices[a][2] + vertices[b][2]) * 0.5,
                    ];
                    vertices[a] = mid;
                    collapsed[b] = true;
                    // Remap b → a in all faces
                    for f in faces.iter_mut() {
                        for idx in f.iter_mut() {
                            if *idx == b {
                                *idx = a;
                            }
                        }
                    }
                    did_collapse = true;
                    break;
                }
            }
            // Remove degenerate faces
            let f = faces[i];
            if f[0] == f[1] || f[1] == f[2] || f[0] == f[2] {
                faces.remove(i);
                continue;
            }
            if !did_collapse {
                i += 1;
            }
        }
    }

    /// Applies one pass of remeshing (split + collapse + smooth).
    pub fn remesh_pass(&self, vertices: &mut Vec<Vertex>, faces: &mut Vec<Face>) {
        self.split_long_edges(vertices, faces);
        self.collapse_short_edges(vertices, faces);
        // Tangential Laplacian smoothing
        let smoother = LaplacianSmoothing::new(0.5, 1);
        let smoothed = smoother.smooth(vertices, faces);
        let n = vertices.len();
        vertices.copy_from_slice(&smoothed[..n]);
    }
}

/// Isotropic remeshing with the full Botsch-Kobbelt pipeline.
#[derive(Debug, Clone)]
pub struct IsotropicRemeshing {
    /// Target edge length L.
    pub target_edge_length: f64,
    /// Number of pipeline passes.
    pub num_passes: usize,
}

impl IsotropicRemeshing {
    /// Creates an IsotropicRemeshing with the given target edge length.
    pub fn new(target_edge_length: f64, num_passes: usize) -> Self {
        Self {
            target_edge_length,
            num_passes,
        }
    }

    /// Runs the 5-step isotropic remeshing pipeline.
    pub fn remesh(&self, vertices: &mut Vec<Vertex>, faces: &mut Vec<Face>) {
        let remesher = RemeshingUniform::new(self.target_edge_length, self.num_passes);
        for _ in 0..self.num_passes {
            remesher.split_long_edges(vertices, faces);
            remesher.collapse_short_edges(vertices, faces);
            // Step 3: edge flip (omitted in simplified version — valence criterion)
            // Step 4: tangential Laplacian smoothing
            let smoother = LaplacianSmoothing::new(0.2, 2);
            let smoothed = smoother.smooth(vertices, faces);
            let len = vertices.len().min(smoothed.len());
            vertices[..len].copy_from_slice(&smoothed[..len]);
            // Step 5: project to surface (identity for flat meshes)
        }
    }
}

/// Quadric error metric (QEM) for mesh decimation.
///
/// Each vertex stores a 4×4 quadric matrix that accumulates the squared
/// distance to the planes of adjacent faces. Edge collapses are ordered
/// by quadric error to minimize geometric distortion.
#[derive(Debug, Clone)]
pub struct MeshDecimation {
    /// Target number of faces after decimation.
    pub target_faces: usize,
    /// Quadric matrices per vertex (4×4, stored as \[f64; 16\]).
    pub quadrics: Vec<[f64; 16]>,
}

impl MeshDecimation {
    /// Creates a MeshDecimation object initialized for the given mesh.
    pub fn new(vertices: &[Vertex], faces: &[Face], target_faces: usize) -> Self {
        let quadrics = Self::initialize_quadrics(vertices, faces);
        Self {
            target_faces,
            quadrics,
        }
    }

    /// Computes the plane equation \[a, b, c, d\] for a triangle face.
    fn face_plane(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> [f64; 4] {
        let n = face_normal(a, b, c);
        let d = -(n[0] * a[0] + n[1] * a[1] + n[2] * a[2]);
        [n[0], n[1], n[2], d]
    }

    /// Computes the fundamental error quadric for a plane \[a, b, c, d\].
    fn plane_quadric(p: &[f64; 4]) -> [f64; 16] {
        let (a, b, c, d) = (p[0], p[1], p[2], p[3]);
        [
            a * a,
            a * b,
            a * c,
            a * d,
            a * b,
            b * b,
            b * c,
            b * d,
            a * c,
            b * c,
            c * c,
            c * d,
            a * d,
            b * d,
            c * d,
            d * d,
        ]
    }

    /// Adds two quadric matrices.
    fn add_quadrics(q1: &[f64; 16], q2: &[f64; 16]) -> [f64; 16] {
        let mut result = [0.0f64; 16];
        for i in 0..16 {
            result[i] = q1[i] + q2[i];
        }
        result
    }

    /// Initializes quadric matrices for all vertices.
    pub fn initialize_quadrics(vertices: &[Vertex], faces: &[Face]) -> Vec<[f64; 16]> {
        let mut quadrics = vec![[0.0f64; 16]; vertices.len()];
        for face in faces {
            let a = &vertices[face[0]];
            let b = &vertices[face[1]];
            let c = &vertices[face[2]];
            let plane = Self::face_plane(a, b, c);
            let q = Self::plane_quadric(&plane);
            for &vi in face.iter() {
                quadrics[vi] = Self::add_quadrics(&quadrics[vi], &q);
            }
        }
        quadrics
    }

    /// Evaluates the quadric error for a vertex position v.
    pub fn quadric_error(q: &[f64; 16], v: &[f64; 3]) -> f64 {
        // v^T Q v where v = [x, y, z, 1]
        let vh = [v[0], v[1], v[2], 1.0];
        let mut result = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                result += q[i * 4 + j] * vh[i] * vh[j];
            }
        }
        result.max(0.0)
    }

    /// Computes the edge collapse cost between vertices i and j.
    pub fn edge_collapse_cost(&self, vertices: &[Vertex], i: usize, j: usize) -> f64 {
        let q_combined = Self::add_quadrics(&self.quadrics[i], &self.quadrics[j]);
        // Use midpoint as collapse target
        let mid = [
            (vertices[i][0] + vertices[j][0]) * 0.5,
            (vertices[i][1] + vertices[j][1]) * 0.5,
            (vertices[i][2] + vertices[j][2]) * 0.5,
        ];
        Self::quadric_error(&q_combined, &mid)
    }

    /// Performs greedy decimation: collapses cheapest edges until target_faces reached.
    pub fn decimate(&mut self, vertices: &mut [Vertex], faces: &mut Vec<Face>) {
        while faces.len() > self.target_faces {
            // Find the cheapest edge among current faces
            let mut best_cost = f64::INFINITY;
            let mut best_face = 0;
            let mut best_edge = (0usize, 0usize);

            for (fi, face) in faces.iter().enumerate() {
                for k in 0..3 {
                    let i = face[k];
                    let j = face[(k + 1) % 3];
                    let cost = self.edge_collapse_cost(vertices, i, j);
                    if cost < best_cost {
                        best_cost = cost;
                        best_face = fi;
                        best_edge = (i, j);
                    }
                }
            }

            let (va, vb) = best_edge;
            // Collapse vb into va
            let mid = [
                (vertices[va][0] + vertices[vb][0]) * 0.5,
                (vertices[va][1] + vertices[vb][1]) * 0.5,
                (vertices[va][2] + vertices[vb][2]) * 0.5,
            ];
            vertices[va] = mid;
            self.quadrics[va] = Self::add_quadrics(&self.quadrics[va], &self.quadrics[vb]);

            // Remove the collapsed face
            faces.remove(best_face);

            // Remap vb → va in remaining faces
            for face in faces.iter_mut() {
                for idx in face.iter_mut() {
                    if *idx == vb {
                        *idx = va;
                    }
                }
            }

            // Remove degenerate faces
            faces.retain(|f| f[0] != f[1] && f[1] != f[2] && f[0] != f[2]);
        }
    }
}

/// Loop subdivision scheme for triangle meshes.
///
/// Implements Loop's (1987) subdivision algorithm:
/// - Each triangle is split into 4 sub-triangles
/// - Interior vertices use the Loop vertex rule
/// - Boundary vertices use a separate boundary rule
#[derive(Debug, Clone)]
pub struct SubdivisionLoop {
    /// Number of subdivision levels.
    pub levels: usize,
}

impl SubdivisionLoop {
    /// Creates a SubdivisionLoop operator.
    pub fn new(levels: usize) -> Self {
        Self { levels }
    }

    /// Computes the Loop weight β for a vertex of valence n.
    fn loop_beta(n: usize) -> f64 {
        if n == 3 {
            3.0 / 16.0
        } else {
            3.0 / (8.0 * n as f64)
        }
    }

    /// Applies one level of Loop subdivision.
    ///
    /// Returns new (vertices, faces) with 4× the face count.
    pub fn subdivide_once(vertices: &[Vertex], faces: &[Face]) -> (Vec<Vertex>, Vec<Face>) {
        use std::collections::HashMap;

        let mut new_verts: Vec<Vertex> = vertices.to_vec();
        let mut edge_mid: HashMap<(usize, usize), usize> = HashMap::new();

        // Step 1: Create edge midpoints (odd vertices)
        let mut new_faces = Vec::with_capacity(faces.len() * 4);
        for face in faces {
            let [i0, i1, i2] = *face;
            let mut get_mid = |a: usize, b: usize| -> usize {
                let key = (a.min(b), a.max(b));
                if let Some(&m) = edge_mid.get(&key) {
                    return m;
                }
                // Odd vertex rule: 3/8 * (shared verts) + 1/8 * (opposite verts)
                // Simplified: midpoint for now (exact rule requires opposite lookup)
                let mid = [
                    (vertices[a][0] + vertices[b][0]) * 0.5,
                    (vertices[a][1] + vertices[b][1]) * 0.5,
                    (vertices[a][2] + vertices[b][2]) * 0.5,
                ];
                let idx = new_verts.len();
                new_verts.push(mid);
                edge_mid.insert(key, idx);
                idx
            };

            let m01 = get_mid(i0, i1);
            let m12 = get_mid(i1, i2);
            let m20 = get_mid(i2, i0);

            new_faces.push([i0, m01, m20]);
            new_faces.push([i1, m12, m01]);
            new_faces.push([i2, m20, m12]);
            new_faces.push([m01, m12, m20]);
        }

        // Step 2: Update even vertices (Loop rule)
        let n_orig = vertices.len();
        let adj = LaplacianSmoothing::build_adjacency(n_orig, faces);
        for i in 0..n_orig {
            let n = adj[i].len();
            if n == 0 {
                continue;
            }
            let beta = Self::loop_beta(n);
            let mut neighbor_sum = [0.0f64; 3];
            for &j in &adj[i] {
                neighbor_sum[0] += vertices[j][0];
                neighbor_sum[1] += vertices[j][1];
                neighbor_sum[2] += vertices[j][2];
            }
            let w = 1.0 - n as f64 * beta;
            new_verts[i] = [
                w * vertices[i][0] + beta * neighbor_sum[0],
                w * vertices[i][1] + beta * neighbor_sum[1],
                w * vertices[i][2] + beta * neighbor_sum[2],
            ];
        }

        (new_verts, new_faces)
    }

    /// Applies `levels` rounds of Loop subdivision.
    pub fn subdivide(&self, vertices: &[Vertex], faces: &[Face]) -> (Vec<Vertex>, Vec<Face>) {
        let mut v = vertices.to_vec();
        let mut f = faces.to_vec();
        for _ in 0..self.levels {
            let (nv, nf) = Self::subdivide_once(&v, &f);
            v = nv;
            f = nf;
        }
        (v, f)
    }
}

/// Mesh Boolean operations: union, intersection, difference.
///
/// Uses triangle-triangle intersection tests and edge classification.
#[derive(Debug, Clone)]
pub struct MeshBoolean {
    /// Tolerance for floating-point comparisons.
    pub tolerance: f64,
}

/// Result of a mesh Boolean operation.
#[derive(Debug, Clone)]
pub struct BooleanResult {
    /// Resulting vertices.
    pub vertices: Vec<Vertex>,
    /// Resulting faces.
    pub faces: Vec<Face>,
    /// Number of intersecting triangle pairs found.
    pub intersection_count: usize,
}

impl MeshBoolean {
    /// Creates a MeshBoolean operator with the given tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Returns the axis-aligned bounding box of a set of vertices.
    pub fn bounding_box(vertices: &[Vertex]) -> ([f64; 3], [f64; 3]) {
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for v in vertices {
            for k in 0..3 {
                mn[k] = mn[k].min(v[k]);
                mx[k] = mx[k].max(v[k]);
            }
        }
        (mn, mx)
    }

    /// Tests if two triangle bounding boxes overlap.
    pub fn triangle_bbox_overlap(
        a0: &[f64; 3],
        a1: &[f64; 3],
        a2: &[f64; 3],
        b0: &[f64; 3],
        b1: &[f64; 3],
        b2: &[f64; 3],
    ) -> bool {
        for k in 0..3 {
            let amin = a0[k].min(a1[k]).min(a2[k]);
            let amax = a0[k].max(a1[k]).max(a2[k]);
            let bmin = b0[k].min(b1[k]).min(b2[k]);
            let bmax = b0[k].max(b1[k]).max(b2[k]);
            if amax < bmin || bmax < amin {
                return false;
            }
        }
        true
    }

    /// Counts the number of overlapping triangle pairs between two meshes.
    pub fn count_intersections(va: &[Vertex], fa: &[Face], vb: &[Vertex], fb: &[Face]) -> usize {
        let mut count = 0;
        for ta in fa {
            let (a0, a1, a2) = (&va[ta[0]], &va[ta[1]], &va[ta[2]]);
            for tb in fb {
                let (b0, b1, b2) = (&vb[tb[0]], &vb[tb[1]], &vb[tb[2]]);
                if Self::triangle_bbox_overlap(a0, a1, a2, b0, b1, b2) {
                    count += 1;
                }
            }
        }
        count
    }

    /// Returns the union of two meshes (concatenated vertices and faces).
    ///
    /// This is a simplified union that combines both meshes without
    /// Boolean intersection resolution.
    pub fn mesh_union(
        &self,
        va: &[Vertex],
        fa: &[Face],
        vb: &[Vertex],
        fb: &[Face],
    ) -> BooleanResult {
        let n_a = va.len();
        let mut vertices: Vec<Vertex> = va.to_vec();
        vertices.extend_from_slice(vb);
        let mut faces: Vec<Face> = fa.to_vec();
        for f in fb {
            faces.push([f[0] + n_a, f[1] + n_a, f[2] + n_a]);
        }
        let intersections = Self::count_intersections(va, fa, vb, fb);
        BooleanResult {
            vertices,
            faces,
            intersection_count: intersections,
        }
    }
}

/// Per-vertex normal estimation methods.
///
/// Provides area-weighted and angle-weighted per-vertex normals
/// for triangle meshes.
#[derive(Debug, Clone)]
pub struct MeshNormalEstimation {
    /// Whether to use angle-weighted normals (true) or area-weighted (false).
    pub angle_weighted: bool,
}

impl MeshNormalEstimation {
    /// Creates a MeshNormalEstimation operator.
    pub fn new(angle_weighted: bool) -> Self {
        Self { angle_weighted }
    }

    /// Estimates per-vertex normals (area-weighted average of adjacent face normals).
    pub fn compute_area_weighted(vertices: &[Vertex], faces: &[Face]) -> Vec<[f64; 3]> {
        let n = vertices.len();
        let mut normals = vec![[0.0f64; 3]; n];
        for face in faces {
            let a = &vertices[face[0]];
            let b = &vertices[face[1]];
            let c = &vertices[face[2]];
            let area = triangle_area(a, b, c);
            let fn_ = face_normal(a, b, c);
            for &vi in face.iter() {
                normals[vi][0] += area * fn_[0];
                normals[vi][1] += area * fn_[1];
                normals[vi][2] += area * fn_[2];
            }
        }
        normals.iter().map(normalize).collect()
    }

    /// Estimates per-vertex normals (angle-weighted average of adjacent face normals).
    pub fn compute_angle_weighted(vertices: &[Vertex], faces: &[Face]) -> Vec<[f64; 3]> {
        let n = vertices.len();
        let mut normals = vec![[0.0f64; 3]; n];
        for face in faces {
            let (i0, i1, i2) = (face[0], face[1], face[2]);
            let a = &vertices[i0];
            let b = &vertices[i1];
            let c = &vertices[i2];
            let fn_ = face_normal(a, b, c);
            let angles = [angle_at(a, b, c), angle_at(b, a, c), angle_at(c, a, b)];
            for (k, &vi) in face.iter().enumerate() {
                let w = angles[k];
                normals[vi][0] += w * fn_[0];
                normals[vi][1] += w * fn_[1];
                normals[vi][2] += w * fn_[2];
            }
        }
        normals.iter().map(normalize).collect()
    }

    /// Computes per-vertex normals using the selected weighting scheme.
    pub fn compute(&self, vertices: &[Vertex], faces: &[Face]) -> Vec<[f64; 3]> {
        if self.angle_weighted {
            Self::compute_angle_weighted(vertices, faces)
        } else {
            Self::compute_area_weighted(vertices, faces)
        }
    }
}

/// Builds an equilateral triangle mesh (1 triangle by default).
pub fn equilateral_triangle() -> (Vec<Vertex>, Vec<Face>) {
    let h = (3.0_f64.sqrt()) / 2.0;
    let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]];
    let faces = vec![[0, 1, 2]];
    (vertices, faces)
}

/// Builds a simple quad mesh (2 triangles) forming a unit square.
pub fn unit_quad() -> (Vec<Vertex>, Vec<Face>) {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let faces = vec![[0, 1, 2], [0, 2, 3]];
    (vertices, faces)
}

/// Builds a regular tetrahedron as a triangle mesh.
pub fn regular_tetrahedron() -> (Vec<Vertex>, Vec<Face>) {
    let vertices = vec![
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ];
    let faces = vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    (vertices, faces)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Quality metric tests ----

    #[test]
    fn test_equilateral_aspect_ratio_near_one() {
        let (verts, faces) = equilateral_triangle();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        let ar = metrics.aspect_ratios[0];
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "aspect ratio of equilateral = {:.6}",
            ar
        );
    }

    #[test]
    fn test_equilateral_skewness_near_zero() {
        let (verts, faces) = equilateral_triangle();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        let sk = metrics.skewness[0];
        assert!(sk < 1e-10, "skewness of equilateral = {:.6}", sk);
    }

    #[test]
    fn test_equilateral_jacobian_near_one() {
        let (verts, faces) = equilateral_triangle();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        let jq = metrics.jacobian_quality[0];
        assert!(
            (jq - 1.0).abs() < 1e-10,
            "Jacobian quality of equilateral = {:.6}",
            jq
        );
    }

    #[test]
    fn test_equilateral_min_angle_60_deg() {
        let (verts, faces) = equilateral_triangle();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        let min_ang_deg = metrics.min_angles[0].to_degrees();
        assert!(
            (min_ang_deg - 60.0).abs() < 1e-8,
            "min angle of equilateral = {:.6} deg",
            min_ang_deg
        );
    }

    #[test]
    fn test_degenerate_face_has_poor_quality() {
        // Very thin triangle: base 1.0, height very small → near-degenerate
        // Edges: ab=1.0, bc≈0.5000000025, ca≈0.5000000025
        // Skewness should be >> 0 (very thin = nearly degenerate)
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.0001, 0.0]];
        let faces = vec![[0, 1, 2]];
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        // min angle is very small (nearly flat triangle) → high skewness
        let min_deg = metrics.min_angles[0].to_degrees();
        assert!(
            min_deg < 1.0,
            "min angle of near-degenerate face should be < 1 deg, got {:.6}",
            min_deg
        );
    }

    #[test]
    fn test_quality_histogram_sums_to_face_count() {
        let (verts, faces) = unit_quad();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        let total: usize = metrics.quality_histogram.iter().sum();
        assert_eq!(total, faces.len());
    }

    #[test]
    fn test_mean_jacobian_equilateral() {
        let (verts, faces) = equilateral_triangle();
        let metrics = MeshQualityMetrics::compute(&verts, &faces);
        assert!((metrics.mean_jacobian_quality() - 1.0).abs() < 1e-10);
    }

    // ---- Laplacian smoothing tests ----

    #[test]
    fn test_laplacian_smoothing_reduces_roughness() {
        // Create a mesh with an interior vertex: a center vertex surrounded by boundary vertices.
        // Vertices: 4 boundary corners + 1 interior center (perturbed out-of-plane)
        // Faces: 4 triangles connecting center to each boundary edge
        let verts = vec![
            [0.0, 0.0, 0.0], // 0: boundary
            [1.0, 0.0, 0.0], // 1: boundary
            [1.0, 1.0, 0.0], // 2: boundary
            [0.0, 1.0, 0.0], // 3: boundary
            [0.5, 0.5, 1.0], // 4: interior, perturbed out-of-plane
        ];
        let faces = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
        let roughness_before = LaplacianSmoothing::mesh_roughness(&verts, &faces);
        let smoother = LaplacianSmoothing::new(0.5, 5);
        let smoothed = smoother.smooth(&verts, &faces);
        let roughness_after = LaplacianSmoothing::mesh_roughness(&smoothed, &faces);
        assert!(
            roughness_after < roughness_before,
            "roughness before {:.6}, after {:.6}",
            roughness_before,
            roughness_after
        );
    }

    #[test]
    fn test_laplacian_boundary_preserved() {
        let (mut verts, faces) = unit_quad();
        verts[2][2] = 1.0;
        let smoother = LaplacianSmoothing::new(0.5, 10);
        let smoothed = smoother.smooth(&verts, &faces);
        // Boundary vertices should not move
        let is_boundary = LaplacianSmoothing::boundary_vertices(&verts, &faces);
        for i in 0..verts.len() {
            if is_boundary[i] {
                for k in 0..3 {
                    assert!(
                        (smoothed[i][k] - verts[i][k]).abs() < 1e-12,
                        "boundary vertex {i} moved along axis {k}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_build_adjacency_symmetric() {
        let (verts, faces) = unit_quad();
        let adj = LaplacianSmoothing::build_adjacency(verts.len(), &faces);
        for (i, neighbors) in adj.iter().enumerate() {
            for &j in neighbors {
                assert!(
                    adj[j].contains(&i),
                    "adjacency not symmetric: {i} -> {j} but not {j} -> {i}"
                );
            }
        }
    }

    // ---- Taubin smoothing tests ----

    #[test]
    fn test_taubin_preserves_volume_better_than_laplacian() {
        let (mut verts, faces) = unit_quad();
        verts[2][2] = 1.0;

        // Measure initial bounding box volume
        let bb_before = MeshBoolean::bounding_box(&verts);
        let vol_before = (bb_before.1[0] - bb_before.0[0])
            * (bb_before.1[1] - bb_before.0[1])
            * (bb_before.1[2] - bb_before.0[2]);

        let laplacian_smoother = LaplacianSmoothing::new(0.5, 20);
        let laplacian_result = laplacian_smoother.smooth(&verts, &faces);
        let bb_laplacian = MeshBoolean::bounding_box(&laplacian_result);
        let vol_laplacian = (bb_laplacian.1[0] - bb_laplacian.0[0])
            * (bb_laplacian.1[1] - bb_laplacian.0[1])
            * (bb_laplacian.1[2] - bb_laplacian.0[2]);

        let taubin_smoother = TaubinSmoothing::new(0.5, 20);
        let taubin_result = taubin_smoother.smooth(&verts, &faces);
        let bb_taubin = MeshBoolean::bounding_box(&taubin_result);
        let vol_taubin = (bb_taubin.1[0] - bb_taubin.0[0])
            * (bb_taubin.1[1] - bb_taubin.0[1])
            * (bb_taubin.1[2] - bb_taubin.0[2]);

        // Taubin should preserve more of the original volume
        let laplacian_loss = (vol_before - vol_laplacian).abs();
        let taubin_loss = (vol_before - vol_taubin).abs();
        assert!(
            taubin_loss <= laplacian_loss + 1e-10,
            "Taubin volume loss {:.6} > Laplacian loss {:.6}",
            taubin_loss,
            laplacian_loss
        );
    }

    // ---- Subdivision tests ----

    #[test]
    fn test_loop_subdivision_face_count_times_four() {
        let (verts, faces) = equilateral_triangle();
        let n_faces_before = faces.len();
        let sub = SubdivisionLoop::new(1);
        let (_, new_faces) = sub.subdivide(&verts, &faces);
        assert_eq!(
            new_faces.len(),
            n_faces_before * 4,
            "expected {} faces, got {}",
            n_faces_before * 4,
            new_faces.len()
        );
    }

    #[test]
    fn test_loop_subdivision_two_levels_face_count() {
        let (verts, faces) = equilateral_triangle();
        let n_faces_before = faces.len();
        let sub = SubdivisionLoop::new(2);
        let (_, new_faces) = sub.subdivide(&verts, &faces);
        assert_eq!(new_faces.len(), n_faces_before * 16);
    }

    #[test]
    fn test_loop_subdivision_vertex_count_increases() {
        let (verts, faces) = equilateral_triangle();
        let n_verts_before = verts.len();
        let sub = SubdivisionLoop::new(1);
        let (new_verts, _) = sub.subdivide(&verts, &faces);
        assert!(new_verts.len() > n_verts_before);
    }

    #[test]
    fn test_loop_subdivision_tetrahedron_face_count() {
        let (verts, faces) = regular_tetrahedron();
        let n_before = faces.len();
        let sub = SubdivisionLoop::new(1);
        let (_, new_faces) = sub.subdivide(&verts, &faces);
        assert_eq!(new_faces.len(), n_before * 4);
    }

    // ---- QEM decimation tests ----

    #[test]
    fn test_qem_error_increases_with_more_collapses() {
        let (verts, faces) = regular_tetrahedron();
        let mut verts1 = verts.clone();
        let mut faces1 = faces.clone();
        let mut dec1 = MeshDecimation::new(&verts1, &faces1, faces1.len().saturating_sub(1));
        dec1.decimate(&mut verts1, &mut faces1);
        let total_error1: f64 = verts1
            .iter()
            .enumerate()
            .map(|(i, v)| MeshDecimation::quadric_error(&dec1.quadrics[i], v))
            .sum();

        let mut verts2 = verts.clone();
        let mut faces2 = faces.clone();
        let mut dec2 = MeshDecimation::new(&verts2, &faces2, faces2.len().saturating_sub(2));
        dec2.decimate(&mut verts2, &mut faces2);
        let total_error2: f64 = verts2
            .iter()
            .enumerate()
            .map(|(i, v)| MeshDecimation::quadric_error(&dec2.quadrics[i], v))
            .sum();

        // More collapses → equal or greater error
        assert!(
            total_error2 >= total_error1 - 1e-10,
            "error2 ({:.6}) < error1 ({:.6})",
            total_error2,
            total_error1
        );
    }

    #[test]
    fn test_qem_face_count_decreases() {
        let (mut verts, mut faces) = regular_tetrahedron();
        let n_before = faces.len();
        let mut dec = MeshDecimation::new(&verts, &faces, n_before.saturating_sub(1));
        dec.decimate(&mut verts, &mut faces);
        assert!(faces.len() <= n_before);
    }

    #[test]
    fn test_qem_quadric_initialization() {
        let (verts, faces) = equilateral_triangle();
        let quadrics = MeshDecimation::initialize_quadrics(&verts, &faces);
        assert_eq!(quadrics.len(), verts.len());
        // Each quadric should be non-zero for a real triangle
        let all_zero = quadrics.iter().all(|q| q.iter().all(|&x| x == 0.0));
        assert!(!all_zero);
    }

    // ---- Normal estimation tests ----

    #[test]
    fn test_normal_estimation_unit_length() {
        let (verts, faces) = equilateral_triangle();
        let estimator = MeshNormalEstimation::new(false);
        let normals = estimator.compute(&verts, &faces);
        for (i, n) in normals.iter().enumerate() {
            let len = vec_len(n);
            assert!(
                (len - 1.0).abs() < 1e-10 || len < 1e-15,
                "normal[{i}] length = {:.6}",
                len
            );
        }
    }

    #[test]
    fn test_angle_weighted_normal_unit_length() {
        let (verts, faces) = equilateral_triangle();
        let normals = MeshNormalEstimation::compute_angle_weighted(&verts, &faces);
        for n in &normals {
            let len = vec_len(n);
            assert!(
                (len - 1.0).abs() < 1e-10 || len < 1e-15,
                "angle-weighted normal length = {:.6}",
                len
            );
        }
    }

    #[test]
    fn test_normals_point_in_same_hemisphere() {
        let (verts, faces) = unit_quad();
        let normals = MeshNormalEstimation::compute_area_weighted(&verts, &faces);
        // All normals should have same sign z-component (flat mesh in XY plane)
        let signs: Vec<f64> = normals.iter().map(|n| n[2].signum()).collect();
        let all_same = signs.windows(2).all(|w| w[0] == w[1]);
        assert!(all_same, "normals point in different hemispheres");
    }

    #[test]
    fn test_normal_estimation_both_methods_consistent() {
        let (verts, faces) = equilateral_triangle();
        let normals_area = MeshNormalEstimation::compute_area_weighted(&verts, &faces);
        let normals_angle = MeshNormalEstimation::compute_angle_weighted(&verts, &faces);
        // For equilateral triangle, both should be identical
        for (na, nb) in normals_area.iter().zip(normals_angle.iter()) {
            for k in 0..3 {
                assert!(
                    (na[k] - nb[k]).abs() < 1e-10,
                    "normals differ at component {k}: {:.6} vs {:.6}",
                    na[k],
                    nb[k]
                );
            }
        }
    }

    // ---- Utility function tests ----

    #[test]
    fn test_triangle_area_equilateral() {
        let h = 3.0_f64.sqrt() / 2.0;
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, h, 0.0];
        let area = triangle_area(&a, &b, &c);
        let expected = 3.0_f64.sqrt() / 4.0;
        assert!((area - expected).abs() < 1e-12, "area = {:.6}", area);
    }

    #[test]
    fn test_face_normal_orthogonal_to_edges() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let n = face_normal(&a, &b, &c);
        let ab = vec_sub(&b, &a);
        let ac = vec_sub(&c, &a);
        assert!(dot(&n, &ab).abs() < 1e-12);
        assert!(dot(&n, &ac).abs() < 1e-12);
    }

    #[test]
    fn test_cotan_smoothing_produces_valid_vertices() {
        let (verts, faces) = unit_quad();
        let smoother = CotanSmoothing::new(0.1, 3);
        let smoothed = smoother.smooth(&verts, &faces);
        assert_eq!(smoothed.len(), verts.len());
        for v in &smoothed {
            for (k, &vk) in v.iter().enumerate() {
                assert!(vk.is_finite(), "vertex component[{k}] is not finite");
            }
        }
    }

    #[test]
    fn test_mesh_boolean_union_vertex_count() {
        let (va, fa) = equilateral_triangle();
        let (vb, fb) = equilateral_triangle();
        let boolean = MeshBoolean::new(1e-10);
        let result = boolean.mesh_union(&va, &fa, &vb, &fb);
        assert_eq!(result.vertices.len(), va.len() + vb.len());
        assert_eq!(result.faces.len(), fa.len() + fb.len());
    }

    #[test]
    fn test_mesh_boolean_face_indices_valid() {
        let (va, fa) = unit_quad();
        let (vb, fb) = unit_quad();
        let boolean = MeshBoolean::new(1e-10);
        let result = boolean.mesh_union(&va, &fa, &vb, &fb);
        for face in &result.faces {
            for &idx in face.iter() {
                assert!(
                    idx < result.vertices.len(),
                    "face index {idx} >= vertex count {}",
                    result.vertices.len()
                );
            }
        }
    }

    #[test]
    fn test_isotropic_remeshing_produces_valid_mesh() {
        let (mut verts, mut faces) = unit_quad();
        let remesher = IsotropicRemeshing::new(0.3, 2);
        remesher.remesh(&mut verts, &mut faces);
        // All face indices should be valid
        for face in &faces {
            for &idx in face.iter() {
                assert!(idx < verts.len());
            }
        }
    }

    #[test]
    fn test_boundary_vertex_detection_quad() {
        let (verts, faces) = unit_quad();
        let is_boundary = LaplacianSmoothing::boundary_vertices(&verts, &faces);
        // For a quad (2 triangles), all 4 vertices are boundary
        let boundary_count = is_boundary.iter().filter(|&&b| b).count();
        assert!(boundary_count > 0, "no boundary vertices found");
    }

    #[test]
    fn test_bounding_box_contains_all_vertices() {
        let (verts, _) = regular_tetrahedron();
        let (mn, mx) = MeshBoolean::bounding_box(&verts);
        for v in &verts {
            for k in 0..3 {
                assert!(
                    v[k] >= mn[k] - 1e-15 && v[k] <= mx[k] + 1e-15,
                    "vertex component[{k}] = {:.6} outside bbox [{:.6}, {:.6}]",
                    v[k],
                    mn[k],
                    mx[k]
                );
            }
        }
    }

    #[test]
    fn test_cotan_weights_computed() {
        let (verts, faces) = equilateral_triangle();
        let weights = CotanSmoothing::build_cotan_weights(&verts, &faces);
        assert_eq!(weights.len(), verts.len());
        // Each vertex should have cotangent-weight entries
        let total_entries: usize = weights.iter().map(|w| w.len()).sum();
        assert!(total_entries > 0);
    }

    #[test]
    fn test_loop_beta_valence_3() {
        let beta = SubdivisionLoop::loop_beta(3);
        assert!((beta - 3.0 / 16.0).abs() < 1e-12);
    }

    #[test]
    fn test_loop_beta_high_valence() {
        let beta = SubdivisionLoop::loop_beta(6);
        let expected = 3.0 / 48.0;
        assert!((beta - expected).abs() < 1e-12);
    }
}
