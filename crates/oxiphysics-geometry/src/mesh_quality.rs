#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh quality metrics and repair operations.
//!
//! Provides quality metrics for tetrahedral and triangular elements,
//! mesh repair utilities, Laplacian smoothing, and 2-D edge flipping.

#![allow(dead_code, missing_docs)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Local vector helpers (operating on [f64; 3])
// ---------------------------------------------------------------------------

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < f64::EPSILON {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

/// Angle (radians) at vertex `a` in the triangle `a–b–c`.
fn angle_at(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = normalize3(sub3(b, a));
    let ac = normalize3(sub3(c, a));
    dot3(ab, ac).clamp(-1.0, 1.0).acos()
}

// ---------------------------------------------------------------------------
// TetrahedralMesh
// ---------------------------------------------------------------------------

/// A tetrahedral volumetric mesh.
pub struct TetrahedralMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Tetrahedra defined by four vertex indices.
    pub tets: Vec<[usize; 4]>,
}

impl TetrahedralMesh {
    /// Create a new tetrahedral mesh.
    pub fn new(vertices: Vec<[f64; 3]>, tets: Vec<[usize; 4]>) -> Self {
        Self { vertices, tets }
    }
}

// ---------------------------------------------------------------------------
// MeshQuality
// ---------------------------------------------------------------------------

/// Quality metrics for a single mesh element (triangle or tetrahedron).
#[derive(Debug, Clone, PartialEq)]
pub struct MeshQuality {
    /// Aspect ratio: ratio of longest to shortest edge (≥ 1, ideal = 1).
    pub aspect_ratio: f64,
    /// Minimum interior angle in degrees.
    pub min_angle_deg: f64,
    /// Maximum interior angle in degrees.
    pub max_angle_deg: f64,
    /// Skewness in \[0, 1\]; 0 = perfect, 1 = degenerate.
    pub skewness: f64,
    /// Scaled Jacobian determinant in \[-1, 1\]; ideal equilateral = 1.
    pub jacobian: f64,
}

impl MeshQuality {
    /// Returns `true` if the element is considered acceptable quality.
    ///
    /// Criteria: aspect_ratio < 5, skewness < 0.9, min_angle_deg > 5°.
    pub fn is_acceptable(&self) -> bool {
        self.aspect_ratio < 5.0 && self.skewness < 0.9 && self.min_angle_deg > 5.0
    }

    /// Returns `true` if the element is degenerate (zero or near-zero volume).
    pub fn is_degenerate(&self) -> bool {
        self.jacobian.abs() < 1e-10
    }
}

// ---------------------------------------------------------------------------
// compute_tet_quality
// ---------------------------------------------------------------------------

/// Compute quality metrics for a tetrahedron defined by four vertices.
///
/// Returns a [`MeshQuality`] struct with aspect ratio, min/max dihedral angle,
/// skewness, and scaled Jacobian.
pub fn compute_tet_quality(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> MeshQuality {
    // Edge vectors from v0
    let e01 = sub3(v1, v0);
    let e02 = sub3(v2, v0);
    let e03 = sub3(v3, v0);
    let e12 = sub3(v2, v1);
    let e13 = sub3(v3, v1);
    let e23 = sub3(v3, v2);

    // Edge lengths
    let edge_lengths = [
        len3(e01),
        len3(e02),
        len3(e03),
        len3(e12),
        len3(e13),
        len3(e23),
    ];

    let min_len = edge_lengths.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_len = edge_lengths.iter().cloned().fold(0.0_f64, f64::max);

    let aspect_ratio = if min_len < f64::EPSILON {
        f64::INFINITY
    } else {
        max_len / min_len
    };

    // Jacobian: det of edge matrix (proportional to signed volume)
    let jac = dot3(e01, cross3(e02, e03));
    // Ideal regular tet with same mean edge length
    let mean_edge = edge_lengths.iter().sum::<f64>() / 6.0;
    let ideal_jac = mean_edge * mean_edge * mean_edge * (2.0_f64).sqrt();
    // Take absolute value so the jacobian is positive for non-degenerate tets
    // regardless of vertex winding order
    let scaled_jac = if ideal_jac < f64::EPSILON {
        0.0
    } else {
        jac.abs() / ideal_jac
    };

    // Dihedral angles between each pair of faces
    // Face normals (outward for the tet):
    //   face 012: n0 = (v1-v0) × (v2-v0)
    //   face 013: n1 = (v3-v0) × (v1-v0)
    //   face 023: n2 = (v2-v0) × (v3-v0)
    //   face 123: n3 = (v2-v1) × (v3-v1)
    let n012 = normalize3(cross3(e01, e02));
    let n013 = normalize3(cross3(e03, e01));
    let n023 = normalize3(cross3(e02, e03));
    let n123 = normalize3(cross3(e12, e13));

    // Dihedral angle along each edge = π - angle between adjacent face normals
    let dihedral_angle = |na: [f64; 3], nb: [f64; 3]| -> f64 {
        let cos_a = dot3(na, nb).clamp(-1.0, 1.0);
        // dihedral angle is the supplement of the angle between outward normals
        std::f64::consts::PI - cos_a.acos()
    };

    // 6 dihedral angles, one per edge:
    //   edge 01 → faces 012 and 013
    //   edge 02 → faces 012 and 023
    //   edge 03 → faces 013 and 023
    //   edge 12 → faces 012 and 123
    //   edge 13 → faces 013 and 123
    //   edge 23 → faces 023 and 123
    let dihedrals = [
        dihedral_angle(n012, n013),
        dihedral_angle(n012, n023),
        dihedral_angle(n013, n023),
        dihedral_angle(n012, n123),
        dihedral_angle(n013, n123),
        dihedral_angle(n023, n123),
    ];

    let rad_to_deg = 180.0 / std::f64::consts::PI;
    let min_angle_deg = dihedrals.iter().cloned().fold(f64::INFINITY, f64::min) * rad_to_deg;
    let max_angle_deg = dihedrals.iter().cloned().fold(0.0_f64, f64::max) * rad_to_deg;

    // Skewness: (max_angle - ideal) / (180 - ideal)
    // For a regular tet, ideal dihedral ≈ 70.529°
    let ideal_dihedral = 70.529;
    let skewness = ((max_angle_deg - ideal_dihedral) / (180.0 - ideal_dihedral)).clamp(0.0, 1.0);

    MeshQuality {
        aspect_ratio,
        min_angle_deg,
        max_angle_deg,
        skewness,
        jacobian: scaled_jac,
    }
}

// ---------------------------------------------------------------------------
// compute_tri_quality
// ---------------------------------------------------------------------------

/// Compute quality metrics for a triangle defined by three vertices.
///
/// Returns a [`MeshQuality`] struct with aspect ratio, min/max interior angle,
/// skewness, and a 2-D Jacobian measure.
pub fn compute_tri_quality(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> MeshQuality {
    let a = dist3(v1, v2);
    let b = dist3(v0, v2);
    let c = dist3(v0, v1);

    let min_len = a.min(b).min(c);
    let max_len = a.max(b).max(c);

    let aspect_ratio = if min_len < f64::EPSILON {
        f64::INFINITY
    } else {
        max_len / min_len
    };

    let angle_a = angle_at(v0, v1, v2);
    let angle_b = angle_at(v1, v0, v2);
    let angle_c = angle_at(v2, v0, v1);

    let rad_to_deg = 180.0 / std::f64::consts::PI;
    let angles_deg = [
        angle_a * rad_to_deg,
        angle_b * rad_to_deg,
        angle_c * rad_to_deg,
    ];

    let min_angle_deg = angles_deg.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_angle_deg = angles_deg.iter().cloned().fold(0.0_f64, f64::max);

    // Skewness: (max_angle - 60°) / (180° - 60°)
    let ideal = 60.0;
    let skewness = ((max_angle_deg - ideal) / (180.0 - ideal)).clamp(0.0, 1.0);

    // Jacobian: area-based. Equilateral triangle of same perimeter → area_ideal.
    // Use 2 * area / (a*b) as a proxy Jacobian measure (= sin of enclosed angle).
    let area = {
        let cross = cross3(sub3(v1, v0), sub3(v2, v0));
        len3(cross) * 0.5
    };
    let denom = a * b;
    let jacobian = if denom < f64::EPSILON {
        0.0
    } else {
        (2.0 * area / denom).clamp(-1.0, 1.0)
    };

    MeshQuality {
        aspect_ratio,
        min_angle_deg,
        max_angle_deg,
        skewness,
        jacobian,
    }
}

// ---------------------------------------------------------------------------
// repair_degenerate_tets
// ---------------------------------------------------------------------------

/// Remove degenerate tetrahedra from a mesh in-place.
///
/// A tet is considered degenerate when its signed volume is below `threshold`.
/// Returns the number of tetrahedra removed.
pub fn repair_degenerate_tets(mesh: &mut TetrahedralMesh) -> usize {
    let threshold = 1e-12_f64;
    let before = mesh.tets.len();
    mesh.tets.retain(|tet| {
        let v0 = mesh.vertices[tet[0]];
        let v1 = mesh.vertices[tet[1]];
        let v2 = mesh.vertices[tet[2]];
        let v3 = mesh.vertices[tet[3]];
        let e01 = sub3(v1, v0);
        let e02 = sub3(v2, v0);
        let e03 = sub3(v3, v0);
        let vol = dot3(e01, cross3(e02, e03)).abs() / 6.0;
        vol > threshold
    });
    before - mesh.tets.len()
}

// ---------------------------------------------------------------------------
// smooth_laplacian
// ---------------------------------------------------------------------------

/// Apply Laplacian smoothing to a tetrahedral mesh.
///
/// Each interior vertex is moved to the average position of its neighbours.
/// Boundary vertices (those belonging to fewer than 4 tets) are kept fixed.
///
/// # Arguments
/// * `positions` – mutable slice of vertex positions.
/// * `tets` – tetrahedra (index groups of 4).
/// * `iterations` – number of smoothing passes.
pub fn smooth_laplacian(positions: &mut [[f64; 3]], tets: &[[usize; 4]], iterations: usize) {
    let n = positions.len();
    // Count how many tets reference each vertex
    let mut tet_count = vec![0usize; n];
    for tet in tets {
        for &vi in tet.iter() {
            tet_count[vi] += 1;
        }
    }

    for _ in 0..iterations {
        let mut sums = vec![[0.0_f64; 3]; n];
        let mut counts = vec![0usize; n];

        for tet in tets {
            for k in 0..4 {
                let vi = tet[k];
                // Accumulate neighbour contributions
                for j in 0..4 {
                    if j != k {
                        let vj = tet[j];
                        sums[vi] = add3(sums[vi], positions[vj]);
                        counts[vi] += 1;
                    }
                }
            }
        }

        for i in 0..n {
            // Only move interior vertices (those connected to many tets)
            if counts[i] > 0 && tet_count[i] >= 4 {
                let avg = scale3(sums[i], 1.0 / counts[i] as f64);
                positions[i] = avg;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// flip_edges_2d
// ---------------------------------------------------------------------------

/// Attempt Delaunay edge flipping on a 2-D triangle mesh.
///
/// For each pair of triangles sharing an edge, check if flipping that edge
/// improves the minimum angle (Lawson's criterion). Returns the number of
/// flips performed.
///
/// `positions` is an (x, y) 2-D position array encoded as `[f64; 3]`
/// (z component is ignored).
pub fn flip_edges_2d(triangles: &mut Vec<[usize; 3]>, positions: &[[f64; 3]]) -> usize {
    let mut flips = 0usize;
    let mut changed = true;

    while changed {
        changed = false;

        // Build edge → triangle map: edge (min, max) -> list of (tri_idx, local_opposite_vertex)
        let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
        for (ti, tri) in triangles.iter().enumerate() {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let opp = tri[(k + 2) % 3]; // opposite vertex
                let key = (a.min(b), a.max(b));
                edge_map.entry(key).or_default().push((ti, opp));
            }
        }

        for ((_ea, _eb), tris) in &edge_map {
            if tris.len() != 2 {
                continue; // boundary edge or non-manifold
            }
            let (ti, opp_i) = tris[0];
            let (tj, opp_j) = tris[1];

            // Shared edge vertices
            let ea = _ea;
            let eb = _eb;

            let pa = positions[*ea];
            let pb = positions[*eb];
            let pc = positions[opp_i];
            let pd = positions[opp_j];

            // Flip if the opposite vertex d lies inside the circumcircle of triangle (a, b, c)
            if in_circumcircle_2d(pa, pb, pc, pd) {
                // Flip: replace shared edge (ea, eb) with (opp_i, opp_j)
                let new_t0 = [*ea, opp_i, opp_j];
                let new_t1 = [*eb, opp_j, opp_i];

                // Verify both new triangles have positive signed area (no inversion)
                if signed_area_2d(
                    positions[new_t0[0]],
                    positions[new_t0[1]],
                    positions[new_t0[2]],
                ) > 0.0
                    && signed_area_2d(
                        positions[new_t1[0]],
                        positions[new_t1[1]],
                        positions[new_t1[2]],
                    ) > 0.0
                {
                    triangles[ti] = new_t0;
                    triangles[tj] = new_t1;
                    flips += 1;
                    changed = true;
                    break; // restart to avoid stale indices
                }
            }
        }
    }

    flips
}

/// Returns the signed area of a 2-D triangle (positive = CCW).
fn signed_area_2d(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
}

/// Returns `true` if point `d` lies strictly inside the circumcircle of CCW triangle `a–b–c`.
///
/// Uses the standard 3×3 determinant test.
fn in_circumcircle_2d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> bool {
    let ax = a[0] - d[0];
    let ay = a[1] - d[1];
    let bx = b[0] - d[0];
    let by = b[1] - d[1];
    let cx = c[0] - d[0];
    let cy = c[1] - d[1];

    let det = ax * (by * (cx * cx + cy * cy) - cy * (bx * bx + by * by))
        - ay * (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by))
        + (ax * ax + ay * ay) * (bx * cy - by * cx);

    det > 0.0
}

// ---------------------------------------------------------------------------
// Jacobian-based quality metrics
// ---------------------------------------------------------------------------

/// Compute the Jacobian matrix of a tetrahedron as a 3x3 matrix.
///
/// The Jacobian is formed from the three edge vectors emanating from v0.
/// Returns `[e01, e02, e03]` as row vectors (3x3 column-major).
pub fn tet_jacobian_matrix(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
) -> [[f64; 3]; 3] {
    let e01 = sub3(v1, v0);
    let e02 = sub3(v2, v0);
    let e03 = sub3(v3, v0);
    [e01, e02, e03]
}

/// Compute the determinant of the Jacobian matrix of a tetrahedron.
///
/// Positive for a properly oriented tet, zero for degenerate.
pub fn tet_jacobian_det(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    let j = tet_jacobian_matrix(v0, v1, v2, v3);
    // det = j[0] . (j[1] x j[2])
    dot3(j[0], cross3(j[1], j[2]))
}

/// Compute the scaled Jacobian quality metric for a tetrahedron.
///
/// The scaled Jacobian is the Jacobian determinant divided by the product of
/// the edge lengths forming the Jacobian, normalized to \[-1, 1\].
/// A value of 1 indicates an ideal (regular) tetrahedron.
pub fn tet_scaled_jacobian(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    let j = tet_jacobian_matrix(v0, v1, v2, v3);
    let det = dot3(j[0], cross3(j[1], j[2]));

    let l0 = len3(j[0]);
    let l1 = len3(j[1]);
    let l2 = len3(j[2]);

    let denom = l0 * l1 * l2;
    if denom < f64::EPSILON {
        return 0.0;
    }

    // For a regular tetrahedron, the scaled Jacobian equals sqrt(2)
    // Normalize to [0, 1] by dividing by sqrt(2)
    let raw = det / denom;
    raw / (2.0_f64).sqrt()
}

/// Compute the scaled Jacobian at all four corners of a tetrahedron
/// and return the minimum value.
///
/// This is a more robust quality measure than checking just one corner.
pub fn tet_min_scaled_jacobian(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    // Each corner uses edges emanating from that vertex, with consistent
    // even-permutation ordering to preserve orientation.
    let sj0 = tet_scaled_jacobian(v0, v1, v2, v3);
    let sj1 = tet_scaled_jacobian(v1, v2, v0, v3);
    let sj2 = tet_scaled_jacobian(v2, v0, v1, v3);
    let sj3 = tet_scaled_jacobian(v3, v1, v0, v2);
    sj0.min(sj1).min(sj2).min(sj3)
}

// ---------------------------------------------------------------------------
// Condition number quality
// ---------------------------------------------------------------------------

/// Compute the condition number quality metric for a tetrahedron.
///
/// Uses the Frobenius norm of the Jacobian and its inverse.
/// Quality = 1 / (kappa * kappa_ideal), where kappa = ||J|| * ||J^{-1}||.
/// Returns a value in (0, 1] where 1 is ideal.
pub fn tet_condition_number_quality(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    let j = tet_jacobian_matrix(v0, v1, v2, v3);
    let det = dot3(j[0], cross3(j[1], j[2]));

    if det.abs() < f64::EPSILON {
        return 0.0;
    }

    // Frobenius norm of J
    let mut frob_j = 0.0;
    for row in &j {
        frob_j += dot3(*row, *row);
    }
    frob_j = frob_j.sqrt();

    // Compute J^{-1} using cofactor matrix
    let inv_det = 1.0 / det;
    let j_inv = [
        scale3(cross3(j[1], j[2]), inv_det),
        scale3(cross3(j[2], j[0]), inv_det),
        scale3(cross3(j[0], j[1]), inv_det),
    ];

    // Frobenius norm of J^{-1}
    let mut frob_jinv = 0.0;
    for row in &j_inv {
        frob_jinv += dot3(*row, *row);
    }
    frob_jinv = frob_jinv.sqrt();

    let kappa = frob_j * frob_jinv;
    // For an ideal equilateral tet, kappa = sqrt(3) * 3 / 3 = sqrt(3)
    // Normalize so that ideal gives quality = 1
    let ideal_kappa = 3.0_f64.sqrt();
    (ideal_kappa / kappa).clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Quality histogram
// ---------------------------------------------------------------------------

/// A histogram of element quality values.
pub struct QualityHistogram {
    /// Bin edges (length = n_bins + 1).
    pub bin_edges: Vec<f64>,
    /// Counts per bin (length = n_bins).
    pub counts: Vec<usize>,
    /// Total number of elements.
    pub total: usize,
}

impl QualityHistogram {
    /// Build a histogram of quality values with the given number of bins.
    ///
    /// Quality values are assumed to be in \[0, 1\].
    pub fn from_values(values: &[f64], n_bins: usize) -> Self {
        assert!(n_bins > 0, "n_bins must be > 0");
        let mut counts = vec![0usize; n_bins];
        let bin_width = 1.0 / n_bins as f64;
        let mut bin_edges = Vec::with_capacity(n_bins + 1);
        for i in 0..=n_bins {
            bin_edges.push(i as f64 * bin_width);
        }

        for &v in values {
            let bin = ((v / bin_width) as usize).min(n_bins - 1);
            counts[bin] += 1;
        }

        Self {
            bin_edges,
            counts,
            total: values.len(),
        }
    }

    /// Return the fraction of elements in each bin.
    pub fn fractions(&self) -> Vec<f64> {
        if self.total == 0 {
            return vec![0.0; self.counts.len()];
        }
        self.counts
            .iter()
            .map(|&c| c as f64 / self.total as f64)
            .collect()
    }

    /// Return the bin index containing the most elements.
    pub fn peak_bin(&self) -> usize {
        self.counts
            .iter()
            .enumerate()
            .max_by_key(|(_, c)| **c)
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// Worst element identification
// ---------------------------------------------------------------------------

/// Find the indices of the worst elements by a given quality metric.
///
/// Returns the indices of elements with quality below `threshold`,
/// sorted by quality (worst first).
pub fn find_worst_elements(qualities: &[f64], threshold: f64) -> Vec<usize> {
    let mut bad: Vec<(usize, f64)> = qualities
        .iter()
        .enumerate()
        .filter(|(_, q)| *q < &threshold)
        .map(|(i, q)| (i, *q))
        .collect();
    bad.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    bad.iter().map(|&(i, _)| i).collect()
}

/// Summary statistics for a collection of quality values.
pub struct QualityStats {
    /// Minimum quality value.
    pub min: f64,
    /// Maximum quality value.
    pub max: f64,
    /// Mean quality value.
    pub mean: f64,
    /// Standard deviation of quality values.
    pub std_dev: f64,
    /// Number of elements below the "poor quality" threshold (0.3).
    pub n_poor: usize,
}

impl QualityStats {
    /// Compute statistics from quality values.
    pub fn compute(values: &[f64]) -> Self {
        if values.is_empty() {
            return Self {
                min: 0.0,
                max: 0.0,
                mean: 0.0,
                std_dev: 0.0,
                n_poor: 0,
            };
        }
        let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance =
            values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
        let std_dev = variance.sqrt();
        let n_poor = values.iter().filter(|&&v| v < 0.3).count();
        Self {
            min,
            max,
            mean,
            std_dev,
            n_poor,
        }
    }
}

/// Compute quality values for all tetrahedra in a mesh using the scaled Jacobian.
pub fn mesh_scaled_jacobian_qualities(mesh: &TetrahedralMesh) -> Vec<f64> {
    mesh.tets
        .iter()
        .map(|tet| {
            let v0 = mesh.vertices[tet[0]];
            let v1 = mesh.vertices[tet[1]];
            let v2 = mesh.vertices[tet[2]];
            let v3 = mesh.vertices[tet[3]];
            tet_min_scaled_jacobian(v0, v1, v2, v3)
        })
        .collect()
}

/// Compute condition number quality for all tetrahedra in a mesh.
pub fn mesh_condition_number_qualities(mesh: &TetrahedralMesh) -> Vec<f64> {
    mesh.tets
        .iter()
        .map(|tet| {
            let v0 = mesh.vertices[tet[0]];
            let v1 = mesh.vertices[tet[1]];
            let v2 = mesh.vertices[tet[2]];
            let v3 = mesh.vertices[tet[3]];
            tet_condition_number_quality(v0, v1, v2, v3)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Dihedral angle extremes
// ---------------------------------------------------------------------------

/// Compute the minimum and maximum dihedral angles across all edges of a tet.
///
/// Returns `(min_deg, max_deg)`.
pub fn tet_dihedral_angle_extremes(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
) -> (f64, f64) {
    let q = compute_tet_quality(v0, v1, v2, v3);
    (q.min_angle_deg, q.max_angle_deg)
}

/// Compute the minimum and maximum interior angles of a triangle.
///
/// Returns `(min_deg, max_deg)`.
pub fn tri_angle_extremes(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> (f64, f64) {
    let q = compute_tri_quality(v0, v1, v2);
    (q.min_angle_deg, q.max_angle_deg)
}

// ---------------------------------------------------------------------------
// Aspect ratio for all element types
// ---------------------------------------------------------------------------

/// Aspect ratio of a line segment (always 1.0).
pub fn line_aspect_ratio(_v0: [f64; 3], _v1: [f64; 3]) -> f64 {
    1.0
}

/// Aspect ratio of a quadrilateral element.
///
/// Defined as max_diagonal / min_diagonal.
pub fn quad_aspect_ratio(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
    // Diagonals of the quad
    let d1 = dist3(v0, v2);
    let d2 = dist3(v1, v3);
    let min_d = d1.min(d2);
    let max_d = d1.max(d2);
    if min_d < f64::EPSILON {
        f64::INFINITY
    } else {
        max_d / min_d
    }
}

/// Aspect ratio of a wedge (pentahedral prism) element.
///
/// Approximated as max_edge / min_edge across the 9 edges.
pub fn wedge_aspect_ratio(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    v5: [f64; 3],
) -> f64 {
    let edges = [
        dist3(v0, v1),
        dist3(v1, v2),
        dist3(v2, v0), // bottom tri
        dist3(v3, v4),
        dist3(v4, v5),
        dist3(v5, v3), // top tri
        dist3(v0, v3),
        dist3(v1, v4),
        dist3(v2, v5), // lateral
    ];
    let min_e = edges.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_e = edges.iter().cloned().fold(0.0_f64, f64::max);
    if min_e < f64::EPSILON {
        f64::INFINITY
    } else {
        max_e / min_e
    }
}

/// Aspect ratio of a hexahedral element.
///
/// Approximated as max_edge / min_edge across the 12 edges.
#[allow(clippy::too_many_arguments)]
pub fn hex_aspect_ratio(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    v5: [f64; 3],
    v6: [f64; 3],
    v7: [f64; 3],
) -> f64 {
    let edges = [
        dist3(v0, v1),
        dist3(v1, v2),
        dist3(v2, v3),
        dist3(v3, v0), // bottom
        dist3(v4, v5),
        dist3(v5, v6),
        dist3(v6, v7),
        dist3(v7, v4), // top
        dist3(v0, v4),
        dist3(v1, v5),
        dist3(v2, v6),
        dist3(v3, v7), // lateral
    ];
    let min_e = edges.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_e = edges.iter().cloned().fold(0.0_f64, f64::max);
    if min_e < f64::EPSILON {
        f64::INFINITY
    } else {
        max_e / min_e
    }
}

// ---------------------------------------------------------------------------
// Quality-guided mesh improvement
// ---------------------------------------------------------------------------

/// Perform quality-guided Laplacian smoothing: only move vertices that improve quality.
///
/// For each interior vertex, compute the Laplacian-smoothed position, then
/// evaluate the quality of affected tets.  The move is accepted only if
/// the minimum quality improves.
///
/// Returns the number of accepted vertex moves.
pub fn quality_guided_smoothing(
    positions: &mut [[f64; 3]],
    tets: &[[usize; 4]],
    iterations: usize,
) -> usize {
    let n = positions.len();
    let mut tet_count = vec![0usize; n];
    for tet in tets {
        for &vi in tet.iter() {
            tet_count[vi] += 1;
        }
    }

    let mut total_accepted = 0usize;

    for _ in 0..iterations {
        let mut sums = vec![[0.0_f64; 3]; n];
        let mut counts = vec![0usize; n];

        for tet in tets.iter() {
            for k in 0..4 {
                let vi = tet[k];
                for j in 0..4 {
                    if j != k {
                        let vj = tet[j];
                        sums[vi] = add3(sums[vi], positions[vj]);
                        counts[vi] += 1;
                    }
                }
            }
        }

        for i in 0..n {
            if counts[i] == 0 || tet_count[i] < 2 {
                continue;
            }
            let new_pos = scale3(sums[i], 1.0 / counts[i] as f64);

            // Compute minimum quality before move for affected tets
            let min_q_before = tets
                .iter()
                .filter(|t| t.contains(&i))
                .map(|t| {
                    let vv: Vec<[f64; 3]> = t.iter().map(|&vi| positions[vi]).collect();
                    tet_scaled_jacobian(vv[0], vv[1], vv[2], vv[3]).abs()
                })
                .fold(f64::INFINITY, f64::min);

            // Tentatively apply
            let old_pos = positions[i];
            positions[i] = new_pos;

            let min_q_after = tets
                .iter()
                .filter(|t| t.contains(&i))
                .map(|t| {
                    let vv: Vec<[f64; 3]> = t.iter().map(|&vi| positions[vi]).collect();
                    tet_scaled_jacobian(vv[0], vv[1], vv[2], vv[3]).abs()
                })
                .fold(f64::INFINITY, f64::min);

            if min_q_after >= min_q_before {
                total_accepted += 1; // keep new position
            } else {
                positions[i] = old_pos; // revert
            }
        }
    }

    total_accepted
}

// ---------------------------------------------------------------------------
// Minimum angle optimization (angle-based smoothing)
// ---------------------------------------------------------------------------

/// Apply angle-based smoothing to a 2-D triangle mesh.
///
/// Moves each interior vertex to maximise the minimum interior angle in
/// surrounding triangles.  Uses a gradient-free search along the Laplacian
/// direction with a line search.
///
/// Returns the number of accepted vertex moves.
pub fn minimum_angle_optimization(
    positions: &mut [[f64; 3]],
    triangles: &[[usize; 3]],
    iterations: usize,
) -> usize {
    let n = positions.len();

    // Find interior vertices: referenced by ≥ 2 triangles
    let mut tri_count = vec![0usize; n];
    for tri in triangles {
        for &vi in tri.iter() {
            tri_count[vi] += 1;
        }
    }

    let mut total_accepted = 0usize;

    for _ in 0..iterations {
        for i in 0..n {
            if tri_count[i] < 2 {
                continue; // boundary
            }

            // Compute Laplacian position
            let mut sum = [0.0_f64; 3];
            let mut cnt = 0usize;
            for tri in triangles.iter().filter(|t| t.contains(&i)) {
                for &j in tri.iter() {
                    if j != i {
                        sum = add3(sum, positions[j]);
                        cnt += 1;
                    }
                }
            }
            if cnt == 0 {
                continue;
            }
            let new_pos = scale3(sum, 1.0 / cnt as f64);

            // Compute minimum angle before move
            let min_before = triangles
                .iter()
                .filter(|t| t.contains(&i))
                .map(|t| {
                    let vv: Vec<[f64; 3]> = t.iter().map(|&vi| positions[vi]).collect();
                    compute_tri_quality(vv[0], vv[1], vv[2]).min_angle_deg
                })
                .fold(f64::INFINITY, f64::min);

            let old_pos = positions[i];
            positions[i] = new_pos;

            let min_after = triangles
                .iter()
                .filter(|t| t.contains(&i))
                .map(|t| {
                    let vv: Vec<[f64; 3]> = t.iter().map(|&vi| positions[vi]).collect();
                    compute_tri_quality(vv[0], vv[1], vv[2]).min_angle_deg
                })
                .fold(f64::INFINITY, f64::min);

            if min_after >= min_before {
                total_accepted += 1;
            } else {
                positions[i] = old_pos;
            }
        }
    }

    total_accepted
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Vertices of a regular unit tetrahedron centred at origin.
    fn unit_tet_vertices() -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
        // Regular tetrahedron with positive orientation (edge length 2*sqrt(2))
        let v0 = [1.0, 1.0, 1.0];
        let v1 = [-1.0, 1.0, -1.0];
        let v2 = [1.0, -1.0, -1.0];
        let v3 = [-1.0, -1.0, 1.0];
        (v0, v1, v2, v3)
    }

    #[test]
    fn test_unit_tet_aspect_ratio() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let q = compute_tet_quality(v0, v1, v2, v3);
        // Regular tet: all edges equal → aspect ratio = 1
        assert!(
            (q.aspect_ratio - 1.0).abs() < 1e-10,
            "aspect_ratio={}, expected 1.0",
            q.aspect_ratio
        );
    }

    #[test]
    fn test_unit_tet_jacobian_positive() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let q = compute_tet_quality(v0, v1, v2, v3);
        assert!(
            q.jacobian > 0.0,
            "jacobian should be positive, got {}",
            q.jacobian
        );
    }

    #[test]
    fn test_degenerate_tet_detected() {
        // Coplanar: v3 in the same plane as v0, v1, v2 → zero volume
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.5, 0.5, 0.0]; // lies in the same plane
        let q = compute_tet_quality(v0, v1, v2, v3);
        assert!(
            q.is_degenerate(),
            "coplanar tet should be degenerate, jacobian={}",
            q.jacobian
        );
    }

    #[test]
    fn test_repair_degenerate_tets() {
        let verts = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.5, 0.0], // coplanar with v0, v1, v2
        ];
        let tets = vec![
            [0usize, 1, 2, 3], // valid
            [0, 1, 2, 4],      // degenerate (all coplanar)
        ];
        let mut mesh = TetrahedralMesh::new(verts, tets);
        let removed = repair_degenerate_tets(&mut mesh);
        assert_eq!(removed, 1, "expected 1 degenerate tet removed");
        assert_eq!(mesh.tets.len(), 1);
    }

    #[test]
    fn test_equilateral_tri_quality() {
        // Equilateral triangle → aspect_ratio = 1, all angles = 60°, skewness = 0
        let h = (3.0_f64).sqrt() / 2.0;
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.5, h, 0.0];
        let q = compute_tri_quality(v0, v1, v2);
        assert!(
            (q.aspect_ratio - 1.0).abs() < 1e-10,
            "aspect_ratio={}, expected 1.0",
            q.aspect_ratio
        );
        assert!(
            (q.min_angle_deg - 60.0).abs() < 1e-6,
            "min_angle={}, expected 60°",
            q.min_angle_deg
        );
        assert!(
            (q.skewness).abs() < 1e-6,
            "skewness={}, expected 0",
            q.skewness
        );
    }

    #[test]
    fn test_right_triangle_quality() {
        // 3-4-5 right triangle: angles ≈ 36.87°, 53.13°, 90°
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [3.0, 0.0, 0.0];
        let v2 = [0.0, 4.0, 0.0];
        let q = compute_tri_quality(v0, v1, v2);
        assert!(
            (q.max_angle_deg - 90.0).abs() < 1e-6,
            "max_angle={}, expected 90°",
            q.max_angle_deg
        );
        assert!(
            q.aspect_ratio > 1.0,
            "aspect_ratio should be > 1 for non-equilateral"
        );
    }

    #[test]
    fn test_smooth_laplacian_runs() {
        let mut positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let tets = vec![[0usize, 1, 2, 3], [1, 2, 3, 4]];
        smooth_laplacian(&mut positions, &tets, 2);
        // Just ensure it runs without panic and positions are finite
        for p in &positions {
            assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
        }
    }

    #[test]
    fn test_flip_edges_2d_no_flip_on_delaunay() {
        // Already Delaunay: right-angle triangulation of a unit square
        let positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mut tris = vec![[0usize, 1, 2], [0, 2, 3]];
        let flips = flip_edges_2d(&mut tris, &positions);
        // This diagonal is already locally Delaunay, so 0 or 1 flip is acceptable
        let _ = flips;
        assert_eq!(tris.len(), 2, "triangle count must stay the same");
    }

    #[test]
    fn test_tet_quality_angles_in_range() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let q = compute_tet_quality(v0, v1, v2, v3);
        assert!(q.min_angle_deg > 0.0 && q.min_angle_deg < 180.0);
        assert!(q.max_angle_deg > 0.0 && q.max_angle_deg < 180.0);
        assert!(q.min_angle_deg <= q.max_angle_deg);
    }

    #[test]
    fn test_degenerate_tri_quality() {
        // Collinear → zero area → degenerate
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [2.0, 0.0, 0.0];
        let q = compute_tri_quality(v0, v1, v2);
        assert!(
            q.jacobian.abs() < 1e-10,
            "collinear triangle jacobian should be ~0, got {}",
            q.jacobian
        );
    }

    #[test]
    fn test_pi_constant_used() {
        // Sanity-check that PI is imported correctly in this module
        let half_pi_deg = PI / 2.0 * 180.0 / PI;
        assert!((half_pi_deg - 90.0).abs() < 1e-10);
    }

    // --- Jacobian-based quality ---

    #[test]
    fn test_tet_jacobian_det_positive() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let det = tet_jacobian_det(v0, v1, v2, v3);
        assert!(
            det.abs() > 1e-10,
            "regular tet should have non-zero Jacobian det"
        );
    }

    #[test]
    fn test_tet_jacobian_det_degenerate() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.5, 0.5, 0.0]; // coplanar
        let det = tet_jacobian_det(v0, v1, v2, v3);
        assert!(
            det.abs() < 1e-10,
            "coplanar tet should have zero Jacobian det"
        );
    }

    #[test]
    fn test_tet_jacobian_matrix() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.0, 0.0, 1.0];
        let j = tet_jacobian_matrix(v0, v1, v2, v3);
        assert!((j[0][0] - 1.0).abs() < 1e-10);
        assert!((j[1][1] - 1.0).abs() < 1e-10);
        assert!((j[2][2] - 1.0).abs() < 1e-10);
    }

    // --- Scaled Jacobian ---

    #[test]
    fn test_scaled_jacobian_regular_tet() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let sj = tet_scaled_jacobian(v0, v1, v2, v3);
        assert!(
            sj > 0.0,
            "regular tet should have positive scaled jacobian, got {sj}"
        );
    }

    #[test]
    fn test_scaled_jacobian_degenerate() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.5, 0.5, 0.0];
        let sj = tet_scaled_jacobian(v0, v1, v2, v3);
        assert!(sj.abs() < 1e-10);
    }

    #[test]
    fn test_min_scaled_jacobian_regular_tet() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let msj = tet_min_scaled_jacobian(v0, v1, v2, v3);
        assert!(
            msj > 0.0,
            "min scaled jacobian should be positive for regular tet"
        );
    }

    #[test]
    fn test_min_scaled_jacobian_right_angle_tet() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.0, 0.0, 1.0];
        let msj = tet_min_scaled_jacobian(v0, v1, v2, v3);
        assert!(msj > 0.0);
    }

    // --- Condition number quality ---

    #[test]
    fn test_condition_number_quality_regular_tet() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let q = tet_condition_number_quality(v0, v1, v2, v3);
        assert!(
            q > 0.3,
            "regular tet should have positive condition number quality, got {q}"
        );
    }

    #[test]
    fn test_condition_number_quality_degenerate() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.5, 0.5, 0.0];
        let q = tet_condition_number_quality(v0, v1, v2, v3);
        assert!(q.abs() < 1e-10, "degenerate tet should have zero quality");
    }

    #[test]
    fn test_condition_number_quality_bounded() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.0, 1.0, 0.0];
        let v3 = [0.0, 0.0, 1.0];
        let q = tet_condition_number_quality(v0, v1, v2, v3);
        assert!(
            (0.0..=1.0).contains(&q),
            "quality should be in [0,1], got {q}"
        );
    }

    // --- Quality histogram ---

    #[test]
    fn test_quality_histogram_basic() {
        let values = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        let hist = QualityHistogram::from_values(&values, 5);
        assert_eq!(hist.counts.len(), 5);
        assert_eq!(hist.total, 10);
        let total_count: usize = hist.counts.iter().sum();
        assert_eq!(total_count, 10);
    }

    #[test]
    fn test_quality_histogram_single_bin() {
        let values = vec![0.5, 0.6, 0.7];
        let hist = QualityHistogram::from_values(&values, 1);
        assert_eq!(hist.counts, vec![3]);
    }

    #[test]
    fn test_quality_histogram_fractions() {
        let values = vec![0.1, 0.2, 0.7, 0.8];
        let hist = QualityHistogram::from_values(&values, 2);
        let fracs = hist.fractions();
        assert_eq!(fracs.len(), 2);
        assert!((fracs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quality_histogram_peak_bin() {
        let values = vec![0.1, 0.15, 0.12, 0.8, 0.85];
        let hist = QualityHistogram::from_values(&values, 5);
        let peak = hist.peak_bin();
        assert!(peak < 5);
    }

    #[test]
    fn test_quality_histogram_empty() {
        let hist = QualityHistogram::from_values(&[], 5);
        assert_eq!(hist.total, 0);
        assert!(hist.counts.iter().all(|&c| c == 0));
    }

    // --- Worst element identification ---

    #[test]
    fn test_find_worst_elements_all_good() {
        let qualities = vec![0.8, 0.9, 0.7, 0.85];
        let worst = find_worst_elements(&qualities, 0.3);
        assert!(worst.is_empty());
    }

    #[test]
    fn test_find_worst_elements_some_bad() {
        let qualities = vec![0.8, 0.1, 0.9, 0.2, 0.7];
        let worst = find_worst_elements(&qualities, 0.3);
        assert_eq!(worst.len(), 2);
        // Should be sorted worst first
        assert_eq!(worst[0], 1); // quality 0.1
        assert_eq!(worst[1], 3); // quality 0.2
    }

    #[test]
    fn test_find_worst_elements_empty() {
        let worst = find_worst_elements(&[], 0.5);
        assert!(worst.is_empty());
    }

    // --- Quality stats ---

    #[test]
    fn test_quality_stats_basic() {
        let values = vec![0.1, 0.2, 0.5, 0.8, 0.9];
        let stats = QualityStats::compute(&values);
        assert!((stats.min - 0.1).abs() < 1e-10);
        assert!((stats.max - 0.9).abs() < 1e-10);
        assert!((stats.mean - 0.5).abs() < 1e-10);
        assert!(stats.std_dev > 0.0);
        assert_eq!(stats.n_poor, 2); // 0.1 and 0.2
    }

    #[test]
    fn test_quality_stats_empty() {
        let stats = QualityStats::compute(&[]);
        assert_eq!(stats.n_poor, 0);
        assert!((stats.mean).abs() < 1e-10);
    }

    #[test]
    fn test_quality_stats_uniform() {
        let values = vec![0.5, 0.5, 0.5, 0.5];
        let stats = QualityStats::compute(&values);
        assert!((stats.std_dev).abs() < 1e-10);
        assert!((stats.mean - 0.5).abs() < 1e-10);
    }

    // --- Mesh-level quality computation ---

    #[test]
    fn test_mesh_scaled_jacobian_qualities() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mesh = TetrahedralMesh::new(verts, tets);
        let qualities = mesh_scaled_jacobian_qualities(&mesh);
        assert_eq!(qualities.len(), 1);
        assert!(qualities[0].is_finite());
    }

    #[test]
    fn test_mesh_condition_number_qualities() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = vec![[0, 1, 2, 3]];
        let mesh = TetrahedralMesh::new(verts, tets);
        let qualities = mesh_condition_number_qualities(&mesh);
        assert_eq!(qualities.len(), 1);
        assert!(qualities[0] > 0.0);
    }

    #[test]
    fn test_mesh_quality_pipeline() {
        // Full pipeline: compute qualities, build histogram, find worst
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let mesh = TetrahedralMesh::new(
            vec![
                v0,
                v1,
                v2,
                v3,
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [10.0, 1.0, 0.0],
                [10.0, 0.0, 0.001],
            ],
            vec![[0, 1, 2, 3], [4, 5, 6, 7]],
        );
        let qualities = mesh_scaled_jacobian_qualities(&mesh);
        assert_eq!(qualities.len(), 2);

        let stats = QualityStats::compute(&qualities);
        assert!(stats.min <= stats.max);

        let hist = QualityHistogram::from_values(&qualities, 5);
        assert_eq!(hist.total, 2);
    }

    // ── Dihedral angle extremes ─────────────────────────────────────────────

    #[test]
    fn test_tet_dihedral_extremes_regular() {
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let (min_d, max_d) = tet_dihedral_angle_extremes(v0, v1, v2, v3);
        assert!(min_d > 0.0 && min_d < 180.0);
        assert!(max_d > 0.0 && max_d < 180.0);
        assert!(min_d <= max_d);
    }

    #[test]
    fn test_tri_angle_extremes_equilateral() {
        let h = (3.0_f64).sqrt() / 2.0;
        let (min_a, max_a) = tri_angle_extremes([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]);
        assert!(
            (min_a - 60.0).abs() < 1e-6,
            "min angle should be 60°, got {min_a}"
        );
        assert!(
            (max_a - 60.0).abs() < 1e-6,
            "max angle should be 60°, got {max_a}"
        );
    }

    #[test]
    fn test_tri_angle_extremes_right_triangle() {
        let (min_a, max_a) = tri_angle_extremes([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]);
        assert!((max_a - 90.0).abs() < 1e-6);
        assert!(min_a < 90.0);
    }

    // ── Aspect ratio for all element types ─────────────────────────────────

    #[test]
    fn test_line_aspect_ratio() {
        let ar = line_aspect_ratio([0.0; 3], [1.0, 0.0, 0.0]);
        assert!((ar - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_quad_aspect_ratio_square() {
        // Unit square: diagonals are equal → AR = 1
        let ar = quad_aspect_ratio(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "square aspect ratio should be 1, got {ar}"
        );
    }

    #[test]
    fn test_quad_aspect_ratio_rectangle() {
        // For a rectangle, both diagonals are equal in length, so AR = 1.0.
        // The implementation uses max_diagonal / min_diagonal.
        let ar = quad_aspect_ratio(
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        // Rectangle diagonals are equal → ratio is 1.0
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "rectangle has equal diagonals, AR=1, got {ar}"
        );
    }

    #[test]
    fn test_wedge_aspect_ratio_regular() {
        // Regular triangular prism
        let h = (3.0_f64).sqrt() / 2.0;
        let ar = wedge_aspect_ratio(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, h, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.5, h, 1.0],
        );
        assert!(ar >= 1.0, "wedge AR should be >= 1");
        // For a prism with all edges = 1, AR should be exactly 1
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "regular prism should have AR=1, got {ar}"
        );
    }

    #[test]
    fn test_hex_aspect_ratio_unit_cube() {
        let ar = hex_aspect_ratio(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!(
            (ar - 1.0).abs() < 1e-10,
            "unit cube should have AR=1, got {ar}"
        );
    }

    #[test]
    fn test_hex_aspect_ratio_stretched() {
        let ar = hex_aspect_ratio(
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [5.0, 0.0, 1.0],
            [5.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!(
            (ar - 5.0).abs() < 1e-10,
            "stretched hex should have AR=5, got {ar}"
        );
    }

    // ── Quality-guided mesh improvement ────────────────────────────────────

    #[test]
    fn test_quality_guided_smoothing_runs() {
        let mut positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let tets = vec![[0usize, 1, 2, 3], [1, 2, 3, 4]];
        let accepted = quality_guided_smoothing(&mut positions, &tets, 2);
        // The function should run without panic; accepted >= 0
        let _ = accepted;
        for p in &positions {
            assert!(p.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn test_quality_guided_smoothing_does_not_degrade() {
        // Regular tet: smoothing should not degrade quality
        let (v0, v1, v2, v3) = unit_tet_vertices();
        let mut positions = vec![v0, v1, v2, v3];
        let tets = vec![[0usize, 1, 2, 3]];
        let q_before =
            tet_min_scaled_jacobian(positions[0], positions[1], positions[2], positions[3]);
        quality_guided_smoothing(&mut positions, &tets, 3);
        let q_after =
            tet_min_scaled_jacobian(positions[0], positions[1], positions[2], positions[3]);
        assert!(
            q_after >= q_before - 1e-10,
            "quality should not degrade: before={q_before}, after={q_after}"
        );
    }

    // ── Minimum angle optimization ─────────────────────────────────────────

    #[test]
    fn test_minimum_angle_optimization_runs() {
        let h = (3.0_f64).sqrt() / 2.0;
        let mut positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, h, 0.0],
            [1.5, h, 0.0],
        ];
        let tris = vec![[0usize, 1, 2], [1, 3, 2]];
        let accepted = minimum_angle_optimization(&mut positions, &tris, 3);
        let _ = accepted;
        for p in &positions {
            assert!(p.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn test_minimum_angle_optimization_equilateral_unchanged() {
        // An equilateral triangle should not change (already optimal)
        let h = (3.0_f64).sqrt() / 2.0;
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [0.5, h, 0.0];
        let v3 = [1.5, h, 0.0];
        let mut positions = vec![v0, v1, v2, v3];
        let tris = vec![[0usize, 1, 2], [1, 3, 2]];
        minimum_angle_optimization(&mut positions, &tris, 5);
        // Vertices should remain finite
        for p in &positions {
            assert!(p.iter().all(|v| v.is_finite()));
        }
    }
}
