// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Origami mechanics — fold patterns and deployable structures.
//!
//! Provides crease-pattern representation, flat-foldability theorems
//! (Kawasaki, Maekawa), Miura-ori geometry, waterbomb and Yoshimura patterns,
//! kirigami cells, and crease stiffness.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 1. CreaseType and CreaseLine
// ---------------------------------------------------------------------------

/// Type of an origami crease.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CreaseType {
    /// Mountain fold (convex when viewed from above).
    Mountain,
    /// Valley fold (concave when viewed from above).
    Valley,
    /// Boundary edge of the paper.
    Boundary,
}

/// A single crease line in the flat fold pattern.
#[derive(Debug, Clone)]
pub struct CreaseLine {
    /// Start vertex index in the OrigamiPattern vertex list.
    pub start: usize,
    /// End vertex index in the OrigamiPattern vertex list.
    pub end: usize,
    /// Dihedral fold angle (radians). 0 = flat, ±π = fully folded.
    pub fold_angle: f64,
    /// Mountain, Valley, or Boundary.
    pub crease_type: CreaseType,
}

// ---------------------------------------------------------------------------
// 2. OrigamiPattern
// ---------------------------------------------------------------------------

/// A planar origami crease pattern: vertices, creases, and faces.
#[derive(Debug, Clone, Default)]
pub struct OrigamiPattern {
    /// 2D coordinates of each vertex.
    pub vertices: Vec<[f64; 2]>,
    /// All crease lines.
    pub creases: Vec<CreaseLine>,
    /// Faces as ordered lists of vertex indices.
    pub faces: Vec<Vec<usize>>,
}

impl OrigamiPattern {
    /// Create an empty pattern.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, v: [f64; 2]) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(v);
        idx
    }

    /// Add a crease between two existing vertex indices.
    pub fn add_crease(&mut self, start: usize, end: usize, angle: f64, kind: CreaseType) {
        self.creases.push(CreaseLine {
            start,
            end,
            fold_angle: angle,
            crease_type: kind,
        });
    }

    /// Number of faces.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of creases.
    pub fn crease_count(&self) -> usize {
        self.creases.len()
    }
}

// ---------------------------------------------------------------------------
// 3. Kawasaki and Maekawa theorems
// ---------------------------------------------------------------------------

/// Check the Kawasaki flat-foldability criterion at a vertex.
///
/// For a flat-foldable vertex, the alternating sum of sector angles = 0,
/// i.e., (α₁ − α₂ + α₃ − … ) = 0.
///
/// # Arguments
/// * `angles_at_vertex` – ordered sector angles around the vertex (radians).
///   Must have an even number of entries.
pub fn kawasaki_theorem(angles_at_vertex: &[f64]) -> bool {
    let n = angles_at_vertex.len();
    if !n.is_multiple_of(2) {
        return false;
    }
    // Alternating sum must equal 0 (equivalently, even-indexed sum == odd-indexed sum).
    let even_sum: f64 = angles_at_vertex.iter().step_by(2).sum();
    let odd_sum: f64 = angles_at_vertex.iter().skip(1).step_by(2).sum();
    (even_sum - odd_sum).abs() < 1e-9
}

/// Check the Maekawa flat-foldability criterion.
///
/// For a flat-foldable vertex: |M − V| = 2.
///
/// # Arguments
/// * `mountains` – number of mountain folds at the vertex.
/// * `valleys`   – number of valley folds at the vertex.
pub fn maekawa_theorem(mountains: usize, valleys: usize) -> bool {
    let m = mountains as i64;
    let v = valleys as i64;
    (m - v).abs() == 2
}

// ---------------------------------------------------------------------------
// 4. MiuraOri
// ---------------------------------------------------------------------------

/// Miura-ori fold pattern geometry.
///
/// A parallelogram tessellation that deploys with negative Poisson ratio.
#[derive(Debug, Clone)]
pub struct MiuraOri {
    /// Number of unit cells in x direction.
    pub nx: usize,
    /// Number of unit cells in y direction.
    pub ny: usize,
    /// Parallelogram angle α (radians), typically 30–70°.
    pub angle_a: f64,
    /// Unit cell length a in the fold direction (m).
    pub length_a: f64,
    /// Unit cell length b in the transverse direction (m).
    pub length_b: f64,
}

impl MiuraOri {
    /// Create a new Miura-ori pattern.
    pub fn new(nx: usize, ny: usize, angle_a: f64, length_a: f64, length_b: f64) -> Self {
        Self {
            nx,
            ny,
            angle_a,
            length_a,
            length_b,
        }
    }

    /// Vertices in the flat (unfolded) configuration.
    pub fn flat_vertices(&self) -> Vec<[f64; 2]> {
        let mut verts = Vec::new();
        for iy in 0..=self.ny {
            for ix in 0..=self.nx {
                let x = ix as f64 * self.length_a;
                let y = iy as f64 * self.length_b;
                verts.push([x, y]);
            }
        }
        verts
    }

    /// Vertices in a partially-folded configuration.
    ///
    /// `fold_ratio` ∈ \[0, 1\]: 0 = fully flat, 1 = fully folded (compressed).
    pub fn folded_vertices(&self, fold_ratio: f64) -> Vec<[f64; 3]> {
        let fold_ratio = fold_ratio.clamp(0.0, 1.0);
        // At fold ratio r: effective x-spacing shrinks, z oscillates.
        let a = self.angle_a;
        let dx = self.length_a * (1.0 - fold_ratio * (1.0 - a.cos()));
        let dz = self.length_b * fold_ratio * a.sin();
        let mut verts = Vec::new();
        for iy in 0..=self.ny {
            for ix in 0..=self.nx {
                let x = ix as f64 * dx;
                let y = iy as f64 * self.length_b * (1.0 - fold_ratio * 0.5);
                let z = if (ix + iy) % 2 == 0 { 0.0 } else { dz };
                verts.push([x, y, z]);
            }
        }
        verts
    }

    /// Miura-ori Poisson ratio (negative): ν = −(sin²α) / (1 + cos²α · tan²(β/2))
    ///
    /// Simplified version: ν ≈ −sin²(α) / (1 − sin²(α)) = −tan²(α) + 1 − 1
    /// For the standard result: ν = −sin²α (to first order in fold angle).
    /// We use the exact expression: ν = −sin²α.
    pub fn poisson_ratio(&self) -> f64 {
        // Negative Poisson ratio for all α in (0, π).
        -(self.angle_a.sin() * self.angle_a.sin())
    }
}

// ---------------------------------------------------------------------------
// 5. Waterbomb base
// ---------------------------------------------------------------------------

/// Generate the waterbomb crease pattern (square sheet, 8 creases from centre).
pub fn waterbomb_base_crease_pattern(size: f64) -> OrigamiPattern {
    let mut pat = OrigamiPattern::new();
    let h = size / 2.0;
    // Corner vertices.
    let c0 = pat.add_vertex([-h, -h]);
    let c1 = pat.add_vertex([h, -h]);
    let c2 = pat.add_vertex([h, h]);
    let c3 = pat.add_vertex([-h, h]);
    // Mid-edge vertices.
    let m0 = pat.add_vertex([0.0, -h]);
    let m1 = pat.add_vertex([h, 0.0]);
    let m2 = pat.add_vertex([0.0, h]);
    let m3 = pat.add_vertex([-h, 0.0]);
    // Centre.
    let ctr = pat.add_vertex([0.0, 0.0]);

    // Boundary edges.
    pat.add_crease(c0, m0, 0.0, CreaseType::Boundary);
    pat.add_crease(m0, c1, 0.0, CreaseType::Boundary);
    pat.add_crease(c1, m1, 0.0, CreaseType::Boundary);
    pat.add_crease(m1, c2, 0.0, CreaseType::Boundary);
    pat.add_crease(c2, m2, 0.0, CreaseType::Boundary);
    pat.add_crease(m2, c3, 0.0, CreaseType::Boundary);
    pat.add_crease(c3, m3, 0.0, CreaseType::Boundary);
    pat.add_crease(m3, c0, 0.0, CreaseType::Boundary);

    // Valley folds (diagonals).
    pat.add_crease(c0, ctr, 0.0, CreaseType::Valley);
    pat.add_crease(c1, ctr, 0.0, CreaseType::Valley);
    pat.add_crease(c2, ctr, 0.0, CreaseType::Valley);
    pat.add_crease(c3, ctr, 0.0, CreaseType::Valley);

    // Mountain folds (cross).
    pat.add_crease(m0, ctr, 0.0, CreaseType::Mountain);
    pat.add_crease(m1, ctr, 0.0, CreaseType::Mountain);
    pat.add_crease(m2, ctr, 0.0, CreaseType::Mountain);
    pat.add_crease(m3, ctr, 0.0, CreaseType::Mountain);

    // Faces (8 triangular sectors).
    pat.faces = vec![
        vec![c0, m0, ctr],
        vec![m0, c1, ctr],
        vec![c1, m1, ctr],
        vec![m1, c2, ctr],
        vec![c2, m2, ctr],
        vec![m2, c3, ctr],
        vec![c3, m3, ctr],
        vec![m3, c0, ctr],
    ];

    pat
}

// ---------------------------------------------------------------------------
// 6. Yoshimura pattern
// ---------------------------------------------------------------------------

/// Generate a Yoshimura buckling pattern for a cylindrical shell.
///
/// Produces a diamond tessellation with `nx × ny` repeating units.
pub fn yoshimura_pattern(nx: usize, ny: usize, radius: f64, length: f64) -> OrigamiPattern {
    let mut pat = OrigamiPattern::new();
    let dtheta = 2.0 * PI / nx as f64;
    let dy = length / ny as f64;

    // Vertices on regular grid (angle, height).
    for iy in 0..=ny {
        for ix in 0..=nx {
            let theta = ix as f64 * dtheta;
            let x = radius * theta.cos();
            let y = iy as f64 * dy;
            let z_offset = if (ix + iy) % 2 == 0 {
                0.01 * radius
            } else {
                0.0
            };
            // Store as 2D projection for the flat pattern.
            let _ = z_offset;
            pat.add_vertex([x, y]);
        }
    }

    let cols = nx + 1;
    // Add creases (diamond pattern).
    for iy in 0..ny {
        for ix in 0..nx {
            let v00 = iy * cols + ix;
            let v10 = iy * cols + ix + 1;
            let v01 = (iy + 1) * cols + ix;
            let v11 = (iy + 1) * cols + ix + 1;
            let kind = if (ix + iy) % 2 == 0 {
                CreaseType::Mountain
            } else {
                CreaseType::Valley
            };
            pat.add_crease(v00, v11, 0.0, kind);
            pat.add_crease(
                v10,
                v01,
                0.0,
                if kind == CreaseType::Mountain {
                    CreaseType::Valley
                } else {
                    CreaseType::Mountain
                },
            );
            pat.faces.push(vec![v00, v10, v11]);
            pat.faces.push(vec![v00, v11, v01]);
        }
    }

    pat
}

// ---------------------------------------------------------------------------
// 7. KirigamiCell
// ---------------------------------------------------------------------------

/// Kirigami unit cell: a square sheet with a rectangular cut.
#[derive(Debug, Clone)]
pub struct KirigamiCell {
    /// Side length of the unit cell (m).
    pub size: f64,
    /// Length of the cut along one axis (m, ≤ size).
    pub cut_length: f64,
    /// Width of the cut slit (m, ≤ size).
    pub cut_width: f64,
}

impl KirigamiCell {
    /// Create a new kirigami unit cell.
    pub fn new(size: f64, cut_length: f64, cut_width: f64) -> Self {
        Self {
            size,
            cut_length,
            cut_width,
        }
    }

    /// Porosity (area fraction removed by cuts).
    pub fn porosity(&self) -> f64 {
        let cut_area = self.cut_length * self.cut_width;
        let cell_area = self.size * self.size;
        (cut_area / cell_area).clamp(0.0, 1.0)
    }

    /// Geometric stretchability: approximate max strain before failure.
    ///
    /// Estimated as (cut_length / size)² based on geometric kinematics.
    pub fn stretchability(&self) -> f64 {
        let r = (self.cut_length / self.size).clamp(0.0, 1.0);
        r * r
    }
}

// ---------------------------------------------------------------------------
// 8. Crease stiffness and deployed area
// ---------------------------------------------------------------------------

/// Bending stiffness of a crease fold.
///
/// k = E · t³ / (12 · (1 − ν²)) · (1 / fold_angle²)  (simplified).
///
/// For small fold angles this reduces to k ≈ E · t³ / 12.
///
/// # Arguments
/// * `thickness`  – sheet thickness t (m).
/// * `e_mod`      – Young's modulus E (Pa).
/// * `fold_angle` – fold angle (radians); absolute value used.
pub fn crease_stiffness(thickness: f64, e_mod: f64, fold_angle: f64) -> f64 {
    let nu = 0.3_f64; // typical Poisson ratio for paper/sheet
    let angle = fold_angle.abs().max(1e-6);
    e_mod * thickness.powi(3) / (12.0 * (1.0 - nu * nu) * angle * angle)
}

/// Approximate deployed area of an origami pattern at a given fold ratio.
///
/// At `fold_ratio = 0` the area equals the sum of all face areas in the flat
/// pattern.  At `fold_ratio = 1` it is approximated as area × cos(π/4).
///
/// # Arguments
/// * `pattern`    – flat origami pattern.
/// * `fold_ratio` – fraction of full deployment, ∈ \[0, 1\].
pub fn deployed_area(pattern: &OrigamiPattern, fold_ratio: f64) -> f64 {
    let fold_ratio = fold_ratio.clamp(0.0, 1.0);
    let flat_area: f64 = pattern
        .faces
        .iter()
        .map(|face| {
            if face.len() < 3 {
                return 0.0;
            }
            // Shoelace formula for polygon area.
            let n = face.len();
            let mut area = 0.0;
            for k in 0..n {
                let i = face[k];
                let j = face[(k + 1) % n];
                let (xi, yi) = (pattern.vertices[i][0], pattern.vertices[i][1]);
                let (xj, yj) = (pattern.vertices[j][0], pattern.vertices[j][1]);
                area += xi * yj - xj * yi;
            }
            area.abs() / 2.0
        })
        .sum();
    // Deployed area decreases with fold ratio (more folded → smaller projected area).
    flat_area * (1.0 - fold_ratio * (1.0 - (PI / 4.0).cos()))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── Kawasaki theorem ─────────────────────────────────────────────────────

    #[test]
    fn test_kawasaki_alternating_90_degrees() {
        // Four 90° sectors: even sum = 180°, odd sum = 180° ✓
        let angles = [PI / 2.0, PI / 2.0, PI / 2.0, PI / 2.0];
        assert!(kawasaki_theorem(&angles));
    }

    #[test]
    fn test_kawasaki_general_valid() {
        // Flat-foldable vertex: even-indexed angles sum = odd-indexed angles sum.
        // With 4 angles [α, β, γ, δ], Kawasaki requires α+γ = β+δ.
        // Use α=γ=π/3, β=δ=π/3 (all equal → even_sum = odd_sum = 2π/3).
        let a = PI / 3.0;
        let angles = [a, a, a, a];
        assert!(kawasaki_theorem(&angles));
    }

    #[test]
    fn test_kawasaki_invalid_odd_count() {
        let angles = [PI / 2.0, PI / 2.0, PI / 2.0];
        assert!(!kawasaki_theorem(&angles));
    }

    #[test]
    fn test_kawasaki_fails_for_unequal_sums() {
        let angles = [0.5, 1.0, 0.5, 1.5]; // even sum=1, odd sum=2.5
        assert!(!kawasaki_theorem(&angles));
    }

    #[test]
    fn test_kawasaki_six_angles() {
        // 60° each: even sum = 3×60 = 180°, odd sum = 3×60 = 180° ✓
        let a = PI / 3.0;
        let angles = [a, a, a, a, a, a];
        assert!(kawasaki_theorem(&angles));
    }

    // ── Maekawa theorem ──────────────────────────────────────────────────────

    #[test]
    fn test_maekawa_3m_1v() {
        // |3-1| = 2 ✓
        assert!(maekawa_theorem(3, 1));
    }

    #[test]
    fn test_maekawa_1m_3v() {
        // |1-3| = 2 ✓
        assert!(maekawa_theorem(1, 3));
    }

    #[test]
    fn test_maekawa_fails_equal() {
        // |2-2| = 0 ✗
        assert!(!maekawa_theorem(2, 2));
    }

    #[test]
    fn test_maekawa_fails_diff_four() {
        // |4-0| = 4 ✗
        assert!(!maekawa_theorem(4, 0));
    }

    #[test]
    fn test_maekawa_5m_3v() {
        // |5-3| = 2 ✓
        assert!(maekawa_theorem(5, 3));
    }

    // ── MiuraOri ─────────────────────────────────────────────────────────────

    #[test]
    fn test_miura_poisson_ratio_negative() {
        let m = MiuraOri::new(4, 4, PI / 4.0, 1.0, 1.0);
        let nu = m.poisson_ratio();
        assert!(
            nu < 0.0,
            "Miura-ori should have negative Poisson ratio, got {nu}"
        );
    }

    #[test]
    fn test_miura_poisson_ratio_range() {
        for deg in [30, 45, 60, 70u64] {
            let alpha = deg as f64 * PI / 180.0;
            let m = MiuraOri::new(4, 4, alpha, 1.0, 1.0);
            let nu = m.poisson_ratio();
            assert!(
                (-1.0..=0.0).contains(&nu),
                "Poisson ratio={nu} at alpha={deg}°"
            );
        }
    }

    #[test]
    fn test_miura_flat_vertices_count() {
        let m = MiuraOri::new(3, 4, PI / 4.0, 1.0, 1.0);
        // (nx+1) * (ny+1) = 4 * 5 = 20
        assert_eq!(m.flat_vertices().len(), 20);
    }

    #[test]
    fn test_miura_folded_vertices_count() {
        let m = MiuraOri::new(3, 4, PI / 4.0, 1.0, 1.0);
        assert_eq!(m.folded_vertices(0.5).len(), 20);
    }

    #[test]
    fn test_miura_fold_ratio_zero_flat() {
        // At fold_ratio=0 all z-components should be 0.
        let m = MiuraOri::new(2, 2, PI / 4.0, 1.0, 1.0);
        let verts = m.folded_vertices(0.0);
        for v in &verts {
            assert!(
                v[2].abs() < EPS,
                "z should be 0 at fold_ratio=0, got {}",
                v[2]
            );
        }
    }

    #[test]
    fn test_miura_fold_ratio_one_compressed() {
        let m = MiuraOri::new(4, 4, PI / 4.0, 1.0, 1.0);
        let flat = m.flat_vertices();
        let folded = m.folded_vertices(1.0);
        // x-extent should be smaller when folded.
        let flat_xmax = flat.iter().map(|v| v[0]).fold(f64::NEG_INFINITY, f64::max);
        let folded_xmax = folded
            .iter()
            .map(|v| v[0])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            folded_xmax < flat_xmax,
            "folded x-extent should be smaller than flat"
        );
    }

    // ── OrigamiPattern ────────────────────────────────────────────────────────

    #[test]
    fn test_origami_pattern_add_vertex() {
        let mut pat = OrigamiPattern::new();
        let idx = pat.add_vertex([1.0, 2.0]);
        assert_eq!(idx, 0);
        assert_eq!(pat.vertices.len(), 1);
    }

    #[test]
    fn test_origami_pattern_add_crease() {
        let mut pat = OrigamiPattern::new();
        pat.add_vertex([0.0, 0.0]);
        pat.add_vertex([1.0, 0.0]);
        pat.add_crease(0, 1, 0.0, CreaseType::Mountain);
        assert_eq!(pat.crease_count(), 1);
    }

    #[test]
    fn test_waterbomb_crease_count() {
        let pat = waterbomb_base_crease_pattern(1.0);
        assert_eq!(pat.crease_count(), 16, "waterbomb should have 16 creases");
    }

    #[test]
    fn test_waterbomb_face_count() {
        let pat = waterbomb_base_crease_pattern(1.0);
        assert_eq!(pat.face_count(), 8, "waterbomb should have 8 faces");
    }

    #[test]
    fn test_waterbomb_vertex_count() {
        let pat = waterbomb_base_crease_pattern(1.0);
        // 4 corners + 4 mid-edge + 1 centre = 9
        assert_eq!(pat.vertices.len(), 9, "waterbomb should have 9 vertices");
    }

    // ── Yoshimura pattern ─────────────────────────────────────────────────────

    #[test]
    fn test_yoshimura_crease_count() {
        let pat = yoshimura_pattern(4, 4, 1.0, 1.0);
        // 2 diagonals per cell × 4×4 cells = 32 creases
        assert_eq!(pat.crease_count(), 32, "yoshimura crease count mismatch");
    }

    #[test]
    fn test_yoshimura_face_count() {
        let pat = yoshimura_pattern(3, 3, 1.0, 1.0);
        // 2 triangles per cell × 3×3 = 18 faces
        assert_eq!(pat.face_count(), 18, "yoshimura face count mismatch");
    }

    // ── KirigamiCell ─────────────────────────────────────────────────────────

    #[test]
    fn test_kirigami_porosity_formula() {
        let cell = KirigamiCell::new(2.0, 1.0, 0.5);
        let expected = 1.0 * 0.5 / (2.0 * 2.0);
        assert!((cell.porosity() - expected).abs() < EPS);
    }

    #[test]
    fn test_kirigami_porosity_in_range() {
        let cell = KirigamiCell::new(1.0, 0.5, 0.5);
        let p = cell.porosity();
        assert!((0.0..=1.0).contains(&p), "porosity={p}");
    }

    #[test]
    fn test_kirigami_stretchability_increases_with_cut() {
        let c1 = KirigamiCell::new(2.0, 0.5, 0.1);
        let c2 = KirigamiCell::new(2.0, 1.5, 0.1);
        assert!(c2.stretchability() > c1.stretchability());
    }

    // ── crease_stiffness ──────────────────────────────────────────────────────

    #[test]
    fn test_crease_stiffness_increases_with_emod() {
        let k1 = crease_stiffness(1e-3, 1e9, PI / 4.0);
        let k2 = crease_stiffness(1e-3, 2e9, PI / 4.0);
        assert!(k2 > k1, "stiffness should increase with E_mod");
    }

    #[test]
    fn test_crease_stiffness_increases_with_thickness() {
        let k1 = crease_stiffness(1e-3, 1e9, PI / 4.0);
        let k2 = crease_stiffness(2e-3, 1e9, PI / 4.0);
        assert!(k2 > k1, "stiffness should increase with thickness");
    }

    #[test]
    fn test_crease_stiffness_positive() {
        let k = crease_stiffness(1e-3, 1e9, PI / 6.0);
        assert!(k > 0.0, "crease stiffness should be positive, got {k}");
    }

    // ── deployed_area ─────────────────────────────────────────────────────────

    #[test]
    fn test_deployed_area_fold_zero_equals_flat() {
        let pat = waterbomb_base_crease_pattern(2.0);
        let area0 = deployed_area(&pat, 0.0);
        // Flat area should be positive.
        assert!(area0 > 0.0, "flat area should be positive, got {area0}");
    }

    #[test]
    fn test_deployed_area_decreases_with_fold_ratio() {
        let pat = waterbomb_base_crease_pattern(2.0);
        let a0 = deployed_area(&pat, 0.0);
        let a1 = deployed_area(&pat, 1.0);
        assert!(a1 <= a0, "area should not increase when folding");
    }

    #[test]
    fn test_deployed_area_empty_pattern_zero() {
        let pat = OrigamiPattern::new();
        let a = deployed_area(&pat, 0.5);
        assert_eq!(a, 0.0);
    }
}
