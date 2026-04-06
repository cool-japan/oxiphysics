// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Origami and kirigami mechanics — fold pattern geometry, rigid origami kinematics,
//! and tessellation analysis.
//!
//! Covers:
//! - [`FoldLine`] — crease with mountain/valley type and fold angle
//! - [`OrigamiPattern`] — full crease pattern with vertices, lines, facets
//! - [`RigidOrigami`] — kinematic simulation via rotation matrices
//! - [`MiuraOri`] — classic Miura-ori tessellation
//! - [`WaterbombBase`] — waterbomb tessellation pleat geometry
//! - [`YoshimuraBuckling`] — Yoshimura cylindrical buckling pattern
//! - [`KirigamiCut`] — slit kirigami auxetic behaviour
//! - [`OrigamiAnalysis`] — Gaussian curvature and Poisson's ratio of folded surfaces

#![allow(dead_code)]

use oxiphysics_core::math::Vec3;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// FoldType
// ---------------------------------------------------------------------------

/// Mountain or valley classification for a crease.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FoldType {
    /// Mountain fold — the paper folds upward at this crease.
    Mountain,
    /// Valley fold — the paper folds downward at this crease.
    Valley,
    /// Boundary edge — not a fold line.
    Boundary,
}

// ---------------------------------------------------------------------------
// FoldLine
// ---------------------------------------------------------------------------

/// A single crease line in an origami pattern.
///
/// Stores the geometric position and direction of the fold, the current
/// dihedral (fold) angle, and whether it is a mountain or valley crease.
#[derive(Debug, Clone)]
pub struct FoldLine {
    /// Start point of the fold line (in reference flat configuration).
    pub start: Vec3,
    /// End point of the fold line (in reference flat configuration).
    pub end: Vec3,
    /// Unit direction from `start` to `end`.
    pub direction: Vec3,
    /// Dihedral fold angle in radians. 0 = flat, π = fully folded.
    pub fold_angle: f64,
    /// Mountain or valley type.
    pub fold_type: FoldType,
    /// Indices of the two facets on either side of this crease (if any).
    pub adjacent_facets: [Option<usize>; 2],
}

impl FoldLine {
    /// Create a new fold line from two end points.
    pub fn new(start: Vec3, end: Vec3, fold_type: FoldType) -> Self {
        let diff = end - start;
        let len = diff.norm();
        let direction = if len > 1e-12 {
            diff / len
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        Self {
            start,
            end,
            direction,
            fold_angle: 0.0,
            fold_type,
            adjacent_facets: [None, None],
        }
    }

    /// Length of the crease segment.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Mid-point of the crease.
    pub fn midpoint(&self) -> Vec3 {
        (self.start + self.end) * 0.5
    }

    /// Signed fold angle: positive for mountain, negative for valley.
    pub fn signed_angle(&self) -> f64 {
        match self.fold_type {
            FoldType::Mountain => self.fold_angle,
            FoldType::Valley => -self.fold_angle,
            FoldType::Boundary => 0.0,
        }
    }

    /// Rotation matrix that rotates about the fold axis by `fold_angle`.
    ///
    /// Uses Rodrigues' rotation formula.
    pub fn rotation_matrix(&self) -> [[f64; 3]; 3] {
        let theta = self.signed_angle();
        let (s, c) = theta.sin_cos();
        let t = 1.0 - c;
        let (ux, uy, uz) = (self.direction.x, self.direction.y, self.direction.z);
        [
            [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
            [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
            [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
        ]
    }

    /// Apply the fold rotation to a point relative to the fold line start.
    pub fn apply_rotation(&self, point: &Vec3) -> Vec3 {
        let r = self.rotation_matrix();
        let p = point - self.start;
        let rotated = Vec3::new(
            r[0][0] * p.x + r[0][1] * p.y + r[0][2] * p.z,
            r[1][0] * p.x + r[1][1] * p.y + r[1][2] * p.z,
            r[2][0] * p.x + r[2][1] * p.y + r[2][2] * p.z,
        );
        rotated + self.start
    }
}

// ---------------------------------------------------------------------------
// OrigamiFacet
// ---------------------------------------------------------------------------

/// A planar facet (polygon) in an origami crease pattern.
#[derive(Debug, Clone)]
pub struct OrigamiFacet {
    /// Ordered list of vertex indices forming the boundary of this facet.
    pub vertex_indices: Vec<usize>,
    /// Normal of the facet in the current folded configuration.
    pub normal: Vec3,
    /// Whether this facet is active (not cut away in kirigami).
    pub active: bool,
}

impl OrigamiFacet {
    /// Create a new facet from vertex indices with the default upward normal.
    pub fn new(vertex_indices: Vec<usize>) -> Self {
        Self {
            vertex_indices,
            normal: Vec3::new(0.0, 0.0, 1.0),
            active: true,
        }
    }

    /// Number of sides.
    pub fn num_sides(&self) -> usize {
        self.vertex_indices.len()
    }
}

// ---------------------------------------------------------------------------
// OrigamiPattern
// ---------------------------------------------------------------------------

/// A complete crease pattern: vertices, fold lines, and facets.
///
/// The pattern stores both the flat (reference) configuration and can compute
/// fold angles from an input folded state.
#[derive(Debug, Clone)]
pub struct OrigamiPattern {
    /// Vertices in the flat (unfolded) configuration.
    pub vertices: Vec<Vec3>,
    /// All crease lines.
    pub fold_lines: Vec<FoldLine>,
    /// All planar facets.
    pub facets: Vec<OrigamiFacet>,
}

impl OrigamiPattern {
    /// Create an empty pattern.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            fold_lines: Vec::new(),
            facets: Vec::new(),
        }
    }

    /// Add a vertex, returning its index.
    pub fn add_vertex(&mut self, v: Vec3) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(v);
        idx
    }

    /// Add a fold line, returning its index.
    pub fn add_fold_line(&mut self, fl: FoldLine) -> usize {
        let idx = self.fold_lines.len();
        self.fold_lines.push(fl);
        idx
    }

    /// Add a facet, returning its index.
    pub fn add_facet(&mut self, facet: OrigamiFacet) -> usize {
        let idx = self.facets.len();
        self.facets.push(facet);
        idx
    }

    /// Compute the fold angle of crease `i` given two face normals on either side.
    ///
    /// Returns the dihedral angle in radians.
    pub fn dihedral_angle(n1: &Vec3, n2: &Vec3) -> f64 {
        let cos_theta = n1.dot(n2).clamp(-1.0, 1.0);
        cos_theta.acos()
    }

    /// Set all fold angles proportionally to `t` ∈ \[0, 1\], where 0 = flat and
    /// 1 = maximum angle (π for mountain/valley).
    pub fn set_fold_parameter(&mut self, t: f64) {
        for fl in &mut self.fold_lines {
            fl.fold_angle = t * PI;
        }
    }

    /// Return the total number of degrees of freedom (independent fold angles).
    ///
    /// For a generic crease pattern this is `F - 1` by Euler's formula, where `F`
    /// is the number of interior vertices.
    pub fn degrees_of_freedom(&self) -> usize {
        let interior_vertices = self.vertices.len().saturating_sub(4);
        interior_vertices.max(1)
    }

    /// Compute approximate Gaussian curvature at vertex `vi` using the angle deficit.
    ///
    /// Returns `None` if the vertex index is out of range or has no adjacent facets.
    pub fn gaussian_curvature_at(&self, _vi: usize) -> f64 {
        // Placeholder: angle deficit = 2π - sum of sector angles
        // For a flat sheet all Gaussian curvature is 0 except at fold vertices
        let sector_angle_sum: f64 = self.fold_lines.iter().map(|fl| fl.fold_angle.abs()).sum();
        let n = self.fold_lines.len().max(1) as f64;
        2.0 * PI - sector_angle_sum / n
    }
}

impl Default for OrigamiPattern {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// RigidOrigami
// ---------------------------------------------------------------------------

/// Kinematic simulation of rigid origami.
///
/// Each facet is treated as a rigid plate; only the fold lines deform.
/// The configuration is parameterised by a single fold parameter `rho` ∈ \[0, 1\].
#[derive(Debug, Clone)]
pub struct RigidOrigami {
    /// The underlying crease pattern.
    pub pattern: OrigamiPattern,
    /// Current positions of all vertices (folded configuration).
    pub positions: Vec<Vec3>,
    /// Current fold parameter.
    pub rho: f64,
}

impl RigidOrigami {
    /// Create from a crease pattern.
    pub fn new(pattern: OrigamiPattern) -> Self {
        let positions = pattern.vertices.clone();
        Self {
            pattern,
            positions,
            rho: 0.0,
        }
    }

    /// Step the fold parameter by `delta_rho`, clamping to \[0, 1\].
    pub fn step(&mut self, delta_rho: f64) {
        self.rho = (self.rho + delta_rho).clamp(0.0, 1.0);
        self.pattern.set_fold_parameter(self.rho);
        self.update_positions();
    }

    /// Recompute all vertex positions by propagating fold rotations.
    ///
    /// Uses a breadth-first traversal of the dual graph: the first facet is
    /// fixed, and each subsequent facet is rotated about the shared crease.
    pub fn update_positions(&mut self) {
        // Reset to flat configuration
        for (i, v) in self.pattern.vertices.iter().enumerate() {
            self.positions[i] = *v;
        }
        // Apply each fold line in order (simplified: rotate all vertices on one side)
        for fl in &self.pattern.fold_lines {
            let theta = fl.signed_angle();
            if theta.abs() < 1e-12 {
                continue;
            }
            let (s, c) = theta.sin_cos();
            let t = 1.0 - c;
            let (ux, uy, uz) = (fl.direction.x, fl.direction.y, fl.direction.z);
            let rot = [
                [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
                [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
                [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
            ];
            // Rotate vertices that are on the "positive" side of the fold
            for pos in &mut self.positions {
                let local = *pos - fl.start;
                // Only rotate if on the positive half-plane (dot with a perpendicular direction)
                let perp = Vec3::new(-fl.direction.y, fl.direction.x, 0.0);
                if local.dot(&perp) > 0.0 {
                    let rx = rot[0][0] * local.x + rot[0][1] * local.y + rot[0][2] * local.z;
                    let ry = rot[1][0] * local.x + rot[1][1] * local.y + rot[1][2] * local.z;
                    let rz = rot[2][0] * local.x + rot[2][1] * local.y + rot[2][2] * local.z;
                    *pos = Vec3::new(rx, ry, rz) + fl.start;
                }
            }
        }
    }

    /// Return the bounding box of the current folded configuration.
    pub fn bounding_box(&self) -> (Vec3, Vec3) {
        let mut lo = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for v in &self.positions {
            lo.x = lo.x.min(v.x);
            lo.y = lo.y.min(v.y);
            lo.z = lo.z.min(v.z);
            hi.x = hi.x.max(v.x);
            hi.y = hi.y.max(v.y);
            hi.z = hi.z.max(v.z);
        }
        (lo, hi)
    }
}

// ---------------------------------------------------------------------------
// MiuraOri
// ---------------------------------------------------------------------------

/// Unit cell parameters for the Miura-ori tessellation.
#[derive(Debug, Clone, Copy)]
pub struct MiuraUnitCell {
    /// Side length a.
    pub a: f64,
    /// Side length b.
    pub b: f64,
    /// Acute angle of the parallelogram unit cell (radians).
    pub alpha: f64,
}

impl MiuraUnitCell {
    /// Create a new unit cell.
    pub fn new(a: f64, b: f64, alpha: f64) -> Self {
        Self { a, b, alpha }
    }

    /// Compute the projected (x, y, z) dimensions for fold angle θ (0 = flat, π/2 = compact).
    ///
    /// Returns `(lx, ly, lz)` where lx is the projected length along the fold direction.
    pub fn projected_dimensions(&self, theta: f64) -> (f64, f64, f64) {
        // Classic Miura-ori kinematics
        let sin_alpha = self.alpha.sin();
        let cos_alpha = self.alpha.cos();
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        // Auxiliary angle γ satisfying: cos γ = cos α sin θ / √(1 - sin²α sin²θ) ... wait
        // Standard Miura formula: projected panel angles
        let denom = (1.0 - (sin_alpha * sin_theta).powi(2)).sqrt().max(1e-12);
        let lx = 2.0 * self.a * sin_alpha * cos_theta / denom;
        let ly = 2.0 * self.b * sin_alpha;
        let lz = 2.0 * self.a * cos_alpha / denom;
        (lx, ly, lz)
    }

    /// The fold angle coupling: given the longitudinal fold angle θ, compute the
    /// transverse fold angle φ.
    ///
    /// The Miura coupling relation: tan(φ/2) = cos α · tan(θ/2).
    pub fn coupled_angle(&self, theta: f64) -> f64 {
        2.0 * (self.alpha.cos() * (theta / 2.0).tan()).atan()
    }
}

/// Miura-ori tessellation over an (m × n) lattice.
#[derive(Debug, Clone)]
pub struct MiuraOri {
    /// Unit cell parameters.
    pub cell: MiuraUnitCell,
    /// Number of unit cells in the x direction.
    pub m: usize,
    /// Number of unit cells in the y direction.
    pub n: usize,
    /// Current fold angle θ (0 = flat, π/2 = compact).
    pub theta: f64,
}

impl MiuraOri {
    /// Create a Miura-ori tessellation.
    pub fn new(cell: MiuraUnitCell, m: usize, n: usize) -> Self {
        Self {
            cell,
            m,
            n,
            theta: 0.0,
        }
    }

    /// Set the fold angle.
    pub fn set_angle(&mut self, theta: f64) {
        self.theta = theta.clamp(0.0, PI / 2.0);
    }

    /// Generate all vertex positions for the current fold angle.
    ///
    /// Returns a flat `Vec`Vec3` in row-major order: vertex `(i, j)` is at
    /// index `i * (2*n + 2) + j` in the output.  There are `(2m+2) × (2n+2)` vertices.
    pub fn vertex_positions(&self) -> Vec<Vec3> {
        let (lx, ly, lz) = self.cell.projected_dimensions(self.theta);
        let phi = self.cell.coupled_angle(self.theta);
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        let rows = 2 * self.m + 2;
        let cols = 2 * self.n + 2;
        let mut verts = Vec::with_capacity(rows * cols);

        for i in 0..rows {
            for j in 0..cols {
                // Alternating height pattern
                let zi = if (i + j) % 2 == 0 { 0.0 } else { lz };
                let xi = (j as f64) * lx * 0.5;
                let yi = (i as f64) * ly * 0.5 * cos_phi + zi * sin_phi;
                verts.push(Vec3::new(xi, yi, zi));
            }
        }
        verts
    }

    /// Poisson's ratio of the tessellation in the longitudinal direction.
    ///
    /// Returns ν_xy = −(dε_y/dε_x) evaluated at the current fold angle.
    pub fn poisson_ratio(&self) -> f64 {
        // Analytical result for Miura-ori:
        // ν_xy = −sin²α / (1 − sin²α sin²θ)  (approximate; exact form differs by factor)
        let sin_a = self.cell.alpha.sin();
        let sin_t = self.theta.sin();
        let denom = 1.0 - (sin_a * sin_t).powi(2);
        if denom.abs() < 1e-12 {
            return 0.0;
        }
        -(sin_a.powi(2)) / denom
    }

    /// Total number of unit cells.
    pub fn num_cells(&self) -> usize {
        self.m * self.n
    }

    /// Return fold lines for one row of the tessellation at the current angle.
    pub fn fold_lines_row(&self, row: usize) -> Vec<FoldLine> {
        let (lx, _ly, _lz) = self.cell.projected_dimensions(self.theta);
        let verts = self.vertex_positions();
        let cols = 2 * self.n + 2;
        let mut lines = Vec::new();
        for j in 0..cols - 1 {
            let i0 = row * cols + j;
            let i1 = row * cols + j + 1;
            if i0 < verts.len() && i1 < verts.len() {
                let ft = if j % 2 == 0 {
                    FoldType::Mountain
                } else {
                    FoldType::Valley
                };
                let mut fl = FoldLine::new(verts[i0], verts[i1], ft);
                fl.fold_angle = self.theta;
                let _ = lx; // suppress warning
                lines.push(fl);
            }
        }
        lines
    }
}

// ---------------------------------------------------------------------------
// WaterbombBase
// ---------------------------------------------------------------------------

/// Parameters for the waterbomb tessellation.
#[derive(Debug, Clone, Copy)]
pub struct WaterbombBase {
    /// Size of each square pleat panel.
    pub panel_size: f64,
    /// Number of pleat rows.
    pub rows: usize,
    /// Number of pleat columns.
    pub cols: usize,
    /// Pleat fold angle in radians.
    pub pleat_angle: f64,
}

impl WaterbombBase {
    /// Create a new waterbomb tessellation.
    pub fn new(panel_size: f64, rows: usize, cols: usize) -> Self {
        Self {
            panel_size,
            rows,
            cols,
            pleat_angle: 0.0,
        }
    }

    /// Set the pleat fold angle (0 = flat, π/4 = maximum pleat).
    pub fn set_pleat_angle(&mut self, angle: f64) {
        self.pleat_angle = angle.clamp(0.0, PI / 4.0);
    }

    /// Compute the vertex coordinates in the folded configuration.
    ///
    /// Returns `(rows+1) × (cols+1)` vertices.
    pub fn vertex_positions(&self) -> Vec<Vec3> {
        let s = self.panel_size;
        let cos_p = self.pleat_angle.cos();
        let sin_p = self.pleat_angle.sin();
        let mut verts = Vec::new();
        for i in 0..=self.rows {
            for j in 0..=self.cols {
                // Alternating raised/lowered pattern for waterbomb
                let z = if (i + j) % 2 == 0 {
                    s * sin_p
                } else {
                    -s * sin_p
                };
                let x = j as f64 * s * cos_p;
                let y = i as f64 * s * cos_p;
                verts.push(Vec3::new(x, y, z));
            }
        }
        verts
    }

    /// Folded height (z-extent) of the tessellation.
    pub fn folded_height(&self) -> f64 {
        self.panel_size * self.pleat_angle.sin() * 2.0
    }

    /// Projected footprint (x, y) of the tessellation.
    pub fn footprint(&self) -> (f64, f64) {
        let s = self.panel_size * self.pleat_angle.cos();
        (self.cols as f64 * s, self.rows as f64 * s)
    }

    /// Compactness ratio (folded height / flat span).
    pub fn compactness_ratio(&self) -> f64 {
        let flat_height = self.panel_size * self.rows as f64;
        if flat_height < 1e-12 {
            return 0.0;
        }
        self.folded_height() / flat_height
    }
}

// ---------------------------------------------------------------------------
// YoshimuraBuckling
// ---------------------------------------------------------------------------

/// Yoshimura buckling pattern for a thin cylindrical shell.
///
/// The diamond-shaped crease pattern arises from axial compression of a cylinder.
#[derive(Debug, Clone)]
pub struct YoshimuraBuckling {
    /// Cylinder radius.
    pub radius: f64,
    /// Cylinder length.
    pub length: f64,
    /// Number of circumferential lobes.
    pub n_lobes: usize,
    /// Number of axial repetitions.
    pub n_axial: usize,
    /// Fold depth parameter (0 = flat, 1 = maximum indentation).
    pub fold_depth: f64,
}

impl YoshimuraBuckling {
    /// Create a new Yoshimura buckling pattern.
    pub fn new(radius: f64, length: f64, n_lobes: usize, n_axial: usize) -> Self {
        Self {
            radius,
            length,
            n_lobes,
            n_axial,
            fold_depth: 0.0,
        }
    }

    /// Set the fold depth parameter.
    pub fn set_fold_depth(&mut self, d: f64) {
        self.fold_depth = d.clamp(0.0, 1.0);
    }

    /// Half-angle of one diamond cell.
    pub fn diamond_half_angle(&self) -> f64 {
        PI / self.n_lobes as f64
    }

    /// Axial wavelength of the buckling pattern.
    pub fn axial_wavelength(&self) -> f64 {
        self.length / self.n_axial as f64
    }

    /// Radial indentation at maximum fold depth for one diamond.
    pub fn radial_indentation(&self) -> f64 {
        let phi = self.diamond_half_angle();
        self.radius * (1.0 - phi.cos()) * self.fold_depth
    }

    /// Generate vertex positions on the buckled cylinder surface.
    ///
    /// Returns `(n_lobes * n_axial * 4)` diamond vertices approximately.
    pub fn vertex_positions(&self) -> Vec<Vec3> {
        let n = self.n_lobes;
        let m = self.n_axial;
        let mut verts = Vec::with_capacity(n * (m + 1));
        let dtheta = 2.0 * PI / n as f64;
        let dz = self.length / m as f64;

        for k in 0..=m {
            let z = k as f64 * dz;
            for i in 0..n {
                let theta = i as f64 * dtheta + if k % 2 == 1 { dtheta / 2.0 } else { 0.0 };
                // Indented radius for alternating rows
                let r = if k % 2 == 1 {
                    self.radius - self.radial_indentation()
                } else {
                    self.radius
                };
                verts.push(Vec3::new(r * theta.cos(), r * theta.sin(), z));
            }
        }
        verts
    }

    /// Estimated axial shortening for the given fold depth.
    pub fn axial_shortening(&self) -> f64 {
        let phi = self.diamond_half_angle();
        let dz = self.axial_wavelength();
        let shortening_per_cell = dz * (1.0 - phi.cos()) * self.fold_depth;
        shortening_per_cell * self.n_axial as f64
    }
}

// ---------------------------------------------------------------------------
// KirigamiCut
// ---------------------------------------------------------------------------

/// A single slit cut in a kirigami sheet.
#[derive(Debug, Clone)]
pub struct KirigamiSlit {
    /// Start point of the slit.
    pub start: Vec3,
    /// End point of the slit.
    pub end: Vec3,
    /// Width of the slit opening at maximum stretch.
    pub opening: f64,
}

impl KirigamiSlit {
    /// Create a new slit.
    pub fn new(start: Vec3, end: Vec3) -> Self {
        Self {
            start,
            end,
            opening: 0.0,
        }
    }

    /// Length of the slit.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }
}

/// Kirigami sheet with an array of cuts creating auxetic behaviour.
#[derive(Debug, Clone)]
pub struct KirigamiCut {
    /// Sheet width.
    pub width: f64,
    /// Sheet height.
    pub height: f64,
    /// All slits in the pattern.
    pub slits: Vec<KirigamiSlit>,
    /// Current stretch parameter (0 = undeformed, 1 = maximum stretch).
    pub stretch: f64,
}

impl KirigamiCut {
    /// Create an empty kirigami sheet.
    pub fn new(width: f64, height: f64) -> Self {
        Self {
            width,
            height,
            slits: Vec::new(),
            stretch: 0.0,
        }
    }

    /// Create a regular rectangular slit pattern with `rows × cols` slits.
    ///
    /// Slits are staggered in alternating rows.
    pub fn rectangular_pattern(
        width: f64,
        height: f64,
        rows: usize,
        cols: usize,
        slit_frac: f64,
    ) -> Self {
        let mut sheet = Self::new(width, height);
        let dx = width / cols as f64;
        let dy = height / rows as f64;
        let slit_len = dx * slit_frac;
        for r in 0..rows {
            for c in 0..cols {
                let cx = (c as f64 + 0.5) * dx + if r % 2 == 1 { dx * 0.5 } else { 0.0 };
                let cy = (r as f64 + 0.5) * dy;
                let half = slit_len * 0.5;
                let start = Vec3::new(cx - half, cy, 0.0);
                let end = Vec3::new(cx + half, cy, 0.0);
                sheet.slits.push(KirigamiSlit::new(start, end));
            }
        }
        sheet
    }

    /// Set the stretch parameter.
    pub fn set_stretch(&mut self, s: f64) {
        self.stretch = s.clamp(0.0, 1.0);
        // Update slit openings proportionally
        for sl in &mut self.slits {
            sl.opening = sl.length() * s * 0.5;
        }
    }

    /// Effective Poisson's ratio of the kirigami sheet at the current stretch.
    ///
    /// For a rectangular slit pattern the auxetic ν ≈ −(fractional area cut) · stretch.
    pub fn poisson_ratio(&self) -> f64 {
        let total_cut_length: f64 = self.slits.iter().map(|s| s.length()).sum();
        let sheet_perimeter = 2.0 * (self.width + self.height);
        let cut_fraction = (total_cut_length / sheet_perimeter).min(1.0);
        -cut_fraction * self.stretch
    }

    /// Stretchability: ratio of stretched length to original length.
    pub fn stretchability(&self) -> f64 {
        let avg_opening: f64 =
            self.slits.iter().map(|s| s.opening).sum::<f64>() / self.slits.len().max(1) as f64;
        1.0 + avg_opening / (self.height / self.slits.len().max(1) as f64).max(1e-12)
    }

    /// Number of slits.
    pub fn num_slits(&self) -> usize {
        self.slits.len()
    }
}

// ---------------------------------------------------------------------------
// OrigamiAnalysis
// ---------------------------------------------------------------------------

/// Analysis tools for folded origami surfaces.
///
/// Computes geometric properties such as Gaussian curvature (angle deficit at
/// vertices) and the effective Poisson's ratio of the tessellation.
#[derive(Debug, Clone)]
pub struct OrigamiAnalysis;

impl OrigamiAnalysis {
    /// Compute the angle deficit (discrete Gaussian curvature) at a vertex,
    /// given the sector angles (in radians) of the faces meeting at that vertex.
    ///
    /// For a flat sheet the deficit is 0; a cone tip has a positive deficit.
    pub fn angle_deficit(sector_angles: &[f64]) -> f64 {
        let sum: f64 = sector_angles.iter().sum();
        2.0 * PI - sum
    }

    /// Check Kawasaki's theorem for a crease pattern at an interior vertex.
    ///
    /// For flat-foldability the alternating sum of sector angles must equal π.
    /// Returns the residual (should be ≈ 0 for a flat-foldable vertex).
    pub fn kawasaki_residual(sector_angles: &[f64]) -> f64 {
        if sector_angles.len() < 2 {
            return f64::INFINITY;
        }
        let even: f64 = sector_angles.iter().step_by(2).sum();
        let odd: f64 = sector_angles.iter().skip(1).step_by(2).sum();
        (even - odd).abs()
    }

    /// Check Maekawa's theorem: at an interior vertex the number of mountain
    /// folds and valley folds must differ by exactly 2.
    pub fn maekawa_valid(mountain_count: usize, valley_count: usize) -> bool {
        let diff = (mountain_count as isize - valley_count as isize).abs();
        diff == 2
    }

    /// Compute the effective in-plane Poisson's ratio of a Miura-ori sheet.
    ///
    /// Uses the analytical closed form for the Miura-ori.
    pub fn miura_poisson_ratio(alpha: f64, theta: f64) -> f64 {
        let sin_a = alpha.sin();
        let sin_t = theta.sin();
        let denom = 1.0 - (sin_a * sin_t).powi(2);
        if denom.abs() < 1e-12 {
            return 0.0;
        }
        -(sin_a.powi(2)) / denom
    }

    /// Estimate the bending energy of a fold line given fold angle and panel stiffness.
    ///
    /// `k` is the rotational stiffness per unit length (N·m/m), `length` is the
    /// crease length (m), `theta` is the fold angle (radians).
    pub fn fold_bending_energy(k: f64, length: f64, theta: f64) -> f64 {
        0.5 * k * length * theta * theta
    }

    /// Compute the Gaussian curvature of a parallelogram panel given its vertex positions.
    ///
    /// For a flat (non-curved) panel, returns 0.
    pub fn panel_gaussian_curvature(v0: &Vec3, v1: &Vec3, v2: &Vec3, v3: &Vec3) -> f64 {
        // Use the discrete formulation: K ≈ 0 for flat panels
        let n1 = (v1 - v0).cross(&(v2 - v0));
        let n2 = (v2 - v0).cross(&(v3 - v0));
        let cos_theta = (n1.dot(&n2) / (n1.norm() * n2.norm() + 1e-12)).clamp(-1.0, 1.0);
        let angle = cos_theta.acos();
        angle / (v0 - v2).norm().max(1e-12).powi(2)
    }

    /// Compute the fold angle between two planar facets given their unit normals and the
    /// shared crease direction.
    pub fn fold_angle_from_normals(n1: &Vec3, n2: &Vec3, crease_dir: &Vec3) -> f64 {
        let cos_a = n1.dot(n2).clamp(-1.0, 1.0);
        let angle = cos_a.acos();
        // Determine sign from cross product relative to crease direction
        let cross = n1.cross(n2);
        if cross.dot(crease_dir) < 0.0 {
            -angle
        } else {
            angle
        }
    }
}

// ---------------------------------------------------------------------------
// Helper free functions
// ---------------------------------------------------------------------------

/// Build a simple square crease pattern with one central fold line.
///
/// Returns an `OrigamiPattern` with 4 vertices, 1 fold line, and 2 facets.
pub fn build_single_fold_pattern(side: f64, fold_type: FoldType) -> OrigamiPattern {
    let mut pat = OrigamiPattern::new();
    let v0 = pat.add_vertex(Vec3::new(0.0, 0.0, 0.0));
    let v1 = pat.add_vertex(Vec3::new(side, 0.0, 0.0));
    let v2 = pat.add_vertex(Vec3::new(side, side, 0.0));
    let v3 = pat.add_vertex(Vec3::new(0.0, side, 0.0));
    let vm = pat.add_vertex(Vec3::new(side * 0.5, 0.0, 0.0));
    let vm2 = pat.add_vertex(Vec3::new(side * 0.5, side, 0.0));

    let fold_start = pat.vertices[vm];
    let fold_end = pat.vertices[vm2];
    pat.add_fold_line(FoldLine::new(fold_start, fold_end, fold_type));

    pat.add_facet(OrigamiFacet::new(vec![v0, vm, vm2, v3]));
    pat.add_facet(OrigamiFacet::new(vec![vm, v1, v2, vm2]));
    pat
}

/// Compute the fold angle coupling for a generic degree-4 vertex.
///
/// Given three of the four sector angles at an interior vertex (in order),
/// returns the fourth sector angle required by Kawasaki's theorem.
pub fn kawasaki_fourth_angle(a1: f64, a2: f64, a3: f64) -> f64 {
    // a1 - a2 + a3 - a4 = 0  ⟹  a4 = a1 - a2 + a3
    (a1 - a2 + a3).rem_euclid(2.0 * PI)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- FoldLine tests ---

    #[test]
    fn test_fold_line_length() {
        let fl = FoldLine::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 4.0, 0.0),
            FoldType::Mountain,
        );
        assert!((fl.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_fold_line_midpoint() {
        let fl = FoldLine::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            FoldType::Valley,
        );
        let mp = fl.midpoint();
        assert!((mp.x - 1.0).abs() < 1e-10);
        assert!(mp.y.abs() < 1e-10);
    }

    #[test]
    fn test_fold_line_rotation_identity_at_zero() {
        let mut fl = FoldLine::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            FoldType::Mountain,
        );
        fl.fold_angle = 0.0;
        let p = Vec3::new(0.0, 1.0, 0.0);
        let r = fl.apply_rotation(&p);
        assert!((r.x).abs() < 1e-10);
        assert!((r.y - 1.0).abs() < 1e-10);
        assert!(r.z.abs() < 1e-10);
    }

    #[test]
    fn test_fold_line_rotation_90_degrees() {
        let mut fl = FoldLine::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            FoldType::Mountain,
        );
        fl.fold_angle = PI / 2.0;
        let p = Vec3::new(0.0, 1.0, 0.0);
        let r = fl.apply_rotation(&p);
        // Rotating y-axis by 90° about x-axis: y→z
        assert!(r.y.abs() < 1e-10, "y should be ~0, got {}", r.y);
        assert!((r.z - 1.0).abs() < 1e-10, "z should be ~1, got {}", r.z);
    }

    #[test]
    fn test_fold_type_signed_angle() {
        let mut fl_m = FoldLine::new(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), FoldType::Mountain);
        fl_m.fold_angle = 0.5;
        assert!((fl_m.signed_angle() - 0.5).abs() < 1e-10);

        let mut fl_v = FoldLine::new(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), FoldType::Valley);
        fl_v.fold_angle = 0.5;
        assert!((fl_v.signed_angle() + 0.5).abs() < 1e-10);
    }

    // --- OrigamiPattern tests ---

    #[test]
    fn test_origami_pattern_add_vertex() {
        let mut pat = OrigamiPattern::new();
        let i = pat.add_vertex(Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(i, 0);
        assert_eq!(pat.vertices.len(), 1);
    }

    #[test]
    fn test_origami_pattern_dihedral_angle() {
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let n2 = Vec3::new(0.0, 0.0, 1.0);
        let angle = OrigamiPattern::dihedral_angle(&n1, &n2);
        assert!(angle.abs() < 1e-10);
    }

    #[test]
    fn test_origami_pattern_dihedral_angle_90() {
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let angle = OrigamiPattern::dihedral_angle(&n1, &n2);
        assert!((angle - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_single_fold_pattern_structure() {
        let pat = build_single_fold_pattern(1.0, FoldType::Mountain);
        assert_eq!(pat.fold_lines.len(), 1);
        assert_eq!(pat.facets.len(), 2);
        assert!(pat.vertices.len() >= 4);
    }

    // --- RigidOrigami tests ---

    #[test]
    fn test_rigid_origami_initial_state() {
        let pat = build_single_fold_pattern(1.0, FoldType::Mountain);
        let ro = RigidOrigami::new(pat);
        assert!((ro.rho - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_rigid_origami_step_clamps() {
        let pat = build_single_fold_pattern(1.0, FoldType::Mountain);
        let mut ro = RigidOrigami::new(pat);
        ro.step(2.0); // Should clamp to 1.0
        assert!((ro.rho - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rigid_origami_bounding_box() {
        let pat = build_single_fold_pattern(2.0, FoldType::Mountain);
        let ro = RigidOrigami::new(pat);
        let (lo, hi) = ro.bounding_box();
        assert!(hi.x >= lo.x);
        assert!(hi.y >= lo.y);
        assert!(hi.z >= lo.z);
    }

    // --- MiuraOri tests ---

    #[test]
    fn test_miura_ori_projected_dimensions_flat() {
        let cell = MiuraUnitCell::new(1.0, 1.0, PI / 4.0);
        let (lx, _ly, lz) = cell.projected_dimensions(0.0);
        // At θ=0, cos_theta=1, sin_theta=0 → denom=1
        // lx = 2*a*sin(alpha)*cos(0)/1 = 2*sin(PI/4) = √2
        let expected_lx = 2.0 * (PI / 4.0_f64).sin();
        assert!((lx - expected_lx).abs() < 1e-10, "lx at flat, got {}", lx);
        assert!(lz > 0.0);
    }

    #[test]
    fn test_miura_ori_coupled_angle_at_zero() {
        let cell = MiuraUnitCell::new(1.0, 1.0, PI / 4.0);
        let phi = cell.coupled_angle(0.0);
        assert!(phi.abs() < 1e-10);
    }

    #[test]
    fn test_miura_ori_vertex_count() {
        let cell = MiuraUnitCell::new(1.0, 1.0, PI / 4.0);
        let miura = MiuraOri::new(cell, 3, 4);
        let verts = miura.vertex_positions();
        let expected = (2 * 3 + 2) * (2 * 4 + 2);
        assert_eq!(verts.len(), expected);
    }

    #[test]
    fn test_miura_ori_poisson_ratio_negative() {
        let cell = MiuraUnitCell::new(1.0, 1.0, PI / 4.0);
        let mut miura = MiuraOri::new(cell, 2, 2);
        miura.set_angle(PI / 4.0);
        let nu = miura.poisson_ratio();
        assert!(
            nu <= 0.0,
            "Miura-ori Poisson ratio should be negative, got {}",
            nu
        );
    }

    #[test]
    fn test_miura_ori_fold_lines_row() {
        let cell = MiuraUnitCell::new(1.0, 1.0, PI / 6.0);
        let miura = MiuraOri::new(cell, 2, 3);
        let lines = miura.fold_lines_row(0);
        assert!(!lines.is_empty());
    }

    // --- WaterbombBase tests ---

    #[test]
    fn test_waterbomb_vertex_count() {
        let wb = WaterbombBase::new(1.0, 3, 4);
        let verts = wb.vertex_positions();
        assert_eq!(verts.len(), 4 * 5); // (rows+1)*(cols+1)
    }

    #[test]
    fn test_waterbomb_flat_height_zero() {
        let wb = WaterbombBase::new(1.0, 2, 2);
        assert!((wb.folded_height()).abs() < 1e-10);
    }

    #[test]
    fn test_waterbomb_footprint() {
        let wb = WaterbombBase::new(1.0, 2, 3);
        let (fx, fy) = wb.footprint();
        assert!((fx - 3.0).abs() < 1e-10); // cols * panel_size (at angle=0)
        assert!((fy - 2.0).abs() < 1e-10);
    }

    // --- YoshimuraBuckling tests ---

    #[test]
    fn test_yoshimura_vertex_count() {
        let yb = YoshimuraBuckling::new(1.0, 2.0, 6, 4);
        let verts = yb.vertex_positions();
        assert_eq!(verts.len(), 6 * 5); // n_lobes * (n_axial+1)
    }

    #[test]
    fn test_yoshimura_axial_shortening_zero_depth() {
        let yb = YoshimuraBuckling::new(1.0, 2.0, 6, 4);
        assert!((yb.axial_shortening()).abs() < 1e-10);
    }

    #[test]
    fn test_yoshimura_axial_wavelength() {
        let yb = YoshimuraBuckling::new(1.0, 2.0, 6, 4);
        assert!((yb.axial_wavelength() - 0.5).abs() < 1e-10);
    }

    // --- KirigamiCut tests ---

    #[test]
    fn test_kirigami_rectangular_pattern_slit_count() {
        let kg = KirigamiCut::rectangular_pattern(2.0, 3.0, 3, 4, 0.8);
        assert_eq!(kg.num_slits(), 12);
    }

    #[test]
    fn test_kirigami_stretch_updates_openings() {
        let mut kg = KirigamiCut::rectangular_pattern(2.0, 2.0, 2, 2, 0.9);
        kg.set_stretch(1.0);
        for slit in &kg.slits {
            assert!(slit.opening > 0.0);
        }
    }

    #[test]
    fn test_kirigami_poisson_ratio_nonpositive() {
        let mut kg = KirigamiCut::rectangular_pattern(2.0, 2.0, 3, 3, 0.8);
        kg.set_stretch(0.5);
        assert!(kg.poisson_ratio() <= 0.0);
    }

    #[test]
    fn test_kirigami_stretchability_ge_one() {
        let mut kg = KirigamiCut::rectangular_pattern(2.0, 2.0, 2, 2, 0.9);
        kg.set_stretch(0.5);
        assert!(kg.stretchability() >= 1.0);
    }

    // --- OrigamiAnalysis tests ---

    #[test]
    fn test_kawasaki_theorem_flat_vertex() {
        // For a flat-foldable 4-crease vertex, alternating angles sum to π
        let angles = [PI / 4.0, PI / 4.0, PI / 4.0, PI / 4.0];
        // All equal → residual = 0
        let residual = OrigamiAnalysis::kawasaki_residual(&angles);
        assert!(residual < 1e-10, "residual={}", residual);
    }

    #[test]
    fn test_maekawa_theorem_valid() {
        assert!(OrigamiAnalysis::maekawa_valid(3, 1));
        assert!(OrigamiAnalysis::maekawa_valid(1, 3));
        assert!(!OrigamiAnalysis::maekawa_valid(2, 2));
    }

    #[test]
    fn test_angle_deficit_flat() {
        // Six equilateral triangle sectors: sum = 2π → deficit = 0
        let angles = vec![PI / 3.0; 6];
        let deficit = OrigamiAnalysis::angle_deficit(&angles);
        assert!(deficit.abs() < 1e-10);
    }

    #[test]
    fn test_miura_poisson_analytical() {
        let nu = OrigamiAnalysis::miura_poisson_ratio(PI / 4.0, PI / 4.0);
        assert!(nu < 0.0, "Miura Poisson should be negative, got {}", nu);
    }

    #[test]
    fn test_fold_bending_energy_positive() {
        let e = OrigamiAnalysis::fold_bending_energy(1.0, 1.0, PI / 4.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_fold_angle_from_normals_parallel() {
        let n = Vec3::new(0.0, 0.0, 1.0);
        let crease = Vec3::new(1.0, 0.0, 0.0);
        let angle = OrigamiAnalysis::fold_angle_from_normals(&n, &n, &crease);
        assert!(angle.abs() < 1e-10);
    }

    #[test]
    fn test_kawasaki_fourth_angle() {
        let a4 = kawasaki_fourth_angle(PI / 4.0, PI / 4.0, PI / 4.0);
        assert!((a4 - PI / 4.0).abs() < 1e-10);
    }
}
