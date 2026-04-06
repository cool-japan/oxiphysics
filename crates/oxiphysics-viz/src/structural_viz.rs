// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Structural visualization: finite element mesh display, deformed-shape
//! overlays, displacement colormaps, strain energy density, mode shapes,
//! stress-tensor glyphs, crack paths, damage fields, fatigue contours, and
//! convergence history plots.
//!
//! # Design notes
//!
//! All types work with plain `f64` arrays (`[f64; 3]` for 3-D positions).
//! No external linear-algebra crate is required.  A downstream renderer
//! consumes the produced [`MeshWireframe`], [`ColoredMesh`], [`GlyphSet`],
//! and [`LinePlot`] data structures.

#![allow(dead_code)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Vec3 arithmetic helpers (local, to avoid any external dep)
// ---------------------------------------------------------------------------

/// Add two 3-D vectors.
#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors.
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector.
#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euclidean length of a 3-D vector.
#[inline]
fn vec3_length(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalise a 3-D vector (returns zero vector if near-zero length).
#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let len = vec3_length(a);
    if len < 1e-15 {
        [0.0; 3]
    } else {
        [a[0] / len, a[1] / len, a[2] / len]
    }
}

/// Dot product.
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// RGBA colour
// ---------------------------------------------------------------------------

/// An RGBA colour with components in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}

impl Rgba {
    /// Construct a colour from RGBA components in `[0, 1]`.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }

    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }

    /// Opaque red.
    pub fn red() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }

    /// Opaque blue.
    pub fn blue() -> Self {
        Self::new(0.0, 0.0, 1.0, 1.0)
    }

    /// Opaque green.
    pub fn green() -> Self {
        Self::new(0.0, 1.0, 0.0, 1.0)
    }

    /// Linear interpolation between two colours.
    pub fn lerp(a: Self, b: Self, t: f32) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self::new(
            a.r + (b.r - a.r) * t,
            a.g + (b.g - a.g) * t,
            a.b + (b.b - a.b) * t,
            a.a + (b.a - a.a) * t,
        )
    }
}

// ---------------------------------------------------------------------------
// Colourmap
// ---------------------------------------------------------------------------

/// Built-in colormaps for scalar field visualisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructColormap {
    /// Blue–cyan–green–yellow–red (classic *jet*).
    Jet,
    /// Perceptually uniform dark-to-bright (viridis-inspired).
    Viridis,
    /// Blue (cold) to red (hot).
    CoolWarm,
    /// Monochrome grey scale.
    Grayscale,
}

/// Map a normalised scalar `t ∈ [0, 1]` to an [`Rgba`] colour.
///
/// Values outside `[0, 1]` are clamped.
pub fn map_scalar_to_color(t: f64, cmap: StructColormap) -> Rgba {
    let t = t.clamp(0.0, 1.0) as f32;
    match cmap {
        StructColormap::Grayscale => Rgba::new(t, t, t, 1.0),
        StructColormap::CoolWarm => Rgba::lerp(Rgba::blue(), Rgba::red(), t),
        StructColormap::Jet => jet_color(t),
        StructColormap::Viridis => viridis_color(t),
    }
}

/// Jet colour at normalised position `t`.
fn jet_color(t: f32) -> Rgba {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    Rgba::new(r, g, b, 1.0)
}

/// Simplified viridis colour at normalised position `t`.
fn viridis_color(t: f32) -> Rgba {
    // Three-stop approximation: dark-blue → teal → yellow
    let stops = [
        (0.0_f32, Rgba::new(0.267_f32, 0.005_f32, 0.329_f32, 1.0)),
        (0.5_f32, Rgba::new(0.128_f32, 0.566_f32, 0.551_f32, 1.0)),
        (1.0_f32, Rgba::new(0.993_f32, 0.906_f32, 0.144_f32, 1.0)),
    ];
    if t <= stops[0].0 {
        return stops[0].1;
    }
    if t >= stops[2].0 {
        return stops[2].1;
    }
    for i in 0..stops.len() - 1 {
        let (t0, c0) = stops[i];
        let (t1, c1) = stops[i + 1];
        if t <= t1 {
            let local = (t - t0) / (t1 - t0);
            return Rgba::lerp(c0, c1, local);
        }
    }
    stops[2].1
}

/// Map raw scalar values to colours by normalising to `[min, max]`.
///
/// If `min == max` the whole array is mapped to `t = 0.5`.
pub fn scalars_to_colors(values: &[f64], cmap: StructColormap) -> Vec<Rgba> {
    if values.is_empty() {
        return vec![];
    }
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max - min;
    values
        .iter()
        .map(|&v| {
            let t = if range.abs() < 1e-15 {
                0.5
            } else {
                (v - min) / range
            };
            map_scalar_to_color(t, cmap)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// FEM mesh types
// ---------------------------------------------------------------------------

/// Type of a finite element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType {
    /// Two-node line element (beam / truss).
    Line2,
    /// Three-node linear triangle.
    Tri3,
    /// Four-node bilinear quadrilateral.
    Quad4,
    /// Four-node linear tetrahedron.
    Tet4,
    /// Eight-node trilinear hexahedron.
    Hex8,
    /// Six-node linear prism (wedge).
    Prism6,
}

impl ElementType {
    /// Number of nodes per element.
    pub fn node_count(self) -> usize {
        match self {
            ElementType::Line2 => 2,
            ElementType::Tri3 => 3,
            ElementType::Quad4 => 4,
            ElementType::Tet4 => 4,
            ElementType::Hex8 => 8,
            ElementType::Prism6 => 6,
        }
    }

    /// Return `true` if this is a 3-D volumetric element.
    pub fn is_volumetric(self) -> bool {
        matches!(
            self,
            ElementType::Tet4 | ElementType::Hex8 | ElementType::Prism6
        )
    }
}

/// A finite element mesh: nodes and element connectivity.
#[derive(Debug, Clone, Default)]
pub struct FemMesh {
    /// Node positions, indexed 0 … N-1.
    pub nodes: Vec<[f64; 3]>,
    /// Element connectivity: each entry is `(element_type, node_indices)`.
    pub elements: Vec<(ElementType, Vec<usize>)>,
    /// Optional element-level integer tags (e.g. material id).
    pub element_tags: Vec<i32>,
    /// Optional node-level group ids.
    pub node_groups: Vec<i32>,
}

impl FemMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node, returning its index.
    pub fn add_node(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(pos);
        idx
    }

    /// Add an element, returning its index.
    pub fn add_element(&mut self, etype: ElementType, nodes: Vec<usize>) -> usize {
        let idx = self.elements.len();
        self.elements.push((etype, nodes));
        idx
    }

    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Total number of elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }

    /// Compute the axis-aligned bounding box of the mesh.
    ///
    /// Returns `None` if the mesh has no nodes.
    pub fn aabb(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut mn = self.nodes[0];
        let mut mx = self.nodes[0];
        for &n in &self.nodes[1..] {
            for k in 0..3 {
                if n[k] < mn[k] {
                    mn[k] = n[k];
                }
                if n[k] > mx[k] {
                    mx[k] = n[k];
                }
            }
        }
        Some((mn, mx))
    }

    /// Centroid of element `idx` (average of its node positions).
    pub fn element_centroid(&self, idx: usize) -> Option<[f64; 3]> {
        let (_, conn) = self.elements.get(idx)?;
        if conn.is_empty() {
            return None;
        }
        let mut sum = [0.0_f64; 3];
        for &ni in conn {
            if let Some(&p) = self.nodes.get(ni) {
                for k in 0..3 {
                    sum[k] += p[k];
                }
            }
        }
        let n = conn.len() as f64;
        Some([sum[0] / n, sum[1] / n, sum[2] / n])
    }
}

// ---------------------------------------------------------------------------
// MeshWireframe
// ---------------------------------------------------------------------------

/// A wireframe representation of a FEM mesh: a set of line segments.
#[derive(Debug, Clone, Default)]
pub struct MeshWireframe {
    /// Line segments as pairs of 3-D endpoints.
    pub segments: Vec<([f64; 3], [f64; 3])>,
    /// Per-segment colour (same length as `segments`).
    pub colors: Vec<Rgba>,
}

impl MeshWireframe {
    /// Create an empty wireframe.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a coloured line segment.
    pub fn add_segment(&mut self, a: [f64; 3], b: [f64; 3], color: Rgba) {
        self.segments.push((a, b));
        self.colors.push(color);
    }

    /// Number of line segments.
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Return `true` if the wireframe has no segments.
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }
}

/// Build a wireframe for a FEM mesh, drawing one edge per element edge.
///
/// Duplicate edges (shared by adjacent elements) are both included.
/// Pass a uniform `color` for all edges.
pub fn fem_mesh_wireframe(mesh: &FemMesh, color: Rgba) -> MeshWireframe {
    let mut wf = MeshWireframe::new();
    for (_, conn) in &mesh.elements {
        let n = conn.len();
        for i in 0..n {
            let j = (i + 1) % n;
            if let (Some(&a), Some(&b)) = (mesh.nodes.get(conn[i]), mesh.nodes.get(conn[j])) {
                wf.add_segment(a, b, color);
            }
        }
    }
    wf
}

// ---------------------------------------------------------------------------
// ColoredMesh — node-coloured mesh for scalar field display
// ---------------------------------------------------------------------------

/// A triangle mesh with per-vertex colours.
#[derive(Debug, Clone, Default)]
pub struct ColoredMesh {
    /// Vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Per-vertex colours.
    pub colors: Vec<Rgba>,
    /// Triangle indices (triplets).
    pub indices: Vec<usize>,
}

impl ColoredMesh {
    /// Create an empty coloured mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
}

/// Build a [`ColoredMesh`] for surface triangles in a FEM mesh, coloured by a
/// per-node scalar field.
///
/// Only `Tri3` elements are tessellated; other element types are skipped.
/// `scalar` must have one value per node in `mesh.nodes`.
#[allow(clippy::too_many_arguments)]
pub fn fem_surface_colored(mesh: &FemMesh, scalar: &[f64], cmap: StructColormap) -> ColoredMesh {
    let colors = scalars_to_colors(scalar, cmap);
    let mut cm = ColoredMesh::new();
    cm.positions = mesh.nodes.clone();
    cm.colors = if colors.len() == mesh.nodes.len() {
        colors
    } else {
        vec![Rgba::white(); mesh.nodes.len()]
    };
    for (etype, conn) in &mesh.elements {
        if *etype == ElementType::Tri3 && conn.len() == 3 {
            cm.indices.extend_from_slice(&[conn[0], conn[1], conn[2]]);
        }
    }
    cm
}

// ---------------------------------------------------------------------------
// Deformed shape overlay
// ---------------------------------------------------------------------------

/// Overlay data containing both original and deformed node positions.
#[derive(Debug, Clone)]
pub struct DeformedOverlay {
    /// Original (reference) node positions.
    pub original: Vec<[f64; 3]>,
    /// Deformed node positions.
    pub deformed: Vec<[f64; 3]>,
    /// Scale factor applied to displacements.
    pub scale: f64,
}

impl DeformedOverlay {
    /// Construct from a mesh and a displacement field.
    ///
    /// `displacements` must have one `[f64; 3]` entry per node.
    /// `scale` amplifies the displacements for clarity.
    pub fn new(mesh: &FemMesh, displacements: &[[f64; 3]], scale: f64) -> Self {
        let original = mesh.nodes.clone();
        let deformed: Vec<[f64; 3]> = original
            .iter()
            .zip(displacements.iter().chain(std::iter::repeat(&[0.0; 3])))
            .map(|(&orig, &disp)| vec3_add(orig, vec3_scale(disp, scale)))
            .collect();
        Self {
            original,
            deformed,
            scale,
        }
    }

    /// Build a wireframe of the **deformed** shape using the element connectivity.
    pub fn deformed_wireframe(&self, mesh: &FemMesh, color: Rgba) -> MeshWireframe {
        let mut wf = MeshWireframe::new();
        for (_, conn) in &mesh.elements {
            let n = conn.len();
            for i in 0..n {
                let j = (i + 1) % n;
                if let (Some(&a), Some(&b)) =
                    (self.deformed.get(conn[i]), self.deformed.get(conn[j]))
                {
                    wf.add_segment(a, b, color);
                }
            }
        }
        wf
    }

    /// Maximum displacement magnitude across all nodes.
    pub fn max_displacement(&self) -> f64 {
        self.original
            .iter()
            .zip(self.deformed.iter())
            .map(|(&o, &d)| vec3_length(vec3_sub(d, o)) / self.scale.abs().max(1e-15))
            .fold(0.0_f64, f64::max)
    }
}

/// Compute per-node displacement magnitudes from a displacement array.
///
/// Returns one `f64` per node.
pub fn displacement_magnitudes(displacements: &[[f64; 3]]) -> Vec<f64> {
    displacements.iter().map(|&d| vec3_length(d)).collect()
}

// ---------------------------------------------------------------------------
// Strain energy density
// ---------------------------------------------------------------------------

/// Compute element strain energy density from a per-element stress and strain
/// Voigt vector (σ, ε each as `[f64; 6]` in order `xx,yy,zz,xy,yz,xz`).
///
/// SED = 0.5 * σ:ε
pub fn strain_energy_density(stress: &[f64; 6], strain: &[f64; 6]) -> f64 {
    // Normal components × 1, shear components × 2 (engineering shear strain)
    0.5 * (stress[0] * strain[0]
        + stress[1] * strain[1]
        + stress[2] * strain[2]
        + 2.0 * stress[3] * strain[3]
        + 2.0 * stress[4] * strain[4]
        + 2.0 * stress[5] * strain[5])
}

/// Map a list of per-element SED values to colours.
pub fn sed_colors(sed_values: &[f64], cmap: StructColormap) -> Vec<Rgba> {
    scalars_to_colors(sed_values, cmap)
}

// ---------------------------------------------------------------------------
// Von Mises stress
// ---------------------------------------------------------------------------

/// Compute von Mises stress from a Voigt stress vector
/// `[σxx, σyy, σzz, σxy, σyz, σxz]`.
pub fn von_mises(stress: &[f64; 6]) -> f64 {
    let [sxx, syy, szz, sxy, syz, sxz] = *stress;
    let v = 0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2))
        + 3.0 * (sxy * sxy + syz * syz + sxz * sxz);
    v.sqrt()
}

/// Map per-element von Mises stress to colours.
pub fn von_mises_colors(stresses: &[[f64; 6]], cmap: StructColormap) -> Vec<Rgba> {
    let vm: Vec<f64> = stresses.iter().map(von_mises).collect();
    scalars_to_colors(&vm, cmap)
}

// ---------------------------------------------------------------------------
// Stress tensor glyphs
// ---------------------------------------------------------------------------

/// A single stress-tensor glyph represented as three principal-axis line
/// segments, colour-coded by sign (tension = red, compression = blue).
#[derive(Debug, Clone)]
pub struct StressGlyph {
    /// Centre position (element centroid).
    pub center: [f64; 3],
    /// Principal stress axes as `(direction, magnitude)`.
    /// Positive magnitude → tensile; negative → compressive.
    pub axes: [([f64; 3], f64); 3],
}

impl StressGlyph {
    /// Return the three line segments (start, end) for rendering.
    /// Each segment spans ±`|magnitude| * scale / 2` along its axis.
    pub fn segments(&self, scale: f64) -> [([f64; 3], [f64; 3]); 3] {
        std::array::from_fn(|i| {
            let (dir, mag) = self.axes[i];
            let half = vec3_scale(dir, mag.abs() * scale * 0.5);
            (vec3_sub(self.center, half), vec3_add(self.center, half))
        })
    }

    /// Return the colour for axis `i`: red if tensile, blue if compressive.
    pub fn axis_color(&self, i: usize) -> Rgba {
        if self.axes[i].1 >= 0.0 {
            Rgba::red()
        } else {
            Rgba::blue()
        }
    }
}

/// Compute eigenvectors and eigenvalues of the 3×3 symmetric stress matrix
/// using the analytical method for symmetric matrices.
///
/// The stress Voigt vector is `[σxx, σyy, σzz, σxy, σyz, σxz]`.
/// Returns three `(eigenvalue, eigenvector)` pairs sorted by descending eigenvalue.
pub fn stress_principal(stress: &[f64; 6]) -> [(f64, [f64; 3]); 3] {
    // Build 3×3 matrix
    let m = [
        [stress[0], stress[3], stress[5]],
        [stress[3], stress[1], stress[4]],
        [stress[5], stress[4], stress[2]],
    ];

    // Cardano / Jacobi iteration (3×3 symmetric, single pass is sufficient
    // for moderate precision)
    let eigenvalues = eigenvalues_sym3(m);
    let mut result = std::array::from_fn(|i| (eigenvalues[i], eigenvector_sym3(m, eigenvalues[i])));
    // Sort descending by eigenvalue
    result.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    result
}

/// Cardano formula for eigenvalues of a 3×3 real symmetric matrix.
fn eigenvalues_sym3(m: [[f64; 3]; 3]) -> [f64; 3] {
    // Standard Cardano method
    let p1 = m[0][1].powi(2) + m[0][2].powi(2) + m[1][2].powi(2);
    if p1 < 1e-30 {
        // Already diagonal
        return [m[0][0], m[1][1], m[2][2]];
    }
    let q = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
    let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    // B = (1/p) * (M - q*I)
    let b: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let diag = if i == j { m[i][j] - q } else { m[i][j] };
            diag / p
        })
    });
    let det_b = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    let r = (det_b / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    use std::f64::consts::PI;
    [
        q + 2.0 * p * phi.cos(),
        q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos(),
        q + 2.0 * p * (phi + 4.0 * PI / 3.0).cos(),
    ]
}

/// Compute an approximate eigenvector for `lambda` using the null-space approach.
fn eigenvector_sym3(m: [[f64; 3]; 3], lambda: f64) -> [f64; 3] {
    // Form (M - λI) and find a null-space vector via cross-products of rows
    let a = [
        [m[0][0] - lambda, m[0][1], m[0][2]],
        [m[1][0], m[1][1] - lambda, m[1][2]],
        [m[2][0], m[2][1], m[2][2] - lambda],
    ];
    // Try cross products of pairs of rows; use the longest
    let candidates = [
        vec3_cross(a[0], a[1]),
        vec3_cross(a[1], a[2]),
        vec3_cross(a[0], a[2]),
    ];
    let best = candidates
        .iter()
        .max_by(|u, v| {
            vec3_length(**u)
                .partial_cmp(&vec3_length(**v))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
        .unwrap_or([1.0, 0.0, 0.0]);
    let n = vec3_normalize(best);
    if vec3_length(n) < 0.5 {
        [1.0, 0.0, 0.0]
    } else {
        n
    }
}

/// Build a [`GlyphSet`] of stress-tensor ellipsoid glyphs for a collection of
/// per-element stress tensors.
///
/// `centroids` must have the same length as `stresses`.
pub fn stress_glyphs(centroids: &[[f64; 3]], stresses: &[[f64; 6]], scale: f64) -> GlyphSet {
    let mut gs = GlyphSet::new();
    for (&center, stress) in centroids.iter().zip(stresses.iter()) {
        let principals = stress_principal(stress);
        let axes = std::array::from_fn(|i| (principals[i].1, principals[i].0));
        let glyph = StressGlyph { center, axes };
        gs.glyphs.push(glyph);
        let _ = scale; // scale is passed to segments() at render time
    }
    gs
}

/// A collection of [`StressGlyph`]s with a shared rendering scale.
#[derive(Debug, Clone, Default)]
pub struct GlyphSet {
    /// Individual stress glyphs.
    pub glyphs: Vec<StressGlyph>,
    /// Scale factor applied to glyph half-lengths.
    pub scale: f64,
}

impl GlyphSet {
    /// Create an empty glyph set with unit scale.
    pub fn new() -> Self {
        Self {
            glyphs: Vec::new(),
            scale: 1.0,
        }
    }

    /// Return line segments for all glyphs.
    pub fn all_segments(&self) -> Vec<([f64; 3], [f64; 3], Rgba)> {
        let mut out = Vec::new();
        for g in &self.glyphs {
            for i in 0..3 {
                let segs = g.segments(self.scale);
                out.push((segs[i].0, segs[i].1, g.axis_color(i)));
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Mode shape animation
// ---------------------------------------------------------------------------

/// Data for animating a single vibration mode shape.
#[derive(Debug, Clone)]
pub struct ModeShape {
    /// Mode number (1-indexed).
    pub mode_number: usize,
    /// Natural frequency in Hz.
    pub frequency_hz: f64,
    /// Per-node displacement vectors (eigenvector).
    pub eigenvector: Vec<[f64; 3]>,
}

impl ModeShape {
    /// Construct a mode shape.
    pub fn new(mode_number: usize, frequency_hz: f64, eigenvector: Vec<[f64; 3]>) -> Self {
        Self {
            mode_number,
            frequency_hz,
            eigenvector,
        }
    }

    /// Sample the mode at time `t` (seconds) with amplitude scale `amp`.
    ///
    /// Returns per-node displacement vectors at the requested time.
    pub fn sample(&self, t: f64, amp: f64) -> Vec<[f64; 3]> {
        use std::f64::consts::TAU;
        let phase = (TAU * self.frequency_hz * t).sin() * amp;
        self.eigenvector
            .iter()
            .map(|&d| vec3_scale(d, phase))
            .collect()
    }

    /// Maximum eigenvector component magnitude (for normalisation).
    pub fn max_amplitude(&self) -> f64 {
        self.eigenvector
            .iter()
            .map(|&d| vec3_length(d))
            .fold(0.0_f64, f64::max)
    }

    /// Normalise the eigenvector to unit maximum amplitude.
    pub fn normalize(&mut self) {
        let amp = self.max_amplitude();
        if amp > 1e-15 {
            for d in &mut self.eigenvector {
                *d = vec3_scale(*d, 1.0 / amp);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Crack path visualisation
// ---------------------------------------------------------------------------

/// A crack path represented as an ordered polyline in 3-D.
#[derive(Debug, Clone, Default)]
pub struct CrackPath {
    /// Ordered crack front / path points.
    pub points: Vec<[f64; 3]>,
    /// Local crack-tip opening displacement (CTOD) at each point.
    pub ctod: Vec<f64>,
}

impl CrackPath {
    /// Create an empty crack path.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a point with an associated CTOD value.
    pub fn push(&mut self, point: [f64; 3], ctod: f64) {
        self.points.push(point);
        self.ctod.push(ctod);
    }

    /// Total crack length (sum of segment lengths).
    pub fn total_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| vec3_length(vec3_sub(w[1], w[0])))
            .sum()
    }

    /// Build a wireframe of the crack path.
    pub fn wireframe(&self, color: Rgba) -> MeshWireframe {
        let mut wf = MeshWireframe::new();
        for w in self.points.windows(2) {
            wf.add_segment(w[0], w[1], color);
        }
        wf
    }

    /// Build a CTOD-coloured wireframe.
    pub fn ctod_colored_wireframe(&self, cmap: StructColormap) -> MeshWireframe {
        let mut wf = MeshWireframe::new();
        if self.ctod.len() < 2 || self.points.len() < 2 {
            return wf;
        }
        let min_ctod = self.ctod.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_ctod = self.ctod.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = (max_ctod - min_ctod).max(1e-15);
        for i in 0..self.points.len().min(self.ctod.len()).saturating_sub(1) {
            let t = (self.ctod[i] - min_ctod) / range;
            let color = map_scalar_to_color(t, cmap);
            wf.add_segment(self.points[i], self.points[i + 1], color);
        }
        wf
    }
}

// ---------------------------------------------------------------------------
// Damage field rendering
// ---------------------------------------------------------------------------

/// Damage field: per-node scalar damage variable d ∈ \[0, 1\].
/// d = 0 → undamaged, d = 1 → fully failed.
#[derive(Debug, Clone)]
pub struct DamageField {
    /// Per-node damage values.
    pub damage: Vec<f64>,
}

impl DamageField {
    /// Construct a damage field with all nodes at zero damage.
    pub fn new(node_count: usize) -> Self {
        Self {
            damage: vec![0.0; node_count],
        }
    }

    /// Set damage at node `i`.
    pub fn set(&mut self, i: usize, value: f64) {
        if let Some(d) = self.damage.get_mut(i) {
            *d = value.clamp(0.0, 1.0);
        }
    }

    /// Return per-node colours: white (undamaged) → black (fully damaged).
    pub fn colors(&self) -> Vec<Rgba> {
        self.damage
            .iter()
            .map(|&d| {
                let v = 1.0 - d as f32;
                Rgba::new(v, v, v, 1.0)
            })
            .collect()
    }

    /// Return the fraction of nodes with damage above `threshold`.
    pub fn damaged_fraction(&self, threshold: f64) -> f64 {
        if self.damage.is_empty() {
            return 0.0;
        }
        let count = self.damage.iter().filter(|&&d| d > threshold).count();
        count as f64 / self.damage.len() as f64
    }

    /// Mean damage across all nodes.
    pub fn mean_damage(&self) -> f64 {
        if self.damage.is_empty() {
            return 0.0;
        }
        self.damage.iter().sum::<f64>() / self.damage.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Fatigue life contours
// ---------------------------------------------------------------------------

/// Fatigue-life field: per-node cycle count to failure.
#[derive(Debug, Clone)]
pub struct FatigueField {
    /// Per-node cycle counts (Nf).
    pub cycles: Vec<f64>,
}

impl FatigueField {
    /// Construct a fatigue field.
    pub fn new(cycles: Vec<f64>) -> Self {
        Self { cycles }
    }

    /// Map cycle counts to colours using a log-scale normalisation.
    ///
    /// Nodes with `Nf ≤ 0` are mapped to `t = 0`.
    pub fn log_colors(&self, cmap: StructColormap) -> Vec<Rgba> {
        let logs: Vec<f64> = self
            .cycles
            .iter()
            .map(|&n| if n > 0.0 { n.log10() } else { 0.0 })
            .collect();
        scalars_to_colors(&logs, cmap)
    }

    /// Generate iso-contour levels.
    ///
    /// Returns `count` evenly spaced values between `min(cycles)` and `max(cycles)`.
    pub fn iso_levels(&self, count: usize) -> Vec<f64> {
        if self.cycles.is_empty() || count == 0 {
            return vec![];
        }
        let min = self.cycles.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = self
            .cycles
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        if (max - min).abs() < 1e-15 {
            return vec![min];
        }
        (0..count)
            .map(|i| min + (max - min) * i as f64 / (count - 1).max(1) as f64)
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Convergence history
// ---------------------------------------------------------------------------

/// A single convergence history: residual norm versus iteration number.
#[derive(Debug, Clone, Default)]
pub struct ConvergenceHistory {
    /// Name of this convergence metric (e.g. `"displacement_residual"`).
    pub name: String,
    /// Iteration indices.
    pub iterations: Vec<usize>,
    /// Residual (or objective) values at each iteration.
    pub residuals: Vec<f64>,
    /// Optional convergence tolerance.
    pub tolerance: Option<f64>,
}

impl ConvergenceHistory {
    /// Create an empty history with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Self::default()
        }
    }

    /// Record a residual at the given iteration.
    pub fn push(&mut self, iteration: usize, residual: f64) {
        self.iterations.push(iteration);
        self.residuals.push(residual);
    }

    /// Return `true` if the last residual is below the tolerance.
    pub fn has_converged(&self) -> bool {
        if let (Some(&last), Some(tol)) = (self.residuals.last(), self.tolerance) {
            last < tol
        } else {
            false
        }
    }

    /// Log-transform all residuals (base-10).
    ///
    /// Values ≤ 0 are replaced with `f64::NAN`.
    pub fn log_residuals(&self) -> Vec<f64> {
        self.residuals
            .iter()
            .map(|&r| if r > 0.0 { r.log10() } else { f64::NAN })
            .collect()
    }

    /// Reduction ratio: first residual divided by last (or 1.0 if ≤ 1 point).
    pub fn reduction_ratio(&self) -> f64 {
        if self.residuals.len() < 2 {
            return 1.0;
        }
        let first = self.residuals[0];
        let last = *self
            .residuals
            .last()
            .expect("collection should not be empty");
        if last.abs() < 1e-300 {
            f64::INFINITY
        } else {
            first / last
        }
    }

    /// Number of recorded iterations.
    pub fn len(&self) -> usize {
        self.iterations.len()
    }

    /// Return `true` if no iterations have been recorded.
    pub fn is_empty(&self) -> bool {
        self.iterations.is_empty()
    }
}

// ---------------------------------------------------------------------------
// LinePlot — convergence plot data
// ---------------------------------------------------------------------------

/// A 2-D line plot data structure for convergence history display.
#[derive(Debug, Clone, Default)]
pub struct LinePlot {
    /// Plot title.
    pub title: String,
    /// X-axis label.
    pub x_label: String,
    /// Y-axis label.
    pub y_label: String,
    /// Named series.
    pub series: Vec<PlotSeries>,
}

impl LinePlot {
    /// Create a new empty plot.
    pub fn new(
        title: impl Into<String>,
        x_label: impl Into<String>,
        y_label: impl Into<String>,
    ) -> Self {
        Self {
            title: title.into(),
            x_label: x_label.into(),
            y_label: y_label.into(),
            series: Vec::new(),
        }
    }

    /// Add a series from a convergence history (log scale).
    pub fn add_convergence(&mut self, history: &ConvergenceHistory, color: Rgba) {
        let xs: Vec<f64> = history.iterations.iter().map(|&i| i as f64).collect();
        let ys = history.log_residuals();
        self.series.push(PlotSeries {
            name: history.name.clone(),
            xs,
            ys,
            color,
        });
    }

    /// Total number of data points across all series.
    pub fn total_points(&self) -> usize {
        self.series.iter().map(|s| s.xs.len()).sum()
    }
}

/// A single named series in a [`LinePlot`].
#[derive(Debug, Clone)]
pub struct PlotSeries {
    /// Series name (displayed in legend).
    pub name: String,
    /// X values.
    pub xs: Vec<f64>,
    /// Y values.
    pub ys: Vec<f64>,
    /// Line colour.
    pub color: Rgba,
}

// ---------------------------------------------------------------------------
// StructuralScene — high-level scene aggregator
// ---------------------------------------------------------------------------

/// A complete structural visualisation scene containing all rendering data.
#[derive(Debug, Default)]
pub struct StructuralScene {
    /// Wireframes (undeformed / deformed mesh edges, crack paths, etc.).
    pub wireframes: Vec<MeshWireframe>,
    /// Coloured meshes (surface fields).
    pub colored_meshes: Vec<ColoredMesh>,
    /// Stress / strain glyph sets.
    pub glyph_sets: Vec<GlyphSet>,
    /// Convergence history plots.
    pub plots: Vec<LinePlot>,
    /// Named metadata fields (e.g. `"step"`, `"time"`, `"load_factor"`).
    pub metadata: HashMap<String, f64>,
}

impl StructuralScene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a wireframe to the scene.
    pub fn add_wireframe(&mut self, wf: MeshWireframe) {
        self.wireframes.push(wf);
    }

    /// Add a coloured mesh.
    pub fn add_colored_mesh(&mut self, cm: ColoredMesh) {
        self.colored_meshes.push(cm);
    }

    /// Add a glyph set.
    pub fn add_glyphs(&mut self, gs: GlyphSet) {
        self.glyph_sets.push(gs);
    }

    /// Add a convergence plot.
    pub fn add_plot(&mut self, plot: LinePlot) {
        self.plots.push(plot);
    }

    /// Total number of rendering primitives (rough count).
    pub fn primitive_count(&self) -> usize {
        let wf: usize = self.wireframes.iter().map(|w| w.len()).sum();
        let cm: usize = self.colored_meshes.iter().map(|m| m.triangle_count()).sum();
        let gl: usize = self.glyph_sets.iter().map(|g| g.glyphs.len() * 3).sum();
        wf + cm + gl
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Vec3 helpers ---

    #[test]
    fn test_vec3_add() {
        let r = vec3_add([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((r[0] - 5.0).abs() < 1e-10);
        assert!((r[1] - 7.0).abs() < 1e-10);
        assert!((r[2] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_length_unit() {
        let len = vec3_length([1.0, 0.0, 0.0]);
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_normalize() {
        let n = vec3_normalize([3.0, 0.0, 0.0]);
        assert!((n[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_normalize_zero() {
        let n = vec3_normalize([0.0, 0.0, 0.0]);
        assert!(vec3_length(n) < 1e-10);
    }

    #[test]
    fn test_vec3_dot() {
        let d = vec3_dot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10);
        let d2 = vec3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(d2.abs() < 1e-10);
    }

    #[test]
    fn test_vec3_cross() {
        let c = vec3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-10);
        assert!(c[0].abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
    }

    // --- Rgba ---

    #[test]
    fn test_rgba_lerp_midpoint() {
        let c = Rgba::lerp(Rgba::black(), Rgba::white(), 0.5);
        assert!((c.r - 0.5).abs() < 1e-6);
        assert!((c.g - 0.5).abs() < 1e-6);
        assert!((c.b - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_rgba_lerp_clamp() {
        let c = Rgba::lerp(Rgba::black(), Rgba::white(), 2.0);
        assert!((c.r - 1.0).abs() < 1e-6);
    }

    // --- Colormap ---

    #[test]
    fn test_jet_endpoints() {
        let c0 = map_scalar_to_color(0.0, StructColormap::Jet);
        let c1 = map_scalar_to_color(1.0, StructColormap::Jet);
        // t=0 → blue end (b > r); t=1 → red end (r > b)
        assert!(c0.b > c0.r);
        assert!(c1.r > c1.b);
    }

    #[test]
    fn test_coolwarm_endpoints() {
        let c0 = map_scalar_to_color(0.0, StructColormap::CoolWarm);
        let c1 = map_scalar_to_color(1.0, StructColormap::CoolWarm);
        assert!(c0.b > c0.r);
        assert!(c1.r > c1.b);
    }

    #[test]
    fn test_grayscale_midpoint() {
        let c = map_scalar_to_color(0.5, StructColormap::Grayscale);
        assert!((c.r - 0.5).abs() < 1e-6);
        assert!((c.g - 0.5).abs() < 1e-6);
        assert!((c.b - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_scalars_to_colors_same_value() {
        let vals = vec![5.0; 4];
        let colors = scalars_to_colors(&vals, StructColormap::Jet);
        assert_eq!(colors.len(), 4);
    }

    // --- ElementType ---

    #[test]
    fn test_element_type_node_count() {
        assert_eq!(ElementType::Tri3.node_count(), 3);
        assert_eq!(ElementType::Hex8.node_count(), 8);
        assert_eq!(ElementType::Tet4.node_count(), 4);
    }

    #[test]
    fn test_element_type_is_volumetric() {
        assert!(ElementType::Tet4.is_volumetric());
        assert!(!ElementType::Tri3.is_volumetric());
    }

    // --- FemMesh ---

    #[test]
    fn test_fem_mesh_add_nodes() {
        let mut mesh = FemMesh::new();
        mesh.add_node([0.0, 0.0, 0.0]);
        mesh.add_node([1.0, 0.0, 0.0]);
        assert_eq!(mesh.node_count(), 2);
    }

    #[test]
    fn test_fem_mesh_aabb() {
        let mut mesh = FemMesh::new();
        mesh.add_node([-1.0, -2.0, -3.0]);
        mesh.add_node([4.0, 5.0, 6.0]);
        let (mn, mx) = mesh.aabb().unwrap();
        assert!((mn[0] - (-1.0)).abs() < 1e-9);
        assert!((mx[2] - 6.0).abs() < 1e-9);
    }

    #[test]
    fn test_fem_mesh_aabb_empty() {
        let mesh = FemMesh::new();
        assert!(mesh.aabb().is_none());
    }

    #[test]
    fn test_fem_mesh_element_centroid() {
        let mut mesh = FemMesh::new();
        mesh.add_node([0.0, 0.0, 0.0]);
        mesh.add_node([3.0, 0.0, 0.0]);
        mesh.add_node([0.0, 3.0, 0.0]);
        mesh.add_element(ElementType::Tri3, vec![0, 1, 2]);
        let c = mesh.element_centroid(0).unwrap();
        assert!((c[0] - 1.0).abs() < 1e-9);
        assert!((c[1] - 1.0).abs() < 1e-9);
    }

    // --- MeshWireframe ---

    #[test]
    fn test_wireframe_from_fem_mesh() {
        let mut mesh = FemMesh::new();
        for i in 0..3 {
            mesh.add_node([i as f64, 0.0, 0.0]);
        }
        mesh.add_element(ElementType::Tri3, vec![0, 1, 2]);
        let wf = fem_mesh_wireframe(&mesh, Rgba::white());
        // Tri3 has 3 edges
        assert_eq!(wf.len(), 3);
    }

    // --- DeformedOverlay ---

    #[test]
    fn test_deformed_overlay_positions() {
        let mut mesh = FemMesh::new();
        mesh.add_node([0.0, 0.0, 0.0]);
        mesh.add_node([1.0, 0.0, 0.0]);
        let disps = vec![[0.0, 1.0, 0.0]; 2];
        let overlay = DeformedOverlay::new(&mesh, &disps, 2.0);
        // First node: [0,0,0] + [0,2,0] = [0,2,0]
        assert!((overlay.deformed[0][1] - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_deformed_overlay_max_displacement() {
        let mut mesh = FemMesh::new();
        mesh.add_node([0.0, 0.0, 0.0]);
        let disps = vec![[3.0, 4.0, 0.0]]; // magnitude = 5
        let overlay = DeformedOverlay::new(&mesh, &disps, 1.0);
        let max_d = overlay.max_displacement();
        assert!((max_d - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_displacement_magnitudes() {
        let disps = vec![[3.0, 4.0, 0.0], [0.0, 0.0, 1.0]];
        let mags = displacement_magnitudes(&disps);
        assert!((mags[0] - 5.0).abs() < 1e-9);
        assert!((mags[1] - 1.0).abs() < 1e-9);
    }

    // --- Strain energy density ---

    #[test]
    fn test_sed_uniaxial() {
        let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let strain = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let sed = strain_energy_density(&stress, &strain);
        // 0.5 * 100 * 0.001 = 0.05
        assert!((sed - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_sed_shear() {
        let stress = [0.0, 0.0, 0.0, 50.0, 0.0, 0.0_f64];
        let strain = [0.0, 0.0, 0.0, 0.002, 0.0, 0.0_f64];
        let sed = strain_energy_density(&stress, &strain);
        // 0.5 * 2 * 50 * 0.002 = 0.1
        assert!((sed - 0.1).abs() < 1e-12);
    }

    // --- Von Mises ---

    #[test]
    fn test_von_mises_uniaxial() {
        let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let vm = von_mises(&stress);
        assert!((vm - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_von_mises_hydrostatic_zero() {
        // Hydrostatic: σ_xx=σ_yy=σ_zz=P, all shear=0 → von Mises = 0
        let s = 200.0_f64;
        let stress = [s, s, s, 0.0, 0.0, 0.0_f64];
        let vm = von_mises(&stress);
        assert!(vm < 1e-9);
    }

    // --- Principal stresses ---

    #[test]
    fn test_stress_principal_diagonal() {
        // Diagonal stress: eigenvalues should be the diagonal entries
        let stress = [3.0, 2.0, 1.0, 0.0, 0.0, 0.0_f64];
        let principals = stress_principal(&stress);
        let eigs: Vec<f64> = principals.iter().map(|p| p.0).collect();
        // Sorted descending: [3, 2, 1]
        assert!((eigs[0] - 3.0).abs() < 1e-6);
        assert!((eigs[1] - 2.0).abs() < 1e-6);
        assert!((eigs[2] - 1.0).abs() < 1e-6);
    }

    // --- ModeShape ---

    #[test]
    fn test_mode_shape_sample_zero_time() {
        let evec = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mode = ModeShape::new(1, 10.0, evec);
        let s = mode.sample(0.0, 1.0);
        // sin(0) = 0 → all zeros
        for d in &s {
            assert!(vec3_length(*d) < 1e-12);
        }
    }

    #[test]
    fn test_mode_shape_max_amplitude() {
        let evec = vec![[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]];
        let mode = ModeShape::new(1, 1.0, evec);
        assert!((mode.max_amplitude() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_mode_shape_normalize() {
        let evec = vec![[3.0, 4.0, 0.0]];
        let mut mode = ModeShape::new(1, 1.0, evec);
        mode.normalize();
        assert!((mode.max_amplitude() - 1.0).abs() < 1e-9);
    }

    // --- CrackPath ---

    #[test]
    fn test_crack_path_total_length() {
        let mut crack = CrackPath::new();
        crack.push([0.0, 0.0, 0.0], 0.0);
        crack.push([1.0, 0.0, 0.0], 0.001);
        crack.push([2.0, 0.0, 0.0], 0.002);
        assert!((crack.total_length() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_crack_path_wireframe() {
        let mut crack = CrackPath::new();
        crack.push([0.0, 0.0, 0.0], 0.0);
        crack.push([1.0, 0.0, 0.0], 0.0);
        crack.push([2.0, 0.0, 0.0], 0.0);
        let wf = crack.wireframe(Rgba::red());
        assert_eq!(wf.len(), 2);
    }

    // --- DamageField ---

    #[test]
    fn test_damage_field_fraction() {
        let mut df = DamageField::new(10);
        for i in 0..5 {
            df.set(i, 0.9);
        }
        let frac = df.damaged_fraction(0.5);
        assert!((frac - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_damage_field_mean() {
        let mut df = DamageField::new(4);
        df.set(0, 0.0);
        df.set(1, 0.5);
        df.set(2, 0.5);
        df.set(3, 1.0);
        assert!((df.mean_damage() - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_damage_field_colors_undamaged_white() {
        let df = DamageField::new(3);
        let colors = df.colors();
        for c in &colors {
            assert!((c.r - 1.0).abs() < 1e-6);
        }
    }

    // --- FatigueField ---

    #[test]
    fn test_fatigue_iso_levels() {
        let ff = FatigueField::new(vec![1e4, 1e5, 1e6]);
        let levels = ff.iso_levels(3);
        assert_eq!(levels.len(), 3);
        assert!((levels[0] - 1e4).abs() / 1e4 < 1e-9);
        assert!((levels[2] - 1e6).abs() / 1e6 < 1e-9);
    }

    #[test]
    fn test_fatigue_log_colors_count() {
        let ff = FatigueField::new(vec![1e3, 1e4, 1e5, 1e6]);
        let colors = ff.log_colors(StructColormap::Jet);
        assert_eq!(colors.len(), 4);
    }

    // --- ConvergenceHistory ---

    #[test]
    fn test_convergence_push_and_len() {
        let mut h = ConvergenceHistory::new("res");
        h.push(0, 1.0);
        h.push(1, 0.1);
        h.push(2, 0.01);
        assert_eq!(h.len(), 3);
    }

    #[test]
    fn test_convergence_has_converged() {
        let mut h = ConvergenceHistory::new("res");
        h.tolerance = Some(1e-6);
        h.push(0, 1.0);
        h.push(1, 1e-8);
        assert!(h.has_converged());
    }

    #[test]
    fn test_convergence_not_converged() {
        let mut h = ConvergenceHistory::new("res");
        h.tolerance = Some(1e-6);
        h.push(0, 0.1);
        assert!(!h.has_converged());
    }

    #[test]
    fn test_convergence_reduction_ratio() {
        let mut h = ConvergenceHistory::new("res");
        h.push(0, 1000.0);
        h.push(1, 1.0);
        let ratio = h.reduction_ratio();
        assert!((ratio - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn test_convergence_log_residuals() {
        let mut h = ConvergenceHistory::new("res");
        h.push(0, 100.0);
        h.push(1, 10.0);
        let lr = h.log_residuals();
        assert!((lr[0] - 2.0).abs() < 1e-9);
        assert!((lr[1] - 1.0).abs() < 1e-9);
    }

    // --- LinePlot ---

    #[test]
    fn test_line_plot_add_convergence() {
        let mut plot = LinePlot::new("residuals", "iteration", "log residual");
        let mut h = ConvergenceHistory::new("disp");
        h.push(0, 1.0);
        h.push(1, 0.1);
        plot.add_convergence(&h, Rgba::red());
        assert_eq!(plot.series.len(), 1);
        assert_eq!(plot.total_points(), 2);
    }

    // --- StructuralScene ---

    #[test]
    fn test_structural_scene_primitive_count() {
        let mut scene = StructuralScene::new();
        let mut wf = MeshWireframe::new();
        wf.add_segment([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], Rgba::white());
        scene.add_wireframe(wf);
        assert_eq!(scene.primitive_count(), 1);
    }

    #[test]
    fn test_structural_scene_add_colored_mesh() {
        let mut scene = StructuralScene::new();
        let mut cm = ColoredMesh::new();
        cm.positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        cm.colors = vec![Rgba::white(); 3];
        cm.indices = vec![0, 1, 2];
        scene.add_colored_mesh(cm);
        assert_eq!(scene.primitive_count(), 1);
    }

    // --- GlyphSet ---

    #[test]
    fn test_glyph_set_segments() {
        let principals = stress_principal(&[100.0, 50.0, 25.0, 0.0, 0.0, 0.0]);
        let axes = std::array::from_fn(|i| (principals[i].1, principals[i].0));
        let glyph = StressGlyph {
            center: [0.0; 3],
            axes,
        };
        let segs = glyph.segments(1.0);
        assert_eq!(segs.len(), 3);
        // All segment midpoints should be at origin
        for (a, b) in &segs {
            let mid = vec3_scale(vec3_add(*a, *b), 0.5);
            assert!(vec3_length(mid) < 1e-9, "midpoint not at origin");
        }
    }

    #[test]
    fn test_fem_surface_colored_tri3() {
        let mut mesh = FemMesh::new();
        for i in 0..3 {
            mesh.add_node([i as f64, 0.0, 0.0]);
        }
        mesh.add_element(ElementType::Tri3, vec![0, 1, 2]);
        let scalar = vec![0.0, 0.5, 1.0];
        let cm = fem_surface_colored(&mesh, &scalar, StructColormap::Jet);
        assert_eq!(cm.triangle_count(), 1);
        assert_eq!(cm.vertex_count(), 3);
    }
}
