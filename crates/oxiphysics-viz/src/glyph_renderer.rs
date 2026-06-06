// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Glyph rendering for vector and tensor field visualization.
//!
//! Glyphs are CPU-side geometry descriptors.  Each function returns a list of
//! [`GlyphInstance`] or vertex data arrays that a downstream renderer can
//! upload to the GPU.  No GPU types are used here.

// ──────────────────────────────────────────────────────────────────────────────
// Types
// ──────────────────────────────────────────────────────────────────────────────

/// Shape of a single visualization glyph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlyphType {
    /// Arrow (shaft + cone head).
    Arrow,
    /// Solid cone.
    Cone,
    /// Sphere.
    Sphere,
    /// Axis-aligned cube.
    Cube,
    /// Single line segment.
    Line,
    /// Flat circular disc.
    Disc,
}

/// A single renderable glyph instance placed in 3-D space.
#[derive(Debug, Clone, PartialEq)]
pub struct GlyphInstance {
    /// World-space position of the glyph origin (x, y, z).
    pub position: [f64; 3],
    /// Unit direction vector along the glyph's principal axis.
    pub direction: [f64; 3],
    /// Uniform scale factor.
    pub scale: f64,
    /// RGBA colour (components in \[0, 1\]).
    pub color: [f32; 4],
    /// Shape variant.
    pub glyph_type: GlyphType,
}

impl GlyphInstance {
    /// Construct a new glyph instance.
    pub fn new(
        position: [f64; 3],
        direction: [f64; 3],
        scale: f64,
        color: [f32; 4],
        glyph_type: GlyphType,
    ) -> Self {
        Self {
            position,
            direction,
            scale,
            color,
            glyph_type,
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Arrow glyph
// ──────────────────────────────────────────────────────────────────────────────

/// Vertex of an arrow glyph triangle mesh (position only, as \[f64; 3\]).
#[derive(Debug, Clone, PartialEq)]
pub struct GlyphVertex {
    /// Position in glyph-local space.
    pub pos: [f64; 3],
    /// Outward-facing normal.
    pub normal: [f64; 3],
}

/// Triangle index triple.
pub type TriIdx = [u32; 3];

/// Simple triangle mesh used for glyph geometry.
#[derive(Debug, Clone)]
pub struct GlyphMesh {
    /// Vertex buffer.
    pub vertices: Vec<GlyphVertex>,
    /// Index buffer (each triple is one triangle).
    pub indices: Vec<TriIdx>,
}

/// Build a local-space arrow glyph mesh.
///
/// The arrow points in the +Z direction in local space, with total length 1.
/// The shaft occupies the first `1 - head_fraction` of the length; the
/// cone head the remainder.
///
/// # Parameters
/// - `shaft_radius` — radius of the cylindrical shaft.
/// - `head_radius` — radius of the cone base.
/// - `head_fraction` — fraction of the total length allocated to the cone head.
/// - `segments` — number of polygonal facets around the circumference.
///
/// # Returns
/// A [`GlyphMesh`] in local space (arrow points along +Z, origin at base).
pub fn arrow_glyph(
    shaft_radius: f64,
    head_radius: f64,
    head_fraction: f64,
    segments: u32,
) -> GlyphMesh {
    let segs = segments.max(3) as usize;
    let shaft_len = 1.0 - head_fraction.clamp(0.1, 0.9);
    let _head_len = 1.0 - shaft_len;

    let mut vertices: Vec<GlyphVertex> = Vec::new();
    let mut indices: Vec<TriIdx> = Vec::new();

    let two_pi = std::f64::consts::TAU;

    // Helper: push a ring of vertices at height z with given radius and outward normal scale.
    let push_ring = |verts: &mut Vec<GlyphVertex>, z: f64, r: f64, nz: f64| {
        for k in 0..segs {
            let theta = (k as f64 / segs as f64) * two_pi;
            let (s, c) = theta.sin_cos();
            verts.push(GlyphVertex {
                pos: [r * c, r * s, z],
                normal: [c, s, nz],
            });
        }
    };

    // Shaft bottom ring (z=0)
    let shaft_bot_start = vertices.len() as u32;
    push_ring(&mut vertices, 0.0, shaft_radius, 0.0);

    // Shaft top ring (z=shaft_len)
    let shaft_top_start = vertices.len() as u32;
    push_ring(&mut vertices, shaft_len, shaft_radius, 0.0);

    // Shaft side triangles
    let s = segs as u32;
    for k in 0..s {
        let k1 = (k + 1) % s;
        let b0 = shaft_bot_start + k;
        let b1 = shaft_bot_start + k1;
        let t0 = shaft_top_start + k;
        let t1 = shaft_top_start + k1;
        indices.push([b0, b1, t0]);
        indices.push([b1, t1, t0]);
    }

    // Cone base ring (z=shaft_len)
    let cone_base_start = vertices.len() as u32;
    push_ring(&mut vertices, shaft_len, head_radius, -1.0);

    // Cone tip vertex
    let tip_idx = vertices.len() as u32;
    vertices.push(GlyphVertex {
        pos: [0.0, 0.0, 1.0],
        normal: [0.0, 0.0, 1.0],
    });

    // Cone side triangles
    for k in 0..s {
        let k1 = (k + 1) % s;
        indices.push([cone_base_start + k, cone_base_start + k1, tip_idx]);
    }

    // Shaft bottom cap
    let cap_center_idx = vertices.len() as u32;
    vertices.push(GlyphVertex {
        pos: [0.0, 0.0, 0.0],
        normal: [0.0, 0.0, -1.0],
    });
    for k in 0..s {
        let k1 = (k + 1) % s;
        indices.push([cap_center_idx, shaft_bot_start + k1, shaft_bot_start + k]);
    }

    GlyphMesh { vertices, indices }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tensor ellipsoid
// ──────────────────────────────────────────────────────────────────────────────

/// Build an ellipsoid mesh from eigenvalues and eigenvectors of a symmetric tensor.
///
/// The ellipsoid semi-axes are `|eigenvalue[i]|`; each axis is aligned with the
/// corresponding eigenvector.
///
/// # Parameters
/// - `eigenvalues` — three eigenvalues `[λ0, λ1, λ2]`.
/// - `eigenvectors` — corresponding unit eigenvectors as rows of a 3×3 array
///   `[[e0x,e0y,e0z\], [e1x,e1y,e1z], [e2x,e2y,e2z]]`.
/// - `center` — world-space center of the ellipsoid.
/// - `stacks` — number of latitude divisions.
/// - `slices` — number of longitude divisions.
///
/// Returns a list of world-space vertex positions `[f64; 3]`.
pub fn tensor_ellipsoid(
    eigenvalues: [f64; 3],
    eigenvectors: [[f64; 3]; 3],
    center: [f64; 3],
    stacks: u32,
    slices: u32,
) -> Vec<[f64; 3]> {
    let stk = stacks.max(2) as usize;
    let slc = slices.max(3) as usize;
    let pi = std::f64::consts::PI;
    let two_pi = std::f64::consts::TAU;

    let semi = [
        eigenvalues[0].abs(),
        eigenvalues[1].abs(),
        eigenvalues[2].abs(),
    ];
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity((stk + 1) * (slc + 1));

    for i in 0..=stk {
        let phi = (i as f64 / stk as f64) * pi; // [0, π]
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        for j in 0..=slc {
            let theta = (j as f64 / slc as f64) * two_pi;
            let sin_th = theta.sin();
            let cos_th = theta.cos();

            // Sphere point scaled by semi-axes
            let local = [
                semi[0] * sin_phi * cos_th,
                semi[1] * sin_phi * sin_th,
                semi[2] * cos_phi,
            ];

            // Rotate into eigenvector frame (e_i are rows → columns of rotation matrix)
            let world = [
                center[0]
                    + eigenvectors[0][0] * local[0]
                    + eigenvectors[1][0] * local[1]
                    + eigenvectors[2][0] * local[2],
                center[1]
                    + eigenvectors[0][1] * local[0]
                    + eigenvectors[1][1] * local[1]
                    + eigenvectors[2][1] * local[2],
                center[2]
                    + eigenvectors[0][2] * local[0]
                    + eigenvectors[1][2] * local[1]
                    + eigenvectors[2][2] * local[2],
            ];
            pts.push(world);
        }
    }
    pts
}

// ──────────────────────────────────────────────────────────────────────────────
// Hedgehog field
// ──────────────────────────────────────────────────────────────────────────────

/// Place arrow glyphs on a regular 3-D grid sampling a vector field.
///
/// # Parameters
/// - `origin` — grid origin `[x0, y0, z0]`.
/// - `spacing` — grid spacing `[dx, dy, dz]`.
/// - `dims` — grid dimensions `[nx, ny, nz]` (number of nodes per axis).
/// - `field` — closure mapping a grid position to a direction vector.
/// - `scale` — uniform scale applied to each glyph.
/// - `color` — RGBA colour for all glyphs.
///
/// Returns one [`GlyphInstance`] per grid node.
pub fn hedgehog_field<F>(
    origin: [f64; 3],
    spacing: [f64; 3],
    dims: [usize; 3],
    field: F,
    scale: f64,
    color: [f32; 4],
) -> Vec<GlyphInstance>
where
    F: Fn([f64; 3]) -> [f64; 3],
{
    let [nx, ny, nz] = dims;
    let mut glyphs = Vec::with_capacity(nx * ny * nz);

    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let pos = [
                    origin[0] + ix as f64 * spacing[0],
                    origin[1] + iy as f64 * spacing[1],
                    origin[2] + iz as f64 * spacing[2],
                ];
                let dir_raw = field(pos);
                let dir = normalize3(dir_raw);
                glyphs.push(GlyphInstance {
                    position: pos,
                    direction: dir,
                    scale,
                    color,
                    glyph_type: GlyphType::Arrow,
                });
            }
        }
    }
    glyphs
}

// ──────────────────────────────────────────────────────────────────────────────
// Streamribbon
// ──────────────────────────────────────────────────────────────────────────────

/// A ribbon quad strip following a streamline, with optional twist.
#[derive(Debug, Clone)]
pub struct StreamRibbon {
    /// Left-edge vertices of the ribbon.
    pub left: Vec<[f64; 3]>,
    /// Right-edge vertices of the ribbon.
    pub right: Vec<[f64; 3]>,
}

/// Build a ribbon along a streamline by sweeping a perpendicular segment.
///
/// At each point along `path`, the ribbon half-width is `width / 2` and the
/// perpendicular is computed from the local tangent and an up vector.  An
/// optional `twist_per_step` (radians) rotates the perpendicular each step.
///
/// # Parameters
/// - `path` — ordered list of 3-D positions along the streamline.
/// - `width` — ribbon width.
/// - `up` — reference up vector for initial frame orientation.
/// - `twist_per_step` — extra twist per step (radians, 0 for no twist).
///
/// Returns a [`StreamRibbon`] with the same number of points as `path`.
pub fn streamribbon(
    path: &[[f64; 3]],
    width: f64,
    up: [f64; 3],
    twist_per_step: f64,
) -> StreamRibbon {
    if path.is_empty() {
        return StreamRibbon {
            left: vec![],
            right: vec![],
        };
    }

    let hw = width * 0.5;
    let mut left = Vec::with_capacity(path.len());
    let mut right = Vec::with_capacity(path.len());
    let mut angle = 0.0_f64;

    for i in 0..path.len() {
        // Tangent: forward difference, or backward at last point
        let tangent = if i + 1 < path.len() {
            normalize3(sub3(path[i + 1], path[i]))
        } else if i > 0 {
            normalize3(sub3(path[i], path[i - 1]))
        } else {
            [0.0, 0.0, 1.0]
        };

        // Perpendicular via cross product with up
        let perp_raw = cross3(tangent, up);
        let perp = if dot3(perp_raw, perp_raw) < 1e-20 {
            // Tangent is parallel to up — use a fallback
            let alt = [up[1], up[2], up[0]]; // cyclic permutation
            normalize3(cross3(tangent, alt))
        } else {
            normalize3(perp_raw)
        };

        // Apply accumulated twist
        angle += twist_per_step;
        let (sa, ca) = angle.sin_cos();
        let twisted = [
            ca * perp[0] - sa * tangent[0],
            ca * perp[1] - sa * tangent[1],
            ca * perp[2] - sa * tangent[2],
        ];
        let t_norm = normalize3(twisted);

        left.push([
            path[i][0] - hw * t_norm[0],
            path[i][1] - hw * t_norm[1],
            path[i][2] - hw * t_norm[2],
        ]);
        right.push([
            path[i][0] + hw * t_norm[0],
            path[i][1] + hw * t_norm[1],
            path[i][2] + hw * t_norm[2],
        ]);
    }

    StreamRibbon { left, right }
}

// ──────────────────────────────────────────────────────────────────────────────
// Oriented discs
// ──────────────────────────────────────────────────────────────────────────────

/// A single oriented disc (a flat circle placed at a contact/surface normal).
#[derive(Debug, Clone)]
pub struct OrientedDisc {
    /// World-space center.
    pub center: [f64; 3],
    /// Unit normal (disc is perpendicular to this).
    pub normal: [f64; 3],
    /// Disc radius.
    pub radius: f64,
    /// RGBA color.
    pub color: [f32; 4],
}

/// Build a list of oriented discs aligned with surface normals.
///
/// Each entry in `contacts` is `(position, normal, radius)`.  The output
/// discs are useful for visualising contact points or surface orientations.
pub fn oriented_discs(
    contacts: &[([f64; 3], [f64; 3], f64)],
    color: [f32; 4],
) -> Vec<OrientedDisc> {
    contacts
        .iter()
        .map(|(pos, nor, r)| OrientedDisc {
            center: *pos,
            normal: normalize3(*nor),
            radius: *r,
            color,
        })
        .collect()
}

// ──────────────────────────────────────────────────────────────────────────────
// Math helpers (plain f64 — no nalgebra in viz)
// ──────────────────────────────────────────────────────────────────────────────

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-20 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    fn white() -> [f32; 4] {
        [1.0, 1.0, 1.0, 1.0]
    }

    // ── GlyphInstance ────────────────────────────────────────────────────────

    #[test]
    fn test_glyph_instance_new() {
        let g = GlyphInstance::new(
            [1.0, 2.0, 3.0],
            [0.0, 1.0, 0.0],
            2.5,
            white(),
            GlyphType::Arrow,
        );
        assert!((g.scale - 2.5).abs() < EPS);
        assert_eq!(g.glyph_type, GlyphType::Arrow);
    }

    #[test]
    fn test_glyph_type_variants() {
        let types = [
            GlyphType::Arrow,
            GlyphType::Cone,
            GlyphType::Sphere,
            GlyphType::Cube,
            GlyphType::Line,
            GlyphType::Disc,
        ];
        assert_eq!(types.len(), 6);
    }

    // ── arrow_glyph ──────────────────────────────────────────────────────────

    #[test]
    fn test_arrow_glyph_has_vertices() {
        let mesh = arrow_glyph(0.05, 0.1, 0.2, 8);
        assert!(!mesh.vertices.is_empty());
        assert!(!mesh.indices.is_empty());
    }

    #[test]
    fn test_arrow_glyph_tip_at_unit_z() {
        let mesh = arrow_glyph(0.05, 0.1, 0.3, 6);
        // The tip vertex has z = 1.0
        let tip = mesh.vertices.iter().find(|v| (v.pos[2] - 1.0).abs() < 1e-9);
        assert!(tip.is_some(), "arrow tip at z=1 must exist");
    }

    #[test]
    fn test_arrow_glyph_base_at_zero() {
        let mesh = arrow_glyph(0.05, 0.1, 0.2, 8);
        // Cap center vertex should be at z=0
        let cap = mesh
            .vertices
            .iter()
            .find(|v| v.pos[0].abs() < 1e-12 && v.pos[1].abs() < 1e-12 && v.pos[2].abs() < 1e-12);
        assert!(cap.is_some(), "cap center at origin must exist");
    }

    #[test]
    fn test_arrow_glyph_min_segments() {
        // segments=1 should be clamped to 3
        let mesh = arrow_glyph(0.05, 0.1, 0.2, 1);
        assert!(!mesh.vertices.is_empty());
    }

    #[test]
    fn test_arrow_glyph_all_indices_in_bounds() {
        let mesh = arrow_glyph(0.05, 0.1, 0.25, 10);
        let n = mesh.vertices.len() as u32;
        for tri in &mesh.indices {
            for &idx in tri.iter() {
                assert!(idx < n, "index {} out of bounds (n={})", idx, n);
            }
        }
    }

    #[test]
    fn test_arrow_glyph_more_segments_more_verts() {
        let m8 = arrow_glyph(0.05, 0.1, 0.2, 8);
        let m16 = arrow_glyph(0.05, 0.1, 0.2, 16);
        assert!(m16.vertices.len() > m8.vertices.len());
    }

    // ── tensor_ellipsoid ─────────────────────────────────────────────────────

    #[test]
    fn test_ellipsoid_vertex_count() {
        let evals = [1.0, 2.0, 3.0];
        let evecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pts = tensor_ellipsoid(evals, evecs, [0.0, 0.0, 0.0], 4, 6);
        assert_eq!(pts.len(), 5 * 7); // (stacks+1)*(slices+1)
    }

    #[test]
    fn test_ellipsoid_identity_frame_unit_sphere() {
        // With all eigenvalues = 1 and identity eigenvectors, points on unit sphere
        let evals = [1.0, 1.0, 1.0];
        let evecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pts = tensor_ellipsoid(evals, evecs, [0.0, 0.0, 0.0], 8, 16);
        for p in &pts {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 1.0).abs() < 1e-10, "point not on unit sphere: r={r}");
        }
    }

    #[test]
    fn test_ellipsoid_negative_eigenvalues() {
        // Negative eigenvalues → semi-axes = |λ|
        let evals = [-1.0, -2.0, -3.0];
        let evecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let pts = tensor_ellipsoid(evals, evecs, [0.0, 0.0, 0.0], 4, 6);
        assert!(!pts.is_empty());
    }

    #[test]
    fn test_ellipsoid_center_translation() {
        let evals = [1.0, 1.0, 1.0];
        let evecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let center = [5.0, 3.0, 1.0];
        let pts = tensor_ellipsoid(evals, evecs, center, 4, 4);
        // All points should be within distance ~1 from center
        for p in &pts {
            let d = ((p[0] - center[0]).powi(2)
                + (p[1] - center[1]).powi(2)
                + (p[2] - center[2]).powi(2))
            .sqrt();
            assert!(d < 1.0 + 1e-9, "point too far from center: d={d}");
        }
    }

    // ── hedgehog_field ───────────────────────────────────────────────────────

    #[test]
    fn test_hedgehog_count() {
        let glyphs = hedgehog_field(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [3, 3, 3],
            |_p| [1.0, 0.0, 0.0],
            1.0,
            white(),
        );
        assert_eq!(glyphs.len(), 27);
    }

    #[test]
    fn test_hedgehog_glyph_type_arrow() {
        let glyphs = hedgehog_field(
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [2, 2, 1],
            |_p| [0.0, 1.0, 0.0],
            1.0,
            white(),
        );
        for g in &glyphs {
            assert_eq!(g.glyph_type, GlyphType::Arrow);
        }
    }

    #[test]
    fn test_hedgehog_directions_normalized() {
        let glyphs = hedgehog_field(
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2, 2, 2],
            |p| [p[0] + 1.0, p[1], p[2]], // unnormalized
            1.0,
            white(),
        );
        for g in &glyphs {
            let d = g.direction;
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 1e-10,
                "direction not normalized: len={len}"
            );
        }
    }

    #[test]
    fn test_hedgehog_grid_positions() {
        let origin = [1.0, 2.0, 3.0];
        let spacing = [0.5, 0.5, 0.5];
        let dims = [2, 2, 2];
        let glyphs = hedgehog_field(origin, spacing, dims, |_| [1.0, 0.0, 0.0], 1.0, white());
        // First glyph should be at origin
        assert!((glyphs[0].position[0] - 1.0).abs() < EPS);
        assert!((glyphs[0].position[1] - 2.0).abs() < EPS);
        assert!((glyphs[0].position[2] - 3.0).abs() < EPS);
    }

    #[test]
    fn test_hedgehog_empty_dims() {
        let glyphs = hedgehog_field(
            [0.0; 3],
            [1.0; 3],
            [0, 1, 1],
            |_| [1.0, 0.0, 0.0],
            1.0,
            white(),
        );
        assert!(glyphs.is_empty());
    }

    // ── streamribbon ─────────────────────────────────────────────────────────

    #[test]
    fn test_streamribbon_empty_path() {
        let ribbon = streamribbon(&[], 0.1, [0.0, 1.0, 0.0], 0.0);
        assert!(ribbon.left.is_empty());
        assert!(ribbon.right.is_empty());
    }

    #[test]
    fn test_streamribbon_same_count_as_path() {
        let path: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let ribbon = streamribbon(&path, 0.2, [0.0, 1.0, 0.0], 0.0);
        assert_eq!(ribbon.left.len(), 10);
        assert_eq!(ribbon.right.len(), 10);
    }

    #[test]
    fn test_streamribbon_width_correct() {
        let path = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let width = 0.4;
        let ribbon = streamribbon(&path, width, [0.0, 1.0, 0.0], 0.0);
        // Each left/right pair should be `width` apart
        for i in 0..ribbon.left.len() {
            let l = ribbon.left[i];
            let r = ribbon.right[i];
            let d = ((l[0] - r[0]).powi(2) + (l[1] - r[1]).powi(2) + (l[2] - r[2]).powi(2)).sqrt();
            assert!(
                (d - width).abs() < 1e-9,
                "ribbon width at {i}: {d:.9} ≠ {width}"
            );
        }
    }

    #[test]
    fn test_streamribbon_single_point() {
        let path = vec![[1.0_f64, 2.0, 3.0]];
        let ribbon = streamribbon(&path, 0.5, [0.0, 1.0, 0.0], 0.0);
        assert_eq!(ribbon.left.len(), 1);
    }

    // ── oriented_discs ───────────────────────────────────────────────────────

    #[test]
    fn test_oriented_discs_count() {
        let contacts = vec![
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5),
            ([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.3),
        ];
        let discs = oriented_discs(&contacts, white());
        assert_eq!(discs.len(), 2);
    }

    #[test]
    fn test_oriented_discs_normals_normalized() {
        let contacts = vec![([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 1.0)];
        let discs = oriented_discs(&contacts, white());
        let n = discs[0].normal;
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_oriented_discs_empty() {
        let discs = oriented_discs(&[], white());
        assert!(discs.is_empty());
    }

    #[test]
    fn test_oriented_discs_radius_preserved() {
        let contacts = vec![([0.0; 3], [0.0, 1.0, 0.0], 2.5)];
        let discs = oriented_discs(&contacts, white());
        assert!((discs[0].radius - 2.5).abs() < EPS);
    }

    // ── math helpers ─────────────────────────────────────────────────────────

    #[test]
    fn test_normalize3_unit_vector() {
        let v = normalize3([3.0, 4.0, 0.0]);
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_zero_returns_fallback() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert!((v[2] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_cross3_orthogonal() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < EPS);
        assert!(c[0].abs() < EPS && c[1].abs() < EPS);
    }

    #[test]
    fn test_dot3_orthogonal() {
        let d = dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(d.abs() < EPS);
    }

    #[test]
    fn test_sub3_basic() {
        let v = sub3([3.0, 5.0, 7.0], [1.0, 2.0, 3.0]);
        assert!((v[0] - 2.0).abs() < EPS && (v[1] - 3.0).abs() < EPS && (v[2] - 4.0).abs() < EPS);
    }
}
