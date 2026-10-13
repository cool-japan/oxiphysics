// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! 2D cross-section slice renderer for 3D fields.
//!
//! Provides an oriented cutting plane (`SlicePlane`), a 2D grid on that plane
//! (`SliceGrid`), scalar and vector field sampling on the grid, an RGBA image
//! buffer (`SliceImage`), colormap rendering, vector-arrow generation, and
//! marching-squares 2D contour extraction.

// ─────────────────────────────────────────────────────────────────────────────
// SlicePlane
// ─────────────────────────────────────────────────────────────────────────────

/// An oriented cutting plane in 3D space, defined by an origin and two
/// orthonormal tangent vectors (`u_axis` and `v_axis`).
#[derive(Debug, Clone, PartialEq)]
pub struct SlicePlane {
    /// Unit normal of the plane.
    pub normal: [f64; 3],
    /// A point on the plane (the 2D origin maps here).
    pub origin: [f64; 3],
    /// First tangent axis (u-direction on the slice).
    pub u_axis: [f64; 3],
    /// Second tangent axis (v-direction on the slice).
    pub v_axis: [f64; 3],
}

impl SlicePlane {
    /// Construct a horizontal XY-plane at height `z`.
    pub fn xy_plane(z: f64) -> Self {
        Self {
            normal: [0.0, 0.0, 1.0],
            origin: [0.0, 0.0, z],
            u_axis: [1.0, 0.0, 0.0],
            v_axis: [0.0, 1.0, 0.0],
        }
    }

    /// Construct an XZ-plane at `y`.
    pub fn xz_plane(y: f64) -> Self {
        Self {
            normal: [0.0, 1.0, 0.0],
            origin: [0.0, y, 0.0],
            u_axis: [1.0, 0.0, 0.0],
            v_axis: [0.0, 0.0, 1.0],
        }
    }

    /// Construct a YZ-plane at `x`.
    pub fn yz_plane(x: f64) -> Self {
        Self {
            normal: [1.0, 0.0, 0.0],
            origin: [x, 0.0, 0.0],
            u_axis: [0.0, 1.0, 0.0],
            v_axis: [0.0, 0.0, 1.0],
        }
    }

    /// Construct an arbitrary plane from a unit `normal` and an `origin`.
    ///
    /// The `u_axis` and `v_axis` are computed via Gram–Schmidt.
    pub fn custom(normal: [f64; 3], origin: [f64; 3]) -> Self {
        // Build an orthonormal basis from the normal
        let n = normalize3(normal);
        // Choose a reference vector not parallel to n
        let ref_vec = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let u = normalize3(cross3(ref_vec, n));
        let v = cross3(n, u);
        Self {
            normal: n,
            origin,
            u_axis: u,
            v_axis: v,
        }
    }

    /// Project a world-space point `p` onto this plane, returning 2D (u, v)
    /// coordinates relative to `origin`.
    pub fn project_point(&self, p: [f64; 3]) -> [f64; 2] {
        let d = sub3(p, self.origin);
        [dot3(d, self.u_axis), dot3(d, self.v_axis)]
    }

    /// Convert 2D slice coordinates `uv` back to a world-space point.
    pub fn world_point(&self, uv: [f64; 2]) -> [f64; 3] {
        [
            self.origin[0] + uv[0] * self.u_axis[0] + uv[1] * self.v_axis[0],
            self.origin[1] + uv[0] * self.u_axis[1] + uv[1] * self.v_axis[1],
            self.origin[2] + uv[0] * self.u_axis[2] + uv[1] * self.v_axis[2],
        ]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SliceGrid
// ─────────────────────────────────────────────────────────────────────────────

/// A regular 2D grid lying on a `SlicePlane`.
///
/// The grid has `nu × nv` cells with spacing `du × dv`.
#[derive(Debug, Clone)]
pub struct SliceGrid {
    /// The cutting plane on which this grid lies.
    pub plane: SlicePlane,
    /// Number of grid points along the u-axis.
    pub nu: usize,
    /// Number of grid points along the v-axis.
    pub nv: usize,
    /// Spacing between grid points in the u-direction.
    pub du: f64,
    /// Spacing between grid points in the v-direction.
    pub dv: f64,
}

impl SliceGrid {
    /// Create a new `SliceGrid`.
    pub fn new(plane: SlicePlane, nu: usize, nv: usize, du: f64, dv: f64) -> Self {
        Self {
            plane,
            nu,
            nv,
            du,
            dv,
        }
    }

    /// Total number of grid cells (`nu * nv`).
    pub fn cell_count(&self) -> usize {
        self.nu * self.nv
    }

    /// World-space position of the grid point at integer indices `(u_idx, v_idx)`.
    pub fn uv_to_world(&self, u_idx: usize, v_idx: usize) -> [f64; 3] {
        let u = u_idx as f64 * self.du;
        let v = v_idx as f64 * self.dv;
        self.plane.world_point([u, v])
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Field sampling
// ─────────────────────────────────────────────────────────────────────────────

/// Sample a 3D scalar field on all grid points of a `SliceGrid`.
///
/// Returns a flat `Vec`f64` of length `nu * nv`, stored in row-major order
/// (u varies fastest).
pub fn sample_scalar_field_on_slice(grid: &SliceGrid, field: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut out = Vec::with_capacity(grid.nu * grid.nv);
    for v_idx in 0..grid.nv {
        for u_idx in 0..grid.nu {
            let p = grid.uv_to_world(u_idx, v_idx);
            out.push(field(p));
        }
    }
    out
}

/// Sample a 3D vector field on all grid points of a `SliceGrid`, projecting
/// each 3D vector onto the slice's 2D tangent plane.
///
/// Returns a flat `Vec<\[f64; 2\]>` of length `nu * nv`.
pub fn sample_vector_field_on_slice(
    grid: &SliceGrid,
    field: impl Fn([f64; 3]) -> [f64; 3],
) -> Vec<[f64; 2]> {
    let mut out = Vec::with_capacity(grid.nu * grid.nv);
    for v_idx in 0..grid.nv {
        for u_idx in 0..grid.nu {
            let p = grid.uv_to_world(u_idx, v_idx);
            let v3 = field(p);
            let u_comp = dot3(v3, grid.plane.u_axis);
            let v_comp = dot3(v3, grid.plane.v_axis);
            out.push([u_comp, v_comp]);
        }
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// SliceImage
// ─────────────────────────────────────────────────────────────────────────────

/// An RGBA image buffer for storing a colorized slice.
#[derive(Debug, Clone)]
pub struct SliceImage {
    /// Pixel width.
    pub width: usize,
    /// Pixel height.
    pub height: usize,
    /// Pixel data in row-major order; each element is `\[R, G, B, A\]` in 0–255.
    pub pixels: Vec<[u8; 4]>,
}

impl SliceImage {
    /// Create a new `SliceImage` filled with opaque black pixels.
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![[0, 0, 0, 255]; width * height],
        }
    }

    /// Set the pixel at `(x, y)` to `rgba`.  `x` is the column, `y` the row.
    pub fn set_pixel(&mut self, x: usize, y: usize, rgba: [u8; 4]) {
        self.pixels[y * self.width + x] = rgba;
    }

    /// Get the pixel at `(x, y)`.
    pub fn get_pixel(&self, x: usize, y: usize) -> [u8; 4] {
        self.pixels[y * self.width + x]
    }

    /// Total number of pixels (`width * height`).
    pub fn pixel_count(&self) -> usize {
        self.width * self.height
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Colormap rendering
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a scalar field to a colormap image.
///
/// `values` holds the `nu * nv` samples of a slice in row-major order (u
/// varies fastest), as produced by [`sample_scalar_field_on_slice`] for a
/// [`SliceGrid`] with `nu` × `nv` points.  The result is a [`SliceImage`] with
/// `width = nu` and `height = nv` whose pixel `(x, y)` shows `values[y * nu + x]`.
/// Values are normalized to `[0, 1]` by the minimum and maximum of `values`
/// (a constant field maps to `0.5`) before being passed to `colormap_fn`.
///
/// # Errors
///
/// Returns [`Error::InvalidParameter`](crate::error::Error::InvalidParameter)
/// when `values.len()` differs from `nu * nv` (or `nu * nv` overflows).
pub fn scalar_to_colormap_slice(
    values: &[f64],
    nu: usize,
    nv: usize,
    colormap_fn: impl Fn(f64) -> [f32; 4],
) -> crate::error::Result<SliceImage> {
    let expected = nu.checked_mul(nv).ok_or_else(|| {
        crate::error::Error::invalid_parameter(
            "nu, nv",
            format!("grid size {nu} x {nv} overflows usize"),
        )
    })?;
    if values.len() != expected {
        return Err(crate::error::Error::invalid_parameter(
            "values",
            format!(
                "expected nu * nv = {nu} * {nv} = {expected} samples, got {}",
                values.len()
            ),
        ));
    }

    let min_v = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_v = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = max_v - min_v;

    // Pixel (x, y) lives at y * width + x, the same row-major, u-fastest order
    // as `values`, so the samples map onto the pixels one to one.
    let mut img = SliceImage::new(nu, nv);
    for (pixel, &v) in img.pixels.iter_mut().zip(values.iter()) {
        let t = if range < 1e-30 {
            0.5
        } else {
            ((v - min_v) / range).clamp(0.0, 1.0)
        };
        let c = colormap_fn(t);
        *pixel = [
            (c[0].clamp(0.0, 1.0) * 255.0) as u8,
            (c[1].clamp(0.0, 1.0) * 255.0) as u8,
            (c[2].clamp(0.0, 1.0) * 255.0) as u8,
            (c[3].clamp(0.0, 1.0) * 255.0) as u8,
        ];
    }
    Ok(img)
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector arrows
// ─────────────────────────────────────────────────────────────────────────────

/// Generate arrow line segments for a 2D vector field.
///
/// For each vector `vectors\[i\]` at position `positions\[i\]`, returns two
/// segments: the main shaft and the arrowhead.  `scale` controls the arrow
/// length.
pub fn vector_arrows_2d(
    vectors: &[[f64; 2]],
    positions: &[[f64; 2]],
    scale: f64,
) -> Vec<[[f64; 2]; 2]> {
    assert_eq!(vectors.len(), positions.len());
    let mut segs = Vec::with_capacity(vectors.len() * 3);
    for (&v, &p) in vectors.iter().zip(positions.iter()) {
        let tip = [p[0] + v[0] * scale, p[1] + v[1] * scale];
        // Main shaft
        segs.push([p, tip]);
        // Arrowhead: two short lines at ~30° from tip
        let len = (v[0] * v[0] + v[1] * v[1]).sqrt();
        if len > 1e-30 {
            let head_len = scale * 0.2;
            let ux = -v[0] / len;
            let uy = -v[1] / len;
            let cos30 = 0.866_f64;
            let sin30 = 0.5_f64;
            let w1 = [
                tip[0] + head_len * (ux * cos30 - uy * sin30),
                tip[1] + head_len * (ux * sin30 + uy * cos30),
            ];
            let w2 = [
                tip[0] + head_len * (ux * cos30 + uy * sin30),
                tip[1] + head_len * (-ux * sin30 + uy * cos30),
            ];
            segs.push([tip, w1]);
            segs.push([tip, w2]);
        }
    }
    segs
}

// ─────────────────────────────────────────────────────────────────────────────
// 2D contour lines (marching squares)
// ─────────────────────────────────────────────────────────────────────────────

/// Extract iso-contour line segments from a 2D scalar field using marching
/// squares.
///
/// `values` has length `nx * ny` stored in row-major order (x varies fastest).
/// Returns a list of `\[\[x0,y0\\],\[x1,y1\]]` segments.
pub fn contour_lines_2d(values: &[f64], nx: usize, ny: usize, level: f64) -> Vec<[[f64; 2]; 2]> {
    if nx < 2 || ny < 2 {
        return vec![];
    }
    let mut segs = Vec::new();
    let idx = |xi: usize, yi: usize| yi * nx + xi;

    for xi in 0..(nx - 1) {
        for yi in 0..(ny - 1) {
            let v00 = values[idx(xi, yi)];
            let v10 = values[idx(xi + 1, yi)];
            let v01 = values[idx(xi, yi + 1)];
            let v11 = values[idx(xi + 1, yi + 1)];

            let corners = [
                [xi as f64, yi as f64],
                [(xi + 1) as f64, yi as f64],
                [(xi + 1) as f64, (yi + 1) as f64],
                [xi as f64, (yi + 1) as f64],
            ];
            let vals = [v00, v10, v11, v01];
            let above = [v00 >= level, v10 >= level, v11 >= level, v01 >= level];

            let edge_pt = |a: usize, b: usize| -> [f64; 2] {
                let t = if (vals[b] - vals[a]).abs() < 1e-30 {
                    0.5
                } else {
                    (level - vals[a]) / (vals[b] - vals[a])
                };
                [
                    corners[a][0] + t * (corners[b][0] - corners[a][0]),
                    corners[a][1] + t * (corners[b][1] - corners[a][1]),
                ]
            };

            let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];
            let mut pts = Vec::new();
            for &(ea, eb) in &edges {
                if above[ea] != above[eb] {
                    pts.push(edge_pt(ea, eb));
                }
            }
            let mut k = 0;
            while k + 1 < pts.len() {
                segs.push([pts[k], pts[k + 1]]);
                k += 2;
            }
        }
    }
    segs
}

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers
// ─────────────────────────────────────────────────────────────────────────────

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
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-30 {
        v
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── SlicePlane ────────────────────────────────────────────────────────────

    #[test]
    fn test_xy_plane_normal() {
        let p = SlicePlane::xy_plane(3.0);
        assert!((p.normal[2] - 1.0).abs() < 1e-12);
        assert!((p.origin[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_xz_plane_normal() {
        let p = SlicePlane::xz_plane(2.0);
        assert!((p.normal[1] - 1.0).abs() < 1e-12);
        assert!((p.origin[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_yz_plane_normal() {
        let p = SlicePlane::yz_plane(1.0);
        assert!((p.normal[0] - 1.0).abs() < 1e-12);
        assert!((p.origin[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_project_world_roundtrip_xy() {
        let plane = SlicePlane::xy_plane(5.0);
        let world = [2.0, 3.0, 5.0];
        let uv = plane.project_point(world);
        let back = plane.world_point(uv);
        for i in 0..3 {
            assert!(
                (back[i] - world[i]).abs() < 1e-10,
                "component {i}: {} vs {}",
                back[i],
                world[i]
            );
        }
    }

    #[test]
    fn test_project_world_roundtrip_xz() {
        let plane = SlicePlane::xz_plane(0.0);
        let world = [1.0, 0.0, 4.0];
        let uv = plane.project_point(world);
        let back = plane.world_point(uv);
        for i in 0..3 {
            assert!((back[i] - world[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_project_world_roundtrip_yz() {
        let plane = SlicePlane::yz_plane(-1.0);
        let world = [-1.0, 2.0, 3.0];
        let uv = plane.project_point(world);
        let back = plane.world_point(uv);
        for i in 0..3 {
            assert!((back[i] - world[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_custom_plane_orthonormal() {
        let plane = SlicePlane::custom([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let u = plane.u_axis;
        let v = plane.v_axis;
        let len_u = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let len_v = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((len_u - 1.0).abs() < 1e-10);
        assert!((len_v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_project_origin_gives_zero_uv() {
        let plane = SlicePlane::xy_plane(2.0);
        let uv = plane.project_point(plane.origin);
        assert!(uv[0].abs() < 1e-12 && uv[1].abs() < 1e-12);
    }

    #[test]
    fn test_xy_plane_project_x_axis() {
        let plane = SlicePlane::xy_plane(0.0);
        let uv = plane.project_point([3.0, 0.0, 0.0]);
        assert!((uv[0] - 3.0).abs() < 1e-12);
        assert!(uv[1].abs() < 1e-12);
    }

    // ── SliceGrid ─────────────────────────────────────────────────────────────

    #[test]
    fn test_cell_count() {
        let g = SliceGrid::new(SlicePlane::xy_plane(0.0), 5, 4, 1.0, 1.0);
        assert_eq!(g.cell_count(), 20);
    }

    #[test]
    fn test_uv_to_world_origin() {
        let plane = SlicePlane::xy_plane(2.0);
        let g = SliceGrid::new(plane.clone(), 4, 4, 1.0, 1.0);
        let w = g.uv_to_world(0, 0);
        for (&wi, &oi) in w.iter().zip(plane.origin.iter()) {
            assert!((wi - oi).abs() < 1e-10);
        }
    }

    #[test]
    fn test_uv_to_world_step() {
        let plane = SlicePlane::xy_plane(0.0);
        let g = SliceGrid::new(plane, 4, 4, 2.0, 3.0);
        let w = g.uv_to_world(1, 0);
        assert!((w[0] - 2.0).abs() < 1e-10);
        let w2 = g.uv_to_world(0, 1);
        assert!((w2[1] - 3.0).abs() < 1e-10);
    }

    // ── sample_scalar_field_on_slice ──────────────────────────────────────────

    #[test]
    fn test_sample_scalar_count() {
        let g = SliceGrid::new(SlicePlane::xy_plane(0.0), 3, 4, 1.0, 1.0);
        let vals = sample_scalar_field_on_slice(&g, |_p| 1.0);
        assert_eq!(vals.len(), 12);
    }

    #[test]
    fn test_sample_scalar_constant_field() {
        let g = SliceGrid::new(SlicePlane::xy_plane(0.0), 3, 3, 1.0, 1.0);
        let vals = sample_scalar_field_on_slice(&g, |_p| 42.0);
        for v in &vals {
            assert!((*v - 42.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_sample_scalar_x_coordinate() {
        let plane = SlicePlane::xy_plane(0.0);
        let g = SliceGrid::new(plane, 4, 1, 1.0, 1.0);
        let vals = sample_scalar_field_on_slice(&g, |p| p[0]);
        // x-coordinates should be 0, 1, 2, 3
        for (k, v) in vals.iter().enumerate() {
            assert!((v - k as f64).abs() < 1e-10);
        }
    }

    // ── sample_vector_field_on_slice ──────────────────────────────────────────

    #[test]
    fn test_sample_vector_count() {
        let g = SliceGrid::new(SlicePlane::xy_plane(0.0), 3, 4, 1.0, 1.0);
        let vecs = sample_vector_field_on_slice(&g, |_p| [1.0, 0.0, 0.0]);
        assert_eq!(vecs.len(), 12);
    }

    #[test]
    fn test_sample_vector_projects_correctly() {
        let plane = SlicePlane::xy_plane(0.0);
        // Field pointing in x: should project to u=1, v=0
        let g = SliceGrid::new(plane, 2, 2, 1.0, 1.0);
        let vecs = sample_vector_field_on_slice(&g, |_p| [1.0, 0.0, 0.0]);
        for v in &vecs {
            assert!((v[0] - 1.0).abs() < 1e-10);
            assert!(v[1].abs() < 1e-10);
        }
    }

    // ── SliceImage ────────────────────────────────────────────────────────────

    #[test]
    fn test_slice_image_pixel_count() {
        let img = SliceImage::new(5, 4);
        assert_eq!(img.pixel_count(), 20);
    }

    #[test]
    fn test_slice_image_default_pixel() {
        let img = SliceImage::new(3, 3);
        let p = img.get_pixel(0, 0);
        assert_eq!(p[3], 255);
    }

    #[test]
    fn test_slice_image_set_get_roundtrip() {
        let mut img = SliceImage::new(4, 4);
        img.set_pixel(2, 3, [10, 20, 30, 255]);
        assert_eq!(img.get_pixel(2, 3), [10, 20, 30, 255]);
    }

    #[test]
    fn test_slice_image_multiple_pixels() {
        let mut img = SliceImage::new(4, 4);
        img.set_pixel(0, 0, [255, 0, 0, 255]);
        img.set_pixel(1, 1, [0, 255, 0, 255]);
        assert_eq!(img.get_pixel(0, 0), [255, 0, 0, 255]);
        assert_eq!(img.get_pixel(1, 1), [0, 255, 0, 255]);
    }

    // ── scalar_to_colormap_slice ──────────────────────────────────────────────

    #[test]
    fn test_scalar_to_colormap_empty() {
        let img = scalar_to_colormap_slice(&[], 0, 0, |_t| [0.5, 0.5, 0.5, 1.0]).unwrap();
        assert_eq!(img.pixel_count(), 0);
    }

    #[test]
    fn test_scalar_to_colormap_pixel_count() {
        let vals: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let img = scalar_to_colormap_slice(&vals, 4, 4, |t| [t as f32, 0.0, 0.0, 1.0]).unwrap();
        assert_eq!((img.width, img.height), (4, 4));
        assert_eq!(img.pixel_count(), vals.len());
    }

    #[test]
    fn test_scalar_to_colormap_alpha_preserved() {
        let vals = [0.0, 0.5, 1.0];
        let img = scalar_to_colormap_slice(&vals, 3, 1, |_t| [0.5, 0.5, 0.5, 0.8]).unwrap();
        // All pixels should have A ≈ 204 (0.8 * 255)
        assert_eq!(img.pixels.len(), vals.len());
        for p in &img.pixels {
            assert_eq!(p[3], 204);
        }
    }

    #[test]
    fn test_scalar_to_colormap_uses_the_requested_layout() {
        // A 4 x 2 ramp 0..8, u fastest: the image must be 4 wide and 2 tall,
        // not the sqrt-shaped square the old code reshaped everything into.
        let vals: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let img = scalar_to_colormap_slice(&vals, 4, 2, |t| [t as f32, 0.0, 0.0, 1.0]).unwrap();
        assert_eq!((img.width, img.height), (4, 2));
        for y in 0..2 {
            for x in 0..4 {
                let k = y * 4 + x;
                let expected = ((k as f64 / 7.0) as f32 * 255.0) as u8;
                assert_eq!(img.get_pixel(x, y)[0], expected, "pixel ({x}, {y})");
            }
        }
        // Brightest pixel of row 0 is its last column.
        let row0: Vec<u8> = (0..4).map(|x| img.get_pixel(x, 0)[0]).collect();
        assert_eq!(row0.iter().max(), Some(&row0[3]));
        assert!(
            row0.windows(2).all(|w| w[0] < w[1]),
            "row 0 not increasing: {row0:?}"
        );
        // Row 1 continues the ramp.
        assert_eq!(img.get_pixel(3, 1)[0], 255);
    }

    #[test]
    fn test_scalar_to_colormap_non_square_slice_round_trips_sampling() {
        // 100 x 10 slice sampled from a field that depends on u only: every
        // image row must be the same horizontal ramp.
        let grid = SliceGrid::new(SlicePlane::xy_plane(0.0), 100, 10, 0.5, 0.5);
        let vals = sample_scalar_field_on_slice(&grid, |p| p[0]);
        let img = scalar_to_colormap_slice(&vals, grid.nu, grid.nv, |t| [t as f32, 0.0, 0.0, 1.0])
            .unwrap();
        assert_eq!((img.width, img.height), (100, 10));
        for y in 0..10 {
            assert_eq!(img.get_pixel(0, y)[0], 0);
            assert_eq!(img.get_pixel(99, y)[0], 255);
            assert_eq!(img.get_pixel(50, y)[0], img.get_pixel(50, 0)[0]);
        }
        assert_eq!(img.get_pixel(50, 3)[0], 128);
    }

    #[test]
    fn test_scalar_to_colormap_constant_field_maps_to_midpoint() {
        let vals = [2.0; 6];
        let img =
            scalar_to_colormap_slice(&vals, 3, 2, |t| [t as f32, t as f32, t as f32, 1.0]).unwrap();
        assert!(img.pixels.iter().all(|p| *p == [127, 127, 127, 255]));
    }

    #[test]
    fn test_scalar_to_colormap_rejects_mismatched_dimensions() {
        use crate::error::Error;
        let cmap = |_t: f64| [0.0_f32, 0.0, 0.0, 1.0];
        let err = scalar_to_colormap_slice(&[1.0; 5], 2, 2, cmap).unwrap_err();
        assert!(matches!(err, Error::InvalidParameter { .. }), "{err:?}");
        assert!(err.to_string().contains("expected nu * nv"), "{err}");
        // Too few samples, samples for a zero-sized grid, and overflow.
        assert!(scalar_to_colormap_slice(&[1.0; 3], 2, 2, cmap).is_err());
        assert!(scalar_to_colormap_slice(&[1.0; 3], 0, 3, cmap).is_err());
        assert!(matches!(
            scalar_to_colormap_slice(&[1.0], usize::MAX, 2, cmap),
            Err(Error::InvalidParameter { .. })
        ));
        // A zero-sized grid with no samples is a valid empty image.
        let empty = scalar_to_colormap_slice(&[], 0, 5, cmap).unwrap();
        assert_eq!((empty.width, empty.height, empty.pixel_count()), (0, 5, 0));
    }

    // ── contour_lines_2d ──────────────────────────────────────────────────────

    #[test]
    fn test_contour_lines_no_panic() {
        let vals: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let _segs = contour_lines_2d(&vals, 4, 4, 8.0);
    }

    #[test]
    fn test_contour_lines_degenerate() {
        // nx < 2: no segments
        let vals = [1.0, 2.0, 3.0];
        assert!(contour_lines_2d(&vals, 1, 3, 1.5).is_empty());
    }

    #[test]
    fn test_contour_lines_all_above() {
        // Level below all values: no crossings
        let vals = [5.0, 6.0, 7.0, 8.0];
        let segs = contour_lines_2d(&vals, 2, 2, 0.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_contour_lines_all_below() {
        // Level above all values: no crossings
        let vals = [1.0, 2.0, 3.0, 4.0];
        let segs = contour_lines_2d(&vals, 2, 2, 10.0);
        assert!(segs.is_empty());
    }

    // ── vector_arrows_2d ─────────────────────────────────────────────────────

    #[test]
    fn test_vector_arrows_empty() {
        let segs = vector_arrows_2d(&[], &[], 1.0);
        assert!(segs.is_empty());
    }

    #[test]
    fn test_vector_arrows_one() {
        let vecs = [[1.0, 0.0]];
        let pos = [[0.0, 0.0]];
        let segs = vector_arrows_2d(&vecs, &pos, 1.0);
        // 1 shaft + 2 arrowhead = 3 segments
        assert_eq!(segs.len(), 3);
    }

    #[test]
    fn test_vector_arrows_shaft_tip() {
        let vecs = [[2.0, 0.0]];
        let pos = [[1.0, 1.0]];
        let segs = vector_arrows_2d(&vecs, &pos, 1.0);
        // First segment tail should be at origin
        assert!((segs[0][0][0] - 1.0).abs() < 1e-10);
        // Tip x = 1 + 2*1 = 3
        assert!((segs[0][1][0] - 3.0).abs() < 1e-10);
    }
}
