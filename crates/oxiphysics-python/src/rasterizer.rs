// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Software rasterizer with Phong shading.
//!
//! Provides orthographic projection, scanline edge-function rasterization, and
//! per-vertex Phong lighting for the [`super::viz_api::PyMeshRenderer`] software
//! rendering path.

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Normalize a 3-component vector. Returns `[0,0,1]` on zero-length input.
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Reflect `incident` around `normal` (all unit vectors).
///
/// Formula: `r = incident - 2 * dot(incident, normal) * normal`
pub fn reflect_vec(incident: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = incident[0] * normal[0] + incident[1] * normal[1] + incident[2] * normal[2];
    [
        incident[0] - 2.0 * dot * normal[0],
        incident[1] - 2.0 * dot * normal[1],
        incident[2] - 2.0 * dot * normal[2],
    ]
}

// ---------------------------------------------------------------------------
// Orthographic projection
// ---------------------------------------------------------------------------

/// Project a world-space position via a fitted orthographic transform.
///
/// `scale` maps the [-1,1] world range to the central 80% of the viewport so
/// geometry does not clip at the screen edge.  The Y axis is flipped because
/// screen Y increases downward.
///
/// Returns `(screen_x, screen_y, depth_z)`.
pub fn project_ortho(pos: [f64; 3], width: u32, height: u32, scale: f64) -> ([f64; 2], f64) {
    let w = width as f64;
    let h = height as f64;
    // Map [-1/scale, +1/scale] → [0.1w, 0.9w] / [0.1h, 0.9h]
    let sx = (pos[0] * scale * 0.4 + 0.5) * w;
    let sy = (1.0 - (pos[1] * scale * 0.4 + 0.5)) * h; // flip Y
    ([sx, sy], pos[2])
}

/// Compute a scale factor that auto-fits `vertices` (flat xyz array) to fill
/// 80% of the viewport.  Falls back to `1.0` for empty meshes.
pub fn auto_fit_scale(vertices: &[f64]) -> f64 {
    if vertices.len() < 3 {
        return 1.0;
    }
    let mut max_abs = 0.0_f64;
    for chunk in vertices.chunks(3) {
        if let (Some(&x), Some(&y)) = (chunk.first(), chunk.get(1)) {
            max_abs = max_abs.max(x.abs()).max(y.abs());
        }
    }
    if max_abs < 1e-15 { 1.0 } else { 1.0 / max_abs }
}

// ---------------------------------------------------------------------------
// Phong shading
// ---------------------------------------------------------------------------

/// Compute Phong-shaded color for a single vertex.
///
/// `normal` and `view_dir` should be unit vectors in world space.
/// `light_dir` should be the **direction toward** the light, normalized.
///
/// Returns `[r, g, b]` in `[0, 1]` each.
pub fn phong_shade(
    normal: [f64; 3],
    view_dir: [f64; 3],
    light_dir: [f64; 3],
    diffuse_color: [f64; 3],
    ambient: f64,
    specular_strength: f64,
    shininess: f64,
) -> [f32; 3] {
    let n = normalize3(normal);
    let n_dot_l = (n[0] * light_dir[0] + n[1] * light_dir[1] + n[2] * light_dir[2]).max(0.0);
    let diffuse = n_dot_l;
    // Reflect light direction around normal
    let reflect = reflect_vec(
        [-light_dir[0], -light_dir[1], -light_dir[2]], // incident = -light_dir
        n,
    );
    let r_dot_v =
        (reflect[0] * view_dir[0] + reflect[1] * view_dir[1] + reflect[2] * view_dir[2]).max(0.0);
    let specular = specular_strength * r_dot_v.powf(shininess);
    let intensity = (ambient + diffuse + specular).clamp(0.0, 1.0);
    [
        (diffuse_color[0] * intensity) as f32,
        (diffuse_color[1] * intensity) as f32,
        (diffuse_color[2] * intensity) as f32,
    ]
}

// ---------------------------------------------------------------------------
// Triangle rasterizer
// ---------------------------------------------------------------------------

/// A vertex as seen by the rasterizer: screen-space position, depth, RGBA color.
#[derive(Copy, Clone, Debug)]
pub struct RastVertex {
    /// Screen-space position `[x, y]` in pixels.
    pub pos: [f64; 2],
    /// Depth (smaller = closer to viewer).
    pub depth: f64,
    /// RGBA color in `[0, 1]` per channel.
    pub color: [f32; 4],
}

/// Rasterize one triangle into `pixels` (flat RGBA u8) respecting `z_buffer`.
///
/// Both `pixels` and `z_buffer` are indexed as `y * width + x`.
/// A pixel is written only if its interpolated depth is **less than** the
/// current z-buffer value (front-wins convention).
pub fn rasterize_triangle(
    pixels: &mut [u8],
    z_buffer: &mut [f64],
    width: u32,
    height: u32,
    v0: RastVertex,
    v1: RastVertex,
    v2: RastVertex,
) {
    if width == 0 || height == 0 {
        return;
    }

    // ── Bounding box clamped to viewport ────────────────────────────────────
    let min_x = v0.pos[0].min(v1.pos[0]).min(v2.pos[0]).max(0.0).floor() as u32;
    let max_x =
        (v0.pos[0].max(v1.pos[0]).max(v2.pos[0]).ceil() as u32).min(width.saturating_sub(1));
    let min_y = v0.pos[1].min(v1.pos[1]).min(v2.pos[1]).max(0.0).floor() as u32;
    let max_y =
        (v0.pos[1].max(v1.pos[1]).max(v2.pos[1]).ceil() as u32).min(height.saturating_sub(1));

    if min_x > max_x || min_y > max_y {
        return;
    }

    // Edge function: positive inside for CCW; we support both windings.
    let edge = |a: [f64; 2], b: [f64; 2], p: [f64; 2]| -> f64 {
        (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    };

    let area = edge(v0.pos, v1.pos, v2.pos);
    if area.abs() < 1e-8 {
        return; // degenerate triangle
    }

    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let p = [x as f64 + 0.5, y as f64 + 0.5];
            let w0 = edge(v1.pos, v2.pos, p);
            let w1 = edge(v2.pos, v0.pos, p);
            let w2 = edge(v0.pos, v1.pos, p);

            // Accept both CCW (all >= 0) and CW (all <= 0) winding.
            let inside =
                (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0) || (w0 <= 0.0 && w1 <= 0.0 && w2 <= 0.0);
            if !inside {
                continue;
            }

            // Barycentric coordinates (sign of area cancels).
            let bc0 = w0 / area;
            let bc1 = w1 / area;
            let bc2 = w2 / area;

            let depth = bc0 * v0.depth + bc1 * v1.depth + bc2 * v2.depth;

            let idx = (y * width + x) as usize;
            if idx >= z_buffer.len() {
                continue;
            }
            if depth < z_buffer[idx] {
                z_buffer[idx] = depth;
                let r = (bc0 * v0.color[0] as f64
                    + bc1 * v1.color[0] as f64
                    + bc2 * v2.color[0] as f64)
                    .clamp(0.0, 1.0);
                let g = (bc0 * v0.color[1] as f64
                    + bc1 * v1.color[1] as f64
                    + bc2 * v2.color[1] as f64)
                    .clamp(0.0, 1.0);
                let b = (bc0 * v0.color[2] as f64
                    + bc1 * v1.color[2] as f64
                    + bc2 * v2.color[2] as f64)
                    .clamp(0.0, 1.0);
                let base = idx * 4;
                if base + 3 < pixels.len() {
                    pixels[base] = (r * 255.0) as u8;
                    pixels[base + 1] = (g * 255.0) as u8;
                    pixels[base + 2] = (b * 255.0) as u8;
                    pixels[base + 3] = 255;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// High-level render entry point
// ---------------------------------------------------------------------------

/// Geometry and colour data for a single mesh render call.
pub struct MeshRenderData<'a> {
    pub vertices: &'a [f64],
    pub normals: &'a [f64],
    pub indices: &'a [u32],
    pub base_color: [f64; 3],
    pub scalars: Option<&'a [f64]>,
    pub scalar_map_fn: Option<&'a dyn Fn(f64) -> [f64; 3]>,
}

/// Render a single [`super::viz_api::PyMeshRenderer`] into a flat RGBA buffer.
///
/// Uses an orthographic projection auto-fitted to the mesh bounding box,
/// with a single directional light from `[1, 1, 1]` and Phong shading.
///
/// `background` is written to every pixel before rendering.
pub fn render_mesh(
    data: MeshRenderData<'_>,
    width: u32,
    height: u32,
    background: [u8; 4],
) -> Vec<u8> {
    let MeshRenderData {
        vertices,
        normals,
        indices,
        base_color,
        scalars,
        scalar_map_fn,
    } = data;
    let pixel_count = (width * height) as usize;
    let mut pixels = vec![0u8; pixel_count * 4];

    // Fill background
    for i in 0..pixel_count {
        let base = i * 4;
        pixels[base] = background[0];
        pixels[base + 1] = background[1];
        pixels[base + 2] = background[2];
        pixels[base + 3] = background[3];
    }

    if width == 0 || height == 0 || indices.len() < 3 || vertices.len() < 3 {
        return pixels;
    }

    let mut z_buffer = vec![f64::MAX; pixel_count];

    // Auto-fit orthographic scale to mesh extents.
    let scale = auto_fit_scale(vertices);

    // Fixed light direction (toward viewer-upper-right, world space).
    let light_dir = normalize3([1.0, 1.0, 1.0]);
    // View direction is +Z (ortho camera looks down -Z axis).
    let view_dir = [0.0, 0.0, 1.0_f64];

    let num_vertices = vertices.len() / 3;

    // Process every triangle.
    for tri in indices.chunks(3) {
        let (i0, i1, i2) = match (tri.first(), tri.get(1), tri.get(2)) {
            (Some(&a), Some(&b), Some(&c)) => (a as usize, b as usize, c as usize),
            _ => continue,
        };
        if i0 >= num_vertices || i1 >= num_vertices || i2 >= num_vertices {
            continue;
        }

        let make_vertex = |vi: usize| -> RastVertex {
            let pos3 = [vertices[3 * vi], vertices[3 * vi + 1], vertices[3 * vi + 2]];
            let (screen, depth) = project_ortho(pos3, width, height, scale);

            // Get normal for this vertex.
            let normal = if 3 * vi + 2 < normals.len() {
                [normals[3 * vi], normals[3 * vi + 1], normals[3 * vi + 2]]
            } else {
                [0.0, 0.0, 1.0]
            };

            // Determine diffuse color: scalar map takes priority over base_color.
            let diffuse = if let (Some(sc), Some(map_fn)) = (scalars, scalar_map_fn) {
                let sv = sc.get(vi).copied().unwrap_or(0.0);
                map_fn(sv)
            } else {
                base_color
            };

            let rgb = phong_shade(
                normal, view_dir, light_dir, diffuse, 0.15, // ambient
                0.4,  // specular strength
                16.0, // shininess
            );

            RastVertex {
                pos: screen,
                depth,
                color: [rgb[0], rgb[1], rgb[2], 1.0],
            }
        };

        let rv0 = make_vertex(i0);
        let rv1 = make_vertex(i1);
        let rv2 = make_vertex(i2);

        rasterize_triangle(&mut pixels, &mut z_buffer, width, height, rv0, rv1, rv2);
    }

    pixels
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helper: build a flat triangle in the XY plane pointing toward +Z ────

    fn flat_triangle_vertices() -> Vec<f64> {
        // Three vertices forming a CCW triangle in XY plane, centered near origin.
        vec![-0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.0, 0.5, 0.0]
    }

    fn flat_triangle_normals() -> Vec<f64> {
        // All normals point in +Z.
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    }

    fn flat_triangle_indices() -> Vec<u32> {
        vec![0, 1, 2]
    }

    // ── 1. render returns correct size ──────────────────────────────────────

    #[test]
    fn test_render_returns_correct_size() {
        let (w, h) = (64u32, 48u32);
        let buf = render_mesh(
            MeshRenderData {
                vertices: &flat_triangle_vertices(),
                normals: &flat_triangle_normals(),
                indices: &flat_triangle_indices(),
                base_color: [0.8, 0.2, 0.2],
                scalars: None,
                scalar_map_fn: None,
            },
            w,
            h,
            [0, 0, 0, 255],
        );
        assert_eq!(buf.len(), (w * h * 4) as usize);
    }

    // ── 2. render differs from empty background after adding mesh ───────────

    #[test]
    fn test_render_not_all_background() {
        let (w, h) = (64u32, 64u32);
        let bg = [20u8, 20, 20, 255];
        let buf = render_mesh(
            MeshRenderData {
                vertices: &flat_triangle_vertices(),
                normals: &flat_triangle_normals(),
                indices: &flat_triangle_indices(),
                base_color: [0.8, 0.5, 0.1],
                scalars: None,
                scalar_map_fn: None,
            },
            w,
            h,
            bg,
        );
        // At least one pixel should differ from the background.
        let differs = buf
            .chunks(4)
            .any(|px| px[0] != bg[0] || px[1] != bg[1] || px[2] != bg[2]);
        assert!(
            differs,
            "render produced no pixels different from background"
        );
    }

    // ── 3. depth test: front triangle wins ──────────────────────────────────

    #[test]
    fn test_render_depth_test_correct() {
        // Two overlapping triangles: one red at z=0 (front), one blue at z=1.
        let (w, h) = (64u32, 64u32);

        let verts_front = vec![
            -0.5, -0.5, 0.0, // z=0  (closer)
            0.5, -0.5, 0.0, 0.0, 0.5, 0.0,
        ];
        let verts_back = vec![
            -0.5, -0.5, 1.0, // z=1  (farther)
            0.5, -0.5, 1.0, 0.0, 0.5, 1.0,
        ];
        let normals = flat_triangle_normals();
        let indices = flat_triangle_indices();
        let bg = [0u8, 0, 0, 255];

        // Render front (red) first.
        let buf_front = render_mesh(
            MeshRenderData {
                vertices: &verts_front,
                normals: &normals,
                indices: &indices,
                base_color: [1.0, 0.0, 0.0],
                scalars: None,
                scalar_map_fn: None,
            },
            w,
            h,
            bg,
        );

        // Render back (blue) first, then front (red) on top.
        // We do this by compositing: render back into separate buffer and manually
        // check the z-buffer logic.
        let buf_back = render_mesh(
            MeshRenderData {
                vertices: &verts_back,
                normals: &normals,
                indices: &indices,
                base_color: [0.0, 0.0, 1.0],
                scalars: None,
                scalar_map_fn: None,
            },
            w,
            h,
            bg,
        );

        // Find a pixel that the front triangle covers.
        let pixel_count = (w * h) as usize;
        let mut found_front_color = false;
        for i in 0..pixel_count {
            let fr = buf_front[4 * i];
            let fg = buf_front[4 * i + 1];
            // Front triangle is red → red channel should be dominant.
            if fr > 50 && fg < 50 {
                found_front_color = true;
                break;
            }
        }
        assert!(found_front_color, "front triangle pixels not found");

        // Back triangle at z=1 → blue channel should be dominant somewhere.
        let mut found_back_color = false;
        for i in 0..pixel_count {
            let br = buf_back[4 * i];
            let bb = buf_back[4 * i + 2];
            if bb > 50 && br < 50 {
                found_back_color = true;
                break;
            }
        }
        assert!(found_back_color, "back triangle pixels not found");

        // The combined render (front wins): simulate two-pass via the z-buffer.
        let mut pixels = vec![0u8; pixel_count * 4];
        let mut z_buf = vec![f64::MAX; pixel_count];
        for i in 0..pixel_count {
            pixels[4 * i] = bg[0];
            pixels[4 * i + 1] = bg[1];
            pixels[4 * i + 2] = bg[2];
            pixels[4 * i + 3] = bg[3];
        }
        // Draw back (blue) first, then front (red) — z-buffer should keep red.
        let scale = auto_fit_scale(&verts_back);
        let light_dir = normalize3([1.0, 1.0, 1.0]);
        let view_dir = [0.0, 0.0, 1.0_f64];

        let to_rv_back = |vi: usize| -> RastVertex {
            let p = [
                verts_back[3 * vi],
                verts_back[3 * vi + 1],
                verts_back[3 * vi + 2],
            ];
            let (screen, depth) = project_ortho(p, w, h, scale);
            let rgb = phong_shade(
                [0.0, 0.0, 1.0],
                view_dir,
                light_dir,
                [0.0, 0.0, 1.0],
                0.15,
                0.4,
                16.0,
            );
            RastVertex {
                pos: screen,
                depth,
                color: [rgb[0], rgb[1], rgb[2], 1.0],
            }
        };
        let to_rv_front = |vi: usize| -> RastVertex {
            let p = [
                verts_front[3 * vi],
                verts_front[3 * vi + 1],
                verts_front[3 * vi + 2],
            ];
            let (screen, depth) = project_ortho(p, w, h, scale);
            let rgb = phong_shade(
                [0.0, 0.0, 1.0],
                view_dir,
                light_dir,
                [1.0, 0.0, 0.0],
                0.15,
                0.4,
                16.0,
            );
            RastVertex {
                pos: screen,
                depth,
                color: [rgb[0], rgb[1], rgb[2], 1.0],
            }
        };

        rasterize_triangle(
            &mut pixels,
            &mut z_buf,
            w,
            h,
            to_rv_back(0),
            to_rv_back(1),
            to_rv_back(2),
        );
        rasterize_triangle(
            &mut pixels,
            &mut z_buf,
            w,
            h,
            to_rv_front(0),
            to_rv_front(1),
            to_rv_front(2),
        );

        // Any pixel that was covered should have red > blue (front wins).
        let mut checked = false;
        for i in 0..pixel_count {
            let r = pixels[4 * i];
            let b = pixels[4 * i + 2];
            if r > 50 || b > 50 {
                assert!(
                    r >= b,
                    "depth test failed: blue (back) pixel won at index {i}, r={r} b={b}"
                );
                checked = true;
            }
        }
        assert!(checked, "no non-background pixels found in depth test");
    }

    // ── 4. Phong diffuse increases with normal aligned to light ─────────────

    #[test]
    fn test_render_phong_shade_diffuse_increases_with_alignment() {
        let light_dir = normalize3([0.0, 0.0, 1.0]);
        let view_dir = [0.0, 0.0, 1.0_f64];
        let diffuse_color = [1.0, 1.0, 1.0];

        // Normal fully aligned to light → maximum diffuse.
        let aligned = phong_shade(
            [0.0, 0.0, 1.0],
            view_dir,
            light_dir,
            diffuse_color,
            0.0,
            0.0,
            1.0,
        );
        // Normal perpendicular to light → zero diffuse contribution.
        let perp = phong_shade(
            [1.0, 0.0, 0.0],
            view_dir,
            light_dir,
            diffuse_color,
            0.0,
            0.0,
            1.0,
        );

        assert!(
            aligned[0] > perp[0],
            "aligned normal ({}) should produce more light than perpendicular ({})",
            aligned[0],
            perp[0]
        );
    }

    // ── 5. z_buffer initialized to MAX before rendering ─────────────────────

    #[test]
    fn test_render_z_buffer_initialized() {
        // Indirectly test by verifying a back triangle is visible in an
        // otherwise empty render (no front triangle to block it).
        let (w, h) = (32u32, 32u32);
        let back_verts = vec![-0.5, -0.5, 5.0, 0.5, -0.5, 5.0, 0.0, 0.5, 5.0];
        let bg = [0u8, 0, 0, 255];
        let buf = render_mesh(
            MeshRenderData {
                vertices: &back_verts,
                normals: &flat_triangle_normals(),
                indices: &flat_triangle_indices(),
                base_color: [0.0, 1.0, 0.0],
                scalars: None,
                scalar_map_fn: None,
            },
            w,
            h,
            bg,
        );
        // If z_buffer started at MAX (not 0), the green triangle is visible.
        let has_green = buf
            .chunks(4)
            .any(|px| px[1] > 50 && px[0] < 50 && px[2] < 50);
        assert!(
            has_green,
            "back triangle not visible – z_buffer may not have been initialized to MAX"
        );
    }

    // ── auto_fit_scale edge cases ─────────────────────────────────────────

    #[test]
    fn test_auto_fit_scale_empty() {
        assert_eq!(auto_fit_scale(&[]), 1.0);
    }

    #[test]
    fn test_auto_fit_scale_unit_box() {
        // Vertices at ±1 → scale = 1.0
        let verts = vec![-1.0f64, -1.0, -1.0, 1.0, 1.0, 1.0];
        let s = auto_fit_scale(&verts);
        assert!((s - 1.0).abs() < 1e-9, "scale={s}");
    }

    #[test]
    fn test_normalize3_zero_vec() {
        let n = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 1.0]);
    }
}
