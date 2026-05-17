// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly canvas rendering path.
//!
//! This module provides a CPU-side rasterizer that writes directly into a flat
//! `Vec<u8>` of RGBA bytes compatible with the HTML5 Canvas 2D API
//! [`ImageData`](https://developer.mozilla.org/en-US/docs/Web/API/ImageData).
//!
//! In a Wasm build the caller can hand this slice to JavaScript via
//! `Uint8ClampedArray.set(new Uint8ClampedArray(memory.buffer, ptr, len))` and
//! then call `ctx.putImageData(imageData, 0, 0)` to display the rendered frame.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_viz::wasm_canvas::{CanvasBuffer, CanvasConfig, CanvasRasterizer};
//!
//! let cfg = CanvasConfig { width: 320, height: 240, ..Default::default() };
//! let mut buf = CanvasBuffer::new(&cfg);
//!
//! let rast = CanvasRasterizer::new(cfg.clone());
//! rast.clear(&mut buf, [0, 0, 0, 255]);
//! rast.draw_filled_rect(&mut buf, 10, 10, 50, 50, [255, 0, 0, 255]);
//!
//! assert_eq!(buf.bytes.len(), 320 * 240 * 4);
//! ```

#![allow(dead_code)]

// ── CanvasConfig ──────────────────────────────────────────────────────────────

/// Configuration for a canvas rendering target.
#[derive(Debug, Clone)]
pub struct CanvasConfig {
    /// Canvas width in pixels.
    pub width: usize,
    /// Canvas height in pixels.
    pub height: usize,
    /// Pixel density multiplier for high-DPI displays (1 = standard, 2 = retina).
    pub device_pixel_ratio: f32,
    /// Whether to pre-multiply alpha values in the output.
    pub premultiply_alpha: bool,
}

impl Default for CanvasConfig {
    fn default() -> Self {
        Self {
            width: 800,
            height: 600,
            device_pixel_ratio: 1.0,
            premultiply_alpha: false,
        }
    }
}

impl CanvasConfig {
    /// Total number of pixels (width × height).
    pub fn pixel_count(&self) -> usize {
        self.width * self.height
    }

    /// Total byte length of the RGBA buffer (pixel_count × 4).
    pub fn byte_len(&self) -> usize {
        self.pixel_count() * 4
    }
}

// ── CanvasBuffer ─────────────────────────────────────────────────────────────

/// A flat RGBA byte buffer compatible with HTML5 Canvas `ImageData`.
///
/// Pixel `(x, y)` starts at byte offset `(y * width + x) * 4`.
#[derive(Debug, Clone)]
pub struct CanvasBuffer {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Raw RGBA bytes (length = width × height × 4).
    pub bytes: Vec<u8>,
}

impl CanvasBuffer {
    /// Create a buffer cleared to opaque black.
    pub fn new(cfg: &CanvasConfig) -> Self {
        let len = cfg.byte_len();
        let mut bytes = vec![0u8; len];
        // Set alpha to 255 (fully opaque)
        for i in 0..cfg.pixel_count() {
            bytes[i * 4 + 3] = 255;
        }
        Self {
            width: cfg.width,
            height: cfg.height,
            bytes,
        }
    }

    /// Return the byte offset for pixel `(x, y)`.
    #[inline]
    pub fn offset(&self, x: usize, y: usize) -> usize {
        (y * self.width + x) * 4
    }

    /// Read pixel `(x, y)` as RGBA bytes.
    #[inline]
    pub fn get_pixel(&self, x: usize, y: usize) -> [u8; 4] {
        let o = self.offset(x, y);
        [
            self.bytes[o],
            self.bytes[o + 1],
            self.bytes[o + 2],
            self.bytes[o + 3],
        ]
    }

    /// Write pixel `(x, y)` with RGBA bytes `rgba`.
    #[inline]
    pub fn set_pixel(&mut self, x: usize, y: usize, rgba: [u8; 4]) {
        let o = self.offset(x, y);
        self.bytes[o] = rgba[0];
        self.bytes[o + 1] = rgba[1];
        self.bytes[o + 2] = rgba[2];
        self.bytes[o + 3] = rgba[3];
    }

    /// Alpha-blend `src` over the existing pixel at `(x, y)`.
    ///
    /// Uses standard Porter-Duff "over" compositing.
    #[inline]
    pub fn blend_pixel(&mut self, x: usize, y: usize, src: [u8; 4]) {
        let dst = self.get_pixel(x, y);
        let sa = src[3] as f32 / 255.0;
        let da = dst[3] as f32 / 255.0;
        let out_a = sa + da * (1.0 - sa);
        if out_a < 1e-6 {
            self.set_pixel(x, y, [0, 0, 0, 0]);
            return;
        }
        let blend = |sc: u8, dc: u8| -> u8 {
            let out = (sc as f32 * sa + dc as f32 * da * (1.0 - sa)) / out_a;
            out.clamp(0.0, 255.0) as u8
        };
        self.set_pixel(
            x,
            y,
            [
                blend(src[0], dst[0]),
                blend(src[1], dst[1]),
                blend(src[2], dst[2]),
                (out_a * 255.0 + 0.5) as u8,
            ],
        );
    }

    /// Return a raw pointer to the byte buffer for Wasm FFI.
    ///
    /// # Safety
    ///
    /// The returned pointer is valid as long as `self` is alive and the buffer
    /// is not reallocated.  The caller must not write past `byte_len()` bytes.
    pub fn as_ptr(&self) -> *const u8 {
        self.bytes.as_ptr()
    }

    /// Return the total byte length of the buffer.
    pub fn byte_len(&self) -> usize {
        self.bytes.len()
    }
}

// ── CanvasRasterizer ─────────────────────────────────────────────────────────

/// CPU rasterizer that writes RGBA output into a [`CanvasBuffer`].
///
/// All coordinates are in pixel space with the origin at the top-left corner,
/// matching HTML5 Canvas conventions.
#[derive(Debug, Clone)]
pub struct CanvasRasterizer {
    /// Canvas configuration.
    pub config: CanvasConfig,
}

impl CanvasRasterizer {
    /// Create a rasterizer for the given canvas configuration.
    pub fn new(config: CanvasConfig) -> Self {
        Self { config }
    }

    // ── Clearing ───────────────────────────────────────────────────────────

    /// Clear the entire buffer to `rgba`.
    pub fn clear(&self, buf: &mut CanvasBuffer, rgba: [u8; 4]) {
        for i in 0..buf.width * buf.height {
            let o = i * 4;
            buf.bytes[o] = rgba[0];
            buf.bytes[o + 1] = rgba[1];
            buf.bytes[o + 2] = rgba[2];
            buf.bytes[o + 3] = rgba[3];
        }
    }

    // ── Filled primitives ─────────────────────────────────────────────────

    /// Draw a filled axis-aligned rectangle.
    pub fn draw_filled_rect(
        &self,
        buf: &mut CanvasBuffer,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        rgba: [u8; 4],
    ) {
        let x1 = x1.min(buf.width);
        let y1 = y1.min(buf.height);
        for y in y0..y1 {
            for x in x0..x1 {
                buf.set_pixel(x, y, rgba);
            }
        }
    }

    /// Draw a filled circle using the midpoint algorithm.
    pub fn draw_filled_circle(
        &self,
        buf: &mut CanvasBuffer,
        cx: i32,
        cy: i32,
        radius: i32,
        rgba: [u8; 4],
    ) {
        let r2 = radius * radius;
        for dy in -radius..=radius {
            let y = cy + dy;
            if y < 0 || y >= buf.height as i32 {
                continue;
            }
            let xspan = ((r2 - dy * dy) as f32).sqrt() as i32;
            let x0 = (cx - xspan).max(0) as usize;
            let x1 = (cx + xspan).min(buf.width as i32 - 1) as usize;
            for x in x0..=x1 {
                buf.set_pixel(x, y as usize, rgba);
            }
        }
    }

    /// Draw a filled triangle with flat-top/flat-bottom decomposition.
    pub fn draw_filled_triangle(
        &self,
        buf: &mut CanvasBuffer,
        v0: [i32; 2],
        v1: [i32; 2],
        v2: [i32; 2],
        rgba: [u8; 4],
    ) {
        // Sort vertices by Y
        let mut verts = [v0, v1, v2];
        verts.sort_by_key(|v| v[1]);
        let [a, b, c] = verts;

        // Rasterize via horizontal scan lines
        let y_start = a[1].max(0);
        let y_end = c[1].min(buf.height as i32 - 1);
        for y in y_start..=y_end {
            let t = if c[1] != a[1] {
                (y - a[1]) as f32 / (c[1] - a[1]) as f32
            } else {
                0.0
            };
            let x_left = a[0] as f32 + t * (c[0] - a[0]) as f32;
            let x_right = if y <= b[1] {
                let dt = if b[1] != a[1] {
                    (y - a[1]) as f32 / (b[1] - a[1]) as f32
                } else {
                    1.0
                };
                a[0] as f32 + dt * (b[0] - a[0]) as f32
            } else {
                let dt = if c[1] != b[1] {
                    (y - b[1]) as f32 / (c[1] - b[1]) as f32
                } else {
                    1.0
                };
                b[0] as f32 + dt * (c[0] - b[0]) as f32
            };
            let (xl, xr) = if x_left < x_right {
                (x_left as i32, x_right as i32)
            } else {
                (x_right as i32, x_left as i32)
            };
            for x in xl.max(0)..=xr.min(buf.width as i32 - 1) {
                buf.set_pixel(x as usize, y as usize, rgba);
            }
        }
    }

    // ── Outlined primitives ────────────────────────────────────────────────

    /// Draw the outline of an axis-aligned rectangle.
    pub fn draw_rect_outline(
        &self,
        buf: &mut CanvasBuffer,
        x0: usize,
        y0: usize,
        x1: usize,
        y1: usize,
        rgba: [u8; 4],
    ) {
        let x1 = x1.min(buf.width - 1);
        let y1 = y1.min(buf.height - 1);
        for x in x0..=x1 {
            buf.set_pixel(x, y0, rgba);
        }
        for x in x0..=x1 {
            buf.set_pixel(x, y1, rgba);
        }
        for y in y0..=y1 {
            buf.set_pixel(x0, y, rgba);
        }
        for y in y0..=y1 {
            buf.set_pixel(x1, y, rgba);
        }
    }

    /// Draw a line using Bresenham's algorithm.
    pub fn draw_line(
        &self,
        buf: &mut CanvasBuffer,
        mut x0: i32,
        mut y0: i32,
        x1: i32,
        y1: i32,
        rgba: [u8; 4],
    ) {
        let dx = (x1 - x0).abs();
        let dy = -(y1 - y0).abs();
        let sx: i32 = if x0 < x1 { 1 } else { -1 };
        let sy: i32 = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            if x0 >= 0 && x0 < buf.width as i32 && y0 >= 0 && y0 < buf.height as i32 {
                buf.set_pixel(x0 as usize, y0 as usize, rgba);
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }

    // ── Text rendering ────────────────────────────────────────────────────────

    /// Render a string of text at pixel position `(x, y)` using a [`crate::font_rendering::BitmapFont`].
    ///
    /// Each character occupies an 8×8 pixel cell.  Characters outside the
    /// printable ASCII range (0x20–0x7E) that have no glyph are silently
    /// skipped.  Pixels that would fall outside the canvas bounds are clipped.
    ///
    /// `rgba` is the foreground colour applied to every set bit.
    pub fn draw_text(
        &self,
        buf: &mut CanvasBuffer,
        x: usize,
        y: usize,
        text: &str,
        font: &crate::font_rendering::BitmapFont,
        rgba: [u8; 4],
    ) {
        let mut char_x = x;
        for ch in text.chars() {
            if let Some(glyph) = font.get_glyph(ch) {
                let w = glyph.width;
                let h = glyph.height;
                for row in 0..h {
                    let py = y + row;
                    if py >= buf.height {
                        break;
                    }
                    for col in 0..w {
                        let px = char_x + col;
                        if px >= buf.width {
                            break;
                        }
                        if glyph.get_pixel(col, row) != 0 {
                            buf.set_pixel(px, py, rgba);
                        }
                    }
                }
                // Advance by the glyph's advance_x (rounded to nearest pixel)
                char_x += glyph.advance_x.round() as usize;
            } else {
                // Unknown character: advance by half a cell width
                char_x += 4;
            }
            if char_x >= buf.width {
                break;
            }
        }
    }

    /// Draw a simple scale indicator / progress bar.
    ///
    /// `rect` is `[x, y, width, height]`, `value` is clamped to `[0, 1]`.
    pub fn draw_progress_bar(
        &self,
        buf: &mut CanvasBuffer,
        rect: [usize; 4],
        value: f32,
        fill_rgba: [u8; 4],
        bg_rgba: [u8; 4],
    ) {
        let [x, y, width, height] = rect;
        // Background
        self.draw_filled_rect(buf, x, y, x + width, y + height, bg_rgba);
        // Fill
        let fill_w = ((value.clamp(0.0, 1.0) * width as f32) as usize).min(width);
        if fill_w > 0 {
            self.draw_filled_rect(buf, x, y, x + fill_w, y + height, fill_rgba);
        }
    }

    // ── Physics debug overlays ────────────────────────────────────────────────

    /// Draw a particle as a filled circle at the given screen coordinate.
    pub fn draw_particle(
        &self,
        buf: &mut CanvasBuffer,
        x: i32,
        y: i32,
        radius: i32,
        rgba: [u8; 4],
    ) {
        self.draw_filled_circle(buf, x, y, radius, rgba);
    }

    /// Draw a velocity arrow: a line from `pos` in direction `vel` scaled by `scale`.
    ///
    /// `pos` is `[x, y]` (screen coords), `vel` is `[vx, vy]`.
    pub fn draw_velocity_arrow(
        &self,
        buf: &mut CanvasBuffer,
        pos: [i32; 2],
        vel: [f32; 2],
        scale: f32,
        rgba: [u8; 4],
    ) {
        let [x, y] = pos;
        let [vx, vy] = vel;
        let tx = x + (vx * scale) as i32;
        let ty = y + (vy * scale) as i32;
        self.draw_line(buf, x, y, tx, ty, rgba);
        // Arrowhead (3-pixel cross)
        self.draw_line(buf, tx - 2, ty - 2, tx + 2, ty + 2, rgba);
        self.draw_line(buf, tx + 2, ty - 2, tx - 2, ty + 2, rgba);
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canvas_buffer_correct_size() {
        let cfg = CanvasConfig {
            width: 320,
            height: 240,
            ..Default::default()
        };
        let buf = CanvasBuffer::new(&cfg);
        assert_eq!(buf.bytes.len(), 320 * 240 * 4);
    }

    #[test]
    fn set_and_get_pixel() {
        let cfg = CanvasConfig {
            width: 10,
            height: 10,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        buf.set_pixel(3, 4, [255, 128, 0, 200]);
        assert_eq!(buf.get_pixel(3, 4), [255, 128, 0, 200]);
    }

    #[test]
    fn clear_fills_all_pixels() {
        let cfg = CanvasConfig {
            width: 8,
            height: 8,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        let rast = CanvasRasterizer::new(cfg);
        rast.clear(&mut buf, [42, 43, 44, 255]);
        for i in 0..8 * 8 {
            let o = i * 4;
            assert_eq!(&buf.bytes[o..o + 4], &[42, 43, 44, 255]);
        }
    }

    #[test]
    fn draw_filled_rect_pixels_correct() {
        let cfg = CanvasConfig {
            width: 20,
            height: 20,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        let rast = CanvasRasterizer::new(cfg);
        rast.clear(&mut buf, [0, 0, 0, 255]);
        rast.draw_filled_rect(&mut buf, 5, 5, 10, 10, [255, 0, 0, 255]);
        // Check inside rect
        assert_eq!(buf.get_pixel(5, 5), [255, 0, 0, 255]);
        assert_eq!(buf.get_pixel(9, 9), [255, 0, 0, 255]);
        // Check outside
        assert_eq!(buf.get_pixel(4, 5), [0, 0, 0, 255]);
    }

    #[test]
    fn draw_line_bresenham() {
        let cfg = CanvasConfig {
            width: 20,
            height: 20,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        let rast = CanvasRasterizer::new(cfg);
        rast.clear(&mut buf, [0, 0, 0, 255]);
        rast.draw_line(&mut buf, 0, 0, 10, 10, [255, 255, 255, 255]);
        // Diagonal should have pixels set
        assert_eq!(buf.get_pixel(0, 0), [255, 255, 255, 255]);
        assert_eq!(buf.get_pixel(10, 10), [255, 255, 255, 255]);
    }

    #[test]
    fn blend_pixel_transparent_overlay() {
        let cfg = CanvasConfig {
            width: 4,
            height: 4,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        buf.set_pixel(1, 1, [0, 0, 0, 255]); // opaque black
        buf.blend_pixel(1, 1, [255, 255, 255, 128]); // semi-transparent white
        let px = buf.get_pixel(1, 1);
        // Result should be grey-ish
        assert!(px[0] > 50, "expected grey after blend, got {:?}", px);
    }

    // ── draw_text (F5) ────────────────────────────────────────────────────────

    #[test]
    fn draw_text_abc_produces_pixels() {
        use crate::font_rendering::BitmapFont;

        let cfg = CanvasConfig {
            width: 200,
            height: 20,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        let rast = CanvasRasterizer::new(cfg);
        rast.clear(&mut buf, [0, 0, 0, 255]);

        let font = BitmapFont::ascii_default();
        rast.draw_text(&mut buf, 0, 0, "ABC", &font, [255, 255, 255, 255]);

        let white_pixels = (0..200usize)
            .flat_map(|x| (0..20usize).map(move |y| (x, y)))
            .filter(|&(x, y)| buf.get_pixel(x, y)[0] == 255)
            .count();
        assert!(
            white_pixels >= 3,
            "draw_text 'ABC' should produce at least 3 foreground pixels, got {white_pixels}"
        );
    }

    #[test]
    fn draw_text_three_distinct_regions() {
        use crate::font_rendering::BitmapFont;

        let cfg = CanvasConfig {
            width: 300,
            height: 16,
            ..Default::default()
        };
        let mut buf = CanvasBuffer::new(&cfg);
        let rast = CanvasRasterizer::new(cfg);
        rast.clear(&mut buf, [0, 0, 0, 255]);

        let font = BitmapFont::ascii_default();
        rast.draw_text(&mut buf, 0, 0, "A", &font, [255, 0, 0, 255]);
        rast.draw_text(&mut buf, 100, 0, "B", &font, [0, 255, 0, 255]);
        rast.draw_text(&mut buf, 200, 0, "C", &font, [0, 0, 255, 255]);

        let has_red = (0..8usize).any(|x| {
            (0..16usize).any(|y| {
                let px = buf.get_pixel(x, y);
                px[0] == 255 && px[1] == 0 && px[2] == 0
            })
        });
        let has_green = (100..108usize).any(|x| {
            (0..16usize).any(|y| {
                let px = buf.get_pixel(x, y);
                px[0] == 0 && px[1] == 255 && px[2] == 0
            })
        });
        let has_blue = (200..208usize).any(|x| {
            (0..16usize).any(|y| {
                let px = buf.get_pixel(x, y);
                px[0] == 0 && px[1] == 0 && px[2] == 255
            })
        });

        assert!(has_red, "first character region should have red pixels");
        assert!(
            has_green,
            "second character region should have green pixels"
        );
        assert!(has_blue, "third character region should have blue pixels");
    }
}
