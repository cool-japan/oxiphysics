// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Text and label rendering for scientific visualization overlays.
//!
//! Provides a CPU-side 8×8 bitmap font rasterizer, text measurement utilities,
//! 3-D-to-2-D label projection, axis label generation and flexible number
//! formatting — without any GPU or windowing dependencies.

// ──────────────────────────────────────────────────────────────────────────────
// FontMetrics
// ──────────────────────────────────────────────────────────────────────────────

/// Metric data for a bitmap font.
#[derive(Debug, Clone)]
pub struct FontMetrics {
    /// Width of each character cell in pixels.
    pub char_width: u32,
    /// Height of each character cell in pixels.
    pub char_height: u32,
    /// Vertical distance between baselines (line height) in pixels.
    pub line_height: u32,
    /// Horizontal gap between characters in pixels.
    pub char_spacing: u32,
}

impl Default for FontMetrics {
    fn default() -> Self {
        Self {
            char_width: 8,
            char_height: 8,
            line_height: 10,
            char_spacing: 1,
        }
    }
}

impl FontMetrics {
    /// Construct metrics for the built-in 8×8 monospace font.
    pub fn monospace_8x8() -> Self {
        Self::default()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// BitmapFont
// ──────────────────────────────────────────────────────────────────────────────

/// An 8×8 monospace bitmap font backed by a glyph table.
///
/// Each glyph is stored as 8 bytes, one per row (MSB = leftmost pixel).
/// Only printable ASCII (0x20–0x7E) is supported; other characters map to a
/// default replacement glyph (a filled rectangle).
#[derive(Debug, Clone)]
pub struct BitmapFont {
    /// Glyph data: 95 glyphs × 8 rows.  Index 0 = space (0x20).
    pub glyphs: Vec<[u8; 8]>,
    /// Metrics for this font.
    pub metrics: FontMetrics,
}

impl Default for BitmapFont {
    fn default() -> Self {
        Self::new()
    }
}

impl BitmapFont {
    /// Construct the built-in 8×8 font.
    ///
    /// Only a small subset of glyphs are hand-coded for illustration; the rest
    /// are filled with a generic block glyph so that every printable ASCII
    /// character maps to _something_ visible.
    pub fn new() -> Self {
        let mut glyphs = vec![[0u8; 8]; 95];

        // Space (0x20 → index 0)
        glyphs[0] = [0x00; 8];

        // '!' (0x21)
        glyphs[0x21 - 0x20] = [0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x18, 0x00];

        // '"' (0x22)
        glyphs[0x22 - 0x20] = [0x66, 0x66, 0x66, 0x00, 0x00, 0x00, 0x00, 0x00];

        // '0' (0x30)
        glyphs[0x30 - 0x20] = [0x3C, 0x66, 0x6E, 0x76, 0x66, 0x66, 0x3C, 0x00];
        // '1'
        glyphs[0x31 - 0x20] = [0x18, 0x38, 0x18, 0x18, 0x18, 0x18, 0x7E, 0x00];
        // '2'
        glyphs[0x32 - 0x20] = [0x3C, 0x66, 0x06, 0x0C, 0x18, 0x30, 0x7E, 0x00];
        // '3'
        glyphs[0x33 - 0x20] = [0x3C, 0x66, 0x06, 0x1C, 0x06, 0x66, 0x3C, 0x00];
        // '4'
        glyphs[0x34 - 0x20] = [0x0C, 0x1C, 0x3C, 0x6C, 0x7E, 0x0C, 0x0C, 0x00];
        // '5'
        glyphs[0x35 - 0x20] = [0x7E, 0x60, 0x7C, 0x06, 0x06, 0x66, 0x3C, 0x00];
        // '6'
        glyphs[0x36 - 0x20] = [0x1C, 0x30, 0x60, 0x7C, 0x66, 0x66, 0x3C, 0x00];
        // '7'
        glyphs[0x37 - 0x20] = [0x7E, 0x06, 0x0C, 0x18, 0x30, 0x30, 0x30, 0x00];
        // '8'
        glyphs[0x38 - 0x20] = [0x3C, 0x66, 0x66, 0x3C, 0x66, 0x66, 0x3C, 0x00];
        // '9'
        glyphs[0x39 - 0x20] = [0x3C, 0x66, 0x66, 0x3E, 0x06, 0x0C, 0x38, 0x00];

        // '.' (0x2E)
        glyphs[0x2E - 0x20] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x00];
        // ',' (0x2C)
        glyphs[0x2C - 0x20] = [0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x18, 0x30];
        // '-' (0x2D)
        glyphs[0x2D - 0x20] = [0x00, 0x00, 0x00, 0x7E, 0x00, 0x00, 0x00, 0x00];
        // '+' (0x2B)
        glyphs[0x2B - 0x20] = [0x00, 0x18, 0x18, 0x7E, 0x18, 0x18, 0x00, 0x00];
        // 'e' / 'E' (0x65 / 0x45)
        glyphs[0x65 - 0x20] = [0x00, 0x00, 0x3C, 0x66, 0x7E, 0x60, 0x3C, 0x00];
        glyphs[0x45 - 0x20] = [0x7E, 0x60, 0x60, 0x7C, 0x60, 0x60, 0x7E, 0x00];

        // Fill remaining slots with a generic block glyph so all ASCII is renderable
        let block: [u8; 8] = [0x7E, 0x42, 0x42, 0x42, 0x42, 0x42, 0x7E, 0x00];
        for g in glyphs.iter_mut() {
            if *g == [0u8; 8] {
                *g = block;
            }
        }
        // Restore space
        glyphs[0] = [0x00; 8];

        Self {
            glyphs,
            metrics: FontMetrics::monospace_8x8(),
        }
    }

    /// Return the 8-byte bitmap for a character, or a replacement block if out of range.
    pub fn glyph(&self, ch: char) -> [u8; 8] {
        let code = ch as u32;
        if (0x20..=0x7E).contains(&code) {
            self.glyphs[(code - 0x20) as usize]
        } else {
            [0xFF, 0x81, 0x81, 0x81, 0x81, 0x81, 0xFF, 0x00]
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// text_width / text_height
// ──────────────────────────────────────────────────────────────────────────────

/// Compute the pixel width of a string using the given font metrics.
///
/// Width = `len(s) * (char_width + char_spacing) − char_spacing` (the last
/// character has no trailing spacing).  Returns 0 for an empty string.
pub fn text_width(s: &str, metrics: &FontMetrics) -> u32 {
    let n = s.chars().count() as u32;
    if n == 0 {
        return 0;
    }
    n * (metrics.char_width + metrics.char_spacing) - metrics.char_spacing
}

/// Compute the pixel height of a (possibly multi-line) string.
///
/// Height = `num_lines * line_height − (line_height − char_height)` so that
/// the last line does not have extra leading below it.
pub fn text_height(s: &str, metrics: &FontMetrics) -> u32 {
    let lines = s.lines().count().max(1) as u32;
    // Total height: (lines - 1) full line_heights + one char_height
    (lines - 1) * metrics.line_height + metrics.char_height
}

// ──────────────────────────────────────────────────────────────────────────────
// render_text
// ──────────────────────────────────────────────────────────────────────────────

/// A simple RGBA pixel buffer produced by the text rasterizer.
#[derive(Debug, Clone)]
pub struct PixelBuffer {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Pixel data in row-major order: each pixel is `[r, g, b, a]` with values in \[0, 255\].
    pub pixels: Vec<[u8; 4]>,
}

impl PixelBuffer {
    /// Construct an empty (transparent black) buffer.
    pub fn new(width: u32, height: u32) -> Self {
        let n = (width * height) as usize;
        Self {
            width,
            height,
            pixels: vec![[0, 0, 0, 0]; n],
        }
    }

    /// Set pixel at `(x, y)` to `color`.  Out-of-bounds writes are silently ignored.
    pub fn set_pixel(&mut self, x: u32, y: u32, color: [u8; 4]) {
        if x < self.width && y < self.height {
            self.pixels[(y * self.width + x) as usize] = color;
        }
    }

    /// Get pixel at `(x, y)`.  Returns transparent black for out-of-bounds.
    pub fn get_pixel(&self, x: u32, y: u32) -> [u8; 4] {
        if x < self.width && y < self.height {
            self.pixels[(y * self.width + x) as usize]
        } else {
            [0, 0, 0, 0]
        }
    }
}

/// Rasterize a string into a pixel buffer at position `(ox, oy)`.
///
/// Only printable ASCII glyphs are rendered.  The buffer must have been
/// pre-allocated by the caller; rendering clips to the buffer bounds.
///
/// # Parameters
/// - `buf` — destination pixel buffer.
/// - `text` — string to render.
/// - `font` — bitmap font to use.
/// - `ox`, `oy` — pixel offset of the top-left corner of the first character.
/// - `fg` — foreground colour `[r, g, b, a]`.
pub fn render_text(
    buf: &mut PixelBuffer,
    text: &str,
    font: &BitmapFont,
    ox: u32,
    oy: u32,
    fg: [u8; 4],
) {
    let cw = font.metrics.char_width;
    let lh = font.metrics.line_height;
    let sp = font.metrics.char_spacing;
    let mut cx = ox;
    let mut cy = oy;

    for ch in text.chars() {
        if ch == '\n' {
            cx = ox;
            cy += lh;
            continue;
        }
        let bitmap = font.glyph(ch);
        for row in 0..8u32 {
            let byte = bitmap[row as usize];
            for col in 0..8u32 {
                let bit = (byte >> (7 - col)) & 1;
                if bit != 0 {
                    buf.set_pixel(cx + col, cy + row, fg);
                }
            }
        }
        cx += cw + sp;
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// label_3d_to_2d
// ──────────────────────────────────────────────────────────────────────────────

/// A simple perspective projection matrix stored as a flat 4×4 row-major array.
pub type ProjMatrix = [[f64; 4]; 4];

/// Project a 3-D world position to 2-D screen coordinates using a
/// perspective (or orthographic) projection matrix and viewport dimensions.
///
/// Returns `None` if the point is behind the camera (negative w or z).
///
/// # Parameters
/// - `world_pos` — 3-D position `[x, y, z]`.
/// - `mvp` — model-view-projection matrix (row-major, 4×4).
/// - `viewport_w`, `viewport_h` — screen dimensions in pixels.
///
/// Returns `(screen_x, screen_y)` in pixel coordinates with origin at
/// the top-left.
pub fn label_3d_to_2d(
    world_pos: [f64; 3],
    mvp: &ProjMatrix,
    viewport_w: u32,
    viewport_h: u32,
) -> Option<(f64, f64)> {
    let [x, y, z] = world_pos;

    // Clip-space position: c = MVP * [x, y, z, 1]ᵀ
    let cx = mvp[0][0] * x + mvp[0][1] * y + mvp[0][2] * z + mvp[0][3];
    let cy = mvp[1][0] * x + mvp[1][1] * y + mvp[1][2] * z + mvp[1][3];
    let _cz = mvp[2][0] * x + mvp[2][1] * y + mvp[2][2] * z + mvp[2][3];
    let cw = mvp[3][0] * x + mvp[3][1] * y + mvp[3][2] * z + mvp[3][3];

    if cw <= 0.0 {
        return None; // behind the camera
    }

    // NDC: [-1, 1] × [-1, 1]
    let ndx = cx / cw;
    let ndy = cy / cw;

    // Screen coordinates (y flipped: +y = down)
    let sx = (ndx + 1.0) * 0.5 * viewport_w as f64;
    let sy = (1.0 - ndy) * 0.5 * viewport_h as f64;

    Some((sx, sy))
}

// ──────────────────────────────────────────────────────────────────────────────
// draw_axis_labels
// ──────────────────────────────────────────────────────────────────────────────

/// Axis identifier for [`draw_axis_labels`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// X axis.
    X,
    /// Y axis.
    Y,
    /// Z axis.
    Z,
}

/// A text label placed at a screen position.
#[derive(Debug, Clone)]
pub struct ScreenLabel {
    /// Label text.
    pub text: String,
    /// Screen x coordinate (pixels, origin top-left).
    pub x: f64,
    /// Screen y coordinate (pixels, origin top-left).
    pub y: f64,
}

/// Generate tick-mark labels for one axis.
///
/// For each value in `ticks`, the world position along `axis` at that value
/// (other coordinates are `axis_offset`) is projected to screen coordinates.
/// Labels that project outside the viewport are omitted.
///
/// # Parameters
/// - `axis` — which axis to label.
/// - `ticks` — tick values along the axis.
/// - `axis_offset` — value of the other two coordinates (e.g. 0.0 for the origin).
/// - `mvp` — model-view-projection matrix.
/// - `viewport_w`, `viewport_h` — screen dimensions.
/// - `precision` — decimal places for each label.
pub fn draw_axis_labels(
    axis: Axis,
    ticks: &[f64],
    axis_offset: f64,
    mvp: &ProjMatrix,
    viewport_w: u32,
    viewport_h: u32,
    precision: usize,
) -> Vec<ScreenLabel> {
    let o = axis_offset;
    ticks
        .iter()
        .filter_map(|&v| {
            let world = match axis {
                Axis::X => [v, o, o],
                Axis::Y => [o, v, o],
                Axis::Z => [o, o, v],
            };
            let (sx, sy) = label_3d_to_2d(world, mvp, viewport_w, viewport_h)?;
            // Clip check
            if sx < 0.0 || sx > viewport_w as f64 || sy < 0.0 || sy > viewport_h as f64 {
                return None;
            }
            let text = format!("{:.prec$}", v, prec = precision);
            Some(ScreenLabel { text, x: sx, y: sy })
        })
        .collect()
}

// ──────────────────────────────────────────────────────────────────────────────
// number_format
// ──────────────────────────────────────────────────────────────────────────────

/// Notation style for [`number_format`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NumberNotation {
    /// Fixed-point: e.g. `3.14159`.
    Fixed,
    /// Scientific: e.g. `3.14e0`.
    Scientific,
    /// Automatic: use fixed when `|v| ∈ [0.001, 9999]`, otherwise scientific.
    Auto,
}

/// Format a floating-point number to a string with the given notation and precision.
///
/// # Parameters
/// - `value` — value to format.
/// - `precision` — number of digits after the decimal point (fixed) or significant figures (scientific).
/// - `notation` — [`NumberNotation`] variant.
///
/// # Examples
/// ```no_run
/// use oxiphysics_viz::text_renderer::{number_format, NumberNotation};
/// let s = number_format(3.14159, 2, NumberNotation::Fixed);
/// assert_eq!(s, "3.14");
/// ```
pub fn number_format(value: f64, precision: usize, notation: NumberNotation) -> String {
    match notation {
        NumberNotation::Fixed => format!("{:.prec$}", value, prec = precision),
        NumberNotation::Scientific => format!("{:.prec$e}", value, prec = precision),
        NumberNotation::Auto => {
            let abs = value.abs();
            if abs == 0.0 || (0.001..10_000.0).contains(&abs) {
                format!("{:.prec$}", value, prec = precision)
            } else {
                format!("{:.prec$e}", value, prec = precision)
            }
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FontMetrics ──────────────────────────────────────────────────────────

    #[test]
    fn test_font_metrics_default_8x8() {
        let m = FontMetrics::default();
        assert_eq!(m.char_width, 8);
        assert_eq!(m.char_height, 8);
    }

    #[test]
    fn test_font_metrics_monospace_constructor() {
        let m = FontMetrics::monospace_8x8();
        assert_eq!(m.char_width, 8);
        assert_eq!(m.char_height, 8);
        assert_eq!(m.line_height, 10);
    }

    // ── BitmapFont ───────────────────────────────────────────────────────────

    #[test]
    fn test_bitmap_font_space_glyph_is_empty() {
        let font = BitmapFont::new();
        let g = font.glyph(' ');
        assert_eq!(g, [0u8; 8]);
    }

    #[test]
    fn test_bitmap_font_digits_non_zero() {
        let font = BitmapFont::new();
        for ch in '0'..='9' {
            let g = font.glyph(ch);
            let sum: u32 = g.iter().map(|&b| b as u32).sum();
            assert!(sum > 0, "digit '{}' glyph should be non-zero", ch);
        }
    }

    #[test]
    fn test_bitmap_font_out_of_range_returns_block() {
        let font = BitmapFont::new();
        let g = font.glyph('\x01'); // control character
        // Replacement block starts with 0xFF
        assert_eq!(g[0], 0xFF);
    }

    #[test]
    fn test_bitmap_font_glyphs_count() {
        let font = BitmapFont::new();
        assert_eq!(font.glyphs.len(), 95); // 0x20..=0x7E
    }

    #[test]
    fn test_bitmap_font_printable_range() {
        let font = BitmapFont::new();
        // All printable ASCII should return a glyph without panic
        for code in 0x20_u8..=0x7E_u8 {
            let _ = font.glyph(code as char);
        }
    }

    // ── text_width / text_height ─────────────────────────────────────────────

    #[test]
    fn test_text_width_empty_string() {
        let m = FontMetrics::default();
        assert_eq!(text_width("", &m), 0);
    }

    #[test]
    fn test_text_width_single_char() {
        let m = FontMetrics::default();
        // Single char: no trailing spacing
        assert_eq!(text_width("A", &m), m.char_width);
    }

    #[test]
    fn test_text_width_two_chars() {
        let m = FontMetrics::default();
        // Two chars: char_width + spacing + char_width
        assert_eq!(text_width("AB", &m), 2 * m.char_width + m.char_spacing);
    }

    #[test]
    fn test_text_width_ten_chars() {
        let m = FontMetrics::default();
        let s = "0123456789"; // 10 chars
        let expected = 10 * (m.char_width + m.char_spacing) - m.char_spacing;
        assert_eq!(text_width(s, &m), expected);
    }

    #[test]
    fn test_text_height_single_line() {
        let m = FontMetrics::default();
        assert_eq!(text_height("Hello", &m), m.char_height);
    }

    #[test]
    fn test_text_height_two_lines() {
        let m = FontMetrics::default();
        let h = text_height("Line1\nLine2", &m);
        let expected = m.line_height + m.char_height;
        assert_eq!(h, expected);
    }

    #[test]
    fn test_text_height_empty_string() {
        let m = FontMetrics::default();
        // Empty string counts as 1 line
        assert_eq!(text_height("", &m), m.char_height);
    }

    // ── PixelBuffer ──────────────────────────────────────────────────────────

    #[test]
    fn test_pixel_buffer_new_transparent() {
        let buf = PixelBuffer::new(4, 4);
        assert_eq!(buf.pixels.len(), 16);
        for px in &buf.pixels {
            assert_eq!(*px, [0, 0, 0, 0]);
        }
    }

    #[test]
    fn test_pixel_buffer_set_get() {
        let mut buf = PixelBuffer::new(10, 10);
        buf.set_pixel(3, 4, [255, 0, 0, 255]);
        assert_eq!(buf.get_pixel(3, 4), [255, 0, 0, 255]);
    }

    #[test]
    fn test_pixel_buffer_out_of_bounds_write_ignored() {
        let mut buf = PixelBuffer::new(4, 4);
        buf.set_pixel(100, 100, [255, 255, 255, 255]); // should not panic
        assert_eq!(buf.get_pixel(100, 100), [0, 0, 0, 0]);
    }

    // ── render_text ──────────────────────────────────────────────────────────

    #[test]
    fn test_render_text_marks_pixels() {
        let font = BitmapFont::new();
        let mut buf = PixelBuffer::new(80, 10);
        render_text(&mut buf, "1", &font, 0, 0, [255, 255, 255, 255]);
        // At least one pixel should be set
        let set: usize = buf.pixels.iter().filter(|&&p| p[3] > 0).count();
        assert!(set > 0, "rendering '1' should set at least one pixel");
    }

    #[test]
    fn test_render_text_empty_string_no_pixels() {
        let font = BitmapFont::new();
        let mut buf = PixelBuffer::new(80, 10);
        render_text(&mut buf, "", &font, 0, 0, [255, 255, 255, 255]);
        let set: usize = buf.pixels.iter().filter(|&&p| p[3] > 0).count();
        assert_eq!(set, 0);
    }

    #[test]
    fn test_render_text_clips_to_buffer() {
        let font = BitmapFont::new();
        let mut buf = PixelBuffer::new(2, 2); // tiny buffer
        // Should not panic even though text overflows
        render_text(&mut buf, "HELLO WORLD", &font, 0, 0, [255, 255, 255, 255]);
    }

    // ── label_3d_to_2d ───────────────────────────────────────────────────────

    #[test]
    fn test_label_3d_to_2d_identity_at_origin() {
        // Identity MVP: point at (0,0,0) → NDC (0,0), screen (w/2, h/2)
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let result = label_3d_to_2d([0.0, 0.0, 0.0], &mvp, 800, 600);
        let (sx, sy) = result.unwrap();
        assert!((sx - 400.0).abs() < 1e-6, "sx={sx}");
        assert!((sy - 300.0).abs() < 1e-6, "sy={sy}");
    }

    #[test]
    fn test_label_3d_to_2d_behind_camera_returns_none() {
        // w = 0 → behind camera
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0], // w = 0
        ];
        let result = label_3d_to_2d([1.0, 0.0, 0.0], &mvp, 800, 600);
        assert!(result.is_none());
    }

    #[test]
    fn test_label_3d_to_2d_positive_x_maps_right() {
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let (sx, _) = label_3d_to_2d([0.5, 0.0, 0.0], &mvp, 800, 600).unwrap();
        // NDC x = 0.5 → screen x = (0.5+1)/2 * 800 = 600
        assert!((sx - 600.0).abs() < 1e-6);
    }

    // ── draw_axis_labels ─────────────────────────────────────────────────────

    #[test]
    fn test_draw_axis_labels_x_generates_labels() {
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let ticks = [-0.5, 0.0, 0.5];
        let labels = draw_axis_labels(Axis::X, &ticks, 0.0, &mvp, 800, 600, 2);
        // All three ticks are in NDC [-0.5, 0.5] → on screen
        assert_eq!(labels.len(), 3);
    }

    #[test]
    fn test_draw_axis_labels_text_format() {
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let labels = draw_axis_labels(Axis::Y, &[0.0], 0.0, &mvp, 800, 600, 1);
        assert!(!labels.is_empty());
        assert_eq!(labels[0].text, "0.0");
    }

    #[test]
    fn test_draw_axis_labels_offscreen_excluded() {
        let mvp: ProjMatrix = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        // x=5.0 → NDC x=5.0 → screen x = 3.0 * 800 = 2400 → off screen
        let labels = draw_axis_labels(Axis::X, &[5.0], 0.0, &mvp, 800, 600, 2);
        assert!(labels.is_empty());
    }

    // ── number_format ─────────────────────────────────────────────────────────

    #[test]
    fn test_number_format_fixed_two_decimal() {
        let s = number_format(std::f64::consts::PI, 2, NumberNotation::Fixed);
        assert_eq!(s, "3.14");
    }

    #[test]
    fn test_number_format_fixed_zero_decimal() {
        let s = number_format(42.7, 0, NumberNotation::Fixed);
        assert_eq!(s, "43");
    }

    #[test]
    fn test_number_format_scientific() {
        let s = number_format(1234.5, 2, NumberNotation::Scientific);
        assert!(s.contains('e') || s.contains('E'));
    }

    #[test]
    fn test_number_format_auto_small_uses_fixed() {
        let s = number_format(3.125, 3, NumberNotation::Auto);
        assert!(!s.contains('e'));
    }

    #[test]
    fn test_number_format_auto_large_uses_scientific() {
        let s = number_format(1e8, 2, NumberNotation::Auto);
        assert!(s.contains('e') || s.contains('E'));
    }

    #[test]
    fn test_number_format_negative() {
        let s = number_format(-0.001, 5, NumberNotation::Fixed);
        assert!(s.starts_with('-'));
    }

    #[test]
    fn test_number_format_zero_auto() {
        let s = number_format(0.0, 2, NumberNotation::Auto);
        assert_eq!(s, "0.00");
    }
}
