// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Basic font and text rendering data structures.
//!
//! Provides CPU-side representations for bitmap fonts, SDF (signed distance
//! field) fonts, text layout, and 3-D annotation projection.  No GPU or
//! windowing dependency is required.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// GlyphBitmap
// ─────────────────────────────────────────────────────────────────────────────

/// A single character's rasterised glyph as a 8-bit alpha bitmap.
#[derive(Debug, Clone)]
pub struct GlyphBitmap {
    /// Unicode code point this glyph represents.
    pub char_code: char,
    /// Pixel width of the bitmap.
    pub width: usize,
    /// Pixel height of the bitmap.
    pub height: usize,
    /// Raw alpha values, row-major, length `width * height`.
    pub pixels: Vec<u8>,
    /// Horizontal advance (in pixels) after rendering this glyph.
    pub advance_x: f64,
}

impl GlyphBitmap {
    /// Construct a new `GlyphBitmap` with the given dimensions.
    ///
    /// All pixels are initialised to `0` (transparent).
    pub fn new(char_code: char, width: usize, height: usize, advance_x: f64) -> Self {
        Self {
            char_code,
            width,
            height,
            pixels: vec![0u8; width * height],
            advance_x,
        }
    }

    /// Set a single pixel at `(x, y)` (column, row), clamping out-of-bounds
    /// accesses silently.
    pub fn set_pixel(&mut self, x: usize, y: usize, alpha: u8) {
        if x < self.width && y < self.height {
            self.pixels[y * self.width + x] = alpha;
        }
    }

    /// Get the alpha value at `(x, y)`, returning `0` for out-of-bounds.
    pub fn get_pixel(&self, x: usize, y: usize) -> u8 {
        if x < self.width && y < self.height {
            self.pixels[y * self.width + x]
        } else {
            0
        }
    }

    /// Total number of pixels (width × height).
    pub fn pixel_count(&self) -> usize {
        self.width * self.height
    }

    /// Return `true` if the glyph has no pixels (both dimensions are zero).
    pub fn is_empty(&self) -> bool {
        self.width == 0 || self.height == 0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BitmapFont
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of [`GlyphBitmap`]s forming a complete bitmap font.
#[derive(Debug, Clone)]
pub struct BitmapFont {
    /// Map from character to glyph.
    pub glyphs: HashMap<char, GlyphBitmap>,
    /// Font size in pixels (e.g. 12, 16, 24).
    pub size_px: u32,
    /// Distance between consecutive baselines, in pixels.
    pub line_height: f64,
}

impl BitmapFont {
    /// Construct an empty `BitmapFont`.
    pub fn new(size_px: u32, line_height: f64) -> Self {
        Self {
            glyphs: HashMap::new(),
            size_px,
            line_height,
        }
    }

    /// Add a glyph to the font, replacing any existing entry for the same
    /// character.
    pub fn add_glyph(&mut self, glyph: GlyphBitmap) {
        self.glyphs.insert(glyph.char_code, glyph);
    }

    /// Look up a glyph by character.  Returns `None` if not present.
    pub fn get_glyph(&self, c: char) -> Option<&GlyphBitmap> {
        self.glyphs.get(&c)
    }

    /// Return the horizontal advance for `c`, or `size_px as f64 * 0.6` as a
    /// fallback when the glyph is missing.
    pub fn advance_for(&self, c: char) -> f64 {
        self.glyphs
            .get(&c)
            .map(|g| g.advance_x)
            .unwrap_or(self.size_px as f64 * 0.6)
    }

    /// Measure the pixel width of a string using glyph advances.
    pub fn measure_width(&self, text: &str) -> f64 {
        text.chars().map(|c| self.advance_for(c)).sum()
    }

    /// Number of glyphs in the font.
    pub fn glyph_count(&self) -> usize {
        self.glyphs.len()
    }

    /// Build a simple synthetic bitmap font populated with printable ASCII
    /// characters (0x20–0x7E).  Each glyph is a blank rectangle of
    /// `size_px × size_px` pixels with a plausible advance width.
    pub fn ascii_placeholder(size_px: u32) -> Self {
        let mut font = Self::new(size_px, size_px as f64 * 1.2);
        let w = size_px as usize;
        let h = size_px as usize;
        for code in 0x20u32..=0x7Eu32 {
            if let Some(c) = char::from_u32(code) {
                let adv = if c == ' ' {
                    size_px as f64 * 0.3
                } else {
                    size_px as f64 * 0.6
                };
                let mut glyph = GlyphBitmap::new(c, w, h, adv);
                // Fill a simple vertical bar pattern to give non-zero coverage
                if c != ' ' && w > 2 && h > 2 {
                    for row in 1..(h - 1) {
                        glyph.set_pixel(1, row, 255);
                        glyph.set_pixel(w - 2, row, 255);
                    }
                }
                font.add_glyph(glyph);
            }
        }
        font
    }

    /// Build a real 8×8 bitmap font for printable ASCII (0x20–0x7E) based on
    /// the public-domain IBM PC CGA/EGA 8×8 ROM font bit patterns.
    ///
    /// Each character is an 8-row × 8-column bitmap.  One byte per row; the
    /// most-significant bit is the left-most pixel.  Characters with no
    /// set bits (e.g. SPACE) are intentionally blank.
    pub fn ascii_default() -> Self {
        let mut font = Self::new(8, 10.0);
        for (idx, rows) in ASCII_FONT_8X8.chunks_exact(8).enumerate() {
            let code = 0x20u32 + idx as u32;
            if let Some(c) = char::from_u32(code) {
                let adv = if c == ' ' { 4.0 } else { 8.0 };
                let mut glyph = GlyphBitmap::new(c, 8, 8, adv);
                for (row, &byte) in rows.iter().enumerate() {
                    for col in 0..8usize {
                        if byte & (0x80 >> col) != 0 {
                            glyph.set_pixel(col, row, 255);
                        }
                    }
                }
                font.add_glyph(glyph);
            }
        }
        font
    }

    /// Return the packed 8-row bit pattern for character `c`, or `None` if
    /// `c` is not in the printable ASCII range (0x20–0x7E).
    ///
    /// Each returned byte represents one row of the 8×8 glyph; bit 7 is the
    /// left-most pixel.
    pub fn glyph_bits(c: char) -> Option<[u8; 8]> {
        let code = c as u32;
        if !(0x20..=0x7E).contains(&code) {
            return None;
        }
        let idx = (code - 0x20) as usize;
        let start = idx * 8;
        let slice = &ASCII_FONT_8X8[start..start + 8];
        let mut out = [0u8; 8];
        out.copy_from_slice(slice);
        Some(out)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public-domain 8×8 IBM PC CGA/EGA ROM bitmap font data.
//
// 95 printable ASCII characters, 0x20 (SPACE) through 0x7E (~).
// Each character is encoded as 8 bytes — one byte per scanline row, MSB = left.
// Source: public domain IBM PC character ROM data widely reproduced in OSS.
// ─────────────────────────────────────────────────────────────────────────────
#[rustfmt::skip]
const ASCII_FONT_8X8: [u8; 95 * 8] = [
    // 0x20 SPACE
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    // 0x21 !
    0x18, 0x3C, 0x3C, 0x18, 0x18, 0x00, 0x18, 0x00,
    // 0x22 "
    0x36, 0x36, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    // 0x23 #
    0x36, 0x36, 0x7F, 0x36, 0x7F, 0x36, 0x36, 0x00,
    // 0x24 $
    0x0C, 0x3E, 0x03, 0x1E, 0x30, 0x1F, 0x0C, 0x00,
    // 0x25 %
    0x00, 0x63, 0x33, 0x18, 0x0C, 0x66, 0x63, 0x00,
    // 0x26 &
    0x1C, 0x36, 0x1C, 0x6E, 0x3B, 0x33, 0x6E, 0x00,
    // 0x27 '
    0x06, 0x06, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    // 0x28 (
    0x18, 0x0C, 0x06, 0x06, 0x06, 0x0C, 0x18, 0x00,
    // 0x29 )
    0x06, 0x0C, 0x18, 0x18, 0x18, 0x0C, 0x06, 0x00,
    // 0x2A *
    0x00, 0x66, 0x3C, 0xFF, 0x3C, 0x66, 0x00, 0x00,
    // 0x2B +
    0x00, 0x0C, 0x0C, 0x3F, 0x0C, 0x0C, 0x00, 0x00,
    // 0x2C ,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C, 0x06,
    // 0x2D -
    0x00, 0x00, 0x00, 0x3F, 0x00, 0x00, 0x00, 0x00,
    // 0x2E .
    0x00, 0x00, 0x00, 0x00, 0x00, 0x0C, 0x0C, 0x00,
    // 0x2F /
    0x60, 0x30, 0x18, 0x0C, 0x06, 0x03, 0x01, 0x00,
    // 0x30 0
    0x3E, 0x63, 0x73, 0x7B, 0x6F, 0x67, 0x3E, 0x00,
    // 0x31 1
    0x0C, 0x0E, 0x0C, 0x0C, 0x0C, 0x0C, 0x3F, 0x00,
    // 0x32 2
    0x1E, 0x33, 0x30, 0x1C, 0x06, 0x33, 0x3F, 0x00,
    // 0x33 3
    0x1E, 0x33, 0x30, 0x1C, 0x30, 0x33, 0x1E, 0x00,
    // 0x34 4
    0x38, 0x3C, 0x36, 0x33, 0x7F, 0x30, 0x78, 0x00,
    // 0x35 5
    0x3F, 0x03, 0x1F, 0x30, 0x30, 0x33, 0x1E, 0x00,
    // 0x36 6
    0x1C, 0x06, 0x03, 0x1F, 0x33, 0x33, 0x1E, 0x00,
    // 0x37 7
    0x3F, 0x33, 0x30, 0x18, 0x0C, 0x0C, 0x0C, 0x00,
    // 0x38 8
    0x1E, 0x33, 0x33, 0x1E, 0x33, 0x33, 0x1E, 0x00,
    // 0x39 9
    0x1E, 0x33, 0x33, 0x3E, 0x30, 0x18, 0x0E, 0x00,
    // 0x3A :
    0x00, 0x0C, 0x0C, 0x00, 0x00, 0x0C, 0x0C, 0x00,
    // 0x3B ;
    0x00, 0x0C, 0x0C, 0x00, 0x00, 0x0C, 0x0C, 0x06,
    // 0x3C <
    0x18, 0x0C, 0x06, 0x03, 0x06, 0x0C, 0x18, 0x00,
    // 0x3D =
    0x00, 0x00, 0x3F, 0x00, 0x00, 0x3F, 0x00, 0x00,
    // 0x3E >
    0x06, 0x0C, 0x18, 0x30, 0x18, 0x0C, 0x06, 0x00,
    // 0x3F ?
    0x1E, 0x33, 0x30, 0x18, 0x0C, 0x00, 0x0C, 0x00,
    // 0x40 @
    0x3E, 0x63, 0x7B, 0x7B, 0x7B, 0x03, 0x1E, 0x00,
    // 0x41 A
    0x0C, 0x1E, 0x33, 0x33, 0x3F, 0x33, 0x33, 0x00,
    // 0x42 B
    0x3F, 0x66, 0x66, 0x3E, 0x66, 0x66, 0x3F, 0x00,
    // 0x43 C
    0x3C, 0x66, 0x03, 0x03, 0x03, 0x66, 0x3C, 0x00,
    // 0x44 D
    0x1F, 0x36, 0x66, 0x66, 0x66, 0x36, 0x1F, 0x00,
    // 0x45 E
    0x7F, 0x46, 0x16, 0x1E, 0x16, 0x46, 0x7F, 0x00,
    // 0x46 F
    0x7F, 0x46, 0x16, 0x1E, 0x16, 0x06, 0x0F, 0x00,
    // 0x47 G
    0x3C, 0x66, 0x03, 0x03, 0x73, 0x66, 0x7C, 0x00,
    // 0x48 H
    0x33, 0x33, 0x33, 0x3F, 0x33, 0x33, 0x33, 0x00,
    // 0x49 I
    0x1E, 0x0C, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00,
    // 0x4A J
    0x78, 0x30, 0x30, 0x30, 0x33, 0x33, 0x1E, 0x00,
    // 0x4B K
    0x67, 0x66, 0x36, 0x1E, 0x36, 0x66, 0x67, 0x00,
    // 0x4C L
    0x0F, 0x06, 0x06, 0x06, 0x46, 0x66, 0x7F, 0x00,
    // 0x4D M
    0x63, 0x77, 0x7F, 0x7F, 0x6B, 0x63, 0x63, 0x00,
    // 0x4E N
    0x63, 0x67, 0x6F, 0x7B, 0x73, 0x63, 0x63, 0x00,
    // 0x4F O
    0x1C, 0x36, 0x63, 0x63, 0x63, 0x36, 0x1C, 0x00,
    // 0x50 P
    0x3F, 0x66, 0x66, 0x3E, 0x06, 0x06, 0x0F, 0x00,
    // 0x51 Q
    0x1E, 0x33, 0x33, 0x33, 0x3B, 0x1E, 0x38, 0x00,
    // 0x52 R
    0x3F, 0x66, 0x66, 0x3E, 0x36, 0x66, 0x67, 0x00,
    // 0x53 S
    0x1E, 0x33, 0x07, 0x0E, 0x38, 0x33, 0x1E, 0x00,
    // 0x54 T
    0x3F, 0x2D, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00,
    // 0x55 U
    0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x3F, 0x00,
    // 0x56 V
    0x33, 0x33, 0x33, 0x33, 0x33, 0x1E, 0x0C, 0x00,
    // 0x57 W
    0x63, 0x63, 0x63, 0x6B, 0x7F, 0x77, 0x63, 0x00,
    // 0x58 X
    0x63, 0x63, 0x36, 0x1C, 0x1C, 0x36, 0x63, 0x00,
    // 0x59 Y
    0x33, 0x33, 0x33, 0x1E, 0x0C, 0x0C, 0x1E, 0x00,
    // 0x5A Z
    0x7F, 0x63, 0x31, 0x18, 0x4C, 0x66, 0x7F, 0x00,
    // 0x5B [
    0x1E, 0x06, 0x06, 0x06, 0x06, 0x06, 0x1E, 0x00,
    // 0x5C backslash
    0x03, 0x06, 0x0C, 0x18, 0x30, 0x60, 0x40, 0x00,
    // 0x5D ]
    0x1E, 0x18, 0x18, 0x18, 0x18, 0x18, 0x1E, 0x00,
    // 0x5E ^
    0x08, 0x1C, 0x36, 0x63, 0x00, 0x00, 0x00, 0x00,
    // 0x5F _
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF,
    // 0x60 `
    0x0C, 0x0C, 0x18, 0x00, 0x00, 0x00, 0x00, 0x00,
    // 0x61 a
    0x00, 0x00, 0x1E, 0x30, 0x3E, 0x33, 0x6E, 0x00,
    // 0x62 b
    0x07, 0x06, 0x06, 0x3E, 0x66, 0x66, 0x3B, 0x00,
    // 0x63 c
    0x00, 0x00, 0x1E, 0x33, 0x03, 0x33, 0x1E, 0x00,
    // 0x64 d
    0x38, 0x30, 0x30, 0x3E, 0x33, 0x33, 0x6E, 0x00,
    // 0x65 e
    0x00, 0x00, 0x1E, 0x33, 0x3F, 0x03, 0x1E, 0x00,
    // 0x66 f
    0x1C, 0x36, 0x06, 0x0F, 0x06, 0x06, 0x0F, 0x00,
    // 0x67 g
    0x00, 0x00, 0x6E, 0x33, 0x33, 0x3E, 0x30, 0x1F,
    // 0x68 h
    0x07, 0x06, 0x36, 0x6E, 0x66, 0x66, 0x67, 0x00,
    // 0x69 i
    0x0C, 0x00, 0x0E, 0x0C, 0x0C, 0x0C, 0x1E, 0x00,
    // 0x6A j
    0x30, 0x00, 0x30, 0x30, 0x30, 0x33, 0x33, 0x1E,
    // 0x6B k
    0x07, 0x06, 0x66, 0x36, 0x1E, 0x36, 0x67, 0x00,
    // 0x6C l
    0x0E, 0x0C, 0x0C, 0x0C, 0x0C, 0x0C, 0x1E, 0x00,
    // 0x6D m
    0x00, 0x00, 0x33, 0x7F, 0x7F, 0x6B, 0x63, 0x00,
    // 0x6E n
    0x00, 0x00, 0x1F, 0x33, 0x33, 0x33, 0x33, 0x00,
    // 0x6F o
    0x00, 0x00, 0x1E, 0x33, 0x33, 0x33, 0x1E, 0x00,
    // 0x70 p
    0x00, 0x00, 0x3B, 0x66, 0x66, 0x3E, 0x06, 0x0F,
    // 0x71 q
    0x00, 0x00, 0x6E, 0x33, 0x33, 0x3E, 0x30, 0x78,
    // 0x72 r
    0x00, 0x00, 0x3B, 0x6E, 0x66, 0x06, 0x0F, 0x00,
    // 0x73 s
    0x00, 0x00, 0x3E, 0x03, 0x1E, 0x30, 0x1F, 0x00,
    // 0x74 t
    0x08, 0x0C, 0x3E, 0x0C, 0x0C, 0x2C, 0x18, 0x00,
    // 0x75 u
    0x00, 0x00, 0x33, 0x33, 0x33, 0x33, 0x6E, 0x00,
    // 0x76 v
    0x00, 0x00, 0x33, 0x33, 0x33, 0x1E, 0x0C, 0x00,
    // 0x77 w
    0x00, 0x00, 0x63, 0x6B, 0x7F, 0x7F, 0x36, 0x00,
    // 0x78 x
    0x00, 0x00, 0x63, 0x36, 0x1C, 0x36, 0x63, 0x00,
    // 0x79 y
    0x00, 0x00, 0x33, 0x33, 0x33, 0x3E, 0x30, 0x1F,
    // 0x7A z
    0x00, 0x00, 0x3F, 0x19, 0x0C, 0x26, 0x3F, 0x00,
    // 0x7B {
    0x38, 0x0C, 0x0C, 0x07, 0x0C, 0x0C, 0x38, 0x00,
    // 0x7C |
    0x18, 0x18, 0x18, 0x00, 0x18, 0x18, 0x18, 0x00,
    // 0x7D }
    0x07, 0x0C, 0x0C, 0x38, 0x0C, 0x0C, 0x07, 0x00,
    // 0x7E ~
    0x6E, 0x3B, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

// ─────────────────────────────────────────────────────────────────────────────
// TextLayout
// ─────────────────────────────────────────────────────────────────────────────

/// Lays out a string of characters into 2-D positions.
///
/// Produces a list of `(x, y, char)` triples.  When `wrapping` is enabled,
/// lines are broken at word boundaries or at `line_width` pixels.
#[derive(Debug, Clone)]
pub struct TextLayout {
    /// The text to lay out.
    pub text: String,
    /// Font size in pixels.
    pub font_size: f64,
    /// Maximum line width in pixels before wrapping.
    pub line_width: f64,
    /// Whether to wrap at `line_width`.
    pub wrapping: bool,
    /// Reference font for measuring glyph advances.
    pub font: BitmapFont,
}

impl TextLayout {
    /// Construct a `TextLayout`.
    pub fn new(
        text: impl Into<String>,
        font: BitmapFont,
        font_size: f64,
        line_width: f64,
        wrapping: bool,
    ) -> Self {
        Self {
            text: text.into(),
            font_size,
            line_width,
            wrapping,
            font,
        }
    }

    /// Lay out the text and return `(x, y, char)` for every character.
    ///
    /// Origin is `(0, 0)` at the top-left.  `y` increases downward by
    /// `font.line_height` for each new line.
    pub fn layout(&self) -> Vec<(f64, f64, char)> {
        let scale = self.font_size / self.font.size_px as f64;
        let line_h = self.font.line_height * scale;

        let mut out = Vec::new();
        let mut x = 0.0_f64;
        let mut y = 0.0_f64;

        for c in self.text.chars() {
            if c == '\n' {
                x = 0.0;
                y += line_h;
                continue;
            }
            // Wrap at word boundary if needed (simple: wrap at any character)
            let adv = self.font.advance_for(c) * scale;
            if self.wrapping && x + adv > self.line_width && x > 0.0 {
                x = 0.0;
                y += line_h;
            }
            out.push((x, y, c));
            x += adv;
        }
        out
    }

    /// Return only the characters that fall within the rectangle
    /// `[0, max_x] × [0, max_y]`.
    pub fn layout_clipped(&self, max_x: f64, max_y: f64) -> Vec<(f64, f64, char)> {
        self.layout()
            .into_iter()
            .filter(|&(x, y, _)| x <= max_x && y <= max_y)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SdfFont
// ─────────────────────────────────────────────────────────────────────────────

/// A Signed Distance Field font representation.
///
/// Each glyph is stored as a flat `Vec`f32` of signed distance values
/// (positive outside, negative inside), arranged in a square of side
/// `sdf_size`.
#[derive(Debug, Clone)]
pub struct SdfFont {
    /// Map from character to SDF data.
    pub glyphs: HashMap<char, Vec<f32>>,
    /// Side length of each SDF tile in texels.
    pub sdf_size: usize,
    /// Spread / padding in pixels for the signed distance field.
    pub spread: f64,
}

impl SdfFont {
    /// Construct an empty `SdfFont`.
    pub fn new(sdf_size: usize, spread: f64) -> Self {
        Self {
            glyphs: HashMap::new(),
            sdf_size,
            spread,
        }
    }

    /// Insert a pre-computed SDF for character `c`.
    pub fn insert(&mut self, c: char, sdf: Vec<f32>) {
        self.glyphs.insert(c, sdf);
    }

    /// Render glyph `c` at the given `scale` factor by thresholding the SDF.
    ///
    /// Returns a rasterised `u8` bitmap of size `(sdf_size * scale) × (sdf_size * scale)`.
    /// Pixels where `sdf > 0` are opaque (`255`), others are transparent (`0`).
    pub fn render_glyph(&self, c: char, scale: f64) -> Vec<u8> {
        let sdf = match self.glyphs.get(&c) {
            Some(s) => s,
            None => return vec![],
        };
        let src = self.sdf_size;
        let dst = (src as f64 * scale).round() as usize;
        if dst == 0 {
            return vec![];
        }
        let mut out = vec![0u8; dst * dst];
        for ry in 0..dst {
            for rx in 0..dst {
                // Nearest-neighbour sample
                let sx = ((rx as f64 / dst as f64) * src as f64) as usize;
                let sy = ((ry as f64 / dst as f64) * src as f64) as usize;
                let sx = sx.min(src - 1);
                let sy = sy.min(src - 1);
                let val = sdf[sy * src + sx];
                out[ry * dst + rx] = if val > 0.0 { 255 } else { 0 };
            }
        }
        out
    }

    /// Generate a synthetic circular disc SDF for character `c` and store it.
    ///
    /// This creates a simple filled-circle glyph, useful for testing.
    pub fn generate_disc_glyph(&mut self, c: char) {
        let n = self.sdf_size;
        let radius = n as f32 / 2.5;
        let cx = n as f32 / 2.0;
        let cy = n as f32 / 2.0;
        let sdf: Vec<f32> = (0..n * n)
            .map(|idx| {
                let px = (idx % n) as f32;
                let py = (idx / n) as f32;
                let dist = ((px - cx).powi(2) + (py - cy).powi(2)).sqrt();
                radius - dist // positive inside
            })
            .collect();
        self.insert(c, sdf);
    }

    /// Return `true` if `c` is present in the font.
    pub fn has_glyph(&self, c: char) -> bool {
        self.glyphs.contains_key(&c)
    }

    /// Number of glyphs in the font.
    pub fn glyph_count(&self) -> usize {
        self.glyphs.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TextRenderer
// ─────────────────────────────────────────────────────────────────────────────

/// Renders a string of text to a list of `(x, y, alpha)` pixel contributions.
///
/// Each glyph is looked up in the font and its bitmap pixels are placed at the
/// appropriate 2-D position.
#[derive(Debug, Clone)]
pub struct TextRenderer {
    /// The bitmap font used for rendering.
    pub font: BitmapFont,
}

impl TextRenderer {
    /// Create a `TextRenderer` using `font`.
    pub fn new(font: BitmapFont) -> Self {
        Self { font }
    }

    /// Render `text` starting at screen position `(x, y)`.
    ///
    /// Returns a list of `(px, py, alpha)` where `(px, py)` is the pixel
    /// coordinate and `alpha` is the coverage value `\[0, 255\]`.
    pub fn render_string(&self, text: &str, x: f64, y: f64) -> Vec<(f64, f64, u8)> {
        let mut out = Vec::new();
        let mut cursor_x = x;

        for c in text.chars() {
            if c == '\n' {
                cursor_x = x;
                // Advance y handled externally; skip newline characters here
                continue;
            }
            if let Some(glyph) = self.font.get_glyph(c) {
                for gy in 0..glyph.height {
                    for gx in 0..glyph.width {
                        let alpha = glyph.get_pixel(gx, gy);
                        if alpha > 0 {
                            out.push((cursor_x + gx as f64, y + gy as f64, alpha));
                        }
                    }
                }
                cursor_x += glyph.advance_x;
            } else {
                // Unknown glyph: advance by fallback width
                cursor_x += self.font.size_px as f64 * 0.6;
            }
        }
        out
    }

    /// Measure the pixel width of `text` using glyph advances.
    pub fn measure_width(&self, text: &str) -> f64 {
        self.font.measure_width(text)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AnnotationLayer
// ─────────────────────────────────────────────────────────────────────────────

/// A collection of text annotations attached to 3-D world positions.
///
/// Annotations can be projected onto a 2-D screen for display using a simple
/// perspective projection.
#[derive(Debug, Clone)]
pub struct AnnotationLayer {
    /// List of `(label, world_position)` pairs.
    pub annotations: Vec<(String, [f64; 3])>,
}

impl AnnotationLayer {
    /// Create an empty `AnnotationLayer`.
    pub fn new() -> Self {
        Self {
            annotations: Vec::new(),
        }
    }

    /// Add an annotation with `label` at world position `pos`.
    pub fn add(&mut self, label: impl Into<String>, pos: [f64; 3]) {
        self.annotations.push((label.into(), pos));
    }

    /// Project all annotations from world space to screen space.
    ///
    /// Uses a simple pinhole camera model:
    /// - The camera is at `camera_pos` looking toward the origin.
    /// - `fov` is the horizontal field of view in radians.
    /// - Returns `(label, screen_xy)` for annotations in front of the camera.
    ///   Annotations behind the camera (negative depth) are omitted.
    pub fn project_to_screen(&self, camera_pos: [f64; 3], fov: f64) -> Vec<(String, [f64; 2])> {
        // Build a simple view matrix: forward = -camera_pos normalised
        let cx = camera_pos[0];
        let cy = camera_pos[1];
        let cz = camera_pos[2];
        let cam_len = (cx * cx + cy * cy + cz * cz).sqrt().max(1e-10);

        // Forward (toward origin)
        let fwd = [-cx / cam_len, -cy / cam_len, -cz / cam_len];
        // Right = fwd × (0,1,0) (simplified; handle degenerate up later)
        let up = [0.0_f64, 1.0, 0.0];
        let right = cross3(fwd, up);
        let right_len = (right[0].powi(2) + right[1].powi(2) + right[2].powi(2))
            .sqrt()
            .max(1e-10);
        let right = [
            right[0] / right_len,
            right[1] / right_len,
            right[2] / right_len,
        ];
        let up2 = cross3(right, fwd);

        let half_fov_tan = (fov * 0.5).tan();

        let mut result = Vec::new();
        for (label, world_pos) in &self.annotations {
            // Vector from camera to annotation
            let rel = [world_pos[0] - cx, world_pos[1] - cy, world_pos[2] - cz];
            let depth = dot3(rel, fwd);
            if depth <= 0.0 {
                // Behind the camera
                continue;
            }
            // Project onto image plane at depth = 1
            let proj_x = dot3(rel, right) / (depth * half_fov_tan);
            let proj_y = dot3(rel, up2) / (depth * half_fov_tan);
            // Convert from [-1,1] to screen pixels (NDC)
            result.push((label.clone(), [proj_x, proj_y]));
        }
        result
    }

    /// Return all annotations as a slice.
    pub fn all(&self) -> &[(String, [f64; 3])] {
        &self.annotations
    }

    /// Remove all annotations.
    pub fn clear(&mut self) {
        self.annotations.clear();
    }

    /// Number of annotations.
    pub fn len(&self) -> usize {
        self.annotations.len()
    }

    /// Return `true` if there are no annotations.
    pub fn is_empty(&self) -> bool {
        self.annotations.is_empty()
    }
}

impl Default for AnnotationLayer {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (local — no nalgebra in this crate)
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── GlyphBitmap ───────────────────────────────────────────────────────────

    #[test]
    fn test_glyph_bitmap_new_all_zero() {
        let g = GlyphBitmap::new('A', 8, 12, 6.0);
        assert_eq!(g.width, 8);
        assert_eq!(g.height, 12);
        assert!(g.pixels.iter().all(|&p| p == 0));
        assert_eq!(g.advance_x, 6.0);
    }

    #[test]
    fn test_glyph_bitmap_set_get_pixel() {
        let mut g = GlyphBitmap::new('B', 4, 4, 3.0);
        g.set_pixel(1, 2, 200);
        assert_eq!(g.get_pixel(1, 2), 200);
        assert_eq!(g.get_pixel(0, 0), 0);
    }

    #[test]
    fn test_glyph_bitmap_oob_set_is_noop() {
        let mut g = GlyphBitmap::new('C', 4, 4, 3.0);
        g.set_pixel(10, 10, 255); // out of bounds — should not panic
        assert!(g.pixels.iter().all(|&p| p == 0));
    }

    #[test]
    fn test_glyph_bitmap_oob_get_returns_zero() {
        let g = GlyphBitmap::new('D', 4, 4, 3.0);
        assert_eq!(g.get_pixel(100, 100), 0);
    }

    #[test]
    fn test_glyph_bitmap_pixel_count() {
        let g = GlyphBitmap::new('E', 5, 7, 4.0);
        assert_eq!(g.pixel_count(), 35);
    }

    #[test]
    fn test_glyph_bitmap_is_empty_zero_dimensions() {
        let g = GlyphBitmap::new('F', 0, 8, 0.0);
        assert!(g.is_empty());
    }

    #[test]
    fn test_glyph_bitmap_is_not_empty() {
        let g = GlyphBitmap::new('G', 8, 8, 5.0);
        assert!(!g.is_empty());
    }

    // ── BitmapFont ────────────────────────────────────────────────────────────

    #[test]
    fn test_bitmap_font_add_and_get_glyph() {
        let mut font = BitmapFont::new(12, 14.0);
        let g = GlyphBitmap::new('H', 8, 12, 7.0);
        font.add_glyph(g);
        assert!(font.get_glyph('H').is_some());
        assert!(font.get_glyph('Z').is_none());
    }

    #[test]
    fn test_bitmap_font_advance_for_known_char() {
        let mut font = BitmapFont::new(12, 14.0);
        let g = GlyphBitmap::new('I', 8, 12, 7.5);
        font.add_glyph(g);
        assert!((font.advance_for('I') - 7.5).abs() < 1e-10);
    }

    #[test]
    fn test_bitmap_font_advance_for_missing_char_fallback() {
        let font = BitmapFont::new(12, 14.0);
        let adv = font.advance_for('?');
        assert!((adv - 12.0 * 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_bitmap_font_measure_width() {
        let mut font = BitmapFont::new(12, 14.0);
        font.add_glyph(GlyphBitmap::new('a', 8, 12, 6.0));
        font.add_glyph(GlyphBitmap::new('b', 8, 12, 6.0));
        assert!((font.measure_width("ab") - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_bitmap_font_glyph_count() {
        let mut font = BitmapFont::new(16, 20.0);
        assert_eq!(font.glyph_count(), 0);
        font.add_glyph(GlyphBitmap::new('X', 12, 16, 8.0));
        assert_eq!(font.glyph_count(), 1);
    }

    #[test]
    fn test_bitmap_font_ascii_placeholder() {
        let font = BitmapFont::ascii_placeholder(16);
        // Space and printable ASCII
        assert!(font.get_glyph(' ').is_some());
        assert!(font.get_glyph('A').is_some());
        assert!(font.get_glyph('z').is_some());
        // 0x7E = '~'
        assert!(font.get_glyph('~').is_some());
        // Non-ASCII not included
        assert!(font.get_glyph('€').is_none());
    }

    #[test]
    fn test_bitmap_font_replace_glyph() {
        let mut font = BitmapFont::new(8, 10.0);
        font.add_glyph(GlyphBitmap::new('Y', 8, 8, 5.0));
        font.add_glyph(GlyphBitmap::new('Y', 8, 8, 7.0)); // replace
        assert!((font.advance_for('Y') - 7.0).abs() < 1e-10);
    }

    // ── TextLayout ────────────────────────────────────────────────────────────

    #[test]
    fn test_text_layout_basic() {
        let font = BitmapFont::ascii_placeholder(12);
        let layout = TextLayout::new("Hi", font, 12.0, 1000.0, false);
        let positions = layout.layout();
        assert_eq!(positions.len(), 2);
        // Characters should be placed left to right
        assert!(
            positions[0].0 < positions[1].0,
            "second char should be to the right"
        );
    }

    #[test]
    fn test_text_layout_all_chars_present() {
        let font = BitmapFont::ascii_placeholder(12);
        let text = "Hello";
        let layout = TextLayout::new(text, font, 12.0, 1000.0, false);
        let positions = layout.layout();
        assert_eq!(positions.len(), text.len());
        for ((_, _, c), expected) in positions.iter().zip(text.chars()) {
            assert_eq!(*c, expected);
        }
    }

    #[test]
    fn test_text_layout_newline_resets_x() {
        let font = BitmapFont::ascii_placeholder(12);
        let layout = TextLayout::new("A\nB", font, 12.0, 1000.0, false);
        let positions = layout.layout();
        // 'A' at (x, 0), 'B' at (0, line_height)
        assert_eq!(positions.len(), 2);
        let (xa, ya, ca) = positions[0];
        let (xb, yb, cb) = positions[1];
        assert_eq!(ca, 'A');
        assert_eq!(cb, 'B');
        assert!(xa >= 0.0);
        assert!(
            (xb).abs() < 1e-6,
            "B should start at x=0 after newline, got {xb}"
        );
        assert!(yb > ya, "B should be on a lower line");
    }

    #[test]
    fn test_text_layout_wrapping() {
        let font = BitmapFont::ascii_placeholder(12);
        // 10 chars each with advance ~7.2 → 72px total; line_width = 40 → should wrap
        let text = "AAAAAAAAAA";
        let layout = TextLayout::new(text, font, 12.0, 40.0, true);
        let positions = layout.layout();
        // Some characters should be on the second line
        let second_line: Vec<_> = positions.iter().filter(|&&(_, y, _)| y > 0.0).collect();
        assert!(
            !second_line.is_empty(),
            "expected wrapping to produce a second line"
        );
    }

    #[test]
    fn test_text_layout_no_wrapping_single_line() {
        let font = BitmapFont::ascii_placeholder(12);
        let text = "AAAA";
        let layout = TextLayout::new(text, font, 12.0, 10.0, false); // narrow but no wrap
        let positions = layout.layout();
        // Without wrapping everything stays on y=0
        for (_, y, _) in &positions {
            assert!(
                (*y).abs() < 1e-6,
                "all characters should be on y=0 without wrapping"
            );
        }
    }

    #[test]
    fn test_text_layout_clipped() {
        let font = BitmapFont::ascii_placeholder(12);
        let text = "ABC";
        let layout = TextLayout::new(text, font, 12.0, 1000.0, false);
        // Clip to a very small area — only 'A' at x=0 should survive
        let clipped = layout.layout_clipped(5.0, 100.0);
        assert!(!clipped.is_empty());
        for (x, _, _) in &clipped {
            assert!(*x <= 5.0);
        }
    }

    #[test]
    fn test_text_layout_empty_string() {
        let font = BitmapFont::ascii_placeholder(12);
        let layout = TextLayout::new("", font, 12.0, 1000.0, false);
        assert!(layout.layout().is_empty());
    }

    // ── SdfFont ───────────────────────────────────────────────────────────────

    #[test]
    fn test_sdf_font_generate_and_has_glyph() {
        let mut sdf = SdfFont::new(32, 4.0);
        sdf.generate_disc_glyph('O');
        assert!(sdf.has_glyph('O'));
        assert!(!sdf.has_glyph('X'));
    }

    #[test]
    fn test_sdf_font_render_glyph_size() {
        let mut sdf = SdfFont::new(16, 2.0);
        sdf.generate_disc_glyph('Q');
        let bitmap = sdf.render_glyph('Q', 2.0);
        assert_eq!(bitmap.len(), 32 * 32); // 16*2 = 32
    }

    #[test]
    fn test_sdf_font_render_missing_glyph() {
        let sdf = SdfFont::new(16, 2.0);
        let bitmap = sdf.render_glyph('?', 1.0);
        assert!(bitmap.is_empty());
    }

    #[test]
    fn test_sdf_font_disc_center_is_opaque() {
        let mut sdf = SdfFont::new(32, 4.0);
        sdf.generate_disc_glyph('P');
        let bitmap = sdf.render_glyph('P', 1.0);
        // Centre pixel should be inside the disc → alpha=255
        let cx = 16usize;
        let cy = 16usize;
        assert_eq!(
            bitmap[cy * 32 + cx],
            255,
            "centre of disc glyph should be opaque"
        );
    }

    #[test]
    fn test_sdf_font_corners_are_transparent() {
        let mut sdf = SdfFont::new(32, 4.0);
        sdf.generate_disc_glyph('R');
        let bitmap = sdf.render_glyph('R', 1.0);
        // Corner pixels at (0,0) should be outside the disc
        assert_eq!(bitmap[0], 0, "corner of disc glyph should be transparent");
    }

    #[test]
    fn test_sdf_font_glyph_count() {
        let mut sdf = SdfFont::new(16, 2.0);
        assert_eq!(sdf.glyph_count(), 0);
        sdf.generate_disc_glyph('A');
        sdf.generate_disc_glyph('B');
        assert_eq!(sdf.glyph_count(), 2);
    }

    #[test]
    fn test_sdf_font_insert_custom() {
        let mut sdf = SdfFont::new(4, 1.0);
        // All-positive SDF: everything inside
        let data = vec![1.0f32; 16];
        sdf.insert('Z', data);
        let bitmap = sdf.render_glyph('Z', 1.0);
        assert!(bitmap.iter().all(|&p| p == 255));
    }

    #[test]
    fn test_sdf_font_all_negative_is_transparent() {
        let mut sdf = SdfFont::new(4, 1.0);
        let data = vec![-1.0f32; 16];
        sdf.insert('W', data);
        let bitmap = sdf.render_glyph('W', 1.0);
        assert!(bitmap.iter().all(|&p| p == 0));
    }

    // ── TextRenderer ──────────────────────────────────────────────────────────

    #[test]
    fn test_text_renderer_renders_known_glyph() {
        let font = BitmapFont::ascii_placeholder(8);
        let renderer = TextRenderer::new(font);
        // 'A' should produce some non-zero pixels
        let pixels = renderer.render_string("A", 0.0, 0.0);
        // ascii_placeholder paints vertical bars for non-space chars
        assert!(!pixels.is_empty(), "rendering 'A' should produce pixels");
    }

    #[test]
    fn test_text_renderer_space_has_no_pixels() {
        let font = BitmapFont::ascii_placeholder(8);
        let renderer = TextRenderer::new(font);
        let pixels = renderer.render_string(" ", 0.0, 0.0);
        assert!(pixels.is_empty(), "space should produce no opaque pixels");
    }

    #[test]
    fn test_text_renderer_measure_width() {
        let font = BitmapFont::ascii_placeholder(12);
        let renderer = TextRenderer::new(font.clone());
        let w = renderer.measure_width("AB");
        assert!(w > 0.0);
        // Should equal sum of advance widths
        let expected = font.advance_for('A') + font.advance_for('B');
        assert!((w - expected).abs() < 1e-10);
    }

    #[test]
    fn test_text_renderer_pixel_coords_advance() {
        let font = BitmapFont::ascii_placeholder(8);
        let renderer = TextRenderer::new(font);
        let pixels = renderer.render_string("AB", 10.0, 5.0);
        // All pixel y values should be >= 5.0
        for (_, py, _) in &pixels {
            assert!(*py >= 5.0, "pixel y should be >= origin y");
        }
        // All pixel x values should be >= 10.0
        for (px, _, _) in &pixels {
            assert!(*px >= 10.0, "pixel x should be >= origin x");
        }
    }

    // ── AnnotationLayer ───────────────────────────────────────────────────────

    #[test]
    fn test_annotation_layer_add_and_len() {
        let mut layer = AnnotationLayer::new();
        assert_eq!(layer.len(), 0);
        assert!(layer.is_empty());
        layer.add("hello", [1.0, 2.0, 3.0]);
        assert_eq!(layer.len(), 1);
        assert!(!layer.is_empty());
    }

    #[test]
    fn test_annotation_layer_clear() {
        let mut layer = AnnotationLayer::new();
        layer.add("A", [0.0, 0.0, 1.0]);
        layer.add("B", [1.0, 0.0, 0.0]);
        layer.clear();
        assert!(layer.is_empty());
    }

    #[test]
    fn test_annotation_projection_in_front() {
        let mut layer = AnnotationLayer::new();
        // Place annotation directly in front of the camera
        layer.add("origin", [0.0, 0.0, 0.0]);
        // Camera at (0, 0, 5) looking toward origin
        let projected = layer.project_to_screen([0.0, 0.0, 5.0], std::f64::consts::FRAC_PI_2);
        assert_eq!(
            projected.len(),
            1,
            "annotation in front should be projected"
        );
        let (label, _) = &projected[0];
        assert_eq!(label, "origin");
    }

    #[test]
    fn test_annotation_projection_behind_camera_omitted() {
        let mut layer = AnnotationLayer::new();
        // Camera at (0,0,5) looking toward origin; put annotation behind at (0,0,10)
        layer.add("behind", [0.0, 0.0, 10.0]);
        let projected = layer.project_to_screen([0.0, 0.0, 5.0], std::f64::consts::FRAC_PI_2);
        assert!(
            projected.is_empty(),
            "annotation behind camera should be omitted"
        );
    }

    #[test]
    fn test_annotation_projection_empty() {
        let layer = AnnotationLayer::default();
        let projected = layer.project_to_screen([0.0, 0.0, 5.0], std::f64::consts::FRAC_PI_2);
        assert!(projected.is_empty());
    }

    #[test]
    fn test_annotation_projection_screen_coords_finite() {
        let mut layer = AnnotationLayer::new();
        layer.add("pt1", [1.0, 0.0, 0.0]);
        layer.add("pt2", [0.0, 1.0, 0.0]);
        let projected = layer.project_to_screen([5.0, 5.0, 5.0], 1.0);
        for (_, xy) in &projected {
            assert!(xy[0].is_finite(), "screen x must be finite");
            assert!(xy[1].is_finite(), "screen y must be finite");
        }
    }

    #[test]
    fn test_annotation_all_returns_all() {
        let mut layer = AnnotationLayer::new();
        layer.add("a", [0.0, 0.0, 0.0]);
        layer.add("b", [1.0, 1.0, 1.0]);
        assert_eq!(layer.all().len(), 2);
    }

    #[test]
    fn test_annotation_projection_centred_is_near_zero() {
        let mut layer = AnnotationLayer::new();
        // Annotation at origin, camera along Z: projection should be near (0,0)
        layer.add("centre", [0.0, 0.0, 0.0]);
        let projected = layer.project_to_screen([0.0, 0.0, 5.0], std::f64::consts::FRAC_PI_2);
        assert_eq!(projected.len(), 1);
        let (_, xy) = &projected[0];
        assert!(
            xy[0].abs() < 1e-6,
            "centre annotation x should be ~0, got {}",
            xy[0]
        );
        assert!(
            xy[1].abs() < 1e-6,
            "centre annotation y should be ~0, got {}",
            xy[1]
        );
    }

    #[test]
    fn test_sdf_font_render_zero_scale() {
        let mut sdf = SdfFont::new(16, 2.0);
        sdf.generate_disc_glyph('S');
        let bitmap = sdf.render_glyph('S', 0.0);
        assert!(bitmap.is_empty(), "zero scale should produce empty bitmap");
    }

    // ── BitmapFont::ascii_default (F4) ────────────────────────────────────────

    #[test]
    fn test_ascii_default_glyph_a_nonzero() {
        let font = BitmapFont::ascii_default();
        let glyph = font.get_glyph('A').expect("'A' must be present");
        assert!(
            glyph.pixels.iter().any(|&b| b != 0),
            "glyph 'A' must have at least one set pixel"
        );
    }

    #[test]
    fn test_ascii_default_glyph_bits_a_nonzero() {
        let bits = BitmapFont::glyph_bits('A').expect("'A' must have bits");
        assert!(
            bits.iter().any(|&b| b != 0),
            "glyph_bits('A') must have at least one non-zero row"
        );
    }

    #[test]
    fn test_ascii_default_has_all_printable() {
        let font = BitmapFont::ascii_default();
        for code in 0x20u32..=0x7Eu32 {
            if let Some(c) = char::from_u32(code) {
                assert!(
                    font.get_glyph(c).is_some(),
                    "ascii_default missing glyph for U+{code:04X} '{c}'"
                );
            }
        }
    }

    #[test]
    fn test_ascii_default_space_is_blank() {
        let font = BitmapFont::ascii_default();
        let glyph = font.get_glyph(' ').expect("space must be present");
        assert!(
            glyph.pixels.iter().all(|&b| b == 0),
            "space glyph must be all-zero pixels"
        );
    }

    #[test]
    fn test_glyph_bits_out_of_range_returns_none() {
        assert!(BitmapFont::glyph_bits('\t').is_none());
        assert!(BitmapFont::glyph_bits('\x7F').is_none());
    }
}
