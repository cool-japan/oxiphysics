// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core types for post-processing: PostColor and Image.

// ─── Color ────────────────────────────────────────────────────────────────────

/// RGBA colour with `f32` components, used by post-processing effects.
///
/// This is a *separate* type from `primitives::Color`; it lives in the
/// post-processing pipeline only.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PostColor {
    /// Red channel.
    pub r: f32,
    /// Green channel.
    pub g: f32,
    /// Blue channel.
    pub b: f32,
    /// Alpha channel.
    pub a: f32,
}

impl PostColor {
    /// Create a new colour from RGBA components.
    #[inline]
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }

    /// Perceptual luminance: 0.2126·r + 0.7152·g + 0.0722·b.
    #[inline]
    pub fn luminance(&self) -> f32 {
        0.2126 * self.r + 0.7152 * self.g + 0.0722 * self.b
    }

    /// Clamp all channels to `[0, 1]`.
    #[inline]
    pub fn clamp(&self) -> Self {
        Self {
            r: self.r.clamp(0.0, 1.0),
            g: self.g.clamp(0.0, 1.0),
            b: self.b.clamp(0.0, 1.0),
            a: self.a.clamp(0.0, 1.0),
        }
    }

    /// Linear interpolation between `self` and `other` at parameter `t ∈ [0,1]`.
    #[inline]
    pub fn lerp(&self, other: &PostColor, t: f32) -> PostColor {
        PostColor {
            r: self.r + (other.r - self.r) * t,
            g: self.g + (other.g - self.g) * t,
            b: self.b + (other.b - self.b) * t,
            a: self.a + (other.a - self.a) * t,
        }
    }

    /// Multiply all channels by scalar `s`.
    #[inline]
    pub fn scale(&self, s: f32) -> PostColor {
        PostColor {
            r: self.r * s,
            g: self.g * s,
            b: self.b * s,
            a: self.a * s,
        }
    }

    /// Component-wise addition.
    #[inline]
    pub fn add(&self, other: &PostColor) -> PostColor {
        PostColor {
            r: self.r + other.r,
            g: self.g + other.g,
            b: self.b + other.b,
            a: self.a + other.a,
        }
    }

    /// Transparent black.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }
}

// ─── Image ────────────────────────────────────────────────────────────────────

/// CPU-side framebuffer holding colour and depth data.
///
/// Pixels are stored in row-major order: pixel at `(x, y)` lives at index
/// `y * width + x`.
pub struct Image {
    /// Width in pixels.
    pub width: usize,
    /// Height in pixels.
    pub height: usize,
    /// Colour buffer (row-major).
    pub pixels: Vec<PostColor>,
    /// Depth buffer in `[0, 1]` (row-major).
    pub depth: Vec<f32>,
}

impl Image {
    /// Create a new image filled with transparent black and depth 1.0.
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            pixels: vec![PostColor::zero(); n],
            depth: vec![1.0_f32; n],
        }
    }

    #[inline]
    pub(crate) fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    /// Read the colour at pixel `(x, y)`.
    #[inline]
    pub fn get(&self, x: usize, y: usize) -> PostColor {
        self.pixels[self.idx(x, y)]
    }

    /// Write a colour to pixel `(x, y)`.
    #[inline]
    pub fn set(&mut self, x: usize, y: usize, c: PostColor) {
        let i = self.idx(x, y);
        self.pixels[i] = c;
    }

    /// Read the depth at pixel `(x, y)`.
    #[inline]
    pub fn get_depth(&self, x: usize, y: usize) -> f32 {
        self.depth[self.idx(x, y)]
    }

    /// Write a depth value to pixel `(x, y)`.
    #[inline]
    pub fn set_depth(&mut self, x: usize, y: usize, d: f32) {
        let i = self.idx(x, y);
        self.depth[i] = d;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── PostColor ────────────────────────────────────────────────────────────

    #[test]
    fn test_color_luminance_white() {
        let white = PostColor::new(1.0, 1.0, 1.0, 1.0);
        // 0.2126 + 0.7152 + 0.0722 = 1.0
        let lum = white.luminance();
        assert!(
            (lum - 1.0).abs() < 1e-5,
            "white luminance should be 1.0, got {lum}"
        );
    }

    #[test]
    fn test_color_luminance_black() {
        let black = PostColor::new(0.0, 0.0, 0.0, 1.0);
        assert_eq!(black.luminance(), 0.0);
    }

    #[test]
    fn test_color_clamp() {
        let c = PostColor::new(-0.5, 0.5, 1.5, 2.0);
        let clamped = c.clamp();
        assert_eq!(clamped.r, 0.0);
        assert_eq!(clamped.g, 0.5);
        assert_eq!(clamped.b, 1.0);
        assert_eq!(clamped.a, 1.0);
    }

    #[test]
    fn test_color_lerp_midpoint() {
        let a = PostColor::new(0.0, 0.0, 0.0, 0.0);
        let b = PostColor::new(1.0, 1.0, 1.0, 1.0);
        let mid = a.lerp(&b, 0.5);
        assert!((mid.r - 0.5).abs() < 1e-6);
        assert!((mid.g - 0.5).abs() < 1e-6);
        assert!((mid.b - 0.5).abs() < 1e-6);
        assert!((mid.a - 0.5).abs() < 1e-6);
    }
}
