// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDR (High Dynamic Range) floating-point framebuffer with tone mapping.
//!
//! [`HdrFramebuffer`] stores pixels as 32-bit floating-point RGBA values, allowing
//! the accumulation of lighting contributions that exceed the `[0, 1]` LDR range.
//! The buffer is then converted to 8-bit LDR via a configurable [`ToneMapper`].
//!
//! ## Supported tone-mapping operators
//!
//! | Operator | Description |
//! |---|---|
//! | [`ToneMapMethod::Reinhard`] | Classic photographic operator; soft roll-off |
//! | [`ToneMapMethod::ReinhardExtended`] | Reinhard with configurable white point |
//! | [`ToneMapMethod::Aces`] | Academy Color Encoding System; cinematic look |
//! | [`ToneMapMethod::Filmic`] | Uncharted 2 / Hable filmic curve |
//! | [`ToneMapMethod::AgX`] | AgX-inspired; low-saturation highlights |
//! | [`ToneMapMethod::Linear`] | No tone mapping; simple clamp to `[0, 1]` |
//! | [`ToneMapMethod::Exposure`] | `exposure * color`, then linear clamp |
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_viz::hdr_framebuffer::{HdrFramebuffer, ToneMapMethod, ToneMapper};
//!
//! let mut fb = HdrFramebuffer::new(4, 4);
//! fb.set_pixel(1, 1, [2.0, 1.5, 0.5, 1.0]); // HDR pixel (overexposed)
//! let mapper = ToneMapper::new(ToneMapMethod::Aces, 1.0, 2.2);
//! let ldr = fb.to_ldr(&mapper);
//! assert_eq!(ldr.len(), 4 * 4 * 4); // RGBA bytes
//! ```

#![allow(dead_code)]

// ── HdrPixel ─────────────────────────────────────────────────────────────────

/// An RGBA pixel with f32 components; values may exceed `[0, 1]`.
pub type HdrPixel = [f32; 4];

// ── HdrFramebuffer ────────────────────────────────────────────────────────────

/// A floating-point RGBA framebuffer for HDR rendering.
///
/// Pixels are stored in row-major order: pixel `(x, y)` is at index
/// `y * width + x`.
#[derive(Debug, Clone)]
pub struct HdrFramebuffer {
    /// Framebuffer width in pixels.
    pub width: usize,
    /// Framebuffer height in pixels.
    pub height: usize,
    /// RGBA pixel data (row-major, length `width * height`).
    pub pixels: Vec<HdrPixel>,
}

impl HdrFramebuffer {
    /// Create a new all-black (0, 0, 0, 1) framebuffer.
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![[0.0, 0.0, 0.0, 1.0]; width * height],
        }
    }

    /// Create a framebuffer cleared to `color`.
    pub fn with_clear(width: usize, height: usize, color: HdrPixel) -> Self {
        Self {
            width,
            height,
            pixels: vec![color; width * height],
        }
    }

    /// Return a reference to the pixel at `(x, y)`.
    #[inline]
    pub fn pixel(&self, x: usize, y: usize) -> &HdrPixel {
        &self.pixels[y * self.width + x]
    }

    /// Return a mutable reference to the pixel at `(x, y)`.
    #[inline]
    pub fn pixel_mut(&mut self, x: usize, y: usize) -> &mut HdrPixel {
        &mut self.pixels[y * self.width + x]
    }

    /// Set pixel `(x, y)` to `value`.
    #[inline]
    pub fn set_pixel(&mut self, x: usize, y: usize, value: HdrPixel) {
        self.pixels[y * self.width + x] = value;
    }

    /// Additively blend `value` into pixel `(x, y)`.
    #[inline]
    pub fn add_pixel(&mut self, x: usize, y: usize, value: HdrPixel) {
        let p = self.pixel_mut(x, y);
        p[0] += value[0];
        p[1] += value[1];
        p[2] += value[2];
        p[3] += value[3];
    }

    /// Clear all pixels to `color`.
    pub fn clear(&mut self, color: HdrPixel) {
        for p in &mut self.pixels {
            *p = color;
        }
    }

    /// Apply exposure and convert to 8-bit LDR RGBA bytes.
    ///
    /// The output is a flat `Vec<u8>` of length `width * height * 4` (RGBA),
    /// compatible with HTML5 Canvas `ImageData.data` and PNG export.
    pub fn to_ldr(&self, mapper: &ToneMapper) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.width * self.height * 4);
        for px in &self.pixels {
            let ldr = mapper.map(*px);
            out.push(linear_to_srgb(ldr[0]));
            out.push(linear_to_srgb(ldr[1]));
            out.push(linear_to_srgb(ldr[2]));
            out.push((ldr[3].clamp(0.0, 1.0) * 255.0 + 0.5) as u8);
        }
        out
    }

    /// Blit another framebuffer into this one at offset `(ox, oy)`.
    ///
    /// Pixels outside the destination bounds are silently clipped.
    pub fn blit(&mut self, src: &HdrFramebuffer, ox: usize, oy: usize) {
        for sy in 0..src.height {
            let dy = oy + sy;
            if dy >= self.height {
                break;
            }
            for sx in 0..src.width {
                let dx = ox + sx;
                if dx >= self.width {
                    break;
                }
                self.set_pixel(dx, dy, *src.pixel(sx, sy));
            }
        }
    }

    /// Compute per-channel average luminance (log-average for Reinhard auto-exposure).
    pub fn log_average_luminance(&self) -> f32 {
        let n = self.pixels.len() as f32;
        let sum: f32 = self
            .pixels
            .iter()
            .map(|p| luminance(p).max(1e-6).ln())
            .sum();
        (sum / n).exp()
    }

    /// Return the maximum luminance across all pixels.
    pub fn max_luminance(&self) -> f32 {
        self.pixels.iter().map(luminance).fold(0.0_f32, f32::max)
    }
}

// ── Luminance helper ────────────────────────────────────────────────────────

#[inline]
fn luminance(px: &HdrPixel) -> f32 {
    // Rec. 709 luminance coefficients
    0.2126 * px[0] + 0.7152 * px[1] + 0.0722 * px[2]
}

// ── ToneMapMethod ────────────────────────────────────────────────────────────

/// Tone-mapping operator selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToneMapMethod {
    /// Simple clamp — no compression, HDR values saturate.
    Linear,
    /// Exposure-only: `exposure * color`, then clamp.
    Exposure,
    /// Reinhard global operator: `color / (1 + color)`.
    Reinhard,
    /// Extended Reinhard: maps white point `L_white` to exactly 1.
    ReinhardExtended {
        /// White-point luminance; values above this are mapped to 1.
        l_white: f32,
    },
    /// ACES filmic tone mapping (approximation by Narkowicz).
    Aces,
    /// Uncharted 2 / Hable filmic curve.
    Filmic,
    /// AgX-inspired low-saturation high-range operator.
    AgX,
}

// ── ToneMapper ───────────────────────────────────────────────────────────────

/// Applies tone mapping and optional gamma correction to an [`HdrPixel`].
#[derive(Debug, Clone)]
pub struct ToneMapper {
    /// Tone-mapping method.
    pub method: ToneMapMethod,
    /// Linear exposure multiplier applied before tone mapping.
    pub exposure: f32,
    /// Gamma value for the sRGB conversion step.  2.2 is the standard.
    ///
    /// Note: [`HdrFramebuffer::to_ldr`] calls [`linear_to_srgb`] which
    /// handles gamma; set this to 1.0 to rely solely on [`linear_to_srgb`].
    pub gamma: f32,
}

impl ToneMapper {
    /// Create a new tone mapper.
    pub fn new(method: ToneMapMethod, exposure: f32, gamma: f32) -> Self {
        Self {
            method,
            exposure,
            gamma,
        }
    }

    /// Apply tone mapping to a single HDR pixel, returning an LDR pixel.
    ///
    /// The returned values are in `[0, 1]` (linear, before gamma encoding).
    pub fn map(&self, px: HdrPixel) -> HdrPixel {
        // Apply exposure
        let rgb = [
            px[0] * self.exposure,
            px[1] * self.exposure,
            px[2] * self.exposure,
        ];
        let mapped = match self.method {
            ToneMapMethod::Linear => [
                rgb[0].clamp(0.0, 1.0),
                rgb[1].clamp(0.0, 1.0),
                rgb[2].clamp(0.0, 1.0),
            ],
            ToneMapMethod::Exposure => [
                rgb[0].clamp(0.0, 1.0),
                rgb[1].clamp(0.0, 1.0),
                rgb[2].clamp(0.0, 1.0),
            ],
            ToneMapMethod::Reinhard => reinhard(rgb),
            ToneMapMethod::ReinhardExtended { l_white } => reinhard_extended(rgb, l_white),
            ToneMapMethod::Aces => aces(rgb),
            ToneMapMethod::Filmic => filmic(rgb),
            ToneMapMethod::AgX => agx(rgb),
        };
        [mapped[0], mapped[1], mapped[2], px[3].clamp(0.0, 1.0)]
    }
}

// ── Tone mapping operators ────────────────────────────────────────────────────

fn reinhard(rgb: [f32; 3]) -> [f32; 3] {
    [
        rgb[0] / (1.0 + rgb[0]),
        rgb[1] / (1.0 + rgb[1]),
        rgb[2] / (1.0 + rgb[2]),
    ]
}

fn reinhard_extended(rgb: [f32; 3], l_white: f32) -> [f32; 3] {
    let lw2 = l_white * l_white;
    let f = move |x: f32| (x * (1.0 + x / lw2) / (1.0 + x)).clamp(0.0, 1.0);
    [f(rgb[0]), f(rgb[1]), f(rgb[2])]
}

/// ACES approximation by Krzysztof Narkowicz.
fn aces(rgb: [f32; 3]) -> [f32; 3] {
    let f = |x: f32| -> f32 {
        let a = 2.51_f32;
        let b = 0.03_f32;
        let c = 2.43_f32;
        let d = 0.59_f32;
        let e = 0.14_f32;
        ((x * (a * x + b)) / (x * (c * x + d) + e)).clamp(0.0, 1.0)
    };
    [f(rgb[0]), f(rgb[1]), f(rgb[2])]
}

/// Uncharted 2 filmic curve by John Hable.
fn filmic(rgb: [f32; 3]) -> [f32; 3] {
    let hable = |x: f32| -> f32 {
        let a = 0.15_f32;
        let b = 0.50_f32;
        let c = 0.10_f32;
        let d = 0.20_f32;
        let e = 0.02_f32;
        let f = 0.30_f32;
        ((x * (a * x + c * b) + d * e) / (x * (a * x + b) + d * f)) - e / f
    };
    let white = hable(11.2);
    let f = |x: f32| (hable(x * 2.0) / white).clamp(0.0, 1.0);
    [f(rgb[0]), f(rgb[1]), f(rgb[2])]
}

/// Simplified AgX-inspired tone mapper (desaturates highlights).
fn agx(rgb: [f32; 3]) -> [f32; 3] {
    // AgX applies a matrix transform to move into a log-linear space,
    // then applies a sigmoid-shaped curve.
    // This is a simplified approximation for CPU preview purposes.
    let f = |x: f32| -> f32 {
        // Soft clip in log domain, similar in shape to AgX
        let l = (x + 0.0001).ln();
        let mapped = (l * 0.18 + 0.5).clamp(0.0, 1.0);
        // S-curve
        let t = mapped;
        t * t * (3.0 - 2.0 * t)
    };
    [f(rgb[0]), f(rgb[1]), f(rgb[2])]
}

/// Convert a linear-light value to sRGB gamma (IEC 61966-2-1).
#[inline]
pub fn linear_to_srgb(x: f32) -> u8 {
    let gamma = if x <= 0.0031308 {
        x * 12.92
    } else {
        1.055 * x.powf(1.0 / 2.4) - 0.055
    };
    (gamma.clamp(0.0, 1.0) * 255.0 + 0.5) as u8
}

/// Convert an sRGB byte value to linear-light f32.
#[inline]
pub fn srgb_to_linear(byte: u8) -> f32 {
    let s = byte as f32 / 255.0;
    if s <= 0.04045 {
        s / 12.92
    } else {
        ((s + 0.055) / 1.055).powf(2.4)
    }
}

// ── HdrRasterizer ────────────────────────────────────────────────────────────

/// A simple triangle rasterizer that writes into an [`HdrFramebuffer`].
///
/// Supports per-vertex linear HDR colors and interpolates them barycentrically.
#[derive(Debug, Default)]
pub struct HdrRasterizer {
    // Target framebuffer is passed mutably to each draw call.
}

impl HdrRasterizer {
    /// Create a new rasterizer.
    pub fn new() -> Self {
        Self {}
    }

    /// Rasterize a single triangle into `fb`.
    ///
    /// * `verts` — screen-space positions `(x, y)` for each vertex
    /// * `colors` — per-vertex HDR colors
    pub fn draw_triangle(
        &self,
        fb: &mut HdrFramebuffer,
        verts: [[f32; 2]; 3],
        colors: [HdrPixel; 3],
    ) {
        let [v0, v1, v2] = verts;
        let [c0, c1, c2] = colors;
        // Bounding box
        let min_x = (v0[0].min(v1[0]).min(v2[0]).max(0.0) as usize).min(fb.width.saturating_sub(1));
        let max_x = (v0[0].max(v1[0]).max(v2[0]).ceil() as usize).min(fb.width.saturating_sub(1));
        let min_y =
            (v0[1].min(v1[1]).min(v2[1]).max(0.0) as usize).min(fb.height.saturating_sub(1));
        let max_y = (v0[1].max(v1[1]).max(v2[1]).ceil() as usize).min(fb.height.saturating_sub(1));

        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let px = [x as f32 + 0.5, y as f32 + 0.5];
                if let Some([u, v, w]) = barycentric(v0, v1, v2, px)
                    && u >= 0.0
                    && v >= 0.0
                    && w >= 0.0
                {
                    let color = lerp_color(c0, c1, c2, u, v, w);
                    fb.set_pixel(x, y, color);
                }
            }
        }
    }

    /// Clear the framebuffer and draw a gradient background.
    pub fn draw_sky_gradient(&self, fb: &mut HdrFramebuffer, top: HdrPixel, bottom: HdrPixel) {
        for y in 0..fb.height {
            let t = y as f32 / fb.height.max(1) as f32;
            let color = [
                top[0] * (1.0 - t) + bottom[0] * t,
                top[1] * (1.0 - t) + bottom[1] * t,
                top[2] * (1.0 - t) + bottom[2] * t,
                1.0,
            ];
            for x in 0..fb.width {
                fb.set_pixel(x, y, color);
            }
        }
    }
}

/// Compute barycentric coordinates of `p` with respect to triangle `a, b, c`.
fn barycentric(a: [f32; 2], b: [f32; 2], c: [f32; 2], p: [f32; 2]) -> Option<[f32; 3]> {
    let v0 = [c[0] - a[0], c[1] - a[1]];
    let v1 = [b[0] - a[0], b[1] - a[1]];
    let v2 = [p[0] - a[0], p[1] - a[1]];
    let dot00 = v0[0] * v0[0] + v0[1] * v0[1];
    let dot01 = v0[0] * v1[0] + v0[1] * v1[1];
    let dot02 = v0[0] * v2[0] + v0[1] * v2[1];
    let dot11 = v1[0] * v1[0] + v1[1] * v1[1];
    let dot12 = v1[0] * v2[0] + v1[1] * v2[1];
    let inv = dot00 * dot11 - dot01 * dot01;
    if inv.abs() < 1e-10 {
        return None;
    }
    let inv = 1.0 / inv;
    let u = (dot11 * dot02 - dot01 * dot12) * inv;
    let v = (dot00 * dot12 - dot01 * dot02) * inv;
    Some([1.0 - u - v, v, u])
}

fn lerp_color(c0: HdrPixel, c1: HdrPixel, c2: HdrPixel, u: f32, v: f32, w: f32) -> HdrPixel {
    [
        c0[0] * u + c1[0] * v + c2[0] * w,
        c0[1] * u + c1[1] * v + c2[1] * w,
        c0[2] * u + c1[2] * v + c2[2] * w,
        c0[3] * u + c1[3] * v + c2[3] * w,
    ]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hdr_framebuffer_set_and_get() {
        let mut fb = HdrFramebuffer::new(4, 4);
        fb.set_pixel(1, 2, [3.0, 2.0, 1.0, 1.0]);
        assert_eq!(*fb.pixel(1, 2), [3.0, 2.0, 1.0, 1.0]);
    }

    #[test]
    fn tone_map_aces_clamps_output() {
        let mapper = ToneMapper::new(ToneMapMethod::Aces, 1.0, 2.2);
        let result = mapper.map([100.0, 100.0, 100.0, 1.0]);
        for c in &result[..3] {
            assert!(*c >= 0.0 && *c <= 1.0, "component {c} out of [0, 1]");
        }
    }

    #[test]
    fn tone_map_reinhard_identity_at_zero() {
        let mapper = ToneMapper::new(ToneMapMethod::Reinhard, 1.0, 2.2);
        let result = mapper.map([0.0, 0.0, 0.0, 1.0]);
        assert!(result[0].abs() < 1e-6 && result[1].abs() < 1e-6 && result[2].abs() < 1e-6);
    }

    #[test]
    fn to_ldr_produces_correct_length() {
        let fb = HdrFramebuffer::new(8, 6);
        let mapper = ToneMapper::new(ToneMapMethod::Filmic, 1.0, 2.2);
        let ldr = fb.to_ldr(&mapper);
        assert_eq!(ldr.len(), 8 * 6 * 4);
    }

    #[test]
    fn add_pixel_accumulates() {
        let mut fb = HdrFramebuffer::new(2, 2);
        fb.set_pixel(0, 0, [1.0, 0.5, 0.25, 1.0]);
        fb.add_pixel(0, 0, [1.0, 0.5, 0.25, 0.0]);
        let p = fb.pixel(0, 0);
        assert!((p[0] - 2.0).abs() < 1e-6);
        assert!((p[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn hdr_rasterizer_draws_triangle() {
        let mut fb = HdrFramebuffer::new(10, 10);
        let r = HdrRasterizer::new();
        r.draw_triangle(
            &mut fb,
            [[2.0, 2.0], [8.0, 2.0], [5.0, 8.0]],
            [
                [2.0, 0.0, 0.0, 1.0],
                [0.0, 2.0, 0.0, 1.0],
                [0.0, 0.0, 2.0, 1.0],
            ],
        );
        // Centre pixel should be inside triangle
        let centre = fb.pixel(5, 4);
        let lum = centre[0] + centre[1] + centre[2];
        assert!(lum > 0.01, "centre of triangle should have been rasterized");
    }

    #[test]
    fn linear_to_srgb_roundtrip() {
        for byte in [0u8, 64, 128, 192, 255] {
            let linear = srgb_to_linear(byte);
            let back = linear_to_srgb(linear);
            assert!(
                (back as i32 - byte as i32).abs() <= 1,
                "roundtrip failed for byte {byte}"
            );
        }
    }

    #[test]
    fn all_tone_mappers_clamp_output() {
        let methods = [
            ToneMapMethod::Linear,
            ToneMapMethod::Exposure,
            ToneMapMethod::Reinhard,
            ToneMapMethod::ReinhardExtended { l_white: 4.0 },
            ToneMapMethod::Aces,
            ToneMapMethod::Filmic,
            ToneMapMethod::AgX,
        ];
        for method in methods {
            let mapper = ToneMapper::new(method, 1.0, 2.2);
            let px = [10.0f32, 5.0, 0.5, 1.0];
            let out = mapper.map(px);
            for c in &out[..3] {
                assert!(
                    *c >= 0.0 && *c <= 1.01,
                    "{method:?}: component {c} out of [0, 1]"
                );
            }
        }
    }
}
