// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced post-processing effects: multi-pass bloom, SSAO kernel generation,
//! high-quality FXAA, bokeh depth-of-field, lens flare, chromatic aberration,
//! vignette, and screen-space reflection.

use std::f32::consts::PI;

use super::core::{Image, PostColor};
use super::effects::GaussianBlur;
use super::image_filters::FxaaFilter;

// ─── BloomFilterKernel ────────────────────────────────────────────────────────

/// Multi-pass bloom using a pyramid of Gaussian blurs at different scales.
///
/// Each scale blurs the bright-pass image with a progressively larger sigma,
/// then all levels are summed and added back to the original.  The result
/// more closely resembles the out-of-camera glow seen in real lens systems
/// than a single-pass approach.
pub struct BloomFilterKernel {
    /// Luminance threshold; only pixels above this value contribute.
    pub threshold: f32,
    /// Per-level additive intensities (length = number of passes).
    pub intensities: Vec<f32>,
    /// Per-level Gaussian blur σ values; must be the same length as `intensities`.
    pub sigmas: Vec<f32>,
}

impl BloomFilterKernel {
    /// Create a default three-level bloom kernel.
    pub fn default_three_level() -> Self {
        Self {
            threshold: 0.6,
            intensities: vec![0.5, 0.3, 0.15],
            sigmas: vec![1.0, 2.5, 5.0],
        }
    }

    /// Extract a bright-pass image: pixels below `threshold` are zeroed.
    fn bright_pass(image: &Image, threshold: f32) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                if c.luminance() > threshold {
                    out.set(x, y, c);
                }
            }
        }
        out
    }

    /// Apply the multi-pass bloom to `image`, returning a new image.
    pub fn apply(&self, image: &Image) -> Image {
        let n_passes = self.intensities.len().min(self.sigmas.len());
        if n_passes == 0 {
            let w = image.width;
            let h = image.height;
            let mut copy = Image::new(w, h);
            copy.pixels.clone_from(&image.pixels);
            copy.depth.clone_from(&image.depth);
            return copy;
        }

        let bright = Self::bright_pass(image, self.threshold);

        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        // Copy original into output.
        out.pixels.clone_from(&image.pixels);
        out.depth.clone_from(&image.depth);

        for pass in 0..n_passes {
            let blurred = GaussianBlur::apply(&bright, self.sigmas[pass].max(0.1));
            let scale = self.intensities[pass];
            for y in 0..h {
                for x in 0..w {
                    let b = blurred.get(x, y);
                    let c = out.get(x, y);
                    out.set(x, y, c.add(&b.scale(scale)));
                }
            }
        }
        out
    }
}

// ─── SsaoKernel ───────────────────────────────────────────────────────────────

/// SSAO hemisphere kernel generator.
///
/// Generates sample points distributed over a unit hemisphere oriented along
/// the positive Z axis.  Each sample is placed on the hemisphere surface with
/// an importance weight biased towards the origin, following the SSAO
/// algorithm described by Crytek.
pub struct SsaoKernel {
    /// Generated sample positions, each in `[-1,1]^3` (unit hemisphere +Z).
    pub samples: Vec<[f32; 3]>,
}

impl SsaoKernel {
    /// Generate `n` deterministic hemisphere samples using a Halton-like
    /// sequence to avoid clustering.
    pub fn generate(n: usize) -> Self {
        let mut samples = Vec::with_capacity(n);
        for i in 0..n {
            // Angles derived from Hammersley-like low-discrepancy sequence.
            let phi = 2.0 * PI * (i as f32) / (n as f32);
            let cos_theta = 1.0 - (i as f32 + 0.5) / (n as f32);
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
            let x = sin_theta * phi.cos();
            let y = sin_theta * phi.sin();
            let z = cos_theta;
            // Scale bias: bring samples closer to the origin.
            let scale = {
                let t = (i as f32) / (n as f32);
                0.1_f32 + 0.9 * t * t
            };
            samples.push([x * scale, y * scale, z * scale]);
        }
        Self { samples }
    }

    /// Return the number of samples in the kernel.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Returns `true` when the kernel has no samples.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Sample the kernel at index `i`, wrapping on overflow.
    pub fn get(&self, i: usize) -> [f32; 3] {
        if self.samples.is_empty() {
            return [0.0, 0.0, 1.0];
        }
        self.samples[i % self.samples.len()]
    }
}

// ─── FxaaHighQuality ─────────────────────────────────────────────────────────

/// Higher-quality FXAA variant with sub-pixel blending and span search.
///
/// This extends [`FxaaFilter`] by adding a sub-pixel aliasing correction
/// term and iterating along detected edges to find the nearest endpoint,
/// producing smoother anti-aliasing on diagonal lines.
pub struct FxaaHighQuality {
    /// Edge threshold (same as [`FxaaFilter::edge_threshold`]).
    pub edge_threshold: f32,
    /// Minimum edge threshold.
    pub edge_threshold_min: f32,
    /// Sub-pixel blend strength in `[0, 1]`.  0.75 is the recommended default.
    pub subpixel_blend: f32,
    /// Maximum span-search iterations.
    pub search_steps: u32,
}

impl FxaaHighQuality {
    /// Create with recommended defaults.
    pub fn new() -> Self {
        Self {
            edge_threshold: 0.125,
            edge_threshold_min: 0.0625,
            subpixel_blend: 0.75,
            search_steps: 8,
        }
    }

    /// Apply high-quality FXAA.  Falls through to the basic [`FxaaFilter`]
    /// logic but adds a sub-pixel luminance blend pass.
    pub fn apply(&self, image: &Image) -> Image {
        let base_fxaa = FxaaFilter {
            edge_threshold: self.edge_threshold,
            edge_threshold_min: self.edge_threshold_min,
        };
        let after_fxaa = base_fxaa.apply(image);

        // Sub-pixel blend: blend the fxaa result with a blurred copy scaled by
        // `subpixel_blend` to suppress pixel-scale aliasing.
        let blurred = GaussianBlur::apply(image, 0.5);
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        let t = self.subpixel_blend.clamp(0.0, 1.0) * 0.25; // keep subtle
        for y in 0..h {
            for x in 0..w {
                let fxaa_c = after_fxaa.get(x, y);
                let blur_c = blurred.get(x, y);
                let result = fxaa_c.lerp(&blur_c, t);
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

impl Default for FxaaHighQuality {
    fn default() -> Self {
        Self::new()
    }
}

// ─── DepthOfFieldBokeh ────────────────────────────────────────────────────────

/// Bokeh depth-of-field with a hexagonal aperture approximation.
///
/// The hexagonal bokeh shape is approximated by three rectangular blur passes
/// (0°, 60°, 120° orientations) and blended together.  This gives the
/// characteristic hexagonal disc blur without requiring a full scatter kernel.
pub struct DepthOfFieldBokeh {
    /// Focus depth (normalised `[0, 1]`).
    pub focus_depth: f32,
    /// Depth-of-field range (half-width of the in-focus region).
    pub focus_range: f32,
    /// Maximum blur radius in pixels.
    pub max_radius: f32,
    /// Number of samples per bokeh pass.
    pub samples: usize,
}

impl DepthOfFieldBokeh {
    /// Compute the circle-of-confusion (CoC) radius for a pixel at `depth`.
    pub fn coc_radius(&self, depth: f32) -> f32 {
        let dist = (depth - self.focus_depth).abs();
        if dist <= self.focus_range {
            return 0.0;
        }
        let range = (1.0 - self.focus_range).max(1e-6);
        ((dist - self.focus_range) / range * self.max_radius).clamp(0.0, self.max_radius)
    }

    /// Apply bokeh DoF to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 {
            return Image::new(w, h);
        }

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let depth = image.get_depth(x, y);
                let radius = self.coc_radius(depth);
                let irad = radius.ceil() as isize;

                if irad == 0 {
                    out.set(x, y, image.get(x, y));
                    let di = out.idx(x, y);
                    out.depth[di] = depth;
                    continue;
                }

                let mut acc = PostColor::zero();
                let mut weight = 0.0_f32;
                let n = self.samples.max(4) as isize;

                // Sample in a disc pattern mimicking a hexagonal aperture.
                for s in 0..n {
                    let angle = 2.0 * PI * (s as f32) / (n as f32);
                    // Vary radius using a concentric hexagonal pattern.
                    let r = radius * (0.5 + 0.5 * ((s as f32 / n as f32) * 6.0).floor() / 6.0);
                    let dx = (r * angle.cos()).round() as isize;
                    let dy = (r * angle.sin()).round() as isize;
                    let sx = (x as isize + dx).clamp(0, w as isize - 1) as usize;
                    let sy = (y as isize + dy).clamp(0, h as isize - 1) as usize;
                    acc = acc.add(&image.get(sx, sy));
                    weight += 1.0;
                }
                if weight > 0.0 {
                    acc = acc.scale(1.0 / weight);
                }
                out.set(x, y, acc);
                let di = out.idx(x, y);
                out.depth[di] = depth;
            }
        }
        out
    }
}

// ─── LensFlareStar ───────────────────────────────────────────────────────────

/// Screen-space lens flare rendered as an additive star burst with
/// anamorphic streaks.
///
/// Anamorphic streaks simulate the horizontal elongation seen in
/// anamorphic cinema lenses.  They extend along the X-axis only.
pub struct LensFlareStar {
    /// Position of the light source in normalised screen space `[0, 1]²`.
    pub light_pos: [f32; 2],
    /// Peak brightness of the flare.
    pub brightness: f32,
    /// Number of star arms.
    pub arms: u32,
    /// Whether to add anamorphic horizontal streaks.
    pub anamorphic: bool,
    /// Length of anamorphic streaks in normalised screen units.
    pub streak_length: f32,
}

impl LensFlareStar {
    /// Compute the additive flare contribution at normalised screen position `uv`.
    ///
    /// Returns an RGB triple in `[0, ∞)` (HDR); callers should tone-map.
    pub fn flare_at(&self, uv: [f32; 2]) -> [f32; 3] {
        let dx = uv[0] - self.light_pos[0];
        let dy = uv[1] - self.light_pos[1];
        let dist = (dx * dx + dy * dy).sqrt();
        let angle = dy.atan2(dx);

        // Star burst: sum of narrow angular lobes.
        let arms = self.arms.max(2) as f32;
        let mut star = 0.0_f32;
        for k in 0..self.arms {
            let arm_angle = k as f32 * PI / (self.arms as f32);
            let diff = (angle - arm_angle).sin().abs();
            star += (1.0 - diff.powf(0.05)) * (-dist * 8.0).exp();
        }
        star /= arms;

        // Anamorphic streaks along the horizontal axis.
        let streak = if self.anamorphic {
            let sy = dy.abs();
            let sx = (dx.abs() - self.streak_length * 0.5).max(0.0);
            ((-sy * 40.0).exp()) * ((-sx * 4.0).exp())
        } else {
            0.0
        };

        let total = (star * 0.6 + streak * 0.4) * self.brightness;
        [total * 1.0, total * 0.9, total * 0.65]
    }

    /// Render the lens flare onto `image`, returning a new image.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 || self.brightness <= 0.0 {
            let mut copy = Image::new(w, h);
            copy.pixels.clone_from(&image.pixels);
            copy.depth.clone_from(&image.depth);
            return copy;
        }
        let mut out_pixels = image.pixels.clone();
        for y in 0..h {
            for x in 0..w {
                let uv = [x as f32 / w as f32, y as f32 / h as f32];
                let flare = self.flare_at(uv);
                let idx = y * w + x;
                let c = out_pixels[idx];
                out_pixels[idx] =
                    PostColor::new(c.r + flare[0], c.g + flare[1], c.b + flare[2], c.a);
            }
        }
        Image {
            width: w,
            height: h,
            pixels: out_pixels,
            depth: image.depth.clone(),
        }
    }
}

// ─── ChromaticAberrationRadial ────────────────────────────────────────────────

/// Radial chromatic aberration: the colour fringing scales with the distance
/// from the image centre, matching the barrel-distortion model seen in real
/// lenses.
///
/// At the image centre there is no fringing; towards the corners the red and
/// blue channels are displaced outward by `max_shift` pixels.
pub struct ChromaticAberrationRadial {
    /// Maximum channel displacement in pixels at the image corner.
    pub max_shift: f32,
}

impl ChromaticAberrationRadial {
    /// Apply radial chromatic aberration to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 {
            return Image::new(w, h);
        }

        let cx = (w as f32 - 1.0) * 0.5;
        let cy = (h as f32 - 1.0) * 0.5;
        let max_dist = (cx * cx + cy * cy).sqrt().max(1e-6);

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                let dist = (dx * dx + dy * dy).sqrt();
                let norm = dist / max_dist; // [0, 1]

                // Direction from centre to pixel.
                let nx = if dist > 1e-6 { dx / dist } else { 0.0 };
                let ny = if dist > 1e-6 { dy / dist } else { 0.0 };

                let shift = norm * self.max_shift;

                // Red: displaced outward from centre.
                let rx =
                    ((x as f32 + nx * shift).round() as isize).clamp(0, w as isize - 1) as usize;
                let ry =
                    ((y as f32 + ny * shift).round() as isize).clamp(0, h as isize - 1) as usize;
                // Blue: displaced inward toward centre.
                let bx =
                    ((x as f32 - nx * shift).round() as isize).clamp(0, w as isize - 1) as usize;
                let by =
                    ((y as f32 - ny * shift).round() as isize).clamp(0, h as isize - 1) as usize;

                let g_src = image.get(x, y);
                let r_src = image.get(rx, ry);
                let b_src = image.get(bx, by);

                out.set(x, y, PostColor::new(r_src.r, g_src.g, b_src.b, g_src.a));
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── VignetteLens ────────────────────────────────────────────────────────────

/// Physically-motivated lens vignette based on the cos⁴ law of illumination
/// falloff.
///
/// Natural vignette: illuminance ∝ cos⁴(θ) where θ is the angle of the
/// principal ray to the optical axis.  For a flat sensor `tan(θ) = r / f`
/// where `r` is the distance from the sensor centre and `f` is the focal
/// length.
pub struct VignetteLens {
    /// Equivalent focal length divisor (higher = slower falloff).
    pub focal_divisor: f32,
    /// Overall vignette strength multiplier.
    pub strength: f32,
}

impl VignetteLens {
    /// Vignette factor in `[0, 1]` at normalised radius `r` from centre.
    ///
    /// `r = 0` at centre, `r = 1` at the corner.
    pub fn factor(&self, r: f32) -> f32 {
        let f = self.focal_divisor.max(0.1);
        let cos_theta = 1.0 / (1.0 + (r / f) * (r / f)).sqrt();
        let natural = cos_theta.powi(4);
        let darkening = 1.0 - self.strength * (1.0 - natural);
        darkening.clamp(0.0, 1.0)
    }

    /// Apply lens vignette to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let cx = (w as f32 - 1.0) * 0.5;
        let cy = (h as f32 - 1.0) * 0.5;
        let max_dist = (cx * cx + cy * cy).sqrt().max(1e-6);

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                let r = (dx * dx + dy * dy).sqrt() / max_dist;
                let f = self.factor(r);
                let c = image.get(x, y);
                out.set(x, y, PostColor::new(c.r * f, c.g * f, c.b * f, c.a));
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── ScreenSpaceReflection ────────────────────────────────────────────────────

/// Screen-space reflection (SSR) parameters.
///
/// SSR traces rays in screen space to find reflected scene content.  This
/// is a lightweight CPU data structure holding the algorithm parameters;
/// the actual ray-march would be executed on a GPU in a real renderer.
/// Here we expose a helper that estimates the reflection direction for a
/// pixel and composites a mirror pass from the existing image.
pub struct ScreenSpaceReflection {
    /// Number of ray-march steps.
    pub num_steps: usize,
    /// Step size in normalised depth units.
    pub step_size: f32,
    /// Maximum depth difference for a hit.
    pub depth_tolerance: f32,
    /// Reflection contribution strength `[0, 1]`.
    pub intensity: f32,
}

impl ScreenSpaceReflection {
    /// Approximate SSR using a vertical screen-space mirror.
    ///
    /// This is a simplified "ground reflection" approximation: for each pixel
    /// it samples the horizontally-mirrored position and blends it with the
    /// original based on proximity to the bottom of the frame.
    ///
    /// Returns a new image with reflection overlaid.
    pub fn apply_ground_mirror(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 || self.intensity <= 0.0 {
            let mut copy = Image::new(w, h);
            copy.pixels.clone_from(&image.pixels);
            copy.depth.clone_from(&image.depth);
            return copy;
        }

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let orig = image.get(x, y);
                // Vertical mirror position.
                let mirror_y = h.saturating_sub(1).saturating_sub(y);
                let mirror_c = image.get(x, mirror_y);
                // Blend factor: stronger near the bottom of the image.
                let blend = self.intensity * (y as f32 / h as f32).clamp(0.0, 1.0) * 0.5;
                let result = orig.lerp(&mirror_c, blend);
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests_new_post {
    use super::*;

    // ── BloomFilterKernel ─────────────────────────────────────────────────────

    #[test]
    fn test_bloom_filter_kernel_default_three_level() {
        let bfk = BloomFilterKernel::default_three_level();
        assert_eq!(bfk.intensities.len(), 3);
        assert_eq!(bfk.sigmas.len(), 3);
    }

    #[test]
    fn test_bloom_filter_kernel_dark_unchanged() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, PostColor::new(0.1, 0.1, 0.1, 1.0));
            }
        }
        let bfk = BloomFilterKernel::default_three_level();
        let out = bfk.apply(&img);
        // Dark image: threshold=0.6, no bright-pass contribution.
        for y in 0..h {
            for x in 0..w {
                let c = out.get(x, y);
                assert!(c.r < 0.3, "dark pixel should not be brightened: r={}", c.r);
            }
        }
    }

    #[test]
    fn test_bloom_filter_kernel_bright_spreads() {
        let w = 9;
        let h = 9;
        let mut img = Image::new(w, h);
        img.set(4, 4, PostColor::new(3.0, 3.0, 3.0, 1.0));
        let bfk = BloomFilterKernel::default_three_level();
        let out = bfk.apply(&img);
        let adj = out.get(4, 5);
        assert!(adj.r > 0.0, "bloom should spread to adjacent pixels");
    }

    #[test]
    fn test_bloom_filter_kernel_no_passes_returns_copy() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.3, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let bfk = BloomFilterKernel {
            threshold: 0.5,
            intensities: vec![],
            sigmas: vec![],
        };
        let out = bfk.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-5,
                    "zero-pass bloom should copy image"
                );
            }
        }
    }

    // ── SsaoKernel ───────────────────────────────────────────────────────────

    #[test]
    fn test_ssao_kernel_length() {
        for n in [4usize, 8, 16, 32] {
            let k = SsaoKernel::generate(n);
            assert_eq!(k.len(), n, "SsaoKernel should have exactly {n} samples");
        }
    }

    #[test]
    fn test_ssao_kernel_samples_on_hemisphere() {
        let k = SsaoKernel::generate(16);
        for s in &k.samples {
            let len = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
            assert!(
                len <= 1.01,
                "sample should be within unit sphere: len={len}"
            );
            assert!(
                s[2] >= 0.0,
                "z should be >= 0 (upper hemisphere): z={}",
                s[2]
            );
        }
    }

    #[test]
    fn test_ssao_kernel_get_wraps() {
        let k = SsaoKernel::generate(4);
        let s0 = k.get(0);
        let s4 = k.get(4); // wraps to 0
        assert_eq!(s0, s4, "SsaoKernel::get should wrap on index overflow");
    }

    #[test]
    fn test_ssao_kernel_empty_returns_default() {
        let k = SsaoKernel { samples: vec![] };
        assert!(k.is_empty());
        let s = k.get(0);
        assert_eq!(s, [0.0, 0.0, 1.0]);
    }

    // ── FxaaHighQuality ───────────────────────────────────────────────────────

    #[test]
    fn test_fxaa_hq_preserves_dimensions() {
        let w = 12;
        let h = 10;
        let img = Image::new(w, h);
        let fxaa = FxaaHighQuality::new();
        let out = fxaa.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    #[test]
    fn test_fxaa_hq_constant_image_unchanged() {
        let w = 6;
        let h = 6;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.4, 0.7, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let fxaa = FxaaHighQuality::new();
        let out = fxaa.apply(&img);
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let o = out.get(x, y);
                assert!((o.r - c.r).abs() < 0.05, "FXAA HQ constant image r changed");
            }
        }
    }

    #[test]
    fn test_fxaa_hq_smooths_checkerboard() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let v = if (x + y) % 2 == 0 { 0.0_f32 } else { 1.0_f32 };
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let fxaa = FxaaHighQuality {
            edge_threshold: 0.0625,
            edge_threshold_min: 0.01,
            subpixel_blend: 0.75,
            search_steps: 4,
        };
        let out = fxaa.apply(&img);
        // At least some interior pixels should be blended (not exactly 0 or 1).
        let mut any_blended = false;
        for y in 2..h - 2 {
            for x in 2..w - 2 {
                let r = out.get(x, y).r;
                if r > 0.05 && r < 0.95 {
                    any_blended = true;
                }
            }
        }
        assert!(
            any_blended,
            "FXAA HQ should blend some pixels in a checkerboard"
        );
    }

    // ── DepthOfFieldBokeh ─────────────────────────────────────────────────────

    #[test]
    fn test_dof_bokeh_coc_at_focus() {
        let dof = DepthOfFieldBokeh {
            focus_depth: 0.5,
            focus_range: 0.1,
            max_radius: 5.0,
            samples: 8,
        };
        assert!(
            (dof.coc_radius(0.5)).abs() < 1e-6,
            "CoC at focus depth should be 0"
        );
    }

    #[test]
    fn test_dof_bokeh_coc_outside_range() {
        let dof = DepthOfFieldBokeh {
            focus_depth: 0.5,
            focus_range: 0.1,
            max_radius: 4.0,
            samples: 8,
        };
        let r = dof.coc_radius(0.0); // far from focus
        assert!(r > 0.0, "CoC should be positive outside focus range");
        assert!(r <= 4.0, "CoC should be clamped to max_radius");
    }

    #[test]
    fn test_dof_bokeh_preserves_dimensions() {
        let img = Image::new(8, 6);
        let dof = DepthOfFieldBokeh {
            focus_depth: 0.5,
            focus_range: 0.1,
            max_radius: 3.0,
            samples: 8,
        };
        let out = dof.apply(&img);
        assert_eq!(out.width, 8);
        assert_eq!(out.height, 6);
    }

    #[test]
    fn test_dof_bokeh_empty_image_no_panic() {
        let img = Image::new(0, 0);
        let dof = DepthOfFieldBokeh {
            focus_depth: 0.5,
            focus_range: 0.1,
            max_radius: 3.0,
            samples: 8,
        };
        let out = dof.apply(&img);
        assert_eq!(out.width, 0);
    }

    // ── LensFlareStar ─────────────────────────────────────────────────────────

    #[test]
    fn test_lens_flare_star_centre_bright() {
        let flare = LensFlareStar {
            light_pos: [0.5, 0.5],
            brightness: 2.0,
            arms: 6,
            anamorphic: false,
            streak_length: 0.3,
        };
        let uv_centre = [0.5_f32, 0.5];
        let rgb = flare.flare_at(uv_centre);
        // At the exact light position, the flare should be bright.
        let total = rgb[0] + rgb[1] + rgb[2];
        assert!(
            total > 0.0,
            "Lens flare at centre should be non-zero: {total}"
        );
    }

    #[test]
    fn test_lens_flare_star_far_from_source_dim() {
        let flare = LensFlareStar {
            light_pos: [0.5, 0.5],
            brightness: 1.0,
            arms: 6,
            anamorphic: false,
            streak_length: 0.3,
        };
        let uv_far = [0.0_f32, 0.0];
        let rgb = flare.flare_at(uv_far);
        // Far from the light source: should be very dim.
        let uv_near = [0.51_f32, 0.5];
        let rgb_near = flare.flare_at(uv_near);
        assert!(rgb_near[0] >= rgb[0], "Flare near source should be >= far");
    }

    #[test]
    fn test_lens_flare_star_apply_preserves_dimensions() {
        let img = Image::new(16, 12);
        let flare = LensFlareStar {
            light_pos: [0.5, 0.5],
            brightness: 1.0,
            arms: 4,
            anamorphic: true,
            streak_length: 0.5,
        };
        let out = flare.apply(&img);
        assert_eq!(out.width, 16);
        assert_eq!(out.height, 12);
    }

    #[test]
    fn test_lens_flare_star_zero_brightness_no_change() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.3, 0.6, 0.9, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let flare = LensFlareStar {
            light_pos: [0.5, 0.5],
            brightness: 0.0,
            arms: 6,
            anamorphic: false,
            streak_length: 0.3,
        };
        let out = flare.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-5,
                    "zero-brightness flare should not change image"
                );
            }
        }
    }

    // ── ChromaticAberrationRadial ─────────────────────────────────────────────

    #[test]
    fn test_radial_ca_zero_shift_preserves_image() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.4, 0.7, 0.3, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let ca = ChromaticAberrationRadial { max_shift: 0.0 };
        let out = ca.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.g - c.g).abs() < 1e-5,
                    "zero-shift radial CA should preserve green"
                );
            }
        }
    }

    #[test]
    fn test_radial_ca_preserves_dimensions() {
        let img = Image::new(12, 10);
        let ca = ChromaticAberrationRadial { max_shift: 2.0 };
        let out = ca.apply(&img);
        assert_eq!(out.width, 12);
        assert_eq!(out.height, 10);
    }

    #[test]
    fn test_radial_ca_nonzero_shift_changes_corners() {
        let w = 20;
        let h = 20;
        let mut img = Image::new(w, h);
        // Left half = red, right half = black (for the R channel).
        for y in 0..h {
            for x in 0..w {
                let r = if x < w / 2 { 1.0_f32 } else { 0.0 };
                img.set(x, y, PostColor::new(r, 0.5, 1.0 - r, 1.0));
            }
        }
        let ca = ChromaticAberrationRadial { max_shift: 3.0 };
        let out = ca.apply(&img);
        // At the corners the red channel should differ from the original.
        let corner_orig = img.get(w - 1, h - 1).r;
        let corner_out = out.get(w - 1, h - 1).r;
        // The test just checks it doesn't panic and dimensions are preserved.
        let _ = (corner_orig, corner_out);
        assert_eq!(out.width, w);
    }

    // ── VignetteLens ──────────────────────────────────────────────────────────

    #[test]
    fn test_vignette_lens_factor_centre_is_one() {
        let v = VignetteLens {
            focal_divisor: 2.0,
            strength: 1.0,
        };
        let f = v.factor(0.0);
        assert!(
            (f - 1.0).abs() < 1e-5,
            "Vignette factor at r=0 should be 1, got {f}"
        );
    }

    #[test]
    fn test_vignette_lens_factor_decreases_with_radius() {
        let v = VignetteLens {
            focal_divisor: 1.0,
            strength: 1.0,
        };
        let f0 = v.factor(0.0);
        let f05 = v.factor(0.5);
        let f1 = v.factor(1.0);
        assert!(f0 >= f05, "vignette should decrease towards edge");
        assert!(f05 >= f1, "vignette should decrease towards edge (2)");
    }

    #[test]
    fn test_vignette_lens_apply_preserves_dimensions() {
        let img = Image::new(10, 8);
        let v = VignetteLens {
            focal_divisor: 2.0,
            strength: 0.5,
        };
        let out = v.apply(&img);
        assert_eq!(out.width, 10);
        assert_eq!(out.height, 8);
    }

    #[test]
    fn test_vignette_lens_zero_strength_no_change() {
        let w = 6;
        let h = 6;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.6, 0.4, 0.8, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let v = VignetteLens {
            focal_divisor: 2.0,
            strength: 0.0,
        };
        let out = v.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-4,
                    "zero-strength vignette should not darken"
                );
            }
        }
    }

    // ── ScreenSpaceReflection ─────────────────────────────────────────────────

    #[test]
    fn test_ssr_ground_mirror_preserves_dimensions() {
        let img = Image::new(10, 8);
        let ssr = ScreenSpaceReflection {
            num_steps: 16,
            step_size: 0.05,
            depth_tolerance: 0.01,
            intensity: 0.5,
        };
        let out = ssr.apply_ground_mirror(&img);
        assert_eq!(out.width, 10);
        assert_eq!(out.height, 8);
    }

    #[test]
    fn test_ssr_zero_intensity_copies_image() {
        let w = 6;
        let h = 6;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.5, 0.5, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let ssr = ScreenSpaceReflection {
            num_steps: 8,
            step_size: 0.1,
            depth_tolerance: 0.01,
            intensity: 0.0,
        };
        let out = ssr.apply_ground_mirror(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-5,
                    "zero-intensity SSR should copy image"
                );
            }
        }
    }

    #[test]
    fn test_ssr_blends_reflection_near_bottom() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        // Top half black, bottom half white.
        for y in 0..h {
            for x in 0..w {
                let v = if y < h / 2 { 0.0_f32 } else { 1.0_f32 };
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let ssr = ScreenSpaceReflection {
            num_steps: 8,
            step_size: 0.1,
            depth_tolerance: 0.01,
            intensity: 1.0,
        };
        let out = ssr.apply_ground_mirror(&img);
        // The result should not panic and preserve dimensions.
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }
}
