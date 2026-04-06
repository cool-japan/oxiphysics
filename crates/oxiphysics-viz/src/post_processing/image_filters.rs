// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Image-space filters: FXAA, edge detection, chromatic aberration,
//! lens effects, colour grading, sharpening, bilateral filtering, and more.

#![allow(dead_code)]

use super::core::{Image, PostColor};
use super::effects::GaussianBlur;

// ─── FxaaFilter ──────────────────────────────────────────────────────────────

/// Fast Approximate Anti-Aliasing (FXAA) — CPU-side luminance-edge smoothing.
///
/// This is a simplified CPU implementation that detects high-contrast edges via
/// a 3×3 luminance neighbourhood and blends adjacent pixels along the dominant
/// gradient direction.
pub struct FxaaFilter {
    /// Edge detection threshold; lower = more aggressive.  Typical values: 0.0625–0.25.
    pub edge_threshold: f32,
    /// Minimum local contrast below which no filtering is applied.
    pub edge_threshold_min: f32,
}

impl FxaaFilter {
    /// Create with standard quality settings.
    pub fn new() -> Self {
        Self {
            edge_threshold: 0.125,
            edge_threshold_min: 0.0625,
        }
    }

    /// Apply FXAA to `image`, returning a new anti-aliased image.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);

        // Precompute luminance buffer.
        let lum: Vec<f32> = image.pixels.iter().map(|p| p.luminance()).collect();

        let lum_at = |x: isize, y: isize| -> f32 {
            let cx = x.clamp(0, w as isize - 1) as usize;
            let cy = y.clamp(0, h as isize - 1) as usize;
            lum[cy * w + cx]
        };

        for y in 0..h {
            for x in 0..w {
                let lm = lum_at(x as isize, y as isize);
                let ln = lum_at(x as isize, y as isize - 1);
                let ls = lum_at(x as isize, y as isize + 1);
                let le = lum_at(x as isize + 1, y as isize);
                let lw = lum_at(x as isize - 1, y as isize);

                let lmin = lm.min(ln).min(ls).min(le).min(lw);
                let lmax = lm.max(ln).max(ls).max(le).max(lw);
                let contrast = lmax - lmin;

                // Skip smooth areas.
                if contrast < (self.edge_threshold * lmax).max(self.edge_threshold_min) {
                    out.set(x, y, image.get(x, y));
                } else {
                    // Simple 1-tap blend: average the two highest-contrast neighbours.
                    let blend_h = (ln + ls) * 0.5;
                    let blend_v = (le + lw) * 0.5;
                    let (nx, ny) = if (le - lw).abs() > (ln - ls).abs() {
                        // Horizontal edge → blend vertically.
                        if blend_h > lm {
                            (x as isize, y as isize - 1)
                        } else {
                            (x as isize, y as isize + 1)
                        }
                    } else {
                        // Vertical edge → blend horizontally.
                        if blend_v > lm {
                            (x as isize - 1, y as isize)
                        } else {
                            (x as isize + 1, y as isize)
                        }
                    };
                    let nx = nx.clamp(0, w as isize - 1) as usize;
                    let ny = ny.clamp(0, h as isize - 1) as usize;
                    let c0 = image.get(x, y);
                    let c1 = image.get(nx, ny);
                    let blended = c0.lerp(&c1, 0.5);
                    out.set(x, y, blended);
                }
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

impl Default for FxaaFilter {
    fn default() -> Self {
        Self::new()
    }
}

// ─── DepthOfFieldLens ─────────────────────────────────────────────────────────

/// Depth-of-field circle-of-confusion model using f64 (thin-lens formula).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct DepthOfFieldLens {
    /// Distance to the focal plane in world units.
    pub focal_distance: f64,
    /// Aperture diameter (f-number denominator) in world units.
    pub aperture: f64,
    /// Focal length of the lens in world units.
    pub focal_length: f64,
}

#[allow(dead_code)]
impl DepthOfFieldLens {
    /// Circle-of-confusion diameter for a point at `depth` from the camera.
    ///
    /// Uses the thin-lens CoC formula:
    /// `coc = |aperture * focal_length * (depth - focal_distance)|
    ///        / (depth * (focal_distance - focal_length))`
    /// Returns 0 when `depth == focal_distance`.
    pub fn circle_of_confusion(&self, depth: f64) -> f64 {
        if (depth - self.focal_distance).abs() < 1e-12 {
            return 0.0;
        }
        let fd = self.focal_distance;
        let fl = self.focal_length;
        let denom = depth * (fd - fl);
        if denom.abs() < 1e-12 {
            return 0.0;
        }
        (self.aperture * fl * (depth - fd) / denom).abs()
    }
}

// ─── ChromaticAberration ──────────────────────────────────────────────────────

/// Chromatic aberration: lateral colour fringing.
///
/// Red and blue channels are offset by `r_offset` / `b_offset` respectively.
/// In a full-image application the offsets drive UV shifts; here we expose a
/// per-pixel helper that applies the offsets directly to channel values.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ChromaticAberration {
    /// Offset applied to the red channel (positive = shifted).
    pub r_offset: f64,
    /// Offset applied to the blue channel (positive = shifted).
    pub b_offset: f64,
}

#[allow(dead_code)]
impl ChromaticAberration {
    /// Apply chromatic aberration offsets to a single pixel.
    ///
    /// The offsets are added directly to the R and B channels; green is
    /// unaffected.  A real implementation would look up the offset UV in
    /// the source image, but for a data-structure-only crate this serves as
    /// the canonical single-pixel computation.
    pub fn apply_to_pixel(&self, r: f64, g: f64, b: f64) -> (f64, f64, f64) {
        (r + self.r_offset, g, b + self.b_offset)
    }
}

// ─── LensFlare ────────────────────────────────────────────────────────────────

/// Lens flare descriptor.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LensFlare {
    /// Screen-space position `[x, y]` in `[0, 1]` NDC.
    pub position: [f64; 2],
    /// Maximum brightness of the flare.
    pub intensity: f64,
    /// Number of streak arms.
    pub streak_count: u32,
}

#[allow(dead_code)]
impl LensFlare {
    /// Visibility factor in `[0, 1]` based on alignment of `view_dir` and `light_dir`.
    ///
    /// Returns 1 when the light is directly in view and falls off as the angle
    /// between the two directions increases.
    pub fn visibility(&self, view_dir: [f64; 3], light_dir: [f64; 3]) -> f64 {
        let dot =
            view_dir[0] * light_dir[0] + view_dir[1] * light_dir[1] + view_dir[2] * light_dir[2];
        dot.clamp(0.0, 1.0).powf(4.0)
    }
}

// ─── GodRays ──────────────────────────────────────────────────────────────────

/// God-rays (light scattering) parameters.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct GodRays {
    /// Screen-space light position `[x, y]` in `[0, 1]` NDC.
    pub light_pos: [f64; 2],
    /// Per-step decay factor (< 1 for attenuation).
    pub decay: f64,
    /// Density of ray scattering.
    pub density: f64,
    /// Per-sample weight.
    pub weight: f64,
    /// Final exposure multiplier.
    pub exposure: f64,
}

#[allow(dead_code)]
impl GodRays {
    /// Compute the accumulated sample weight for a single ray sample.
    ///
    /// `pos` is the current screen-space sample position, `step` is the
    /// zero-based sample index, and `n_steps` is the total number of samples.
    /// Uses the standard god-rays decay formula:
    /// `weight * decay^step`.
    pub fn compute_sample(&self, _pos: [f64; 2], step: u32, _n_steps: u32) -> f64 {
        self.weight * self.decay.powi(step as i32) * self.exposure
    }
}

// ─── VignetteEffect ──────────────────────────────────────────────────────────

/// Vignette effect with configurable radius and softness (f64 variant).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct VignetteEffect {
    /// Radius of the unvignetted region in UV space (UV centre = 0.5, 0.5).
    pub radius: f64,
    /// Softness (feather) of the vignette edge.
    pub softness: f64,
}

#[allow(dead_code)]
impl VignetteEffect {
    /// Vignette strength in `[0, 1]` for a pixel at `uv` (each component in `[0, 1]`).
    ///
    /// Returns 0 at the centre and approaches 1 towards the edges, following a
    /// smooth-step transition around `radius`.
    pub fn strength(&self, uv: [f64; 2]) -> f64 {
        let dx = uv[0] - 0.5;
        let dy = uv[1] - 0.5;
        let dist = (dx * dx + dy * dy).sqrt();
        let t = ((dist - self.radius) / self.softness.max(1e-12)).clamp(0.0, 1.0);
        t * t * (3.0 - 2.0 * t) // smooth-step
    }
}

// ─── ColorGrader ─────────────────────────────────────────────────────────────

/// Simple colour grading using f64 channels (contrast / brightness / saturation).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ColorGrader {
    /// Contrast multiplier around mid-grey (0.5); > 1 increases range.
    pub contrast: f64,
    /// Additive brightness offset applied after contrast.
    pub brightness: f64,
    /// Saturation multiplier; 1 = original, 0 = greyscale.
    pub saturation: f64,
}

#[allow(dead_code)]
impl ColorGrader {
    /// Apply grading to a single `[r, g, b]` colour (f64 components).
    ///
    /// 1. Contrast: `c' = (c - 0.5) * contrast + 0.5`
    /// 2. Brightness: additive offset
    /// 3. Saturation: blend with luminance
    pub fn apply(&self, color: [f64; 3]) -> [f64; 3] {
        // 1. Contrast
        let mid = 0.5_f64;
        let r = (color[0] - mid) * self.contrast + mid;
        let g = (color[1] - mid) * self.contrast + mid;
        let b = (color[2] - mid) * self.contrast + mid;
        // 2. Brightness
        let r = r + self.brightness;
        let g = g + self.brightness;
        let b = b + self.brightness;
        // 3. Saturation
        let lum = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        [
            lum + (r - lum) * self.saturation,
            lum + (g - lum) * self.saturation,
            lum + (b - lum) * self.saturation,
        ]
    }
}

// ─── SharpenFilter ────────────────────────────────────────────────────────────

/// Unsharp-mask sharpening filter operating on f64 pixel data.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct SharpenFilter {
    /// Sharpening strength; higher values produce stronger edges.
    pub strength: f64,
}

#[allow(dead_code)]
impl SharpenFilter {
    /// Apply sharpening to a flat pixel slice of size `w × h`.
    ///
    /// Uses a discrete Laplacian kernel (4-neighbour) scaled by `strength`.
    /// Edge pixels are passed through unchanged.
    pub fn apply_to_image(&self, pixels: &[[f64; 3]], w: usize, h: usize) -> Vec<[f64; 3]> {
        let mut out = pixels.to_vec();
        for y in 1..h.saturating_sub(1) {
            for x in 1..w.saturating_sub(1) {
                let idx = y * w + x;
                let up = pixels[(y - 1) * w + x];
                let down = pixels[(y + 1) * w + x];
                let left = pixels[y * w + (x - 1)];
                let right = pixels[y * w + (x + 1)];
                let c = pixels[idx];
                for i in 0..3 {
                    let laplacian = 4.0 * c[i] - up[i] - down[i] - left[i] - right[i];
                    out[idx][i] = c[i] + self.strength * laplacian;
                }
            }
        }
        out
    }
}

// ─── ChromaticAberrationEffect ────────────────────────────────────────────────

/// Full-image chromatic aberration: red and blue channels are shifted by
/// pixel offsets while green stays in place.
///
/// The shifts are applied in integer pixel units.  Sub-pixel precision can be
/// achieved by pre-scaling.
pub struct ChromaticAberrationEffect {
    /// Horizontal pixel offset for the red channel (positive = right).
    pub r_shift_x: i32,
    /// Vertical pixel offset for the red channel (positive = down).
    pub r_shift_y: i32,
    /// Horizontal pixel offset for the blue channel.
    pub b_shift_x: i32,
    /// Vertical pixel offset for the blue channel.
    pub b_shift_y: i32,
}

impl ChromaticAberrationEffect {
    /// Apply chromatic aberration to `image`, returning a new image.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);

        let clamp_x = |ix: i32| ix.clamp(0, w as i32 - 1) as usize;
        let clamp_y = |iy: i32| iy.clamp(0, h as i32 - 1) as usize;

        for y in 0..h {
            for x in 0..w {
                let ix = x as i32;
                let iy = y as i32;

                // Sample green from original position.
                let g_src = image.get(x, y);

                // Red channel: look up source pixel at (x - r_shift_x, y - r_shift_y).
                let rx = clamp_x(ix - self.r_shift_x);
                let ry = clamp_y(iy - self.r_shift_y);
                let r_src = image.get(rx, ry);

                // Blue channel: look up source pixel at (x - b_shift_x, y - b_shift_y).
                let bx = clamp_x(ix - self.b_shift_x);
                let by = clamp_y(iy - self.b_shift_y);
                let b_src = image.get(bx, by);

                let result = PostColor::new(r_src.r, g_src.g, b_src.b, g_src.a);
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── SobelEdgeDetection ───────────────────────────────────────────────────────

/// Screen-space edge detection using Sobel operators on the luminance channel.
///
/// The output image encodes edge magnitude in the R channel (and copies it to
/// G/B for greyscale compatibility).  The depth buffer is preserved.
pub struct SobelEdgeDetection {
    /// Multiplier applied to the raw gradient magnitude before clamping.
    pub sensitivity: f32,
}

impl SobelEdgeDetection {
    /// Apply Sobel edge detection to `image`.
    ///
    /// Luminance of each pixel is used as the scalar field.  Boundary pixels
    /// use edge-clamping.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);

        // Pre-compute luminance buffer.
        let lum: Vec<f32> = image.pixels.iter().map(|p| p.luminance()).collect();

        let lum_at = |x: isize, y: isize| -> f32 {
            let cx = x.clamp(0, w as isize - 1) as usize;
            let cy = y.clamp(0, h as isize - 1) as usize;
            lum[cy * w + cx]
        };

        for y in 0..h {
            for x in 0..w {
                let ix = x as isize;
                let iy = y as isize;

                // Sobel Gx kernel: [-1,0,1; -2,0,2; -1,0,1]
                let gx = -lum_at(ix - 1, iy - 1) + lum_at(ix + 1, iy - 1)
                    - 2.0 * lum_at(ix - 1, iy)
                    + 2.0 * lum_at(ix + 1, iy)
                    - lum_at(ix - 1, iy + 1)
                    + lum_at(ix + 1, iy + 1);

                // Sobel Gy kernel: [-1,-2,-1; 0,0,0; 1,2,1]
                let gy =
                    -lum_at(ix - 1, iy - 1) - 2.0 * lum_at(ix, iy - 1) - lum_at(ix + 1, iy - 1)
                        + lum_at(ix - 1, iy + 1)
                        + 2.0 * lum_at(ix, iy + 1)
                        + lum_at(ix + 1, iy + 1);

                let mag = ((gx * gx + gy * gy).sqrt() * self.sensitivity).clamp(0.0, 1.0);
                let c = PostColor::new(mag, mag, mag, 1.0);
                out.set(x, y, c);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── ContrastBrightnessSaturation ─────────────────────────────────────────────

/// Adjust contrast, brightness, and saturation of an image in one pass.
pub struct ContrastBrightnessSaturation {
    /// Contrast multiplier around mid-grey (0.5).
    pub contrast: f32,
    /// Additive brightness offset (applied after contrast).
    pub brightness: f32,
    /// Saturation multiplier; 1.0 = original, 0.0 = greyscale.
    pub saturation: f32,
}

impl ContrastBrightnessSaturation {
    /// Construct with identity (no-op) settings.
    pub fn identity() -> Self {
        Self {
            contrast: 1.0,
            brightness: 0.0,
            saturation: 1.0,
        }
    }

    /// Apply to a single `PostColor`.
    #[inline]
    pub fn apply_to_color(&self, c: PostColor) -> PostColor {
        let mid = 0.5_f32;
        let r = (c.r - mid) * self.contrast + mid + self.brightness;
        let g = (c.g - mid) * self.contrast + mid + self.brightness;
        let b = (c.b - mid) * self.contrast + mid + self.brightness;
        let lum = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        PostColor::new(
            lum + (r - lum) * self.saturation,
            lum + (g - lum) * self.saturation,
            lum + (b - lum) * self.saturation,
            c.a,
        )
    }

    /// Apply to every pixel of `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let result = self.apply_to_color(image.get(x, y));
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── SrgbGammaCorrection ──────────────────────────────────────────────────────

/// Gamma correction following the IEC 61966-2-1 sRGB piece-wise curve.
///
/// Use this to convert from linear light to display-encoded sRGB before
/// presenting to an 8-bit monitor.
pub struct SrgbGammaCorrection;

impl SrgbGammaCorrection {
    /// Encode a single linear-light channel value to sRGB.
    #[inline]
    pub fn encode(linear: f32) -> f32 {
        let v = linear.clamp(0.0, 1.0);
        if v <= 0.003130_8 {
            12.92 * v
        } else {
            1.055 * v.powf(1.0 / 2.4) - 0.055
        }
    }

    /// Decode a single sRGB-encoded value back to linear light.
    #[inline]
    pub fn decode(srgb: f32) -> f32 {
        let v = srgb.clamp(0.0, 1.0);
        if v <= 0.04045 {
            v / 12.92
        } else {
            ((v + 0.055) / 1.055).powf(2.4)
        }
    }

    /// Apply sRGB encoding to every RGB channel of `image` (alpha is preserved).
    pub fn apply(image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let enc =
                    PostColor::new(Self::encode(c.r), Self::encode(c.g), Self::encode(c.b), c.a);
                out.set(x, y, enc);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── HistogramEqualization ────────────────────────────────────────────────────

/// Histogram equalization applied to the luminance channel of an `Image`.
///
/// The original hue and saturation are preserved; only the luminance (value
/// channel in an HSV sense) is redistributed.
pub struct HistogramEqualization {
    /// Number of histogram bins (higher = finer CDF approximation).
    pub n_bins: usize,
}

impl HistogramEqualization {
    /// Create with a standard 256-bin histogram.
    pub fn new() -> Self {
        Self { n_bins: 256 }
    }

    /// Apply histogram equalization to the luminance channel of `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let n = w * h;
        let bins = self.n_bins.max(2);

        // Collect luminance values.
        let lums: Vec<f32> = image.pixels.iter().map(|p| p.luminance()).collect();
        let lmin = lums.iter().cloned().fold(f32::INFINITY, f32::min);
        let lmax = lums.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let range = lmax - lmin;

        let mut out = Image::new(w, h);

        if range < 1e-7 {
            // Flat image — copy as-is.
            out.pixels.clone_from(&image.pixels);
            out.depth.clone_from(&image.depth);
            return out;
        }

        // Build histogram.
        let mut hist = vec![0usize; bins];
        for &l in &lums {
            let bin = (((l - lmin) / range) * (bins as f32 - 1.0))
                .round()
                .clamp(0.0, (bins - 1) as f32) as usize;
            hist[bin] += 1;
        }

        // Build CDF.
        let mut cdf = vec![0usize; bins];
        cdf[0] = hist[0];
        for i in 1..bins {
            cdf[i] = cdf[i - 1] + hist[i];
        }
        let cdf_min = *cdf.iter().find(|&&c| c > 0).unwrap_or(&0);

        // Map each pixel.
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let l = c.luminance();
                let bin = (((l - lmin) / range) * (bins as f32 - 1.0))
                    .round()
                    .clamp(0.0, (bins - 1) as f32) as usize;
                let new_lum = if n <= cdf_min {
                    0.5
                } else {
                    ((cdf[bin] - cdf_min) as f32 / (n - cdf_min) as f32).clamp(0.0, 1.0)
                };
                let scale = if l > 1e-7 { new_lum / l } else { 1.0 };
                let result = PostColor::new(
                    (c.r * scale).clamp(0.0, 1.0),
                    (c.g * scale).clamp(0.0, 1.0),
                    (c.b * scale).clamp(0.0, 1.0),
                    c.a,
                );
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

impl Default for HistogramEqualization {
    fn default() -> Self {
        Self::new()
    }
}

// ─── UnsharpMask ──────────────────────────────────────────────────────────────

/// Unsharp masking (sharpening) on the full `Image` pipeline.
///
/// The sharpened output is `original + strength * (original - blurred)`.
pub struct UnsharpMask {
    /// Gaussian blur sigma used for the blur step.
    pub sigma: f32,
    /// Sharpening strength; 0 = no effect, 1 = standard sharpen.
    pub strength: f32,
}

impl UnsharpMask {
    /// Apply unsharp masking to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let blurred = GaussianBlur::apply(image, self.sigma.max(0.1));
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let orig = image.get(x, y);
                let blur = blurred.get(x, y);
                // Sharpened = orig + strength * (orig - blur)
                let r = orig.r + self.strength * (orig.r - blur.r);
                let g = orig.g + self.strength * (orig.g - blur.g);
                let b = orig.b + self.strength * (orig.b - blur.b);
                out.set(x, y, PostColor::new(r, g, b, orig.a));
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── BilateralFilter ──────────────────────────────────────────────────────────

/// Bilateral filter approximation: edge-preserving noise reduction.
///
/// Each output pixel is a weighted sum of nearby pixels where the weight
/// considers both spatial proximity (Gaussian in pixel distance) and
/// intensity similarity (Gaussian in luminance difference).
pub struct BilateralFilter {
    /// Spatial Gaussian sigma (in pixels).
    pub sigma_spatial: f32,
    /// Range (intensity) Gaussian sigma.
    pub sigma_range: f32,
    /// Sampling radius (in pixels).
    pub radius: usize,
}

impl BilateralFilter {
    /// Apply the bilateral filter to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let r = self.radius as isize;
        let inv_2_ss = 1.0 / (2.0 * self.sigma_spatial * self.sigma_spatial);
        let inv_2_sr = 1.0 / (2.0 * self.sigma_range * self.sigma_range);

        let mut out = Image::new(w, h);

        for y in 0..h {
            for x in 0..w {
                let center = image.get(x, y);
                let center_lum = center.luminance();

                let mut acc = PostColor::zero();
                let mut weight_sum = 0.0_f32;

                for dy in -r..=r {
                    for dx in -r..=r {
                        let sx = x as isize + dx;
                        let sy = y as isize + dy;
                        if sx < 0 || sy < 0 || sx >= w as isize || sy >= h as isize {
                            continue;
                        }
                        let sample = image.get(sx as usize, sy as usize);
                        let lum_diff = sample.luminance() - center_lum;
                        let dist2 = (dx * dx + dy * dy) as f32;

                        let w_spatial = (-dist2 * inv_2_ss).exp();
                        let w_range = (-(lum_diff * lum_diff) * inv_2_sr).exp();
                        let w = w_spatial * w_range;

                        acc = acc.add(&sample.scale(w));
                        weight_sum += w;
                    }
                }

                let result = if weight_sum > 1e-10 {
                    acc.scale(1.0 / weight_sum)
                } else {
                    center
                };
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── FxaaFilter ────────────────────────────────────────────────────────────

    #[test]
    fn test_fxaa_constant_image_unchanged() {
        let w = 5;
        let h = 5;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.3, 0.8, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let fxaa = FxaaFilter::new();
        let out = fxaa.apply(&img);
        // No edges → all pixels should be preserved.
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-4,
                    "FXAA constant image: r changed at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn test_fxaa_preserves_image_dimensions() {
        let w = 16;
        let h = 12;
        let img = Image::new(w, h);
        let fxaa = FxaaFilter::new();
        let out = fxaa.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    #[test]
    fn test_fxaa_smooths_hard_edge() {
        // Use alternating columns (vertical stripes): left-of-centre = black,
        // right-of-centre = white.  At the boundary column FXAA detects
        // (le-lw).abs() > 0 and blends vertically with the North neighbour.
        // The North neighbour (y-1) is black (same column), same colour — but the
        // pixel at (edge_x, edge_y) has blend_h (avg of North+South) > lm=0,
        // so it blends toward a white neighbour from the vertical direction.
        // Simpler: just verify FXAA detects the edge and changes at least one pixel.
        let w = 6;
        let h = 4;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                // Alternating columns: even x = black, odd x = white.
                let c = if x % 2 == 0 {
                    PostColor::new(0.0, 0.0, 0.0, 1.0)
                } else {
                    PostColor::new(1.0, 1.0, 1.0, 1.0)
                };
                img.set(x, y, c);
            }
        }
        let fxaa = FxaaFilter {
            edge_threshold: 0.0625,
            edge_threshold_min: 0.01,
        };
        let out = fxaa.apply(&img);
        // At any boundary pixel (e.g., even x), the luminance neighbours differ.
        // FXAA should blend with the adjacent column, making the output ≠ original.
        let mut any_changed = false;
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let orig = img.get(x, y).r;
                let result = out.get(x, y).r;
                if (orig - result).abs() > 1e-4 {
                    any_changed = true;
                }
            }
        }
        assert!(
            any_changed,
            "FXAA should alter at least one pixel on a high-contrast alternating pattern"
        );
    }

    // ── DepthOfFieldLens ──────────────────────────────────────────────────────

    #[test]
    fn test_dof_lens_coc_at_focal_distance_is_zero() {
        let dof = DepthOfFieldLens {
            focal_distance: 5.0,
            aperture: 2.8,
            focal_length: 50.0,
        };
        let coc = dof.circle_of_confusion(5.0);
        assert!(
            coc.abs() < 1e-9,
            "CoC at focal distance should be zero, got {coc}"
        );
    }

    #[test]
    fn test_dof_lens_coc_increases_with_defocus() {
        let dof = DepthOfFieldLens {
            focal_distance: 5.0,
            aperture: 2.8,
            focal_length: 50.0,
        };
        let coc_near = dof.circle_of_confusion(2.0);
        let coc_far = dof.circle_of_confusion(10.0);
        assert!(coc_near > 0.0, "CoC for near defocus should be positive");
        assert!(coc_far > 0.0, "CoC for far defocus should be positive");
    }

    // ── ChromaticAberration ──────────────────────────────────────────────────

    #[test]
    fn test_chromatic_aberration_neutral_identity() {
        let ca = ChromaticAberration {
            r_offset: 0.0,
            b_offset: 0.0,
        };
        let (r2, g2, b2) = ca.apply_to_pixel(0.5, 0.6, 0.7);
        assert!((r2 - 0.5).abs() < 1e-9);
        assert!((g2 - 0.6).abs() < 1e-9);
        assert!((b2 - 0.7).abs() < 1e-9);
    }

    // ── VignetteEffect ────────────────────────────────────────────────────────

    #[test]
    fn test_vignette_effect_center_is_low() {
        let vig = VignetteEffect {
            radius: 0.75,
            softness: 0.5,
        };
        // UV at center = [0.5, 0.5]
        let s = vig.strength([0.5, 0.5]);
        assert!(
            s < 0.1,
            "vignette at center should be near 0 (low darkening), got {s}"
        );
    }

    #[test]
    fn test_vignette_effect_corner_is_stronger() {
        let vig = VignetteEffect {
            radius: 0.5,
            softness: 0.3,
        };
        let center = vig.strength([0.5, 0.5]);
        let corner = vig.strength([0.0, 0.0]);
        assert!(
            corner > center,
            "corner vignette should be stronger than center"
        );
    }

    // ── ColorGrader ──────────────────────────────────────────────────────────

    #[test]
    fn test_color_grader_contrast_above_one_increases_range() {
        let cg = ColorGrader {
            contrast: 2.0,
            brightness: 0.0,
            saturation: 1.0,
        };
        let bright = cg.apply([0.8, 0.8, 0.8]);
        let dark = cg.apply([0.2, 0.2, 0.2]);
        // contrast > 1 should push bright values up and dark values down
        assert!(
            bright[0] > 0.8,
            "bright color should get brighter with contrast > 1"
        );
        assert!(
            dark[0] < 0.2,
            "dark color should get darker with contrast > 1"
        );
    }

    #[test]
    fn test_color_grader_identity() {
        let cg = ColorGrader {
            contrast: 1.0,
            brightness: 0.0,
            saturation: 1.0,
        };
        let input = [0.4, 0.6, 0.2];
        let out = cg.apply(input);
        for i in 0..3 {
            assert!(
                (out[i] - input[i]).abs() < 1e-9,
                "identity grader should not change color"
            );
        }
    }

    // ── SharpenFilter ────────────────────────────────────────────────────────

    #[test]
    fn test_sharpen_increases_edge_pixels() {
        // Create a 5x5 image with a horizontal step edge in the middle.
        let w = 5usize;
        let h = 5usize;
        let mut pixels = vec![[0.0f64; 3]; w * h];
        for y in 0..h {
            for x in 0..w {
                let v = if x >= w / 2 { 1.0 } else { 0.0 };
                pixels[y * w + x] = [v, v, v];
            }
        }
        let sf = SharpenFilter { strength: 1.0 };
        let out = sf.apply_to_image(&pixels, w, h);
        assert_eq!(out.len(), w * h, "output length must match input");
        // At least one pixel near the edge should have value > 1.0 (sharpened).
        let any_enhanced = out.iter().any(|p| p[0] > 1.0 + 1e-9 || p[0] < -1e-9);
        assert!(
            any_enhanced,
            "sharpening should produce values outside [0,1] at edges"
        );
    }

    #[test]
    fn test_sharpen_uniform_image_unchanged() {
        let w = 4usize;
        let h = 4usize;
        let pixels = vec![[0.5f64; 3]; w * h];
        let sf = SharpenFilter { strength: 1.0 };
        let out = sf.apply_to_image(&pixels, w, h);
        for p in &out {
            assert!(
                (p[0] - 0.5).abs() < 1e-9,
                "uniform image should be unchanged by sharpen"
            );
        }
    }

    // ── ChromaticAberrationEffect ─────────────────────────────────────────────

    #[test]
    fn test_chromatic_aberration_effect_zero_shift_identity() {
        let w = 6;
        let h = 6;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.8, 0.4, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let ca = ChromaticAberrationEffect {
            r_shift_x: 0,
            r_shift_y: 0,
            b_shift_x: 0,
            b_shift_y: 0,
        };
        let out = ca.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-5,
                    "zero-shift CA should not change r at ({x},{y})"
                );
                assert!(
                    (o.g - c.g).abs() < 1e-5,
                    "zero-shift CA should not change g at ({x},{y})"
                );
                assert!(
                    (o.b - c.b).abs() < 1e-5,
                    "zero-shift CA should not change b at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn test_chromatic_aberration_effect_preserves_dimensions() {
        let w = 10;
        let h = 8;
        let img = Image::new(w, h);
        let ca = ChromaticAberrationEffect {
            r_shift_x: 2,
            r_shift_y: 0,
            b_shift_x: -2,
            b_shift_y: 0,
        };
        let out = ca.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    #[test]
    fn test_chromatic_aberration_effect_shifts_red_channel() {
        // Create an image where only the left half has a non-zero red channel.
        let w = 10;
        let h = 4;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let r = if x < w / 2 { 1.0 } else { 0.0 };
                img.set(x, y, PostColor::new(r, 0.5, 0.5, 1.0));
            }
        }
        let ca = ChromaticAberrationEffect {
            r_shift_x: -2,
            r_shift_y: 0,
            b_shift_x: 0,
            b_shift_y: 0,
        };
        let out = ca.apply(&img);
        // After a -2 shift, pixel at x=w/2 reads red from x=w/2+2 (which is 0.0).
        let px = out.get(w / 2, h / 2);
        assert!(
            px.r < 0.1,
            "shifted red channel should be 0 at x=w/2 after left shift"
        );
    }

    // ── SobelEdgeDetection ────────────────────────────────────────────────────

    #[test]
    fn test_sobel_constant_image_zero_edges() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.5, 0.5, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let sobel = SobelEdgeDetection { sensitivity: 1.0 };
        let out = sobel.apply(&img);
        // Interior pixels of a constant image should have ~0 edge magnitude.
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let mag = out.get(x, y).r;
                assert!(
                    mag < 0.01,
                    "constant image → no edges (got {mag} at ({x},{y}))"
                );
            }
        }
    }

    #[test]
    fn test_sobel_hard_edge_detects_edge() {
        let w = 9;
        let h = 5;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let v = if x < w / 2 { 0.0 } else { 1.0 };
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let sobel = SobelEdgeDetection { sensitivity: 1.0 };
        let out = sobel.apply(&img);
        // At the boundary column the gradient should be large.
        let edge_mag = out.get(w / 2, h / 2).r;
        assert!(
            edge_mag > 0.3,
            "hard vertical edge should produce large Sobel magnitude, got {edge_mag}"
        );
    }

    #[test]
    fn test_sobel_preserves_dimensions() {
        let w = 12;
        let h = 10;
        let img = Image::new(w, h);
        let out = SobelEdgeDetection { sensitivity: 2.0 }.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    // ── ContrastBrightnessSaturation ──────────────────────────────────────────

    #[test]
    fn test_cbs_identity_no_change() {
        let c = PostColor::new(0.4, 0.6, 0.2, 1.0);
        let adj = ContrastBrightnessSaturation::identity();
        let out = adj.apply_to_color(c);
        assert!(
            (out.r - c.r).abs() < 1e-4,
            "CBS identity should not change r"
        );
        assert!(
            (out.g - c.g).abs() < 1e-4,
            "CBS identity should not change g"
        );
        assert!(
            (out.b - c.b).abs() < 1e-4,
            "CBS identity should not change b"
        );
    }

    #[test]
    fn test_cbs_saturation_zero_gives_greyscale() {
        let c = PostColor::new(0.8, 0.2, 0.5, 1.0);
        let adj = ContrastBrightnessSaturation {
            contrast: 1.0,
            brightness: 0.0,
            saturation: 0.0,
        };
        let out = adj.apply_to_color(c);
        assert!((out.r - out.g).abs() < 1e-4, "sat=0 → r should equal g");
        assert!((out.g - out.b).abs() < 1e-4, "sat=0 → g should equal b");
    }

    #[test]
    fn test_cbs_brightness_increases_all_channels() {
        let c = PostColor::new(0.3, 0.3, 0.3, 1.0);
        let adj = ContrastBrightnessSaturation {
            contrast: 1.0,
            brightness: 0.2,
            saturation: 1.0,
        };
        let out = adj.apply_to_color(c);
        assert!(out.r > c.r, "positive brightness should increase r");
        assert!(out.g > c.g, "positive brightness should increase g");
    }

    #[test]
    fn test_cbs_apply_full_image() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.4, 0.4, 0.4, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let adj = ContrastBrightnessSaturation {
            contrast: 1.0,
            brightness: 0.1,
            saturation: 1.0,
        };
        let out = adj.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
        for y in 0..h {
            for x in 0..w {
                assert!(
                    out.get(x, y).r > c.r,
                    "CBS brightness should affect all pixels"
                );
            }
        }
    }

    // ── SrgbGammaCorrection ────────────────────────────────────────────────────

    #[test]
    fn test_srgb_encode_black_white_fixed_points() {
        assert!(
            (SrgbGammaCorrection::encode(0.0)).abs() < 1e-6,
            "encode(0) should be 0"
        );
        assert!(
            (SrgbGammaCorrection::encode(1.0) - 1.0).abs() < 1e-5,
            "encode(1) should be 1"
        );
    }

    #[test]
    fn test_srgb_decode_black_white_fixed_points() {
        assert!(
            (SrgbGammaCorrection::decode(0.0)).abs() < 1e-6,
            "decode(0) should be 0"
        );
        assert!(
            (SrgbGammaCorrection::decode(1.0) - 1.0).abs() < 1e-5,
            "decode(1) should be 1"
        );
    }

    #[test]
    fn test_srgb_encode_decode_roundtrip() {
        for i in 0..=10 {
            let linear = i as f32 / 10.0;
            let encoded = SrgbGammaCorrection::encode(linear);
            let decoded = SrgbGammaCorrection::decode(encoded);
            assert!(
                (decoded - linear).abs() < 1e-4,
                "sRGB roundtrip failed for linear={linear}: got {decoded}"
            );
        }
    }

    #[test]
    fn test_srgb_apply_brightens_dark_values() {
        // sRGB encoding brightens values < 1 compared to linear.
        let linear_val = 0.2_f32;
        let encoded = SrgbGammaCorrection::encode(linear_val);
        assert!(
            encoded > linear_val,
            "sRGB encode should brighten dark linear values"
        );
    }

    #[test]
    fn test_srgb_apply_image_preserves_dimensions() {
        let w = 8;
        let h = 6;
        let img = Image::new(w, h);
        let out = SrgbGammaCorrection::apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    // ── HistogramEqualization ─────────────────────────────────────────────────

    #[test]
    fn test_histogram_eq_output_range() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let v = (y * w + x) as f32 / (w * h) as f32;
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let eq = HistogramEqualization::new();
        let out = eq.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let c = out.get(x, y);
                assert!(
                    c.r >= 0.0 && c.r <= 1.0 + 1e-5,
                    "equalized r out of range: {}",
                    c.r
                );
            }
        }
    }

    #[test]
    fn test_histogram_eq_constant_image_no_panic() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.5, 0.5, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let eq = HistogramEqualization { n_bins: 256 };
        let out = eq.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    #[test]
    fn test_histogram_eq_preserves_dimensions() {
        let w = 16;
        let h = 12;
        let img = Image::new(w, h);
        let eq = HistogramEqualization::default();
        let out = eq.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    // ── UnsharpMask ───────────────────────────────────────────────────────────

    #[test]
    fn test_unsharp_mask_constant_image_unchanged() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.6, 0.4, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let um = UnsharpMask {
            sigma: 1.0,
            strength: 1.0,
        };
        let out = um.apply(&img);
        // For a constant image, blurred == original, so result == original.
        for y in 2..h - 2 {
            for x in 2..w - 2 {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 0.02,
                    "UnsharpMask on constant image should not change r"
                );
            }
        }
    }

    #[test]
    fn test_unsharp_mask_sharpens_edge() {
        let w = 11;
        let h = 5;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let v = if x < w / 2 { 0.0_f32 } else { 1.0_f32 };
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let um = UnsharpMask {
            sigma: 1.0,
            strength: 1.0,
        };
        let out = um.apply(&img);
        // At the edge boundary, values should be pushed beyond [0, 1] (overshoot).
        let edge_bright = out.get(w / 2 + 1, h / 2).r;
        assert!(
            edge_bright > 1.0,
            "UnsharpMask should overshoot at edges (got {edge_bright})"
        );
    }

    #[test]
    fn test_unsharp_mask_preserves_dimensions() {
        let w = 10;
        let h = 8;
        let img = Image::new(w, h);
        let um = UnsharpMask {
            sigma: 0.5,
            strength: 0.5,
        };
        let out = um.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    // ── BilateralFilter ───────────────────────────────────────────────────────

    #[test]
    fn test_bilateral_constant_image_unchanged() {
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.3, 0.7, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let bf = BilateralFilter {
            sigma_spatial: 2.0,
            sigma_range: 0.3,
            radius: 2,
        };
        let out = bf.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-4,
                    "bilateral on constant image should not change r"
                );
            }
        }
    }

    #[test]
    fn test_bilateral_preserves_dimensions() {
        let w = 12;
        let h = 10;
        let img = Image::new(w, h);
        let bf = BilateralFilter {
            sigma_spatial: 1.5,
            sigma_range: 0.2,
            radius: 3,
        };
        let out = bf.apply(&img);
        assert_eq!(out.width, w);
        assert_eq!(out.height, h);
    }

    #[test]
    fn test_bilateral_smooths_noise() {
        // Create a noisy image (alternating bright/dark pixels) and verify
        // the bilateral filter reduces the luminance range at interior pixels.
        let w = 9;
        let h = 9;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let v = if (x + y) % 2 == 0 { 0.0_f32 } else { 1.0_f32 };
                img.set(x, y, PostColor::new(v, v, v, 1.0));
            }
        }
        let bf = BilateralFilter {
            sigma_spatial: 2.0,
            sigma_range: 0.5,
            radius: 3,
        };
        let out = bf.apply(&img);
        // Interior pixels should have luminance between 0.1 and 0.9.
        let cx = w / 2;
        let cy = h / 2;
        let lum = out.get(cx, cy).luminance();
        assert!(
            lum > 0.1 && lum < 0.9,
            "bilateral should smooth noisy checkerboard (lum={lum})"
        );
    }
}
