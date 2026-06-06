// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Basic post-processing effects: Gaussian blur, depth-of-field, bloom,
//! tone mapping, SSAO, vignette, and the composable pipeline.

use std::f32::consts::PI;

use super::core::{Image, PostColor};

// ─── GaussianBlur ─────────────────────────────────────────────────────────────

/// Separable Gaussian blur for [`Image`].
pub struct GaussianBlur;

impl GaussianBlur {
    /// Build a 1-D normalised Gaussian kernel with standard deviation `sigma`
    /// and half-width `radius` (total length = `2*radius + 1`).
    pub fn kernel_1d(sigma: f32, radius: usize) -> Vec<f32> {
        let len = 2 * radius + 1;
        let mut k = Vec::with_capacity(len);
        let inv_denom = 1.0 / (2.0 * sigma * sigma);
        let mut sum = 0.0_f32;
        for i in 0..len {
            let x = i as f32 - radius as f32;
            let v = (-x * x * inv_denom).exp();
            k.push(v);
            sum += v;
        }
        // Normalise so weights sum to 1.
        for v in &mut k {
            *v /= sum;
        }
        k
    }

    /// Apply a 1-D kernel horizontally (along x).
    pub fn apply_horizontal(image: &Image, kernel: &[f32]) -> Image {
        let w = image.width;
        let h = image.height;
        let radius = kernel.len() / 2;
        let mut out = Image::new(w, h);

        for y in 0..h {
            for x in 0..w {
                let mut acc = PostColor::zero();
                let mut weight_sum = 0.0_f32;
                for (ki, &kw) in kernel.iter().enumerate() {
                    let sx = x as isize + ki as isize - radius as isize;
                    if sx >= 0 && sx < w as isize {
                        let src = image.get(sx as usize, y);
                        acc = acc.add(&src.scale(kw));
                        weight_sum += kw;
                    }
                }
                if weight_sum > 0.0 {
                    acc = acc.scale(1.0 / weight_sum);
                }
                out.set(x, y, acc);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }

    /// Apply a 1-D kernel vertically (along y).
    pub fn apply_vertical(image: &Image, kernel: &[f32]) -> Image {
        let w = image.width;
        let h = image.height;
        let radius = kernel.len() / 2;
        let mut out = Image::new(w, h);

        for y in 0..h {
            for x in 0..w {
                let mut acc = PostColor::zero();
                let mut weight_sum = 0.0_f32;
                for (ki, &kw) in kernel.iter().enumerate() {
                    let sy = y as isize + ki as isize - radius as isize;
                    if sy >= 0 && sy < h as isize {
                        let src = image.get(x, sy as usize);
                        acc = acc.add(&src.scale(kw));
                        weight_sum += kw;
                    }
                }
                if weight_sum > 0.0 {
                    acc = acc.scale(1.0 / weight_sum);
                }
                out.set(x, y, acc);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }

    /// Separable 2-D Gaussian blur with standard deviation `sigma`.
    ///
    /// The kernel radius is chosen as `ceil(3*sigma)` (at least 1).
    pub fn apply(image: &Image, sigma: f32) -> Image {
        let radius = ((3.0 * sigma).ceil() as usize).max(1);
        let kernel = Self::kernel_1d(sigma, radius);
        let h_pass = Self::apply_horizontal(image, &kernel);
        Self::apply_vertical(&h_pass, &kernel)
    }
}

// ─── DepthOfField ─────────────────────────────────────────────────────────────

/// Depth-of-field effect based on a simple circle-of-confusion model.
pub struct DepthOfField {
    /// Depth at which the scene is in perfect focus (`[0, 1]`).
    pub focus_distance: f32,
    /// Aperture size — scales the circle of confusion.
    pub aperture: f32,
    /// Maximum blur radius in pixels.
    pub max_blur_radius: f32,
}

impl DepthOfField {
    /// Compute the blur radius (in pixels) for a pixel at the given `depth`.
    ///
    /// `coc = aperture * |depth - focus_distance| / focus_distance`
    /// clamped to `[0, max_blur_radius]`.
    pub fn blur_radius_at_depth(&self, depth: f32) -> f32 {
        let denom = self.focus_distance.max(1e-6);
        let coc = self.aperture * (depth - self.focus_distance).abs() / denom;
        coc.min(self.max_blur_radius)
    }

    /// Apply depth-of-field blur to `image`.
    ///
    /// For each pixel the blur radius is derived from its depth value.
    /// Nearby pixels within that radius are accumulated with Gaussian weights.
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);

        for y in 0..h {
            for x in 0..w {
                let depth = image.get_depth(x, y);
                let radius = self.blur_radius_at_depth(depth);
                let irad = radius.ceil() as isize;

                if irad == 0 {
                    out.set(x, y, image.get(x, y));
                    let di = out.idx(x, y);
                    out.depth[di] = depth;
                    continue;
                }

                let sigma = (radius / 3.0).max(0.1);
                let inv_two_sigma_sq = 1.0 / (2.0 * sigma * sigma);

                let mut acc = PostColor::zero();
                let mut weight_sum = 0.0_f32;

                for dy in -irad..=irad {
                    for dx in -irad..=irad {
                        let sx = x as isize + dx;
                        let sy = y as isize + dy;
                        if sx < 0 || sy < 0 || sx >= w as isize || sy >= h as isize {
                            continue;
                        }
                        let dist2 = (dx * dx + dy * dy) as f32;
                        if dist2 > (irad * irad) as f32 {
                            continue;
                        }
                        let w_g = (-dist2 * inv_two_sigma_sq).exp();
                        let src = image.get(sx as usize, sy as usize);
                        acc = acc.add(&src.scale(w_g));
                        weight_sum += w_g;
                    }
                }

                if weight_sum > 0.0 {
                    acc = acc.scale(1.0 / weight_sum);
                }
                out.set(x, y, acc);
                let di = out.idx(x, y);
                out.depth[di] = depth;
            }
        }
        out
    }
}

// ─── Bloom ────────────────────────────────────────────────────────────────────

/// Screen-space bloom effect.
pub struct Bloom {
    /// Luminance threshold — only pixels brighter than this contribute.
    pub threshold: f32,
    /// Bloom strength (additive intensity).
    pub intensity: f32,
    /// Gaussian blur sigma for the bright pass.
    pub blur_sigma: f32,
}

impl Bloom {
    /// Extract pixels whose luminance exceeds `threshold`; dark pixels become zero.
    pub fn extract_bright(image: &Image, threshold: f32) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                if c.luminance() > threshold {
                    out.set(x, y, c);
                }
                // else leave as zero (default)
            }
        }
        out
    }

    /// Apply bloom: bright-pass → blur → add back to original.
    pub fn apply(&self, image: &Image) -> Image {
        let bright = Self::extract_bright(image, self.threshold);
        let blurred = GaussianBlur::apply(&bright, self.blur_sigma);

        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let orig = image.get(x, y);
                let b = blurred.get(x, y);
                let result = orig.add(&b.scale(self.intensity));
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── ToneMapping ──────────────────────────────────────────────────────────────

/// Collection of tone-mapping operators.
pub struct ToneMapping;

impl ToneMapping {
    /// Reinhard operator: `c / (1 + c)`.
    #[inline]
    pub fn reinhard(c: f32) -> f32 {
        c / (1.0 + c)
    }

    /// Extended Reinhard operator that preserves `max_white`.
    ///
    /// `c * (1 + c / max_white²) / (1 + c)`
    #[inline]
    pub fn reinhard_extended(c: f32, max_white: f32) -> f32 {
        let mw2 = max_white * max_white;
        c * (1.0 + c / mw2) / (1.0 + c)
    }

    /// ACES filmic tone-mapping approximation (Krzysztof Narkowicz, 2015).
    #[inline]
    pub fn aces_filmic(c: f32) -> f32 {
        let a = 2.51_f32;
        let b = 0.03_f32;
        let cc = 2.43_f32;
        let d = 0.59_f32;
        let e = 0.14_f32;
        ((c * (a * c + b)) / (c * (cc * c + d) + e)).clamp(0.0, 1.0)
    }

    /// Apply Reinhard tone mapping channel-wise then gamma-correct (γ = 2.2).
    pub fn apply_reinhard(image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let mapped = PostColor::new(
                    Self::reinhard(c.r),
                    Self::reinhard(c.g),
                    Self::reinhard(c.b),
                    c.a,
                );
                out.set(x, y, mapped);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }

    /// Apply ACES filmic tone mapping channel-wise.
    pub fn apply_aces(image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let mapped = PostColor::new(
                    Self::aces_filmic(c.r),
                    Self::aces_filmic(c.g),
                    Self::aces_filmic(c.b),
                    c.a,
                );
                out.set(x, y, mapped);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }

    /// Apply gamma correction: raise each channel to the power `1 / gamma`.
    pub fn gamma_correct(image: &Image, gamma: f32) -> Image {
        let w = image.width;
        let h = image.height;
        let inv_gamma = 1.0 / gamma;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let corrected = PostColor::new(
                    c.r.max(0.0).powf(inv_gamma),
                    c.g.max(0.0).powf(inv_gamma),
                    c.b.max(0.0).powf(inv_gamma),
                    c.a,
                );
                out.set(x, y, corrected);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── SSAO ─────────────────────────────────────────────────────────────────────

/// Screen-Space Ambient Occlusion (SSAO) approximation.
///
/// Uses a simple depth-based approach: for each pixel, samples a disc of
/// neighbouring pixels and computes how many are "occluded" (i.e. shallower).
pub struct Ssao {
    /// Number of samples in the occlusion kernel.
    pub num_samples: usize,
    /// Sampling radius in pixels.
    pub radius: f32,
    /// Bias to avoid self-occlusion artefacts.
    pub bias: f32,
    /// Overall strength of the AO effect `[0, 1]`.
    pub strength: f32,
}

impl Ssao {
    /// Compute an AO factor for every pixel and return a greyscale occlusion map.
    ///
    /// The returned image has each pixel set to `(ao, ao, ao, 1.0)` where
    /// `ao ∈ [0, 1]`; 1.0 = fully unoccluded, 0.0 = fully occluded.
    pub fn compute_ao_map(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);

        // Precompute a simple disc of sample offsets distributed evenly.
        let mut offsets: Vec<(f32, f32)> = Vec::with_capacity(self.num_samples);
        for i in 0..self.num_samples {
            let angle = 2.0 * PI * (i as f32) / (self.num_samples as f32);
            // Use a golden-ratio-like spiral for coverage.
            let r = self.radius * ((i as f32 + 1.0) / self.num_samples as f32).sqrt();
            offsets.push((r * angle.cos(), r * angle.sin()));
        }

        for y in 0..h {
            for x in 0..w {
                let center_depth = image.get_depth(x, y);
                let mut occlusion = 0.0_f32;

                for &(ox, oy) in &offsets {
                    let sx = x as isize + ox.round() as isize;
                    let sy = y as isize + oy.round() as isize;
                    if sx < 0 || sy < 0 || sx >= w as isize || sy >= h as isize {
                        continue;
                    }
                    let sample_depth = image.get_depth(sx as usize, sy as usize);
                    // A sample occludes if it is shallower (closer to camera) by
                    // more than the bias.
                    if center_depth - sample_depth > self.bias {
                        // Range-check: samples outside a reasonable depth range
                        // don't occlude.
                        let range_check =
                            (1.0 - ((center_depth - sample_depth) / self.radius).abs()).max(0.0);
                        occlusion += range_check;
                    }
                }

                let ao_raw = 1.0 - self.strength * (occlusion / self.num_samples as f32);
                let ao = ao_raw.clamp(0.0, 1.0);
                let pixel = PostColor::new(ao, ao, ao, 1.0);
                out.set(x, y, pixel);
                let di = out.idx(x, y);
                out.depth[di] = center_depth;
            }
        }
        out
    }

    /// Apply SSAO: compute the AO map and multiply each pixel's colour.
    pub fn apply(&self, image: &Image) -> Image {
        let ao_map = self.compute_ao_map(image);
        let w = image.width;
        let h = image.height;
        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let c = image.get(x, y);
                let ao = ao_map.get(x, y).r; // greyscale factor
                let result = PostColor::new(c.r * ao, c.g * ao, c.b * ao, c.a);
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── Vignette ─────────────────────────────────────────────────────────────────

/// Lens vignette darkening towards corners.
pub struct Vignette {
    /// Vignette strength in `[0, 1]`: 0 = no vignette, 1 = full darkening at corners.
    pub strength: f32,
}

impl Vignette {
    /// Apply vignette: multiply each pixel by `1 - strength * r²` where `r` is
    /// the normalised distance from the image centre (`r = 0` at centre,
    /// `r = 1` at corner).
    pub fn apply(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        let cx = (w as f32 - 1.0) * 0.5;
        let cy = (h as f32 - 1.0) * 0.5;
        // Maximum possible distance from centre (corner distance).
        let max_dist = ((cx * cx) + (cy * cy)).sqrt().max(1e-6);

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let dx = x as f32 - cx;
                let dy = y as f32 - cy;
                let r = (dx * dx + dy * dy).sqrt() / max_dist; // normalised [0,1]
                let factor = (1.0 - self.strength * r * r).max(0.0);
                let c = image.get(x, y);
                let result = PostColor::new(c.r * factor, c.g * factor, c.b * factor, c.a);
                out.set(x, y, result);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── PostProcessPipeline ──────────────────────────────────────────────────────

/// Composable post-processing pipeline.
///
/// Effects are applied in order: bloom → depth-of-field → tone mapping → vignette.
pub struct PostProcessPipeline {
    /// Optional bloom effect.
    pub bloom: Option<Bloom>,
    /// Optional depth-of-field effect.
    pub dof: Option<DepthOfField>,
    /// Whether to apply ACES filmic tone mapping.
    pub tone_mapping: bool,
    /// Optional vignette.
    pub vignette: Option<Vignette>,
}

impl PostProcessPipeline {
    /// Create an empty pipeline with all effects disabled.
    pub fn new() -> Self {
        Self {
            bloom: None,
            dof: None,
            tone_mapping: false,
            vignette: None,
        }
    }

    /// Run the pipeline on `image`, returning a new processed image.
    pub fn process(&self, image: &Image) -> Image {
        // We need to chain through owned `Image` values.
        // Using a small helper to avoid repeated borrowing.

        // Step 1: bloom
        let after_bloom: Image = match &self.bloom {
            Some(b) => b.apply(image),
            None => {
                // Clone the image without having Clone on Image — build manually.
                let w = image.width;
                let h = image.height;
                let mut img = Image::new(w, h);
                img.pixels.clone_from(&image.pixels);
                img.depth.clone_from(&image.depth);
                img
            }
        };

        // Step 2: depth of field
        let after_dof: Image = match &self.dof {
            Some(d) => d.apply(&after_bloom),
            None => after_bloom,
        };

        // Step 3: tone mapping
        let after_tone: Image = if self.tone_mapping {
            ToneMapping::apply_aces(&after_dof)
        } else {
            after_dof
        };

        // Step 4: vignette
        match &self.vignette {
            Some(v) => v.apply(&after_tone),
            None => after_tone,
        }
    }
}

impl Default for PostProcessPipeline {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── GaussianBlur ─────────────────────────────────────────────────────────

    #[test]
    fn test_gaussian_kernel_sums_to_one() {
        for &sigma in &[0.5_f32, 1.0, 2.0, 3.0] {
            let radius = ((3.0 * sigma).ceil() as usize).max(1);
            let k = GaussianBlur::kernel_1d(sigma, radius);
            let sum: f32 = k.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-5,
                "kernel sum for sigma={sigma} should be 1.0, got {sum}"
            );
        }
    }

    #[test]
    fn test_gaussian_blur_constant_image() {
        // Blurring a constant image should produce the same constant image.
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.4, 0.3, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let blurred = GaussianBlur::apply(&img, 1.0);
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let out = blurred.get(x, y);
                assert!(
                    (out.r - c.r).abs() < 1e-4,
                    "blur of constant image should be constant (r)"
                );
                assert!(
                    (out.g - c.g).abs() < 1e-4,
                    "blur of constant image should be constant (g)"
                );
                assert!(
                    (out.b - c.b).abs() < 1e-4,
                    "blur of constant image should be constant (b)"
                );
            }
        }
    }

    // ── Bloom ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_bloom_extract_bright_passes_bright_rejects_dark() {
        let mut img = Image::new(4, 4);
        // Bright pixel at (1,1)
        img.set(1, 1, PostColor::new(1.0, 1.0, 1.0, 1.0));
        // Dark pixel at (2,2)
        img.set(2, 2, PostColor::new(0.1, 0.1, 0.1, 1.0));

        let threshold = 0.5;
        let bright = Bloom::extract_bright(&img, threshold);

        // Bright pixel should pass.
        let bp = bright.get(1, 1);
        assert!(
            bp.luminance() > threshold,
            "bright pixel should pass through"
        );

        // Dark pixel should become zero.
        let dp = bright.get(2, 2);
        assert!(
            dp.r.abs() < 1e-6 && dp.g.abs() < 1e-6 && dp.b.abs() < 1e-6,
            "dark pixel should be zeroed by bright-pass extract"
        );
    }

    // ── ToneMapping ───────────────────────────────────────────────────────────

    #[test]
    fn test_reinhard_half_at_one() {
        let result = ToneMapping::reinhard(1.0);
        assert!(
            (result - 0.5).abs() < 1e-6,
            "reinhard(1.0) should equal 0.5, got {result}"
        );
    }

    #[test]
    fn test_reinhard_maps_zero_to_zero() {
        assert_eq!(ToneMapping::reinhard(0.0), 0.0);
    }

    #[test]
    fn test_aces_output_in_range() {
        // ACES output should be in [0,1] for inputs in [0, 100].
        for i in 0..=100 {
            let input = i as f32;
            let output = ToneMapping::aces_filmic(input);
            assert!(
                (0.0..=1.0).contains(&output),
                "ACES output {output} out of [0,1] for input {input}"
            );
        }
    }

    #[test]
    fn test_aces_monotone_increasing() {
        let mut prev = ToneMapping::aces_filmic(0.0);
        for i in 1..=50 {
            let cur = ToneMapping::aces_filmic(i as f32 * 0.1);
            assert!(
                cur >= prev - 1e-6,
                "ACES should be monotone increasing: f({}) = {} < f({}) = {}",
                i as f32 * 0.1,
                cur,
                (i - 1) as f32 * 0.1,
                prev
            );
            prev = cur;
        }
    }

    // ── DepthOfField ──────────────────────────────────────────────────────────

    #[test]
    fn test_dof_blur_radius_zero_at_focus() {
        let dof = DepthOfField {
            focus_distance: 0.5,
            aperture: 2.0,
            max_blur_radius: 10.0,
        };
        let r = dof.blur_radius_at_depth(0.5);
        assert!(
            r.abs() < 1e-6,
            "blur radius should be zero at focus distance, got {r}"
        );
    }

    #[test]
    fn test_dof_blur_radius_increases_away_from_focus() {
        let dof = DepthOfField {
            focus_distance: 0.5,
            aperture: 2.0,
            max_blur_radius: 100.0,
        };
        let r_near = dof.blur_radius_at_depth(0.1);
        let r_far = dof.blur_radius_at_depth(0.9);
        assert!(
            r_near > 0.0,
            "blur radius should be positive away from focus"
        );
        assert!(
            r_far > 0.0,
            "blur radius should be positive away from focus"
        );
    }

    #[test]
    fn test_dof_blur_radius_clamped() {
        let max = 5.0;
        let dof = DepthOfField {
            focus_distance: 0.5,
            aperture: 1000.0,
            max_blur_radius: max,
        };
        let r = dof.blur_radius_at_depth(1.0);
        assert!(
            r <= max + 1e-6,
            "blur radius should be clamped to max_blur_radius, got {r}"
        );
    }

    // ── Vignette ──────────────────────────────────────────────────────────────

    #[test]
    fn test_vignette_darkens_corners_vs_center() {
        let w = 11;
        let h = 11;
        let mut img = Image::new(w, h);
        // Fill with white.
        let white = PostColor::new(1.0, 1.0, 1.0, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, white);
            }
        }

        let vignette = Vignette { strength: 0.8 };
        let out = vignette.apply(&img);

        let cx = w / 2;
        let cy = h / 2;
        let center_lum = out.get(cx, cy).luminance();
        let corner_lum = out.get(0, 0).luminance();

        assert!(
            center_lum > corner_lum,
            "vignette should make corners darker than center: center={center_lum}, corner={corner_lum}"
        );
    }

    #[test]
    fn test_vignette_center_unaffected_with_zero_strength() {
        let w = 5;
        let h = 5;
        let mut img = Image::new(w, h);
        let white = PostColor::new(1.0, 1.0, 1.0, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, white);
            }
        }
        let vignette = Vignette { strength: 0.0 };
        let out = vignette.apply(&img);
        for y in 0..h {
            for x in 0..w {
                let c = out.get(x, y);
                assert!(
                    (c.r - 1.0).abs() < 1e-6,
                    "zero-strength vignette should not change image"
                );
            }
        }
    }

    // ── PostProcessPipeline ───────────────────────────────────────────────────

    #[test]
    fn test_pipeline_empty_pass_through() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.3, 0.1, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let pipeline = PostProcessPipeline::new();
        let out = pipeline.process(&img);
        for y in 0..h {
            for x in 0..w {
                let oc = out.get(x, y);
                assert!(
                    (oc.r - c.r).abs() < 1e-6,
                    "empty pipeline should not alter image"
                );
            }
        }
    }

    #[test]
    fn test_pipeline_with_tone_mapping() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        // HDR values above 1.0.
        let c = PostColor::new(3.0, 2.0, 1.5, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let mut pipeline = PostProcessPipeline::new();
        pipeline.tone_mapping = true;
        let out = pipeline.process(&img);
        for y in 0..h {
            for x in 0..w {
                let oc = out.get(x, y);
                assert!(
                    oc.r <= 1.0 + 1e-5 && oc.g <= 1.0 + 1e-5 && oc.b <= 1.0 + 1e-5,
                    "tone-mapped output should be in [0,1]"
                );
            }
        }
    }
}
