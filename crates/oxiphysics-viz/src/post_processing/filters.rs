// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tone mapping operators, bloom effect, motion blur, ambient occlusion,
//! and colour grading filters.

#![allow(dead_code)]

use super::core::{Image, PostColor};
use super::effects::{Bloom, GaussianBlur, Ssao, ToneMapping};

// ─── ToneMappingOperator ──────────────────────────────────────────────────────

/// Tone-mapping operator selection.
///
/// Each variant maps an HDR linear-light value (potentially > 1.0) to
/// a displayable LDR value in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToneMappingOperator {
    /// Simple Reinhard: `c / (1 + c)`.
    Reinhard,
    /// ACES filmic approximation (Narkowicz 2015).
    Aces,
    /// Uncharted 2 / Hable "filmic" curve.
    Uncharted2,
    /// Clamped linear: `clamp(c, 0, 1)`.
    Filmic,
}

/// Apply a [`ToneMappingOperator`] to every pixel in `image` channel-wise.
///
/// Returns a new [`Image`]; the depth buffer is preserved.
pub fn apply_tone_mapping(image: &Image, operator: ToneMappingOperator) -> Image {
    let w = image.width;
    let h = image.height;
    let mut out = Image::new(w, h);

    /// Uncharted 2 partial function (Hable).
    #[inline]
    fn uc2_partial(x: f32) -> f32 {
        let a = 0.15_f32;
        let b = 0.50_f32;
        let c = 0.10_f32;
        let d = 0.20_f32;
        let e = 0.02_f32;
        let f = 0.30_f32;
        ((x * (a * x + c * b) + d * e) / (x * (a * x + b) + d * f)) - e / f
    }

    #[inline]
    fn uncharted2(x: f32) -> f32 {
        let w = 11.2_f32; // white point
        uc2_partial(x) / uc2_partial(w)
    }

    for y in 0..h {
        for x in 0..w {
            let c = image.get(x, y);
            let r = match operator {
                ToneMappingOperator::Reinhard => ToneMapping::reinhard(c.r),
                ToneMappingOperator::Aces => ToneMapping::aces_filmic(c.r),
                ToneMappingOperator::Uncharted2 => uncharted2(c.r).clamp(0.0, 1.0),
                ToneMappingOperator::Filmic => c.r.clamp(0.0, 1.0),
            };
            let g = match operator {
                ToneMappingOperator::Reinhard => ToneMapping::reinhard(c.g),
                ToneMappingOperator::Aces => ToneMapping::aces_filmic(c.g),
                ToneMappingOperator::Uncharted2 => uncharted2(c.g).clamp(0.0, 1.0),
                ToneMappingOperator::Filmic => c.g.clamp(0.0, 1.0),
            };
            let b = match operator {
                ToneMappingOperator::Reinhard => ToneMapping::reinhard(c.b),
                ToneMappingOperator::Aces => ToneMapping::aces_filmic(c.b),
                ToneMappingOperator::Uncharted2 => uncharted2(c.b).clamp(0.0, 1.0),
                ToneMappingOperator::Filmic => c.b.clamp(0.0, 1.0),
            };
            let mapped = PostColor::new(r, g, b, c.a);
            out.set(x, y, mapped);
            let di = out.idx(x, y);
            out.depth[di] = image.get_depth(x, y);
        }
    }
    out
}

// ─── BloomEffect ─────────────────────────────────────────────────────────────

/// Overexposed-pixel bloom using a configurable bright-pass threshold.
///
/// This is a higher-level alternative to the lower-level [`Bloom`] struct that
/// exposes explicit control over the blur radius independently of sigma.
pub struct BloomEffect {
    /// Luminance threshold; pixels below this value do not contribute.
    pub threshold: f32,
    /// Additive intensity of the bloom pass.
    pub intensity: f32,
    /// Standard deviation (σ) for the Gaussian blur applied to bright pixels.
    pub blur_radius: f32,
}

impl BloomEffect {
    /// Apply bloom: bright-pass → Gaussian blur → additive composite.
    pub fn apply(&self, image: &Image) -> Image {
        let bright = Bloom::extract_bright(image, self.threshold);
        let blurred = GaussianBlur::apply(&bright, self.blur_radius.max(0.1));

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

// ─── MotionBlur ──────────────────────────────────────────────────────────────

/// Camera-motion blur using a simple radial velocity-buffer approximation.
///
/// In a full renderer the velocity buffer would contain per-pixel screen-space
/// motion vectors; here we simulate that by accepting explicit per-pixel
/// velocities or falling back to a uniform screen-centre radial field.
pub struct MotionBlur {
    /// Number of samples along the velocity vector.
    pub samples: usize,
    /// Shutter angle in degrees (0–360).  180° is the cinematic default.
    ///
    /// The effective blur length (in pixels) at a given velocity magnitude `v`
    /// is `v * shutter_angle / 360`.
    pub shutter_angle: f32,
}

impl MotionBlur {
    /// Apply motion blur using a uniform velocity (in pixels) for all pixels.
    ///
    /// Each pixel is blended with `samples` neighbours sampled along `velocity`.
    pub fn apply_uniform(&self, image: &Image, velocity: [f32; 2]) -> Image {
        let w = image.width;
        let h = image.height;
        let n = self.samples.max(1) as f32;
        let scale = self.shutter_angle / 360.0;
        let vx = velocity[0] * scale;
        let vy = velocity[1] * scale;

        let mut out = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                let mut acc = PostColor::zero();
                for s in 0..self.samples.max(1) {
                    let t = s as f32 / n;
                    let sx = x as f32 - vx * t;
                    let sy = y as f32 - vy * t;
                    let ix = sx.round() as isize;
                    let iy = sy.round() as isize;
                    let sample = if ix >= 0 && iy >= 0 && ix < w as isize && iy < h as isize {
                        image.get(ix as usize, iy as usize)
                    } else {
                        PostColor::zero()
                    };
                    acc = acc.add(&sample);
                }
                let blended = acc.scale(1.0 / n);
                out.set(x, y, blended);
                let di = out.idx(x, y);
                out.depth[di] = image.get_depth(x, y);
            }
        }
        out
    }
}

// ─── AmbientOcclusion ─────────────────────────────────────────────────────────

/// SSAO-style ambient occlusion parameters, as a standalone helper struct.
///
/// Differs from [`Ssao`] in that it stores the radius in world-space units and
/// exposes a `bias` to prevent self-occlusion.  Internally delegates to [`Ssao`].
pub struct AmbientOcclusion {
    /// Sampling radius in pixels.
    pub radius: f32,
    /// Number of occlusion samples.
    pub samples: usize,
    /// Depth-space bias to prevent self-occlusion artefacts.
    pub bias: f32,
}

impl AmbientOcclusion {
    /// Compute and apply ambient occlusion to `image`.
    pub fn apply(&self, image: &Image) -> Image {
        let ssao = Ssao {
            num_samples: self.samples.max(1),
            radius: self.radius.max(1e-4),
            bias: self.bias,
            strength: 1.0,
        };
        ssao.apply(image)
    }

    /// Return only the AO factor map (greyscale) without blending with the image.
    pub fn compute_map(&self, image: &Image) -> Image {
        let ssao = Ssao {
            num_samples: self.samples.max(1),
            radius: self.radius.max(1e-4),
            bias: self.bias,
            strength: 1.0,
        };
        ssao.compute_ao_map(image)
    }
}

// ─── ColorGrading ────────────────────────────────────────────────────────────

/// PBR-style colour grading applied in linear light.
pub struct ColorGrading {
    /// EV (stops) exposure offset; 0 = no change, +1 = one stop brighter.
    pub exposure: f32,
    /// Contrast around mid-grey (0.18); > 1 increases contrast, < 1 reduces.
    pub contrast: f32,
    /// Colour saturation; 1 = original, 0 = greyscale, 2 = double-saturated.
    pub saturation: f32,
    /// Hue shift in degrees `[-180, 180]`.
    pub hue_shift: f32,
}

impl ColorGrading {
    /// Construct with neutral (identity) settings.
    pub fn identity() -> Self {
        Self {
            exposure: 0.0,
            contrast: 1.0,
            saturation: 1.0,
            hue_shift: 0.0,
        }
    }

    /// Apply the grading to a single `PostColor`.
    ///
    /// Operations are applied in order: exposure → contrast → saturation → hue shift.
    pub fn grade_color(&self, c: PostColor) -> PostColor {
        // 1. Exposure: multiply by 2^exposure.
        let exp_mul = (2.0_f32).powf(self.exposure);
        let r = c.r * exp_mul;
        let g = c.g * exp_mul;
        let b = c.b * exp_mul;

        // 2. Contrast around 0.18 mid-grey.
        let mid = 0.18_f32;
        let r = mid + (r - mid) * self.contrast;
        let g = mid + (g - mid) * self.contrast;
        let b = mid + (b - mid) * self.contrast;

        // 3. Saturation via luminance.
        let lum = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        let r = lum + (r - lum) * self.saturation;
        let g = lum + (g - lum) * self.saturation;
        let b = lum + (b - lum) * self.saturation;

        // 4. Hue shift: convert RGB → HSV → rotate H → convert back.
        let (h, s, v) = rgb_to_hsv(r, g, b);
        let h2 = (h + self.hue_shift).rem_euclid(360.0);
        let (r2, g2, b2) = hsv_to_rgb(h2, s, v);

        PostColor::new(r2, g2, b2, c.a)
    }
}

/// Apply [`ColorGrading`] to every pixel in `image`.
pub fn apply_color_grading(image: &Image, grading: &ColorGrading) -> Image {
    let w = image.width;
    let h = image.height;
    let mut out = Image::new(w, h);
    for y in 0..h {
        for x in 0..w {
            let c = image.get(x, y);
            out.set(x, y, grading.grade_color(c));
            let di = out.idx(x, y);
            out.depth[di] = image.get_depth(x, y);
        }
    }
    out
}

/// Convert RGB to HSV.  All inputs/outputs in `[0, 1]` except H which is in `[0, 360)`.
fn rgb_to_hsv(r: f32, g: f32, b: f32) -> (f32, f32, f32) {
    let cmax = r.max(g).max(b);
    let cmin = r.min(g).min(b);
    let delta = cmax - cmin;

    let h = if delta < 1e-7 {
        0.0
    } else if cmax == r {
        60.0 * (((g - b) / delta).rem_euclid(6.0))
    } else if cmax == g {
        60.0 * ((b - r) / delta + 2.0)
    } else {
        60.0 * ((r - g) / delta + 4.0)
    };

    let s = if cmax < 1e-7 { 0.0 } else { delta / cmax };
    let v = cmax;
    (h, s, v)
}

/// Convert HSV to RGB.  H in `[0, 360)`, S/V in `[0, 1\]`.
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> (f32, f32, f32) {
    if s < 1e-7 {
        return (v, v, v);
    }
    let hi = (h / 60.0).floor() as u32 % 6;
    let f = h / 60.0 - (h / 60.0).floor();
    let p = v * (1.0 - s);
    let q = v * (1.0 - f * s);
    let t = v * (1.0 - (1.0 - f) * s);
    match hi {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── ToneMappingOperator ───────────────────────────────────────────────────

    #[test]
    fn test_tone_mapping_operator_reinhard_maps_zero_to_zero() {
        let mut img = Image::new(2, 2);
        img.set(0, 0, PostColor::new(0.0, 0.0, 0.0, 1.0));
        let out = apply_tone_mapping(&img, ToneMappingOperator::Reinhard);
        let c = out.get(0, 0);
        assert!(c.r.abs() < 1e-6 && c.g.abs() < 1e-6 && c.b.abs() < 1e-6);
    }

    #[test]
    fn test_tone_mapping_operator_aces_output_in_range() {
        let mut img = Image::new(4, 1);
        for x in 0..4 {
            img.set(x, 0, PostColor::new(x as f32 * 3.0 + 0.1, 1.0, 0.5, 1.0));
        }
        let out = apply_tone_mapping(&img, ToneMappingOperator::Aces);
        for x in 0..4 {
            let c = out.get(x, 0);
            assert!(
                c.r >= 0.0 && c.r <= 1.0,
                "ACES r out of range at x={x}: {}",
                c.r
            );
            assert!(
                c.g >= 0.0 && c.g <= 1.0,
                "ACES g out of range at x={x}: {}",
                c.g
            );
        }
    }

    #[test]
    fn test_tone_mapping_operator_uncharted2_output_in_range() {
        let mut img = Image::new(3, 1);
        img.set(0, 0, PostColor::new(0.0, 0.0, 0.0, 1.0));
        img.set(1, 0, PostColor::new(1.0, 1.0, 1.0, 1.0));
        img.set(2, 0, PostColor::new(5.0, 5.0, 5.0, 1.0));
        let out = apply_tone_mapping(&img, ToneMappingOperator::Uncharted2);
        for x in 0..3 {
            let c = out.get(x, 0);
            assert!(
                c.r >= 0.0 && c.r <= 1.0 + 1e-5,
                "Uncharted2 r out of range: {}",
                c.r
            );
        }
    }

    #[test]
    fn test_tone_mapping_operator_filmic_clamps() {
        let mut img = Image::new(2, 1);
        img.set(0, 0, PostColor::new(2.0, -0.5, 0.5, 1.0));
        img.set(1, 0, PostColor::new(0.3, 0.6, 0.9, 1.0));
        let out = apply_tone_mapping(&img, ToneMappingOperator::Filmic);
        let c0 = out.get(0, 0);
        assert!((c0.r - 1.0).abs() < 1e-6, "filmic should clamp 2.0 to 1.0");
        assert!(c0.g.abs() < 1e-6, "filmic should clamp -0.5 to 0.0");
        let c1 = out.get(1, 0);
        assert!((c1.b - 0.9).abs() < 1e-5, "filmic should preserve 0.9");
    }

    // ── BloomEffect ──────────────────────────────────────────────────────────

    #[test]
    fn test_bloom_effect_dark_image_unchanged() {
        // A completely dark image should produce zero bloom addition.
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let dark = PostColor::new(0.05, 0.05, 0.05, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, dark);
            }
        }
        let bloom = BloomEffect {
            threshold: 0.5,
            intensity: 2.0,
            blur_radius: 1.0,
        };
        let out = bloom.apply(&img);
        // Pixels should stay close to the original (no bright-pass contribution).
        let c = out.get(2, 2);
        assert!(c.r < 0.15, "dark pixel should not be brightened by bloom");
    }

    #[test]
    fn test_bloom_effect_bright_pixel_spreads() {
        let w = 9;
        let h = 9;
        let mut img = Image::new(w, h);
        // Single bright pixel in the centre.
        img.set(4, 4, PostColor::new(2.0, 2.0, 2.0, 1.0));
        let bloom = BloomEffect {
            threshold: 0.5,
            intensity: 1.0,
            blur_radius: 2.0,
        };
        let out = bloom.apply(&img);
        // Adjacent pixel should have picked up some bloom.
        let adj = out.get(4, 5);
        assert!(adj.r > 0.0, "bloom should spread to adjacent pixels");
    }

    // ── MotionBlur ────────────────────────────────────────────────────────────

    #[test]
    fn test_motion_blur_zero_velocity_no_change() {
        let w = 4;
        let h = 4;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.3, 0.1, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let mb = MotionBlur {
            samples: 4,
            shutter_angle: 180.0,
        };
        let out = mb.apply_uniform(&img, [0.0, 0.0]);
        // With zero velocity all samples hit the same pixel, result equals original.
        let o = out.get(2, 2);
        assert!(
            (o.r - c.r).abs() < 1e-4,
            "zero velocity should not change image"
        );
    }

    #[test]
    fn test_motion_blur_constant_field_preserves_color() {
        // A uniform-color image blurred with any velocity should stay uniform.
        let w = 8;
        let h = 8;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.6, 0.4, 0.2, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let mb = MotionBlur {
            samples: 6,
            shutter_angle: 180.0,
        };
        let out = mb.apply_uniform(&img, [2.0, 1.0]);
        // Interior pixels should be unchanged (all samples land within the image).
        for y in 4..6 {
            for x in 4..6 {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-4,
                    "constant field: r mismatch at ({x},{y})"
                );
            }
        }
    }

    // ── AmbientOcclusion ─────────────────────────────────────────────────────

    #[test]
    fn test_ao_map_flat_depth_no_occlusion() {
        // Perfectly flat depth buffer → no pixel occludes another → AO = 1.
        let w = 5;
        let h = 5;
        let img = Image::new(w, h); // depth initialised to 1.0
        let ao = AmbientOcclusion {
            radius: 2.0,
            samples: 8,
            bias: 0.01,
        };
        let map = ao.compute_map(&img);
        for y in 0..h {
            for x in 0..w {
                let v = map.get(x, y).r;
                assert!(v > 0.9, "flat depth → AO factor should be near 1, got {v}");
            }
        }
    }

    #[test]
    fn test_ao_apply_darkens_occluded_region() {
        let w = 7;
        let h = 7;
        let mut img = Image::new(w, h);
        let white = PostColor::new(1.0, 1.0, 1.0, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, white);
            }
        }
        // Create a depth discontinuity: one pixel much closer.
        let centre_idx = (h / 2) * w + (w / 2);
        img.depth[centre_idx] = 0.0;
        let ao = AmbientOcclusion {
            radius: 3.0,
            samples: 16,
            bias: 0.01,
        };
        let out = ao.apply(&img);
        // The centre pixel should not be darkened (it is the closest point).
        // A surrounding pixel should be somewhat darkened.
        let centre = out.get(w / 2, h / 2);
        assert!(centre.r > 0.0, "output should have non-zero luminance");
    }

    // ── ColorGrading ─────────────────────────────────────────────────────────

    #[test]
    fn test_color_grading_identity_no_change() {
        let c = PostColor::new(0.4, 0.6, 0.2, 1.0);
        let g = ColorGrading::identity();
        let out = g.grade_color(c);
        assert!(
            (out.r - c.r).abs() < 1e-4,
            "identity grading should not change r"
        );
        assert!(
            (out.g - c.g).abs() < 1e-4,
            "identity grading should not change g"
        );
        assert!(
            (out.b - c.b).abs() < 1e-4,
            "identity grading should not change b"
        );
    }

    #[test]
    fn test_color_grading_exposure_doubles_brightness() {
        let c = PostColor::new(0.2, 0.2, 0.2, 1.0);
        let g = ColorGrading {
            exposure: 1.0,
            contrast: 1.0,
            saturation: 1.0,
            hue_shift: 0.0,
        };
        let out = g.grade_color(c);
        // +1 EV → multiply by 2.  Contrast also applied around 0.18 midgrey.
        // Verify output is brighter.
        assert!(out.r > c.r, "positive exposure should increase brightness");
    }

    #[test]
    fn test_color_grading_saturation_zero_greyscale() {
        let c = PostColor::new(0.8, 0.2, 0.4, 1.0);
        let g = ColorGrading {
            exposure: 0.0,
            contrast: 1.0,
            saturation: 0.0,
            hue_shift: 0.0,
        };
        let out = g.grade_color(c);
        // With saturation=0 all channels should equal luminance.
        assert!(
            (out.r - out.g).abs() < 1e-4,
            "zero saturation: r should equal g"
        );
        assert!(
            (out.g - out.b).abs() < 1e-4,
            "zero saturation: g should equal b"
        );
    }

    #[test]
    fn test_apply_color_grading_applies_to_all_pixels() {
        let w = 3;
        let h = 3;
        let mut img = Image::new(w, h);
        let c = PostColor::new(0.5, 0.5, 0.5, 1.0);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, c);
            }
        }
        let g = ColorGrading {
            exposure: 1.0,
            contrast: 1.0,
            saturation: 1.0,
            hue_shift: 0.0,
        };
        let out = apply_color_grading(&img, &g);
        // All pixels should be brighter than the original.
        for y in 0..h {
            for x in 0..w {
                assert!(
                    out.get(x, y).r > c.r,
                    "apply_color_grading should affect all pixels"
                );
            }
        }
    }
}
