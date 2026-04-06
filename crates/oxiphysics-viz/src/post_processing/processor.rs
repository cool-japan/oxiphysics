// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! High-level post-processing dispatcher (PostProcessor).

#![allow(dead_code)]

use super::core::{Image, PostColor};
use super::image_filters::ChromaticAberrationEffect;

// ─── PostProcessor ────────────────────────────────────────────────────────────

/// High-level post-processing dispatcher that owns effect parameters and
/// applies them sequentially to an [`Image`].
///
/// Each method consumes an `&Image` and returns a new processed `Image`.  This
/// lets callers build a per-frame pipeline without mutating the original.
#[allow(dead_code)]
pub struct PostProcessor {
    /// Focal distance for depth-of-field (in normalised depth units 0–1).
    pub focal_distance: f32,
    /// Focus range half-width (pixels within `±dof_range` of focal depth are sharp).
    pub dof_range: f32,
    /// Maximum circle-of-confusion radius for depth-of-field blur (in pixels).
    pub dof_max_coc: f32,
    /// Strength of the analytical lens flare (0 = off, 1 = full).
    pub flare_intensity: f32,
    /// Number of flare streak arms.
    pub flare_streak_count: u32,
    /// Chromatic aberration red-channel horizontal offset (in pixels).
    pub ca_r_shift: i32,
    /// Chromatic aberration blue-channel horizontal offset (in pixels).
    pub ca_b_shift: i32,
}

#[allow(dead_code)]
impl PostProcessor {
    /// Construct a default `PostProcessor` with all effects disabled.
    pub fn new() -> Self {
        Self {
            focal_distance: 0.5,
            dof_range: 0.1,
            dof_max_coc: 4.0,
            flare_intensity: 0.0,
            flare_streak_count: 6,
            ca_r_shift: 0,
            ca_b_shift: 0,
        }
    }

    // -----------------------------------------------------------------------
    // Bokeh depth-of-field
    // -----------------------------------------------------------------------

    /// Apply a physically-motivated bokeh depth-of-field (DoF) effect to `image`.
    ///
    /// Pixels are blurred proportionally to their circle-of-confusion (CoC)
    /// radius, which grows linearly with the absolute distance from
    /// `focal_distance`.  The maximum blur kernel radius is clamped to
    /// `dof_max_coc`.  Pixels within `dof_range` of the focal depth are kept
    /// sharp.
    ///
    /// The effect uses a variable-radius box blur as an efficient CoC
    /// approximation (a full scatter-gather bokeh would require a GPU).
    pub fn compute_depth_of_field(&self, image: &Image) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 {
            return Image::new(w, h);
        }
        let mut out = Image::new(w, h);

        for y in 0..h {
            for x in 0..w {
                let depth = image.get_depth(x, y);
                let dist = (depth - self.focal_distance).abs();

                if dist <= self.dof_range {
                    // In-focus — copy directly.
                    out.set(x, y, image.get(x, y));
                    let di = out.idx(x, y);
                    out.depth[di] = depth;
                    continue;
                }

                // Compute CoC radius in pixels.
                let coc =
                    ((dist - self.dof_range) / (1.0 - self.dof_range).max(1e-7)) * self.dof_max_coc;
                let radius = coc.clamp(0.0, self.dof_max_coc) as isize;

                // Variable-radius box blur.
                let mut acc = PostColor::zero();
                let mut weight = 0.0_f32;
                let ix = x as isize;
                let iy = y as isize;
                for dy in -radius..=radius {
                    for dx in -radius..=radius {
                        let sx = (ix + dx).clamp(0, w as isize - 1) as usize;
                        let sy = (iy + dy).clamp(0, h as isize - 1) as usize;
                        acc = acc.add(&image.get(sx, sy).scale(1.0));
                        weight += 1.0;
                    }
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

    // -----------------------------------------------------------------------
    // Analytical lens flare
    // -----------------------------------------------------------------------

    /// Add an analytical lens flare to `image`.
    ///
    /// The flare is rendered as a radial starburst pattern centred at
    /// `flare_pos` (in normalised screen-space `[0, 1]²`).  Each of the
    /// `flare_streak_count` streak arms contributes a linearly-attenuated glow
    /// that decays with distance from the source.
    ///
    /// The returned image is `input + flare`, so it may exceed `[0, 1]` in
    /// bright regions (callers should tone-map afterwards).
    #[allow(clippy::too_many_arguments)]
    pub fn compute_lens_flare(&self, image: &Image, flare_pos: [f32; 2]) -> Image {
        let w = image.width;
        let h = image.height;
        if w == 0 || h == 0 || self.flare_intensity <= 0.0 {
            return Image::new(w, h);
        }
        let mut out = image.pixels.clone();

        let cx = flare_pos[0];
        let cy = flare_pos[1];
        let streaks = self.flare_streak_count.max(1);
        let intensity = self.flare_intensity;

        for y in 0..h {
            for x in 0..w {
                let ux = x as f32 / w as f32 - cx;
                let uy = y as f32 / h as f32 - cy;
                let dist = (ux * ux + uy * uy).sqrt();

                // Halo: radially attenuated ring.
                let halo_r = 0.25_f32;
                let halo_w = 0.05_f32;
                let halo = (1.0 - ((dist - halo_r) / halo_w).abs().min(1.0)) * 0.5;

                // Streaks: angular dependency.
                let angle = uy.atan2(ux);
                let mut streak_sum = 0.0_f32;
                for k in 0..streaks {
                    let streak_angle = k as f32 * std::f32::consts::PI / streaks as f32;
                    let angular_diff = (angle - streak_angle).sin().abs().min(1.0);
                    let streak = (1.0 - angular_diff.powf(0.1)) * (-dist * 6.0).exp();
                    streak_sum += streak;
                }
                streak_sum /= streaks as f32;

                let flare_val = (halo + streak_sum * 0.5) * intensity;
                let idx = y * w + x;
                let c = out[idx];
                out[idx] = PostColor::new(
                    c.r + flare_val * 1.0,
                    c.g + flare_val * 0.85,
                    c.b + flare_val * 0.6,
                    c.a,
                );
            }
        }

        Image {
            width: w,
            height: h,
            pixels: out,
            depth: image.depth.clone(),
        }
    }

    // -----------------------------------------------------------------------
    // Chromatic aberration
    // -----------------------------------------------------------------------

    /// Apply chromatic aberration by independently shifting the red and blue
    /// channels horizontally.
    ///
    /// Uses `self.ca_r_shift` and `self.ca_b_shift` (pixel offsets).  Green
    /// remains at the original position.  A zero shift leaves the channel
    /// unchanged.
    pub fn compute_chromatic_aberration(&self, image: &Image) -> Image {
        let ca = ChromaticAberrationEffect {
            r_shift_x: self.ca_r_shift,
            r_shift_y: 0,
            b_shift_x: self.ca_b_shift,
            b_shift_y: 0,
        };
        ca.apply(image)
    }
}

#[allow(dead_code)]
impl Default for PostProcessor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── PostProcessor ────────────────────────────────────────────────────────

    #[test]
    fn test_post_processor_default_constructs() {
        let pp = PostProcessor::new();
        assert!(pp.focal_distance > 0.0);
        assert_eq!(pp.ca_r_shift, 0);
        assert_eq!(pp.ca_b_shift, 0);
    }

    #[test]
    fn test_dof_in_focus_unchanged() {
        let mut img = Image::new(4, 4);
        let c = PostColor::new(0.8, 0.4, 0.2, 1.0);
        for y in 0..4 {
            for x in 0..4 {
                img.set(x, y, c);
            }
        }
        // Set all depths to focal distance so no blur occurs.
        let pp = PostProcessor {
            focal_distance: 0.5,
            dof_range: 0.3,
            dof_max_coc: 4.0,
            ..PostProcessor::new()
        };
        for i in 0..16 {
            img.depth[i] = 0.5;
        }
        let out = pp.compute_depth_of_field(&img);
        for y in 0..4 {
            for x in 0..4 {
                let o = out.get(x, y);
                assert!(
                    (o.r - c.r).abs() < 1e-5,
                    "In-focus pixel r changed: {} vs {}",
                    o.r,
                    c.r
                );
            }
        }
    }

    #[test]
    fn test_dof_out_of_focus_blurs() {
        let w = 5;
        let h = 5;
        let mut img = Image::new(w, h);
        // All black except centre pixel is bright.
        img.set(2, 2, PostColor::new(1.0, 1.0, 1.0, 1.0));
        // Set all depths far from focus.
        for i in 0..w * h {
            img.depth[i] = 1.0;
        }
        let pp = PostProcessor {
            focal_distance: 0.0,
            dof_range: 0.05,
            dof_max_coc: 2.0,
            ..PostProcessor::new()
        };
        let out = pp.compute_depth_of_field(&img);
        // The bright pixel should have spread to neighbours.
        let total: f32 = (0..w * h)
            .map(|i| {
                let x = i % w;
                let y = i / w;
                out.get(x, y).r
            })
            .sum();
        assert!(
            total > 0.9,
            "DoF blur should spread brightness (total={total})"
        );
    }

    #[test]
    fn test_dof_preserves_dimensions() {
        let img = Image::new(8, 6);
        let pp = PostProcessor::new();
        let out = pp.compute_depth_of_field(&img);
        assert_eq!(out.width, 8);
        assert_eq!(out.height, 6);
    }

    #[test]
    fn test_dof_empty_image() {
        let img = Image::new(0, 0);
        let pp = PostProcessor::new();
        let out = pp.compute_depth_of_field(&img);
        assert_eq!(out.width, 0);
        assert_eq!(out.height, 0);
    }

    #[test]
    fn test_lens_flare_zero_intensity_no_change() {
        let img = Image::new(8, 8);
        let pp = PostProcessor {
            flare_intensity: 0.0,
            ..PostProcessor::new()
        };
        let out = pp.compute_lens_flare(&img, [0.5, 0.5]);
        // Zero intensity → same as fresh black image.
        for y in 0..8 {
            for x in 0..8 {
                assert!((out.get(x, y).r).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_lens_flare_nonzero_intensity_brightens() {
        let img = Image::new(16, 16);
        let pp = PostProcessor {
            flare_intensity: 1.0,
            flare_streak_count: 4,
            ..PostProcessor::new()
        };
        let out = pp.compute_lens_flare(&img, [0.5, 0.5]);
        // At least some pixels near the flare centre should be brighter than 0.
        let centre = out.get(8, 8);
        let bright = centre.r + centre.g + centre.b;
        assert!(bright > 0.0, "Lens flare should add brightness near centre");
    }

    #[test]
    fn test_lens_flare_preserves_dimensions() {
        let img = Image::new(10, 7);
        let pp = PostProcessor {
            flare_intensity: 0.5,
            ..PostProcessor::new()
        };
        let out = pp.compute_lens_flare(&img, [0.5, 0.5]);
        assert_eq!(out.width, 10);
        assert_eq!(out.height, 7);
    }

    #[test]
    fn test_chromatic_aberration_zero_shift_unchanged() {
        let w = 6;
        let h = 4;
        let mut img = Image::new(w, h);
        for y in 0..h {
            for x in 0..w {
                img.set(x, y, PostColor::new(0.3, 0.6, 0.9, 1.0));
            }
        }
        let pp = PostProcessor {
            ca_r_shift: 0,
            ca_b_shift: 0,
            ..PostProcessor::new()
        };
        let out = pp.compute_chromatic_aberration(&img);
        for y in 0..h {
            for x in 0..w {
                let o = out.get(x, y);
                assert!(
                    (o.g - 0.6).abs() < 1e-5,
                    "Green should not change with zero shift"
                );
            }
        }
    }

    #[test]
    fn test_chromatic_aberration_nonzero_shift_changes_channels() {
        let w = 8;
        let h = 4;
        let mut img = Image::new(w, h);
        // Set alternating columns to different R values.
        for y in 0..h {
            for x in 0..w {
                let r = if x % 2 == 0 { 0.0_f32 } else { 1.0_f32 };
                img.set(x, y, PostColor::new(r, 0.5, r, 1.0));
            }
        }
        let pp = PostProcessor {
            ca_r_shift: 1,
            ca_b_shift: -1,
            ..PostProcessor::new()
        };
        let out = pp.compute_chromatic_aberration(&img);
        // Verify some pixel has a different R than the original at its position.
        let orig_r = img.get(2, 2).r;
        let out_r = out.get(2, 2).r;
        // With shift=1 the red channel at x=2 comes from x=1, which is 1.0 != 0.0.
        assert!(
            (orig_r - out_r).abs() > 0.5,
            "CA should shift R channel (orig={orig_r} out={out_r})"
        );
    }
}
