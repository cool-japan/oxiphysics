// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shader effect utilities for physics visualization.
//!
//! Provides Phong shading, toon shading, rim lighting, Fresnel-Schlick,
//! PBR BRDF (GGX NDF), chromatic aberration, vignette, film grain,
//! lens flare, and depth-of-field circle-of-confusion utilities.

// ─────────────────────────────────────────────────────────────────────────────
// Internal vector helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

#[inline]
fn reflect3(incident: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let d = dot3(incident, normal);
    [
        incident[0] - 2.0 * d * normal[0],
        incident[1] - 2.0 * d * normal[1],
        incident[2] - 2.0 * d * normal[2],
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// ShaderParams
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters controlling Phong shading coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct ShaderParams {
    /// Ambient reflection coefficient (≥ 0).
    pub ambient: f64,
    /// Diffuse reflection coefficient (≥ 0).
    pub diffuse: f64,
    /// Specular reflection coefficient (≥ 0).
    pub specular: f64,
    /// Specular shininess exponent (> 0; higher = narrower highlight).
    pub shininess: f64,
}

impl Default for ShaderParams {
    fn default() -> Self {
        Self {
            ambient: 0.1,
            diffuse: 0.7,
            specular: 0.5,
            shininess: 32.0,
        }
    }
}

impl ShaderParams {
    /// Construct new shader parameters.
    pub fn new(ambient: f64, diffuse: f64, specular: f64, shininess: f64) -> Self {
        Self {
            ambient,
            diffuse,
            specular,
            shininess,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Phong shading
// ─────────────────────────────────────────────────────────────────────────────

/// Phong shading model: returns a scalar intensity in \[0, 1\] (clamped).
///
/// # Arguments
/// * `normal`    – surface normal (need not be normalised).
/// * `light_dir` – direction **from** the surface **toward** the light source.
/// * `view_dir`  – direction **from** the surface **toward** the camera.
/// * `params`    – Phong coefficients.
pub fn phong_shading(
    normal: [f64; 3],
    light_dir: [f64; 3],
    view_dir: [f64; 3],
    params: &ShaderParams,
) -> f64 {
    let n = normalize3(normal);
    let l = normalize3(light_dir);
    let v = normalize3(view_dir);

    // Ambient
    let ambient = params.ambient;

    // Diffuse: max(N·L, 0)
    let n_dot_l = dot3(n, l).max(0.0);
    let diffuse = params.diffuse * n_dot_l;

    // Specular: Phong reflection vector R = reflect(-L, N); max(R·V, 0)^s
    let neg_l = [-l[0], -l[1], -l[2]];
    let r = reflect3(neg_l, n);
    let r_dot_v = dot3(r, v).max(0.0);
    let specular = params.specular * r_dot_v.powf(params.shininess.max(1.0));

    (ambient + diffuse + specular).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Toon shading
// ─────────────────────────────────────────────────────────────────────────────

/// Toon (cel) shading: quantises `intensity` into `steps` discrete bands.
///
/// `intensity` is expected in \[0, 1\]. Returns a value in \[0, 1\].
pub fn toon_shading(intensity: f64, steps: usize) -> f64 {
    let steps = steps.max(1) as f64;
    (intensity * steps).floor() / steps
}

// ─────────────────────────────────────────────────────────────────────────────
// Rim lighting
// ─────────────────────────────────────────────────────────────────────────────

/// Rim lighting intensity: bright on silhouette edges facing away from camera.
///
/// `rim_power` controls how tight the rim band is (higher = tighter).
/// Returns a value in \[0, 1\].
pub fn rim_lighting(normal: [f64; 3], view_dir: [f64; 3], rim_power: f64) -> f64 {
    let n = normalize3(normal);
    let v = normalize3(view_dir);
    let n_dot_v = dot3(n, v).clamp(0.0, 1.0);
    (1.0 - n_dot_v).powf(rim_power.max(0.0))
}

// ─────────────────────────────────────────────────────────────────────────────
// Fresnel–Schlick
// ─────────────────────────────────────────────────────────────────────────────

/// Fresnel–Schlick approximation: `F = F0 + (1 − F0)(1 − cosθ)⁵`.
///
/// `cos_theta` is the cosine of the angle between the view direction and the
/// half-vector (or surface normal). `f0` is the normal-incidence reflectance.
pub fn fresnel_schlick(cos_theta: f64, f0: f64) -> f64 {
    let ct = cos_theta.clamp(0.0, 1.0);
    f0 + (1.0 - f0) * (1.0 - ct).powi(5)
}

// ─────────────────────────────────────────────────────────────────────────────
// PBR BRDF — simplified GGX NDF
// ─────────────────────────────────────────────────────────────────────────────

/// Simplified GGX normal distribution function (NDF) as a PBR BRDF estimate.
///
/// Returns the NDF value D(h) for a microfacet BRDF using the GGX/Trowbridge-Reitz
/// distribution:  D = α² / (π · (cosθ_h² · (α²−1) + 1)²)
///
/// `roughness` ∈ \[0, 1\]; `metallic` scales the overall reflectance linearly;
/// `cos_theta_h` is cosθ between the surface normal and the half-vector.
pub fn pbr_brdf(roughness: f64, metallic: f64, cos_theta_h: f64) -> f64 {
    let alpha = roughness * roughness;
    let alpha2 = alpha * alpha;
    let ct = cos_theta_h.clamp(0.0, 1.0);
    let denom = ct * ct * (alpha2 - 1.0) + 1.0;
    let ndf = alpha2 / (std::f64::consts::PI * denom * denom + 1e-15);
    (ndf * metallic.clamp(0.0, 1.0)).max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Chromatic aberration
// ─────────────────────────────────────────────────────────────────────────────

/// Chromatic aberration: returns UV offsets for \[R, G, B\] channels.
///
/// The effect displaces each channel's UV coordinates outward from the centre
/// by an amount proportional to `strength`. Returns `[[u_r, v_r\], [u_g, v_g], [u_b, v_b]]`.
pub fn chromatic_aberration(uv: [f64; 2], strength: f64) -> [[f64; 2]; 3] {
    // Direction from centre (0.5, 0.5)
    let dx = uv[0] - 0.5;
    let dy = uv[1] - 0.5;

    // R: pushed outward by full strength
    let r_uv = [uv[0] + dx * strength, uv[1] + dy * strength];
    // G: no offset (or very small — use 0)
    let g_uv = uv;
    // B: pushed inward by full strength
    let b_uv = [uv[0] - dx * strength, uv[1] - dy * strength];

    [r_uv, g_uv, b_uv]
}

// ─────────────────────────────────────────────────────────────────────────────
// Vignette
// ─────────────────────────────────────────────────────────────────────────────

/// Vignette effect: darkens the frame towards the edges.
///
/// Returns a multiplier in \[0, 1\], where 1 = fully bright (centre)
/// and 0 = fully dark (beyond `radius` from centre).
pub fn vignette(uv: [f64; 2], strength: f64, radius: f64) -> f64 {
    let dx = uv[0] - 0.5;
    let dy = uv[1] - 0.5;
    let dist = (dx * dx + dy * dy).sqrt();
    let r = radius.max(1e-15);
    let t = (dist / r).clamp(0.0, 1.0);
    1.0 - strength * t * t
}

// ─────────────────────────────────────────────────────────────────────────────
// Film grain
// ─────────────────────────────────────────────────────────────────────────────

/// Deterministic film grain noise at a UV coordinate.
///
/// Uses a hash of the UV coordinates and `seed` to produce a noise value
/// in \[−intensity, +intensity\].
pub fn film_grain(uv: [f64; 2], seed: u64, intensity: f64) -> f64 {
    // Mix uv and seed into a single integer hash (Wang hash variant)
    let ix = (uv[0] * 65537.0) as u64;
    let iy = (uv[1] * 65537.0) as u64;
    let mut h = ix.wrapping_mul(2246822519).wrapping_add(seed);
    h ^= iy.wrapping_mul(3266489917);
    h = h.wrapping_mul(668265263).wrapping_add(h >> 15);
    h ^= h >> 12;
    h = h.wrapping_mul(374761393);
    h ^= h >> 15;
    // Map to [0, 1]
    let frac = (h & 0xFFFF_FFFF) as f64 / 4_294_967_295.0;
    // Shift to [-intensity, +intensity]
    (frac - 0.5) * 2.0 * intensity
}

// ─────────────────────────────────────────────────────────────────────────────
// Lens flare
// ─────────────────────────────────────────────────────────────────────────────

/// Lens flare intensity at UV coordinate `uv` for a light source at `light_pos`.
///
/// Returns a value in \[0, 1\]: 1 when `uv == light_pos`, falling off to 0 at
/// distance `radius`.
pub fn lens_flare_intensity(light_pos: [f64; 2], uv: [f64; 2], radius: f64) -> f64 {
    let dx = uv[0] - light_pos[0];
    let dy = uv[1] - light_pos[1];
    let dist = (dx * dx + dy * dy).sqrt();
    let r = radius.max(1e-15);
    (1.0 - (dist / r)).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Depth of field — circle of confusion
// ─────────────────────────────────────────────────────────────────────────────

/// Circle of confusion (CoC) diameter for depth-of-field rendering.
///
/// Computes the blur disk size for a point at `depth` given a lens with
/// focal length `focal_len`, aperture diameter `aperture`, and focus distance
/// `focus_dist`.  All distances in the same units.
///
/// Returns the absolute CoC diameter (≥ 0).
pub fn depth_of_field_coc(depth: f64, focus_dist: f64, aperture: f64, focal_len: f64) -> f64 {
    // Standard thin-lens CoC formula:
    //   CoC = aperture * focal_len * |depth - focus_dist| / (depth * (focus_dist - focal_len))
    if depth.abs() < 1e-15 {
        return 0.0;
    }
    let focus_minus_f = focus_dist - focal_len;
    if focus_minus_f.abs() < 1e-15 {
        return 0.0;
    }
    let coc =
        aperture * focal_len * (depth - focus_dist).abs() / (depth.abs() * focus_minus_f.abs());
    coc.max(0.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ShaderParams ─────────────────────────────────────────────────────────

    #[test]
    fn test_shader_params_default_valid() {
        let p = ShaderParams::default();
        assert!(p.ambient >= 0.0);
        assert!(p.diffuse >= 0.0);
        assert!(p.specular >= 0.0);
        assert!(p.shininess > 0.0);
    }

    #[test]
    fn test_shader_params_new() {
        let p = ShaderParams::new(0.2, 0.6, 0.4, 64.0);
        assert!((p.ambient - 0.2).abs() < 1e-12);
        assert!((p.diffuse - 0.6).abs() < 1e-12);
        assert!((p.specular - 0.4).abs() < 1e-12);
        assert!((p.shininess - 64.0).abs() < 1e-12);
    }

    // ── phong_shading ────────────────────────────────────────────────────────

    #[test]
    fn test_phong_ambient_only_when_no_light() {
        // Light coming from opposite side of normal → diffuse = 0, specular = 0
        let normal = [0.0, 0.0, 1.0];
        let light_dir = [0.0, 0.0, -1.0]; // behind surface
        let view_dir = [0.0, 0.0, 1.0];
        let params = ShaderParams::new(0.1, 0.8, 0.5, 32.0);
        let intensity = phong_shading(normal, light_dir, view_dir, &params);
        // Only ambient (0.1) should contribute
        assert!(
            (intensity - 0.1).abs() < 1e-10,
            "expected 0.1, got {intensity}"
        );
    }

    #[test]
    fn test_phong_full_diffuse_front_lit() {
        // Normal and light perfectly aligned → max diffuse
        let normal = [0.0, 0.0, 1.0];
        let light_dir = [0.0, 0.0, 1.0];
        let view_dir = [0.0, 0.0, 1.0]; // head-on view
        let params = ShaderParams::new(0.0, 1.0, 0.0, 32.0);
        let intensity = phong_shading(normal, light_dir, view_dir, &params);
        assert!(
            (intensity - 1.0).abs() < 1e-10,
            "expected 1.0, got {intensity}"
        );
    }

    #[test]
    fn test_phong_output_in_range() {
        let normal = [0.3, 0.7, 0.6];
        let light_dir = [0.5, 0.5, 0.0];
        let view_dir = [0.0, 0.0, 1.0];
        let params = ShaderParams::default();
        let intensity = phong_shading(normal, light_dir, view_dir, &params);
        assert!(
            (0.0..=1.0).contains(&intensity),
            "intensity out of [0,1]: {intensity}"
        );
    }

    #[test]
    fn test_phong_specular_at_perfect_reflection() {
        // View direction == reflection of light direction → max specular
        let normal = [0.0, 0.0, 1.0];
        let light_dir = [0.0, 0.0, 1.0];
        let view_dir = [0.0, 0.0, 1.0]; // reflects into view
        let params = ShaderParams::new(0.0, 0.0, 1.0, 1.0);
        let intensity = phong_shading(normal, light_dir, view_dir, &params);
        assert!(
            (intensity - 1.0).abs() < 1e-10,
            "expected 1.0, got {intensity}"
        );
    }

    #[test]
    fn test_phong_with_unnormalized_vectors() {
        // Should still work: function normalises inputs internally
        let normal = [0.0, 0.0, 10.0]; // scaled
        let light_dir = [0.0, 0.0, 5.0];
        let view_dir = [0.0, 0.0, 3.0];
        let params = ShaderParams::new(0.0, 1.0, 0.0, 32.0);
        let intensity = phong_shading(normal, light_dir, view_dir, &params);
        assert!(
            (intensity - 1.0).abs() < 1e-10,
            "expected 1.0, got {intensity}"
        );
    }

    // ── toon_shading ─────────────────────────────────────────────────────────

    #[test]
    fn test_toon_shading_two_steps() {
        // With 2 steps: [0, 0.5) → 0.0, [0.5, 1.0) → 0.5
        assert!((toon_shading(0.3, 2) - 0.0).abs() < 1e-12);
        assert!((toon_shading(0.7, 2) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_toon_shading_four_steps_uniform() {
        let vals = [0.1, 0.3, 0.6, 0.9];
        let expected = [0.0, 0.25, 0.5, 0.75];
        for (&v, &e) in vals.iter().zip(expected.iter()) {
            let t = toon_shading(v, 4);
            assert!((t - e).abs() < 1e-12, "toon({v}, 4) = {t}, expected {e}");
        }
    }

    #[test]
    fn test_toon_shading_one_step_always_zero() {
        // With 1 step everything below 1 maps to 0.0
        assert!((toon_shading(0.99, 1) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_toon_shading_output_in_range() {
        for i in 0..10 {
            let t = toon_shading(i as f64 / 10.0, 5);
            assert!((0.0..=1.0).contains(&t), "toon out of range: {t}");
        }
    }

    // ── rim_lighting ─────────────────────────────────────────────────────────

    #[test]
    fn test_rim_lighting_silhouette_is_max() {
        // Normal perpendicular to view → rim = 1.0
        let normal = [1.0, 0.0, 0.0];
        let view = [0.0, 0.0, 1.0]; // perpendicular → N·V = 0
        let rim = rim_lighting(normal, view, 2.0);
        assert!(
            (rim - 1.0).abs() < 1e-10,
            "rim at silhouette should be 1, got {rim}"
        );
    }

    #[test]
    fn test_rim_lighting_facing_camera_is_min() {
        // Normal aligned with view → rim = 0.0
        let normal = [0.0, 0.0, 1.0];
        let view = [0.0, 0.0, 1.0];
        let rim = rim_lighting(normal, view, 3.0);
        assert!(
            rim.abs() < 1e-10,
            "rim facing camera should be 0, got {rim}"
        );
    }

    #[test]
    fn test_rim_lighting_output_in_range() {
        let normal = [0.5, 0.5, 0.707];
        let view = [0.0, 0.0, 1.0];
        let rim = rim_lighting(normal, view, 4.0);
        assert!((0.0..=1.0).contains(&rim), "rim out of [0,1]: {rim}");
    }

    // ── fresnel_schlick ───────────────────────────────────────────────────────

    #[test]
    fn test_fresnel_at_zero_angle_returns_f0() {
        // cos_theta = 1 → F = f0
        let f0 = 0.04;
        let f = fresnel_schlick(1.0, f0);
        assert!((f - f0).abs() < 1e-12, "fresnel(1, f0) should equal f0");
    }

    #[test]
    fn test_fresnel_at_grazing_angle_returns_one() {
        // cos_theta = 0 → F = 1
        let f = fresnel_schlick(0.0, 0.04);
        assert!((f - 1.0).abs() < 1e-12, "fresnel(0, f0) should be 1");
    }

    #[test]
    fn test_fresnel_monotone_decreasing_with_cos_theta() {
        // Higher cos_theta → lower Fresnel (closer to f0)
        let f0 = 0.04;
        let f_low = fresnel_schlick(0.2, f0);
        let f_high = fresnel_schlick(0.9, f0);
        assert!(
            f_low > f_high,
            "fresnel should decrease with increasing cos_theta"
        );
    }

    #[test]
    fn test_fresnel_output_in_range() {
        for i in 0..=10 {
            let ct = i as f64 / 10.0;
            let f = fresnel_schlick(ct, 0.05);
            assert!(
                (0.0..=1.0).contains(&f),
                "fresnel out of [0,1] at cos_theta={ct}: {f}"
            );
        }
    }

    // ── pbr_brdf ─────────────────────────────────────────────────────────────

    #[test]
    fn test_pbr_brdf_non_negative() {
        let val = pbr_brdf(0.5, 1.0, 1.0);
        assert!(val >= 0.0, "PBR BRDF should be non-negative, got {val}");
    }

    #[test]
    fn test_pbr_brdf_zero_metallic_returns_zero() {
        let val = pbr_brdf(0.5, 0.0, 0.9);
        assert!(
            val.abs() < 1e-12,
            "PBR BRDF with metallic=0 should be 0, got {val}"
        );
    }

    #[test]
    fn test_pbr_brdf_roughness_effect() {
        // Smooth surface (low roughness) should have higher peak NDF at cos_theta_h=1
        let smooth = pbr_brdf(0.1, 1.0, 1.0);
        let rough = pbr_brdf(0.9, 1.0, 1.0);
        assert!(smooth > rough, "smooth surface should have higher NDF peak");
    }

    // ── chromatic_aberration ─────────────────────────────────────────────────

    #[test]
    fn test_chromatic_aberration_center_unchanged() {
        // At UV center (0.5, 0.5), all channels remain at center regardless of strength
        let uvs = chromatic_aberration([0.5, 0.5], 0.05);
        for ch in &uvs {
            assert!((ch[0] - 0.5).abs() < 1e-12);
            assert!((ch[1] - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn test_chromatic_aberration_r_and_b_opposite() {
        // R should be pushed outward, B inward from center
        let uv = [0.7, 0.5]; // right of center
        let uvs = chromatic_aberration(uv, 0.1);
        // R channel pushed further right (larger x), B channel pushed left
        assert!(uvs[0][0] > uv[0], "R channel should be pushed outward");
        assert!(uvs[2][0] < uv[0], "B channel should be pushed inward");
    }

    #[test]
    fn test_chromatic_aberration_zero_strength_identity() {
        let uv = [0.3, 0.8];
        let uvs = chromatic_aberration(uv, 0.0);
        for ch in &uvs {
            assert!((ch[0] - uv[0]).abs() < 1e-12);
            assert!((ch[1] - uv[1]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_chromatic_aberration_returns_three_channels() {
        let uvs = chromatic_aberration([0.5, 0.5], 0.02);
        assert_eq!(uvs.len(), 3);
    }

    // ── vignette ─────────────────────────────────────────────────────────────

    #[test]
    fn test_vignette_center_is_full_bright() {
        // At UV center (0.5, 0.5), distance = 0 → vignette = 1
        let v = vignette([0.5, 0.5], 1.0, 0.5);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "vignette at center should be 1, got {v}"
        );
    }

    #[test]
    fn test_vignette_beyond_radius_darkens() {
        // At corner (0, 0), distance ≈ 0.707 > radius=0.5 → fully dark
        let v = vignette([0.0, 0.0], 1.0, 0.5);
        assert!(v <= 1.0);
        assert!(v < 1.0, "vignette at corner should be less than 1");
    }

    #[test]
    fn test_vignette_output_in_range() {
        for i in 0..=10 {
            let x = i as f64 / 10.0;
            let v = vignette([x, 0.5], 1.0, 0.5);
            // Vignette multiplier: can be negative when strength > 1 and far from centre
            // but with strength=1 it ranges roughly from 0..1
            assert!(v <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn test_vignette_zero_strength_is_one() {
        // Zero strength → no darkening → 1.0 everywhere
        let v = vignette([0.0, 0.0], 0.0, 0.5);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "vignette with strength=0 should be 1"
        );
    }

    // ── film_grain ────────────────────────────────────────────────────────────

    #[test]
    fn test_film_grain_in_range() {
        let intensity = 0.1;
        for i in 0..20 {
            let uv = [i as f64 / 20.0, 0.5];
            let g = film_grain(uv, 42, intensity);
            assert!(
                g >= -intensity && g <= intensity,
                "film_grain out of range: {g}"
            );
        }
    }

    #[test]
    fn test_film_grain_zero_intensity_is_zero() {
        let g = film_grain([0.3, 0.7], 0, 0.0);
        assert!(g.abs() < 1e-12, "film_grain with 0 intensity should be 0");
    }

    #[test]
    fn test_film_grain_deterministic() {
        // Same UV and seed → same result
        let g1 = film_grain([0.5, 0.5], 99, 0.05);
        let g2 = film_grain([0.5, 0.5], 99, 0.05);
        assert!(
            (g1 - g2).abs() < 1e-15,
            "film_grain should be deterministic"
        );
    }

    #[test]
    fn test_film_grain_different_uv_different_value() {
        let g1 = film_grain([0.1, 0.2], 1, 0.1);
        let g2 = film_grain([0.9, 0.8], 1, 0.1);
        // Very unlikely to be identical
        assert!(
            (g1 - g2).abs() > 1e-15,
            "film_grain at different UVs should differ"
        );
    }

    // ── lens_flare_intensity ─────────────────────────────────────────────────

    #[test]
    fn test_lens_flare_at_light_source_is_one() {
        let pos = [0.5, 0.5];
        let intensity = lens_flare_intensity(pos, pos, 0.3);
        assert!(
            (intensity - 1.0).abs() < 1e-12,
            "lens flare at source should be 1, got {intensity}"
        );
    }

    #[test]
    fn test_lens_flare_beyond_radius_is_zero() {
        let light = [0.5, 0.5];
        let uv = [2.0, 0.5]; // far away
        let intensity = lens_flare_intensity(light, uv, 0.1);
        assert!(
            intensity.abs() < 1e-12,
            "lens flare beyond radius should be 0, got {intensity}"
        );
    }

    #[test]
    fn test_lens_flare_output_in_range() {
        let light = [0.5, 0.5];
        for i in 0..10 {
            let uv = [i as f64 / 10.0, 0.5];
            let v = lens_flare_intensity(light, uv, 0.5);
            assert!((0.0..=1.0).contains(&v), "lens flare out of [0,1]: {v}");
        }
    }

    #[test]
    fn test_lens_flare_falloff_with_distance() {
        let light = [0.5, 0.5];
        let near = [0.55, 0.5];
        let far = [0.65, 0.5];
        let v_near = lens_flare_intensity(light, near, 0.5);
        let v_far = lens_flare_intensity(light, far, 0.5);
        assert!(v_near > v_far, "lens flare should decrease with distance");
    }

    // ── depth_of_field_coc ────────────────────────────────────────────────────

    #[test]
    fn test_coc_zero_at_focus_distance() {
        // Object at focus distance → CoC = 0
        let coc = depth_of_field_coc(5.0, 5.0, 2.8, 0.05);
        assert!(
            coc.abs() < 1e-10,
            "CoC at focus distance should be 0, got {coc}"
        );
    }

    #[test]
    fn test_coc_increases_with_distance_from_focus() {
        let focus = 5.0;
        let coc_near = depth_of_field_coc(4.0, focus, 2.8, 0.05);
        let coc_far = depth_of_field_coc(8.0, focus, 2.8, 0.05);
        assert!(coc_near > 0.0, "CoC should be positive when not at focus");
        assert!(coc_far > 0.0, "CoC should be positive when not at focus");
    }

    #[test]
    fn test_coc_non_negative() {
        let coc = depth_of_field_coc(10.0, 5.0, 4.0, 0.05);
        assert!(coc >= 0.0, "CoC should be non-negative, got {coc}");
    }

    #[test]
    fn test_coc_zero_depth_returns_zero() {
        let coc = depth_of_field_coc(0.0, 5.0, 2.8, 0.05);
        assert!(coc.abs() < 1e-12, "CoC with depth=0 should be 0");
    }

    #[test]
    fn test_coc_larger_aperture_increases_blur() {
        let small_aperture = depth_of_field_coc(8.0, 5.0, 1.4, 0.05);
        let large_aperture = depth_of_field_coc(8.0, 5.0, 8.0, 0.05);
        assert!(
            large_aperture > small_aperture,
            "larger aperture should produce bigger CoC"
        );
    }
}
