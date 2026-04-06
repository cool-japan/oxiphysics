//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
use std::f64::consts::PI;

use super::types::{
    AreaLight, Bvh, Camera, HitRecord, Material, MaterialType, PathState, PointLight, Ray,
    RenderConfig, Triangle,
};

/// Add two 3-D vectors.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-D vectors.
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-D vector by a scalar.
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Dot product of two 3-D vectors.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-D vectors.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Length of a 3-D vector.
pub fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3-D vector (returns zero vector if near zero).
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = length3(v);
    if len < 1e-15 {
        return [0.0; 3];
    }
    scale3(v, 1.0 / len)
}
/// Reflect direction `d` about normal `n` (both normalized).
pub fn reflect3(d: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    let dn2 = 2.0 * dot3(d, n);
    sub3(d, scale3(n, dn2))
}
/// Component-wise multiply two RGB colors.
pub fn mul_color(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}
/// Clamp each component of a color to \[0,1\].
pub fn clamp_color(c: [f64; 3]) -> [f64; 3] {
    [
        c[0].clamp(0.0, 1.0),
        c[1].clamp(0.0, 1.0),
        c[2].clamp(0.0, 1.0),
    ]
}
/// Phong shading model.
///
/// Computes diffuse + specular contribution from one point light.
pub fn phong_shading(
    hit: &HitRecord,
    light: &PointLight,
    view_dir: [f64; 3],
    mat: &Material,
    shadow: bool,
) -> [f64; 3] {
    if shadow {
        return [0.0; 3];
    }
    let light_vec = sub3(light.position, hit.position);
    let dist = length3(light_vec);
    let light_dir = normalize3(light_vec);
    let n_dot_l = dot3(hit.normal, light_dir).max(0.0);
    let diffuse = scale3(
        mul_color(mat.albedo, light.color),
        n_dot_l * light.intensity * light.attenuate(dist),
    );
    let reflect_dir = reflect3(scale3(light_dir, -1.0), hit.normal);
    let r_dot_v = dot3(reflect_dir, view_dir).max(0.0);
    let spec_factor = r_dot_v.powf(mat.shininess.max(1.0));
    let specular = scale3(
        mul_color(light.color, [1.0; 3]),
        spec_factor * light.intensity * light.attenuate(dist),
    );
    add3(diffuse, specular)
}
/// Schlick Fresnel approximation.
pub fn fresnel_schlick(cos_theta: f64, f0: [f64; 3]) -> [f64; 3] {
    let c = (1.0 - cos_theta).clamp(0.0, 1.0);
    let c5 = c * c * c * c * c;
    [
        f0[0] + (1.0 - f0[0]) * c5,
        f0[1] + (1.0 - f0[1]) * c5,
        f0[2] + (1.0 - f0[2]) * c5,
    ]
}
/// GGX normal distribution function.
pub fn distribution_ggx(n_dot_h: f64, roughness: f64) -> f64 {
    let a = roughness * roughness;
    let a2 = a * a;
    let n_dot_h2 = n_dot_h * n_dot_h;
    let denom = n_dot_h2 * (a2 - 1.0) + 1.0;
    if denom.abs() < 1e-15 {
        return 0.0;
    }
    a2 / (PI * denom * denom)
}
/// Smith's geometry function (GGX).
pub fn geometry_smith(n_dot_v: f64, n_dot_l: f64, roughness: f64) -> f64 {
    let r = roughness + 1.0;
    let k = r * r / 8.0;
    let ggx1 = n_dot_v / (n_dot_v * (1.0 - k) + k);
    let ggx2 = n_dot_l / (n_dot_l * (1.0 - k) + k);
    ggx1 * ggx2
}
/// Cook-Torrance PBR BRDF shading.
pub fn pbr_shading(
    hit: &HitRecord,
    light: &PointLight,
    view_dir: [f64; 3],
    mat: &Material,
    shadow: bool,
) -> [f64; 3] {
    if shadow {
        return [0.0; 3];
    }
    let light_vec = sub3(light.position, hit.position);
    let dist = length3(light_vec);
    let l = normalize3(light_vec);
    let n = hit.normal;
    let v = view_dir;
    let h = normalize3(add3(v, l));
    let n_dot_l = dot3(n, l).max(0.0);
    if n_dot_l < 1e-10 {
        return [0.0; 3];
    }
    let n_dot_v = dot3(n, v).max(1e-4);
    let n_dot_h = dot3(n, h).max(0.0);
    let h_dot_v = dot3(h, v).max(0.0);
    let f0_dielectric = [0.04; 3];
    let f0 = [
        f0_dielectric[0] * (1.0 - mat.metallic) + mat.albedo[0] * mat.metallic,
        f0_dielectric[1] * (1.0 - mat.metallic) + mat.albedo[1] * mat.metallic,
        f0_dielectric[2] * (1.0 - mat.metallic) + mat.albedo[2] * mat.metallic,
    ];
    let f = fresnel_schlick(h_dot_v, f0);
    let d = distribution_ggx(n_dot_h, mat.roughness);
    let g = geometry_smith(n_dot_v, n_dot_l, mat.roughness);
    let denom = 4.0 * n_dot_v * n_dot_l;
    let specular = if denom > 1e-10 {
        scale3(f, d * g / denom)
    } else {
        [0.0; 3]
    };
    let kd = [
        (1.0 - f[0]) * (1.0 - mat.metallic),
        (1.0 - f[1]) * (1.0 - mat.metallic),
        (1.0 - f[2]) * (1.0 - mat.metallic),
    ];
    let diffuse = [
        kd[0] * mat.albedo[0] / PI,
        kd[1] * mat.albedo[1] / PI,
        kd[2] * mat.albedo[2] / PI,
    ];
    let radiance = scale3(light.color, light.intensity * light.attenuate(dist));
    let result = add3(diffuse, specular);
    [
        result[0] * radiance[0] * n_dot_l,
        result[1] * radiance[1] * n_dot_l,
        result[2] * radiance[2] * n_dot_l,
    ]
}
/// Compute soft shadow factor (0=fully shadowed, 1=fully lit) by sampling
/// an area light with `num_samples` jittered samples.
pub fn soft_shadow_factor(
    hit_pos: [f64; 3],
    hit_normal: [f64; 3],
    light: &AreaLight,
    bvh: &Bvh,
    triangles: &[Triangle],
    num_samples: usize,
    samples: &[[f64; 2]],
) -> f64 {
    if num_samples == 0 || samples.is_empty() {
        return 1.0;
    }
    let actual_samples = num_samples.min(samples.len());
    let mut unblocked = 0u32;
    for i in 0..actual_samples {
        let su = samples[i][0] * 2.0 - 1.0;
        let sv = samples[i][1] * 2.0 - 1.0;
        let light_point = light.sample_point(su, sv);
        let to_light = sub3(light_point, hit_pos);
        let dist = length3(to_light);
        if dist < 1e-10 {
            unblocked += 1;
            continue;
        }
        let dir = scale3(to_light, 1.0 / dist);
        if dot3(hit_normal, dir) <= 0.0 {
            continue;
        }
        let mut shadow_ray = Ray::new(add3(hit_pos, scale3(hit_normal, 1e-4)), dir);
        shadow_ray.t_max = dist - 1e-4;
        if !bvh.intersect_any(&shadow_ray, triangles) {
            unblocked += 1;
        }
    }
    unblocked as f64 / actual_samples as f64
}
/// Compute ambient occlusion by casting rays in the hemisphere.
///
/// Returns a value in \[0, 1\] where 0 = fully occluded, 1 = fully open.
pub fn ambient_occlusion(
    hit_pos: [f64; 3],
    hit_normal: [f64; 3],
    bvh: &Bvh,
    triangles: &[Triangle],
    hemisphere_samples: &[[f64; 3]],
    max_dist: f64,
) -> f64 {
    if hemisphere_samples.is_empty() {
        return 1.0;
    }
    let tangent = if hit_normal[0].abs() < 0.9 {
        normalize3(cross3(hit_normal, [1.0, 0.0, 0.0]))
    } else {
        normalize3(cross3(hit_normal, [0.0, 1.0, 0.0]))
    };
    let bitangent = cross3(hit_normal, tangent);
    let mut unoccluded = 0u32;
    let n = hemisphere_samples.len();
    for s in hemisphere_samples {
        let world_dir = normalize3(add3(
            add3(scale3(tangent, s[0]), scale3(bitangent, s[1])),
            scale3(hit_normal, s[2].abs()),
        ));
        if dot3(world_dir, hit_normal) <= 0.0 {
            unoccluded += 1;
            continue;
        }
        let origin = add3(hit_pos, scale3(hit_normal, 1e-4));
        let mut ao_ray = Ray::new(origin, world_dir);
        ao_ray.t_max = max_dist;
        if !bvh.intersect_any(&ao_ray, triangles) {
            unoccluded += 1;
        }
    }
    unoccluded as f64 / n as f64
}
/// Schlick approximation for reflectance at a dielectric interface.
pub fn schlick_reflectance(cos_theta: f64, ior_ratio: f64) -> f64 {
    let r0 = ((1.0 - ior_ratio) / (1.0 + ior_ratio)).powi(2);
    r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
}
/// Compute refraction direction using Snell's law.
///
/// Returns `None` if total internal reflection occurs.
pub fn refract(d: [f64; 3], n: [f64; 3], ior_ratio: f64) -> Option<[f64; 3]> {
    let cos_theta = dot3(scale3(d, -1.0), n).min(1.0);
    let sin_theta_sq = 1.0 - cos_theta * cos_theta;
    if sin_theta_sq * ior_ratio * ior_ratio > 1.0 {
        return None;
    }
    let r_out_perp = scale3(add3(d, scale3(n, cos_theta)), ior_ratio);
    let r_out_parallel = scale3(n, -(1.0 - dot3(r_out_perp, r_out_perp)).abs().sqrt());
    Some(normalize3(add3(r_out_perp, r_out_parallel)))
}
/// Simple hemisphere sampler using cosine-weighted distribution.
///
/// Takes two uniform samples u1, u2 in \[0,1) and returns a direction
/// in the upper hemisphere (y > 0) in local space.
pub fn cosine_sample_hemisphere(u1: f64, u2: f64) -> [f64; 3] {
    let r = u1.sqrt();
    let theta = 2.0 * PI * u2;
    let x = r * theta.cos();
    let z = r * theta.sin();
    let y = (1.0 - u1).max(0.0).sqrt();
    [x, y, z]
}
/// Uniform hemisphere sampler.
pub fn uniform_sample_hemisphere(u1: f64, u2: f64) -> [f64; 3] {
    let cos_theta = u1;
    let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
    let phi = 2.0 * PI * u2;
    [sin_theta * phi.cos(), cos_theta, sin_theta * phi.sin()]
}
/// Trace a single path and accumulate radiance.
///
/// This is a simplified CPU-side path trace step (one bounce).
pub fn path_trace_step(
    state: &mut PathState,
    bvh: &Bvh,
    triangles: &[Triangle],
    materials: &[Material],
    background: [f64; 3],
    u1: f64,
    u2: f64,
) -> bool {
    if !state.should_continue() {
        if state.depth == 0 {
            state.radiance = background;
        }
        return false;
    }
    match bvh.intersect(&state.ray, triangles) {
        None => {
            let contrib = [
                state.throughput[0] * background[0],
                state.throughput[1] * background[1],
                state.throughput[2] * background[2],
            ];
            state.radiance = add3(state.radiance, contrib);
            false
        }
        Some((hit, _tri)) => {
            let mat = if (hit.material_id as usize) < materials.len() {
                &materials[hit.material_id as usize]
            } else {
                &materials[0]
            };
            let emission_contrib = [
                state.throughput[0] * mat.emission[0],
                state.throughput[1] * mat.emission[1],
                state.throughput[2] * mat.emission[2],
            ];
            state.radiance = add3(state.radiance, emission_contrib);
            let local_dir = cosine_sample_hemisphere(u1, u2);
            let tangent = if hit.normal[0].abs() < 0.9 {
                normalize3(cross3(hit.normal, [1.0, 0.0, 0.0]))
            } else {
                normalize3(cross3(hit.normal, [0.0, 1.0, 0.0]))
            };
            let bitangent = cross3(hit.normal, tangent);
            let world_dir = normalize3(add3(
                add3(
                    scale3(tangent, local_dir[0]),
                    scale3(bitangent, local_dir[2]),
                ),
                scale3(hit.normal, local_dir[1]),
            ));
            let cos_theta = dot3(hit.normal, world_dir).max(0.0);
            state.throughput = [
                state.throughput[0] * mat.albedo[0] * cos_theta * 2.0,
                state.throughput[1] * mat.albedo[1] * cos_theta * 2.0,
                state.throughput[2] * mat.albedo[2] * cos_theta * 2.0,
            ];
            state.ray = Ray::new(add3(hit.position, scale3(hit.normal, 1e-4)), world_dir);
            state.depth += 1;
            true
        }
    }
}
/// A-trous wavelet denoising pass for progressive path traced images.
///
/// `color` is a flat buffer of (r, g, b) triples in row-major order.
/// `width` and `height` are the image dimensions.
/// `step_width` is the kernel step (power of two: 1, 2, 4, 8, ...).
pub fn atrous_denoise(
    color: &[[f64; 3]],
    normal: &[[f64; 3]],
    position: &[[f64; 3]],
    width: usize,
    height: usize,
    step_width: usize,
    sigma_color: f64,
    sigma_normal: f64,
    sigma_position: f64,
) -> Vec<[f64; 3]> {
    let kernel = [
        [
            1.0f64 / 256.0,
            1.0 / 64.0,
            3.0 / 128.0,
            1.0 / 64.0,
            1.0 / 256.0,
        ],
        [1.0 / 64.0, 1.0 / 16.0, 3.0 / 32.0, 1.0 / 16.0, 1.0 / 64.0],
        [3.0 / 128.0, 3.0 / 32.0, 9.0 / 64.0, 3.0 / 32.0, 3.0 / 128.0],
        [1.0 / 64.0, 1.0 / 16.0, 3.0 / 32.0, 1.0 / 16.0, 1.0 / 64.0],
        [
            1.0 / 256.0,
            1.0 / 64.0,
            3.0 / 128.0,
            1.0 / 64.0,
            1.0 / 256.0,
        ],
    ];
    let n = width * height;
    let mut output = vec![[0.0f64; 3]; n];
    for py in 0..height {
        for px in 0..width {
            let idx = py * width + px;
            let c_center = color[idx];
            let n_center = normal[idx];
            let p_center = position[idx];
            let mut accum = [0.0f64; 3];
            let mut weight_sum = 0.0f64;
            for ky in 0..5i32 {
                for kx in 0..5i32 {
                    let oy = ky - 2;
                    let ox = kx - 2;
                    let nx = px as i32 + ox * step_width as i32;
                    let ny = py as i32 + oy * step_width as i32;
                    if nx < 0 || ny < 0 || nx >= width as i32 || ny >= height as i32 {
                        continue;
                    }
                    let sidx = ny as usize * width + nx as usize;
                    let c_s = color[sidx];
                    let n_s = normal[sidx];
                    let p_s = position[sidx];
                    let dc = [
                        c_center[0] - c_s[0],
                        c_center[1] - c_s[1],
                        c_center[2] - c_s[2],
                    ];
                    let dist_c = dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2];
                    let w_c = (-dist_c / (sigma_color * sigma_color)).exp();
                    let dn_x = n_center[0] - n_s[0];
                    let dn_y = n_center[1] - n_s[1];
                    let dn_z = n_center[2] - n_s[2];
                    let dist_n = dn_x * dn_x + dn_y * dn_y + dn_z * dn_z;
                    let w_n = (-dist_n / (sigma_normal * sigma_normal)).exp();
                    let dp_x = p_center[0] - p_s[0];
                    let dp_y = p_center[1] - p_s[1];
                    let dp_z = p_center[2] - p_s[2];
                    let dist_p = dp_x * dp_x + dp_y * dp_y + dp_z * dp_z;
                    let w_p = (-dist_p / (sigma_position * sigma_position)).exp();
                    let h_weight = kernel[ky as usize][kx as usize];
                    let w = h_weight * w_c * w_n * w_p;
                    accum[0] += w * c_s[0];
                    accum[1] += w * c_s[1];
                    accum[2] += w * c_s[2];
                    weight_sum += w;
                }
            }
            if weight_sum > 1e-10 {
                output[idx] = [
                    accum[0] / weight_sum,
                    accum[1] / weight_sum,
                    accum[2] / weight_sum,
                ];
            } else {
                output[idx] = c_center;
            }
        }
    }
    output
}
/// Temporal accumulation (exponential moving average) denoising.
///
/// Blends current frame with previous frame using history buffer.
pub fn temporal_accumulate(
    current: &[[f64; 3]],
    history: &[[f64; 3]],
    alpha: f64,
) -> Vec<[f64; 3]> {
    let n = current.len().min(history.len());
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let c = current[i];
        let h = history[i];
        result.push([
            alpha * c[0] + (1.0 - alpha) * h[0],
            alpha * c[1] + (1.0 - alpha) * h[1],
            alpha * c[2] + (1.0 - alpha) * h[2],
        ]);
    }
    result
}
/// Box filter denoising pass (simple spatial blur).
pub fn box_filter(color: &[[f64; 3]], width: usize, height: usize, radius: usize) -> Vec<[f64; 3]> {
    let n = width * height;
    let mut output = vec![[0.0f64; 3]; n];
    for py in 0..height {
        for px in 0..width {
            let mut accum = [0.0f64; 3];
            let mut count = 0u32;
            let y0 = py.saturating_sub(radius);
            let y1 = (py + radius + 1).min(height);
            let x0 = px.saturating_sub(radius);
            let x1 = (px + radius + 1).min(width);
            for sy in y0..y1 {
                for sx in x0..x1 {
                    let sidx = sy * width + sx;
                    accum[0] += color[sidx][0];
                    accum[1] += color[sidx][1];
                    accum[2] += color[sidx][2];
                    count += 1;
                }
            }
            let inv = 1.0 / count as f64;
            output[py * width + px] = [accum[0] * inv, accum[1] * inv, accum[2] * inv];
        }
    }
    output
}
/// Reinhard tone mapping operator.
pub fn tonemap_reinhard(color: [f64; 3]) -> [f64; 3] {
    [
        color[0] / (1.0 + color[0]),
        color[1] / (1.0 + color[1]),
        color[2] / (1.0 + color[2]),
    ]
}
/// Filmic tone mapping (Hejl and Burgess-Dawson approximation).
pub fn tonemap_filmic(color: [f64; 3]) -> [f64; 3] {
    let f = |x: f64| {
        let x = (x - 0.004).max(0.0);
        (x * (6.2 * x + 0.5)) / (x * (6.2 * x + 1.7) + 0.06)
    };
    [f(color[0]), f(color[1]), f(color[2])]
}
/// ACES filmic tone mapping.
pub fn tonemap_aces(color: [f64; 3]) -> [f64; 3] {
    let aces = |x: f64| {
        let a = 2.51;
        let b = 0.03;
        let c = 2.43;
        let d = 0.59;
        let e = 0.14;
        ((x * (a * x + b)) / (x * (c * x + d) + e)).clamp(0.0, 1.0)
    };
    [aces(color[0]), aces(color[1]), aces(color[2])]
}
/// Linear to sRGB gamma correction.
pub fn linear_to_srgb(c: f64) -> f64 {
    if c <= 0.0031308 {
        c * 12.92
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    }
}
/// Apply sRGB gamma to a color.
pub fn gamma_correct(color: [f64; 3]) -> [f64; 3] {
    [
        linear_to_srgb(color[0].clamp(0.0, 1.0)),
        linear_to_srgb(color[1].clamp(0.0, 1.0)),
        linear_to_srgb(color[2].clamp(0.0, 1.0)),
    ]
}
/// Ray trace a simple scene (direct illumination only) on the CPU.
///
/// Returns a flat buffer of (r, g, b) per pixel in row-major order.
pub fn render_direct(
    config: &RenderConfig,
    camera: &Camera,
    bvh: &Bvh,
    triangles: &[Triangle],
    materials: &[Material],
    lights: &[PointLight],
) -> Vec<[f64; 3]> {
    let n = config.width * config.height;
    let mut image = vec![[0.0f64; 3]; n];
    let w = config.width as f64;
    let h = config.height as f64;
    for py in 0..config.height {
        for px in 0..config.width {
            let ray = camera.generate_ray(px as f64, py as f64, w, h);
            let color = trace_direct(&ray, bvh, triangles, materials, lights, config);
            let idx = py * config.width + px;
            image[idx] = match config.tonemap {
                1 => gamma_correct(tonemap_reinhard(color)),
                2 => gamma_correct(tonemap_filmic(color)),
                3 => gamma_correct(tonemap_aces(color)),
                _ => gamma_correct(color),
            };
        }
    }
    image
}
/// Trace a single ray for direct illumination.
pub(super) fn trace_direct(
    ray: &Ray,
    bvh: &Bvh,
    triangles: &[Triangle],
    materials: &[Material],
    lights: &[PointLight],
    config: &RenderConfig,
) -> [f64; 3] {
    match bvh.intersect(ray, triangles) {
        None => config.background,
        Some((hit, _tri)) => {
            let mat = if (hit.material_id as usize) < materials.len() {
                &materials[hit.material_id as usize]
            } else {
                return config.background;
            };
            if mat.mat_type == MaterialType::Emissive {
                return mat.emission;
            }
            let view_dir = normalize3(scale3(ray.direction, -1.0));
            let mut color = config.ambient;
            for light in lights {
                let to_light = sub3(light.position, hit.position);
                let dist = length3(to_light);
                let light_dir = normalize3(to_light);
                let shadow_origin = add3(hit.position, scale3(hit.normal, 1e-4));
                let mut shadow_ray = Ray::new(shadow_origin, light_dir);
                shadow_ray.t_max = dist - 1e-4;
                let in_shadow = bvh.intersect_any(&shadow_ray, triangles);
                let contrib = if mat.mat_type == MaterialType::Pbr {
                    pbr_shading(&hit, light, view_dir, mat, in_shadow)
                } else {
                    phong_shading(&hit, light, view_dir, mat, in_shadow)
                };
                color = add3(color, contrib);
            }
            clamp_color(color)
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::raytracing::Aabb;
    use crate::raytracing::Scene;
    fn make_floor_triangle() -> Triangle {
        Triangle::new(
            [-5.0, 0.0, -5.0],
            [5.0, 0.0, -5.0],
            [0.0, 0.0, 5.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            0,
        )
    }
    #[test]
    fn test_ray_at() {
        let ray = Ray::new([0.0; 3], [0.0, 0.0, 1.0]);
        let p = ray.at(3.0);
        assert!((p[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_triangle_intersect_hit() {
        let tri = make_floor_triangle();
        let ray = Ray::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(tri.intersect(&ray).is_some());
    }
    #[test]
    fn test_triangle_intersect_miss_parallel() {
        let tri = make_floor_triangle();
        let ray = Ray::new([0.0, 5.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(tri.intersect(&ray).is_none());
    }
    #[test]
    fn test_triangle_intersect_miss_outside() {
        let tri = make_floor_triangle();
        let ray = Ray::new([100.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(tri.intersect(&ray).is_none());
    }
    #[test]
    fn test_triangle_geometric_normal() {
        let tri = Triangle::new(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            0,
        );
        let n = tri.geometric_normal();
        assert!((n[2] - 1.0).abs() < 1e-10 || (n[2] + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_triangle_area() {
        let tri = Triangle::new(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            0,
        );
        assert!((tri.area() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_aabb_merge() {
        let a = Aabb::new([0.0; 3], [1.0; 3]);
        let b = Aabb::new([-1.0; 3], [2.0; 3]);
        let merged = a.merge(&b);
        assert!((merged.min[0] + 1.0).abs() < 1e-12);
        assert!((merged.max[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_ray_hit() {
        let aabb = Aabb::new([-1.0; 3], [1.0; 3]);
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(aabb.intersect_ray(&ray).is_some());
    }
    #[test]
    fn test_aabb_ray_miss() {
        let aabb = Aabb::new([-1.0; 3], [1.0; 3]);
        let ray = Ray::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(aabb.intersect_ray(&ray).is_none());
    }
    #[test]
    fn test_aabb_longest_axis() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [3.0, 1.0, 2.0]);
        assert_eq!(aabb.longest_axis(), 0);
    }
    #[test]
    fn test_bvh_build_empty() {
        let bvh = Bvh::build(&[]);
        assert_eq!(bvh.prim_count, 0);
    }
    #[test]
    fn test_bvh_build_single_triangle() {
        let tri = make_floor_triangle();
        let bvh = Bvh::build(&[tri]);
        assert_eq!(bvh.prim_count, 1);
    }
    #[test]
    fn test_bvh_intersect_hit() {
        let tri = make_floor_triangle();
        let triangles = vec![tri];
        let bvh = Bvh::build(&triangles);
        let ray = Ray::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(bvh.intersect(&ray, &triangles).is_some());
    }
    #[test]
    fn test_bvh_intersect_miss() {
        let tri = make_floor_triangle();
        let triangles = vec![tri];
        let bvh = Bvh::build(&triangles);
        let ray = Ray::new([0.0, 5.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(bvh.intersect(&ray, &triangles).is_none());
    }
    #[test]
    fn test_bvh_intersect_any_shadow() {
        let tri = make_floor_triangle();
        let triangles = vec![tri];
        let bvh = Bvh::build(&triangles);
        let ray = Ray::new([0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(bvh.intersect_any(&ray, &triangles));
    }
    #[test]
    fn test_camera_generate_ray() {
        let cam = Camera::look_at(
            [0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            60.0,
            16.0 / 9.0,
            0.0,
            5.0,
        );
        let ray = cam.generate_ray(400.0, 300.0, 800.0, 600.0);
        assert!(ray.direction[2] < 0.0);
    }
    #[test]
    fn test_fresnel_schlick_zero_angle() {
        let f0 = [0.04; 3];
        let f = fresnel_schlick(0.0, f0);
        assert!((f[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fresnel_schlick_one_angle() {
        let f0 = [0.04; 3];
        let f = fresnel_schlick(1.0, f0);
        assert!((f[0] - 0.04).abs() < 1e-10);
    }
    #[test]
    fn test_distribution_ggx_smooth() {
        let d = distribution_ggx(1.0, 0.01);
        assert!(d > 100.0);
    }
    #[test]
    fn test_refract_no_tir() {
        let d = normalize3([0.0, -1.0, 0.0]);
        let n = [0.0, 1.0, 0.0];
        let result = refract(d, n, 1.0 / 1.5);
        assert!(result.is_some());
    }
    #[test]
    fn test_refract_tir() {
        let d = normalize3([0.9, -0.1, 0.0]);
        let n = [0.0, 1.0, 0.0];
        let result = refract(d, n, 1.5);
        assert!(result.is_none());
    }
    #[test]
    fn test_cosine_sample_hemisphere() {
        let s = cosine_sample_hemisphere(0.5, 0.5);
        let len = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
        assert!(s[1] >= 0.0);
    }
    #[test]
    fn test_tonemap_reinhard() {
        let c = tonemap_reinhard([2.0, 1.0, 0.5]);
        assert!((c[0] - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_tonemap_aces_clamp() {
        let c = tonemap_aces([100.0, 100.0, 100.0]);
        assert!(c[0] <= 1.0);
        assert!(c[0] >= 0.0);
    }
    #[test]
    fn test_linear_to_srgb_zero() {
        assert!((linear_to_srgb(0.0)).abs() < 1e-10);
    }
    #[test]
    fn test_linear_to_srgb_one() {
        assert!((linear_to_srgb(1.0) - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_box_filter_single_pixel() {
        let image = vec![[1.0f64, 0.0, 0.0]];
        let result = box_filter(&image, 1, 1, 1);
        assert_eq!(result.len(), 1);
        assert!((result[0][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_temporal_accumulate() {
        let current = vec![[1.0f64; 3]];
        let history = vec![[0.0f64; 3]];
        let result = temporal_accumulate(&current, &history, 0.5);
        assert!((result[0][0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_scene_add_box() {
        let mut scene = Scene::new();
        let mid = scene.add_material(Material::diffuse([0.8, 0.2, 0.2]));
        scene.add_box([0.0; 3], [1.0; 3], mid);
        assert_eq!(scene.triangles.len(), 12);
    }
    #[test]
    fn test_scene_build_and_intersect() {
        let mut scene = Scene::new();
        let mid = scene.add_material(Material::diffuse([0.8, 0.8, 0.8]));
        scene.add_box([0.0; 3], [1.0; 3], mid);
        scene.add_light(PointLight::new([5.0, 5.0, 5.0], [1.0; 3], 1.0));
        scene.build_bvh();
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        assert!(scene.intersect(&ray).is_some());
    }
    #[test]
    fn test_render_direct_small() {
        let mut scene = Scene::new();
        let mid = scene.add_material(Material::diffuse([0.8, 0.8, 0.8]));
        scene.add_box([0.0; 3], [1.0; 3], mid);
        scene.add_light(PointLight::new([5.0, 5.0, 5.0], [1.0; 3], 1.0));
        scene.build_bvh();
        let cam = Camera::look_at(
            [0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            60.0,
            4.0 / 3.0,
            0.0,
            5.0,
        );
        let config = RenderConfig {
            width: 4,
            height: 3,
            ..Default::default()
        };
        let image = render_direct(
            &config,
            &cam,
            scene.bvh.as_ref().unwrap(),
            &scene.triangles,
            &scene.materials,
            &scene.lights,
        );
        assert_eq!(image.len(), 12);
        for pixel in &image {
            for c in pixel {
                assert!(*c >= 0.0 && *c <= 1.0);
            }
        }
    }
    #[test]
    fn test_path_state_should_continue() {
        let ray = Ray::new([0.0; 3], [0.0, 0.0, 1.0]);
        let mut state = PathState::new(ray, 4);
        assert!(state.should_continue());
        state.depth = 4;
        assert!(!state.should_continue());
    }
    #[test]
    fn test_reflect3_down_up_normal() {
        let d = [0.0, -1.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let r = reflect3(d, n);
        assert!((r[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_point_light_attenuation() {
        let light = PointLight::new([0.0; 3], [1.0; 3], 1.0);
        let att = light.attenuate(0.0);
        assert!((att - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pbr_material_creation() {
        let mat = Material::pbr([0.5; 3], 0.0, 0.5, 1.0);
        assert_eq!(mat.mat_type, MaterialType::Pbr);
        assert!((mat.roughness - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_atrous_denoise_flat_image() {
        let n = 4 * 4;
        let color = vec![[0.5f64, 0.5, 0.5]; n];
        let normal = vec![[0.0f64, 1.0, 0.0]; n];
        let position = vec![[0.0f64; 3]; n];
        let result = atrous_denoise(&color, &normal, &position, 4, 4, 1, 0.1, 0.1, 0.1);
        assert_eq!(result.len(), n);
        for p in &result {
            assert!((p[0] - 0.5).abs() < 0.01);
        }
    }
}
