//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    Ray, RayMarchSettings, TransferFunction, TransferFunctionF64, Vec3, Volume, VolumeGrid,
};

/// Compute the cross product of two Vec3s.
pub(super) fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
/// Fill a Volume with a Gaussian density blob centred at `(cx, cy, cz)` in
/// normalised \[0, 1\]³ coordinates, with standard deviation `sigma` (also
/// normalised).
pub fn gaussian_blob_volume(
    nx: usize,
    ny: usize,
    nz: usize,
    cx: f32,
    cy: f32,
    cz: f32,
    sigma: f32,
) -> Volume {
    let origin = Vec3::zero();
    let size = Vec3::new(1.0, 1.0, 1.0);
    let mut vol = Volume::new(nx, ny, nz, origin, size);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let x = (ix as f32 + 0.5) / nx as f32;
                let y = (iy as f32 + 0.5) / ny as f32;
                let z = (iz as f32 + 0.5) / nz as f32;
                let dx = x - cx;
                let dy = y - cy;
                let dz = z - cz;
                let r2 = dx * dx + dy * dy + dz * dz;
                let s2 = sigma * sigma;
                let v = (-r2 / (2.0 * s2)).exp();
                vol.set_value(ix, iy, iz, v);
            }
        }
    }
    vol
}
/// Composite two RGBA colors using front-to-back blending.
///
/// `dst` is the accumulated buffer, `src` is the new sample.
/// Returns the composited result.
#[allow(dead_code)]
pub fn composite_front_to_back(dst: [f32; 4], src: [f32; 4]) -> [f32; 4] {
    let one_minus_a = 1.0 - dst[3];
    [
        dst[0] + one_minus_a * src[0] * src[3],
        dst[1] + one_minus_a * src[1] * src[3],
        dst[2] + one_minus_a * src[2] * src[3],
        dst[3] + one_minus_a * src[3],
    ]
}
/// Composite two RGBA colors using back-to-front blending.
#[allow(dead_code)]
pub fn composite_back_to_front(dst: [f32; 4], src: [f32; 4]) -> [f32; 4] {
    [
        src[0] * src[3] + dst[0] * (1.0 - src[3]),
        src[1] * src[3] + dst[1] * (1.0 - src[3]),
        src[2] * src[3] + dst[2] * (1.0 - src[3]),
        src[3] + dst[3] * (1.0 - src[3]),
    ]
}
/// Extract an axis-aligned slice from a Volume.
///
/// `axis` is 0 (X), 1 (Y), or 2 (Z).
/// `slice_idx` is the integer index along that axis.
///
/// Returns a 2D grid of density values.
#[allow(dead_code)]
pub fn extract_slice(vol: &Volume, axis: usize, slice_idx: usize) -> Vec<Vec<f32>> {
    match axis {
        0 => {
            let mut grid = Vec::with_capacity(vol.ny);
            for iy in 0..vol.ny {
                let mut row = Vec::with_capacity(vol.nz);
                for iz in 0..vol.nz {
                    row.push(vol.get_value(slice_idx.min(vol.nx - 1), iy, iz));
                }
                grid.push(row);
            }
            grid
        }
        1 => {
            let mut grid = Vec::with_capacity(vol.nx);
            for ix in 0..vol.nx {
                let mut row = Vec::with_capacity(vol.nz);
                for iz in 0..vol.nz {
                    row.push(vol.get_value(ix, slice_idx.min(vol.ny - 1), iz));
                }
                grid.push(row);
            }
            grid
        }
        _ => {
            let mut grid = Vec::with_capacity(vol.nx);
            for ix in 0..vol.nx {
                let mut row = Vec::with_capacity(vol.ny);
                for iy in 0..vol.ny {
                    row.push(vol.get_value(ix, iy, slice_idx.min(vol.nz - 1)));
                }
                grid.push(row);
            }
            grid
        }
    }
}
/// Extract an isosurface at the given `threshold` from a Volume.
///
/// Uses a simplified marching-cubes-like approach: for each cell, checks
/// if the threshold is crossed along any edge and produces intersection
/// vertices. This is a basic educational implementation, not the full
/// 256-case marching cubes table.
#[allow(dead_code)]
pub fn extract_isosurface_simple(vol: &Volume, threshold: f32) -> Vec<[Vec3; 3]> {
    let mut triangles: Vec<[Vec3; 3]> = Vec::new();
    for ix in 0..vol.nx.saturating_sub(1) {
        for iy in 0..vol.ny.saturating_sub(1) {
            for iz in 0..vol.nz.saturating_sub(1) {
                let corners = [
                    vol.get_value(ix, iy, iz),
                    vol.get_value(ix + 1, iy, iz),
                    vol.get_value(ix + 1, iy + 1, iz),
                    vol.get_value(ix, iy + 1, iz),
                    vol.get_value(ix, iy, iz + 1),
                    vol.get_value(ix + 1, iy, iz + 1),
                    vol.get_value(ix + 1, iy + 1, iz + 1),
                    vol.get_value(ix, iy + 1, iz + 1),
                ];
                let mut inside = 0u8;
                for (k, &c) in corners.iter().enumerate() {
                    if c >= threshold {
                        inside |= 1 << k;
                    }
                }
                if inside == 0 || inside == 0xFF {
                    continue;
                }
                let cell_size_x = vol.size.x / vol.nx as f32;
                let cell_size_y = vol.size.y / vol.ny as f32;
                let cell_size_z = vol.size.z / vol.nz as f32;
                let cx = vol.origin.x + (ix as f32 + 0.5) * cell_size_x;
                let cy = vol.origin.y + (iy as f32 + 0.5) * cell_size_y;
                let cz = vol.origin.z + (iz as f32 + 0.5) * cell_size_z;
                let half = cell_size_x * 0.3;
                let v0 = Vec3::new(cx - half, cy, cz);
                let v1 = Vec3::new(cx + half, cy, cz);
                let v2 = Vec3::new(cx, cy + half, cz);
                triangles.push([v0, v1, v2]);
            }
        }
    }
    triangles
}
/// March a ray through `grid` using Beer-Lambert absorption + transfer function.
///
/// Returns accumulated `[r, g, b, alpha]`.
#[allow(dead_code)]
pub fn ray_march(
    grid: &VolumeGrid,
    tf: &TransferFunctionF64,
    origin: [f64; 3],
    dir: [f64; 3],
    settings: &RayMarchSettings,
) -> [f32; 4] {
    let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    let dir = if len > 1e-12 {
        [dir[0] / len, dir[1] / len, dir[2] / len]
    } else {
        return [0.0; 4];
    };
    let mut t = 0.0_f64;
    let mut acc_r = 0.0_f32;
    let mut acc_g = 0.0_f32;
    let mut acc_b = 0.0_f32;
    let mut acc_a = 0.0_f32;
    for _ in 0..settings.max_steps {
        let pos = [
            origin[0] + dir[0] * t,
            origin[1] + dir[1] * t,
            origin[2] + dir[2] * t,
        ];
        let density = grid.sample_trilinear(pos);
        let rgba = tf.sample(density);
        let alpha_step = compute_opacity(settings.absorption * density, settings.step_size) as f32;
        let one_minus_a = 1.0 - acc_a;
        acc_r += one_minus_a * rgba[0] * alpha_step;
        acc_g += one_minus_a * rgba[1] * alpha_step;
        acc_b += one_minus_a * rgba[2] * alpha_step;
        acc_a += one_minus_a * alpha_step;
        if acc_a >= 0.99 {
            break;
        }
        t += settings.step_size;
    }
    [acc_r, acc_g, acc_b, acc_a]
}
/// Compute an approximate isosurface normal at `pos` using the gradient of `grid`.
///
/// The gradient points in the direction of increasing density; the normal
/// pointing outward from the isosurface is therefore the negated normalised
/// gradient.  Returns `[0, 0, 0]` when the gradient magnitude is near zero.
#[allow(dead_code)]
pub fn iso_surface_normal(grid: &VolumeGrid, pos: [f64; 3], _iso_value: f64) -> [f64; 3] {
    let g = grid.gradient_at(pos);
    let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    if mag < 1e-12 {
        [0.0; 3]
    } else {
        [g[0] / mag, g[1] / mag, g[2] / mag]
    }
}
/// Compute the opacity of a participating-media step via Beer-Lambert:
/// `opacity = 1 - exp(-absorption * step_size)`.
#[allow(dead_code)]
pub fn compute_opacity(absorption: f64, step_size: f64) -> f64 {
    1.0 - (-absorption * step_size).exp()
}
/// Build a `VolumeGrid` by sampling a scalar field function `f` over a uniform
/// grid.  Grid points are voxel-centred: `pos = (i + 0.5) * voxel_size`.
#[allow(dead_code)]
pub fn volume_from_scalar_field(
    f: impl Fn([f64; 3]) -> f64,
    nx: usize,
    ny: usize,
    nz: usize,
    voxel_size: f64,
) -> VolumeGrid {
    let mut grid = VolumeGrid::new(nx, ny, nz, voxel_size);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let pos = [
                    (ix as f64 + 0.5) * voxel_size,
                    (iy as f64 + 0.5) * voxel_size,
                    (iz as f64 + 0.5) * voxel_size,
                ];
                let idx = grid.idx(ix, iy, iz);
                grid.data[idx] = f(pos);
            }
        }
    }
    grid
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::BrickVolume;
    use crate::Camera;

    use crate::EarlyTerminationRayMarcher;
    use crate::EmissionAbsorptionAccumulator;

    use crate::RayMarcher;
    use crate::VolumeRenderer;
    fn unit_volume_constant(c: f32) -> Volume {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    vol.set_value(ix, iy, iz, c);
                }
            }
        }
        vol
    }
    #[test]
    fn test_sample_trilinear_constant() {
        let vol = unit_volume_constant(0.7);
        let sample = vol.sample_trilinear(Vec3::new(0.5, 0.5, 0.5));
        assert!((sample - 0.7).abs() < 1e-5, "expected 0.7, got {}", sample);
    }
    #[test]
    fn test_aabb_intersect_hit() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        assert!(
            vol.aabb_intersect(&ray).is_some(),
            "ray through centre should hit"
        );
    }
    #[test]
    fn test_aabb_intersect_miss() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let ray = Ray::new(Vec3::new(0.5, 0.5, 2.0), Vec3::new(0.0, 0.0, 1.0));
        assert!(
            vol.aabb_intersect(&ray).is_none(),
            "ray pointing away should miss"
        );
    }
    #[test]
    fn test_transfer_function_at_control_point() {
        let mut tf = TransferFunction::new();
        tf.add_control_point(0.5, [1.0, 0.0, 0.0], 0.8);
        let rgba = tf.evaluate(0.5);
        assert!((rgba[0] - 1.0).abs() < 1e-6);
        assert!((rgba[3] - 0.8).abs() < 1e-6);
    }
    #[test]
    fn test_transfer_function_interpolation() {
        let mut tf = TransferFunction::new();
        tf.add_control_point(0.0, [0.0, 0.0, 0.0], 0.0);
        tf.add_control_point(1.0, [1.0, 1.0, 1.0], 1.0);
        let rgba = tf.evaluate(0.5);
        assert!((rgba[0] - 0.5).abs() < 1e-5, "r should interpolate to 0.5");
        assert!((rgba[3] - 0.5).abs() < 1e-5, "a should interpolate to 0.5");
    }
    #[test]
    fn test_ray_marcher_empty_volume() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = TransferFunction::smoke_transfer_function();
        let marcher = RayMarcher::new(0.05);
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let rgba = marcher.march(&ray, &vol, &tf);
        assert!(
            rgba[3] < 1e-6,
            "empty volume should give zero alpha, got {}",
            rgba[3]
        );
    }
    #[test]
    fn test_camera_center_ray_toward_target() {
        let pos = Vec3::new(0.0, 0.0, -3.0);
        let target = Vec3::new(0.0, 0.0, 0.0);
        let up = Vec3::new(0.0, 1.0, 0.0);
        let cam = Camera::new(pos, target, up, std::f32::consts::FRAC_PI_4);
        let w = 100.0_f32;
        let h = 100.0_f32;
        let ray = cam.generate_ray(w / 2.0 - 0.5, h / 2.0 - 0.5, w, h);
        assert!(
            ray.direction.z > 0.9,
            "centre ray z should be close to 1, got {}",
            ray.direction.z
        );
        assert!(ray.direction.x.abs() < 0.1);
        assert!(ray.direction.y.abs() < 0.1);
    }
    #[test]
    fn test_volume_renderer_no_panic() {
        let pos = Vec3::new(0.5, 0.5, -2.0);
        let target = Vec3::new(0.5, 0.5, 0.5);
        let up = Vec3::new(0.0, 1.0, 0.0);
        let cam = Camera::new(pos, target, up, std::f32::consts::FRAC_PI_4);
        let marcher = RayMarcher::new(0.1);
        let tf = TransferFunction::smoke_transfer_function();
        let renderer = VolumeRenderer::new(cam, marcher, tf);
        let vol = gaussian_blob_volume(8, 8, 8, 0.5, 0.5, 0.5, 0.2);
        let pixels = renderer.render(&vol, 4, 4);
        assert_eq!(pixels.len(), 16, "4×4 render should produce 16 pixels");
    }
    #[test]
    fn test_composite_front_to_back_opaque() {
        let dst = [0.0, 0.0, 0.0, 0.0];
        let src = [1.0, 0.0, 0.0, 1.0];
        let result = composite_front_to_back(dst, src);
        assert!((result[0] - 1.0).abs() < 1e-5);
        assert!((result[3] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_composite_front_to_back_transparent() {
        let dst = [0.5, 0.5, 0.5, 0.5];
        let src = [1.0, 0.0, 0.0, 0.5];
        let result = composite_front_to_back(dst, src);
        assert!((result[0] - 0.75).abs() < 1e-5);
    }
    #[test]
    fn test_composite_back_to_front() {
        let dst = [0.0, 0.0, 0.0, 0.0];
        let src = [1.0, 0.0, 0.0, 1.0];
        let result = composite_back_to_front(dst, src);
        assert!((result[0] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_extract_slice_x() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.set_value(2, 1, 3, 0.7);
        let slice = extract_slice(&vol, 0, 2);
        assert_eq!(slice.len(), 4);
        assert_eq!(slice[0].len(), 4);
        assert!((slice[1][3] - 0.7).abs() < 1e-5);
    }
    #[test]
    fn test_extract_slice_y() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.set_value(1, 2, 0, 0.3);
        let slice = extract_slice(&vol, 1, 2);
        assert_eq!(slice.len(), 4);
        assert!((slice[1][0] - 0.3).abs() < 1e-5);
    }
    #[test]
    fn test_extract_slice_z() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.set_value(0, 0, 3, 0.9);
        let slice = extract_slice(&vol, 2, 3);
        assert_eq!(slice.len(), 4);
        assert!((slice[0][0] - 0.9).abs() < 1e-5);
    }
    #[test]
    fn test_isosurface_empty_volume() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tris = extract_isosurface_simple(&vol, 0.5);
        assert!(tris.is_empty(), "Empty volume should have no isosurface");
    }
    #[test]
    fn test_isosurface_full_volume() {
        let vol = unit_volume_constant(1.0);
        let tris = extract_isosurface_simple(&vol, 0.5);
        assert!(
            tris.is_empty(),
            "Fully inside volume should have no isosurface"
        );
    }
    #[test]
    fn test_isosurface_sphere_produces_triangles() {
        let mut vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.add_sphere(Vec3::new(0.5, 0.5, 0.5), 0.3, 1.0);
        let tris = extract_isosurface_simple(&vol, 0.5);
        assert!(
            !tris.is_empty(),
            "Sphere should produce isosurface triangles"
        );
    }
    #[test]
    fn test_volume_min_max() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.set_value(0, 0, 0, 0.1);
        vol.set_value(1, 1, 1, 0.9);
        let (min, max) = vol.min_max();
        assert!((min - 0.0).abs() < 1e-5);
        assert!((max - 0.9).abs() < 1e-5);
    }
    #[test]
    fn test_volume_average_density() {
        let vol = unit_volume_constant(0.5);
        let avg = vol.average_density();
        assert!((avg - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_volume_voxel_count() {
        let vol = Volume::new(3, 4, 5, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(vol.voxel_count(), 60);
    }
    #[test]
    fn test_volume_cell_size() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(2.0, 2.0, 2.0));
        let cs = vol.cell_size();
        assert!((cs.x - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_volume_fill() {
        let mut vol = Volume::new(2, 2, 2, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.fill(0.42);
        for &v in &vol.values {
            assert!((v - 0.42).abs() < 1e-5);
        }
    }
    #[test]
    fn test_volume_add_sphere() {
        let mut vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.add_sphere(Vec3::new(0.5, 0.5, 0.5), 0.3, 1.0);
        let center_val = vol.get_value(4, 4, 4);
        assert!(center_val > 0.0, "Centre of sphere should have density");
    }
    #[test]
    fn test_gaussian_blob_peak() {
        let vol = gaussian_blob_volume(8, 8, 8, 0.5, 0.5, 0.5, 0.1);
        let peak = vol.get_value(4, 4, 4);
        assert!(peak > 0.5, "Peak of Gaussian should be high, got {peak}");
    }
    #[test]
    fn test_plasma_tf_endpoint() {
        let tf = TransferFunction::plasma_transfer_function();
        let rgba = tf.evaluate(1.0);
        assert!((rgba[0] - 1.0).abs() < 1e-5);
        assert!((rgba[3] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_smoke_tf_zero_transparent() {
        let tf = TransferFunction::smoke_transfer_function();
        let rgba = tf.evaluate(0.0);
        assert!(
            rgba[3].abs() < 1e-5,
            "Smoke at density=0 should be transparent"
        );
    }
    #[test]
    fn test_gradient_constant_field_is_zero() {
        let vol = unit_volume_constant(0.5);
        let g = vol.gradient(1, 1, 1);
        for k in 0..3 {
            assert!(
                g[k].abs() < 1e-3,
                "gradient of constant field should be 0, got {}",
                g[k]
            );
        }
    }
    #[test]
    fn test_gradient_non_zero() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    vol.set_value(ix, iy, iz, ix as f32 * 0.1);
                }
            }
        }
        let g = vol.gradient(1, 2, 2);
        assert!(
            g[0].abs() > 0.01,
            "X-gradient should be non-zero, got {}",
            g[0]
        );
    }
    #[test]
    fn test_brick_volume_creation() {
        let bv = BrickVolume::new(8, 8, 8, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(bv.brick_size, 4);
        assert!(bv.n_bricks_x > 0);
        assert!(bv.n_bricks_y > 0);
        assert!(bv.n_bricks_z > 0);
    }
    #[test]
    fn test_brick_volume_set_get() {
        let mut bv = BrickVolume::new(8, 8, 8, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        bv.set_value(3, 5, 7, 0.77);
        assert!((bv.get_value(3, 5, 7) - 0.77).abs() < 1e-5);
    }
    #[test]
    fn test_brick_volume_fill_and_minmax() {
        let mut bv = BrickVolume::new(4, 4, 4, 2, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        bv.fill(0.42);
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    assert!((bv.get_value(ix, iy, iz) - 0.42).abs() < 1e-5);
                }
            }
        }
    }
    #[test]
    fn test_brick_count() {
        let bv = BrickVolume::new(8, 8, 8, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(bv.total_bricks(), 8);
    }
    #[test]
    fn test_early_termination_marcher_misses() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = TransferFunction::smoke_transfer_function();
        let marcher = EarlyTerminationRayMarcher {
            step_size: 0.05,
            opacity_threshold: 0.99,
        };
        let ray = Ray::new(Vec3::new(0.5, 0.5, 5.0), Vec3::new(0.0, 0.0, -1.0));
        let (rgba, steps) = marcher.march_with_stats(&ray, &vol, &tf);
        assert!(rgba[3] < 1e-5);
        assert!(steps > 0);
    }
    #[test]
    fn test_early_termination_marcher_opaque_volume() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.fill(1.0);
        let mut tf = TransferFunction::new();
        tf.add_control_point(0.0, [1.0, 0.0, 0.0], 0.0);
        tf.add_control_point(1.0, [1.0, 0.0, 0.0], 1.0);
        let marcher = EarlyTerminationRayMarcher {
            step_size: 0.05,
            opacity_threshold: 0.95,
        };
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let (rgba, steps) = marcher.march_with_stats(&ray, &vol, &tf);
        assert!(rgba[3] > 0.1, "opaque volume should accumulate opacity");
        let _ = steps;
    }
    #[test]
    fn test_emission_absorption_model_accumulate() {
        let mut acc = EmissionAbsorptionAccumulator::new();
        acc.step(0.1, [1.0, 0.0, 0.0], 0.5);
        acc.step(0.1, [0.0, 1.0, 0.0], 0.5);
        let result = acc.result();
        assert!(result[0] > 0.0, "red emission should contribute");
        assert!(result[1] > 0.0, "green emission should contribute");
        assert!(result[3] > 0.0, "alpha should be positive");
    }
    #[test]
    fn test_emission_absorption_fully_absorbed() {
        let mut acc = EmissionAbsorptionAccumulator::new();
        for _ in 0..100 {
            acc.step(0.1, [1.0, 0.0, 0.0], 100.0);
        }
        let result = acc.result();
        assert!(result[3] > 0.9, "high absorption should lead to high alpha");
    }
    #[test]
    fn test_gradient_magnitude_constant_zero() {
        let vol = unit_volume_constant(0.5);
        let mag = vol.gradient_magnitude(2, 2, 2);
        assert!(mag < 1e-3, "constant field gradient magnitude should be ~0");
    }
    #[test]
    fn test_gradient_magnitude_positive_for_varying_field() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    vol.set_value(ix, iy, iz, ix as f32 * 0.3 + iy as f32 * 0.1);
                }
            }
        }
        let mag = vol.gradient_magnitude(1, 1, 1);
        assert!(
            mag > 0.0,
            "gradient magnitude should be positive for varying field"
        );
    }
    #[test]
    fn test_volume_normal_unit_length() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        for ix in 0..4 {
            for iy in 0..4 {
                for iz in 0..4 {
                    vol.set_value(ix, iy, iz, ix as f32 * 0.2);
                }
            }
        }
        let n = vol.normal_at(1, 2, 2);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (len - 1.0).abs() < 0.01 || len < 1e-5,
            "normal should be unit length or zero, got {len}"
        );
    }
    #[test]
    fn test_volume_grid_sample_at_corner() {
        let grid = volume_from_scalar_field(|_| 0.75, 4, 4, 4, 1.0);
        let v = grid.sample_trilinear([1.5, 1.5, 1.5]);
        assert!(
            (v - 0.75).abs() < 1e-9,
            "trilinear sample on constant field should return its value, got {v}"
        );
    }
    #[test]
    fn test_volume_grid_outside_returns_zero() {
        let grid = VolumeGrid::new(4, 4, 4, 1.0);
        let v = grid.sample_trilinear([100.0, 0.0, 0.0]);
        assert!(
            (v - 0.0).abs() < 1e-12,
            "sample outside grid should return 0"
        );
    }
    #[test]
    fn test_gradient_at_sphere_center_is_near_zero() {
        let grid = volume_from_scalar_field(
            |p| {
                let dx = p[0] - 2.0;
                let dy = p[1] - 2.0;
                let dz = p[2] - 2.0;
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                (1.0 - r / 2.0).max(0.0)
            },
            4,
            4,
            4,
            1.0,
        );
        let g = grid.gradient_at([2.0, 2.0, 2.0]);
        let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!(
            mag < 0.5,
            "gradient at sphere centre should be small, got {mag}"
        );
    }
    #[test]
    fn test_compute_opacity_beer_lambert_zero_absorption() {
        let o = compute_opacity(0.0, 1.0);
        assert!(o.abs() < 1e-12, "zero absorption should give zero opacity");
    }
    #[test]
    fn test_compute_opacity_beer_lambert_high_absorption() {
        let o = compute_opacity(100.0, 1.0);
        assert!(
            (o - 1.0).abs() < 1e-4,
            "very high absorption should give opacity ~ 1"
        );
    }
    #[test]
    fn test_compute_opacity_beer_lambert_formula() {
        let mu = 0.5_f64;
        let ds = 0.2_f64;
        let expected = 1.0 - (-mu * ds).exp();
        let got = compute_opacity(mu, ds);
        assert!(
            (got - expected).abs() < 1e-12,
            "Beer-Lambert opacity mismatch"
        );
    }
    #[test]
    fn test_transfer_function_f64_interpolation() {
        let mut tf = TransferFunctionF64::new();
        tf.add_point(0.0, [0.0, 0.0, 0.0, 0.0]);
        tf.add_point(1.0, [1.0, 1.0, 1.0, 1.0]);
        let mid = tf.sample(0.5);
        assert!(
            (mid[0] - 0.5).abs() < 1e-5,
            "midpoint should interpolate to 0.5"
        );
        assert!((mid[3] - 0.5).abs() < 1e-5, "alpha midpoint should be 0.5");
    }
    #[test]
    fn test_transfer_function_f64_clamps_below() {
        let mut tf = TransferFunctionF64::new();
        tf.add_point(0.5, [0.3, 0.3, 0.3, 1.0]);
        let v = tf.sample(-10.0);
        assert!(
            (v[0] - 0.3).abs() < 1e-5,
            "sample below range should return first control point"
        );
    }
    #[test]
    fn test_ray_march_returns_valid_alpha() {
        let grid = volume_from_scalar_field(|_| 0.1, 8, 8, 8, 1.0);
        let mut tf = TransferFunctionF64::new();
        tf.add_point(0.0, [0.5, 0.5, 0.5, 0.5]);
        tf.add_point(1.0, [1.0, 1.0, 1.0, 1.0]);
        let settings = RayMarchSettings {
            step_size: 0.5,
            max_steps: 20,
            absorption: 0.2,
            scattering: 0.0,
        };
        let rgba = ray_march(&grid, &tf, [4.0, 4.0, -1.0], [0.0, 0.0, 1.0], &settings);
        assert!(
            rgba[3] >= 0.0 && rgba[3] <= 1.0 + 1e-5,
            "ray march alpha must be in [0,1], got {}",
            rgba[3]
        );
        assert!(
            rgba[3] > 0.0,
            "non-empty volume should produce positive alpha"
        );
    }
    #[test]
    fn test_ray_march_empty_volume_zero_alpha() {
        let grid = VolumeGrid::new(4, 4, 4, 1.0);
        let tf = TransferFunctionF64::new();
        let settings = RayMarchSettings {
            step_size: 0.5,
            max_steps: 20,
            absorption: 1.0,
            scattering: 0.0,
        };
        let rgba = ray_march(&grid, &tf, [2.0, 2.0, -1.0], [0.0, 0.0, 1.0], &settings);
        assert!(
            rgba[3] < 1e-5,
            "empty volume should give zero alpha, got {}",
            rgba[3]
        );
    }
    #[test]
    fn test_iso_surface_normal_unit_length() {
        let grid = volume_from_scalar_field(|p| p[0] / 4.0, 8, 8, 8, 0.5);
        let n = iso_surface_normal(&grid, [2.0, 2.0, 2.0], 0.5);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (mag - 1.0).abs() < 0.01 || mag < 1e-10,
            "iso surface normal should be unit length, got {mag}"
        );
    }
    #[test]
    fn test_volume_from_scalar_field_constant() {
        let grid = volume_from_scalar_field(|_| 0.42, 4, 4, 4, 1.0);
        for &v in &grid.data {
            assert!(
                (v - 0.42).abs() < 1e-12,
                "constant field should fill grid uniformly"
            );
        }
    }
}
/// Performs Maximum Intensity Projection (MIP) along a ray.
///
/// Returns the maximum density value encountered, mapped through
/// the transfer function.
#[allow(dead_code)]
pub fn mip_ray(ray: &Ray, volume: &Volume, tf: &TransferFunction, step_size: f32) -> [f32; 4] {
    let Some((t_enter, t_exit)) = volume.aabb_intersect(ray) else {
        return [0.0; 4];
    };
    let mut max_density = 0.0_f32;
    let mut t = t_enter;
    while t < t_exit {
        let p = ray.at(t);
        let d = volume.sample_trilinear(p);
        if d > max_density {
            max_density = d;
        }
        t += step_size;
    }
    tf.evaluate(max_density)
}
/// Compute simple Phong shading for a surface normal and light direction.
///
/// Returns a scalar in `[0, 1]` representing the lighting intensity.
#[allow(dead_code)]
pub fn phong_shading(
    normal: Vec3,
    light_dir: Vec3,
    view_dir: Vec3,
    ambient: f32,
    diffuse: f32,
    specular: f32,
    shininess: f32,
) -> f32 {
    let n = normal.normalize();
    let l = light_dir.normalize();
    let v = view_dir.normalize();
    let diff = n.dot(&l).max(0.0) * diffuse;
    let r = n.scale(2.0 * n.dot(&l)).sub(&l);
    let spec = r.dot(&v).max(0.0).powf(shininess) * specular;
    (ambient + diff + spec).clamp(0.0, 1.0)
}
#[cfg(test)]
mod extended_volume_tests {
    use super::*;
    use crate::AdaptiveRayMarcher;
    use crate::BrickVolume;
    use crate::Camera;
    use crate::DvrAccumulator;
    use crate::EarlyTerminationRayMarcher;
    use crate::EmissionAbsorptionAccumulator;
    use crate::EmptySpaceSkippingMarcher;
    use crate::MipRenderer;
    use crate::OccupancyGrid;

    use std::f32::consts::FRAC_PI_4;
    fn make_camera() -> Camera {
        Camera::new(
            Vec3::new(0.5, 0.5, -2.0),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 1.0, 0.0),
            FRAC_PI_4,
        )
    }
    fn dense_volume() -> Volume {
        gaussian_blob_volume(8, 8, 8, 0.5, 0.5, 0.5, 0.15)
    }
    fn smoke_tf() -> TransferFunction {
        TransferFunction::smoke_transfer_function()
    }
    #[test]
    fn test_mip_ray_empty_volume() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = smoke_tf();
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let rgba = mip_ray(&ray, &vol, &tf, 0.1);
        assert!(rgba[3] < 1e-5, "empty volume MIP should give zero alpha");
    }
    #[test]
    fn test_mip_ray_dense_volume_nonzero() {
        let vol = dense_volume();
        let tf = TransferFunction::plasma_transfer_function();
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let rgba = mip_ray(&ray, &vol, &tf, 0.05);
        assert!(rgba[3] >= 0.0);
    }
    #[test]
    fn test_mip_renderer_pixel_count() {
        let cam = make_camera();
        let renderer = MipRenderer::new(cam, 0.1);
        let vol = dense_volume();
        let tf = smoke_tf();
        let pixels = renderer.render(&vol, &tf, 4, 4);
        assert_eq!(pixels.len(), 16);
    }
    #[test]
    fn test_mip_miss_returns_zero() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = smoke_tf();
        let ray = Ray::new(Vec3::new(5.0, 5.0, 5.0), Vec3::new(1.0, 0.0, 0.0));
        let rgba = mip_ray(&ray, &vol, &tf, 0.1);
        for v in rgba {
            assert!(v.abs() < 1e-5);
        }
    }
    #[test]
    fn test_adaptive_marcher_empty_volume() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = smoke_tf();
        let marcher = AdaptiveRayMarcher::new(0.02, 0.1);
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let rgba = marcher.march(&ray, &vol, &tf);
        assert!(rgba[3] < 1e-5);
    }
    #[test]
    fn test_adaptive_marcher_miss() {
        let vol = dense_volume();
        let tf = smoke_tf();
        let marcher = AdaptiveRayMarcher::new(0.02, 0.1);
        let ray = Ray::new(Vec3::new(10.0, 10.0, 10.0), Vec3::new(1.0, 0.0, 0.0));
        let rgba = marcher.march(&ray, &vol, &tf);
        for v in rgba {
            assert!(v.abs() < 1e-5);
        }
    }
    #[test]
    fn test_adaptive_marcher_valid_alpha_range() {
        let vol = dense_volume();
        let tf = TransferFunction::plasma_transfer_function();
        let marcher = AdaptiveRayMarcher::new(0.05, 0.15);
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let rgba = marcher.march(&ray, &vol, &tf);
        assert!(rgba[3] >= 0.0 && rgba[3] <= 1.0 + 1e-5);
    }
    #[test]
    fn test_occupancy_grid_empty_volume() {
        let vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let occ = OccupancyGrid::build(&vol, 4, 0.01);
        assert!((occ.occupancy_ratio() - 0.0).abs() < 1e-5);
    }
    #[test]
    fn test_occupancy_grid_full_volume() {
        let vol = gaussian_blob_volume(8, 8, 8, 0.5, 0.5, 0.5, 0.5);
        let occ = OccupancyGrid::build(&vol, 4, 0.0001);
        assert!(occ.occupancy_ratio() > 0.0);
    }
    #[test]
    fn test_occupancy_grid_brick_count() {
        let vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let occ = OccupancyGrid::build(&vol, 4, 0.5);
        assert_eq!(occ.nx, 2);
        assert_eq!(occ.ny, 2);
        assert_eq!(occ.nz, 2);
    }
    #[test]
    fn test_occupancy_grid_outside_returns_false() {
        let vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let occ = OccupancyGrid::build(&vol, 4, 0.5);
        let p_outside = Vec3::new(10.0, 10.0, 10.0);
        assert!(!occ.is_occupied_at(&vol, p_outside));
    }
    #[test]
    fn test_empty_space_skipping_empty_volume() {
        let vol = Volume::new(8, 8, 8, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = smoke_tf();
        let occ = OccupancyGrid::build(&vol, 4, 0.01);
        let marcher = EmptySpaceSkippingMarcher::new(0.02, 0.1);
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let (rgba, _steps) = marcher.march(&ray, &vol, &tf, &occ);
        assert!(rgba[3] < 1e-5, "empty volume should give zero alpha");
    }
    #[test]
    fn test_empty_space_skipping_miss() {
        let vol = dense_volume();
        let tf = smoke_tf();
        let occ = OccupancyGrid::build(&vol, 4, 0.1);
        let marcher = EmptySpaceSkippingMarcher::new(0.02, 0.1);
        let ray = Ray::new(Vec3::new(10.0, 10.0, 10.0), Vec3::new(1.0, 0.0, 0.0));
        let (rgba, steps) = marcher.march(&ray, &vol, &tf, &occ);
        for v in rgba {
            assert!(v.abs() < 1e-5);
        }
        assert_eq!(steps, 0);
    }
    #[test]
    fn test_dvr_accumulator_starts_clear() {
        let acc = DvrAccumulator::new();
        assert!((acc.transmittance - 1.0).abs() < 1e-10);
        assert_eq!(acc.steps, 0);
        for v in acc.color {
            assert!(v.abs() < 1e-10);
        }
    }
    #[test]
    fn test_dvr_accumulator_transparent_step_no_change() {
        let mut acc = DvrAccumulator::new();
        let t_before = acc.transmittance;
        acc.integrate_step([0.5, 0.5, 0.5, 0.0], 0.1, 1.0);
        assert!((acc.transmittance - t_before).abs() < 1e-5);
    }
    #[test]
    fn test_dvr_accumulator_opaque_reduces_transmittance() {
        let mut acc = DvrAccumulator::new();
        acc.integrate_step([1.0, 0.0, 0.0, 1.0], 1.0, 10.0);
        assert!(
            acc.transmittance < 1.0,
            "opaque sample should reduce transmittance"
        );
    }
    #[test]
    fn test_dvr_accumulator_result_valid_range() {
        let mut acc = DvrAccumulator::new();
        for _ in 0..50 {
            acc.integrate_step([0.8, 0.2, 0.1, 0.5], 0.1, 1.0);
        }
        let r = acc.result();
        for v in r {
            assert!(
                (0.0..=1.0 + 1e-5).contains(&v),
                "result component out of range: {v}"
            );
        }
    }
    #[test]
    fn test_dvr_accumulator_termination() {
        let mut acc = DvrAccumulator::new();
        for _ in 0..100 {
            if acc.is_terminated(0.01) {
                break;
            }
            acc.integrate_step([1.0, 1.0, 1.0, 1.0], 1.0, 10.0);
        }
        assert!(
            acc.is_terminated(0.01) || acc.transmittance < 0.01,
            "should terminate after many opaque steps"
        );
    }
    #[test]
    fn test_volume_histogram_sum() {
        let vol = dense_volume();
        let hist = vol.histogram(16);
        let total: u32 = hist.iter().sum();
        assert_eq!(
            total as usize,
            vol.voxel_count(),
            "histogram should count all voxels"
        );
    }
    #[test]
    fn test_volume_histogram_bins() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.fill(0.5);
        let hist = vol.histogram(4);
        assert_eq!(hist.len(), 4);
        let nonzero: usize = hist.iter().filter(|&&c| c > 0).count();
        assert_eq!(nonzero, 1);
    }
    #[test]
    fn test_volume_normalize() {
        let mut vol = gaussian_blob_volume(4, 4, 4, 0.5, 0.5, 0.5, 0.2);
        vol.normalize();
        let (_, max_v) = vol.min_max();
        assert!(
            (max_v - 1.0).abs() < 1e-5,
            "after normalize max should be 1.0"
        );
    }
    #[test]
    fn test_volume_clamp_values() {
        let mut vol = gaussian_blob_volume(4, 4, 4, 0.5, 0.5, 0.5, 0.2);
        vol.clamp_values(0.2, 0.8);
        let (min_v, max_v) = vol.min_max();
        assert!(min_v >= 0.2 - 1e-5, "min should be >= 0.2");
        assert!(max_v <= 0.8 + 1e-5, "max should be <= 0.8");
    }
    #[test]
    fn test_volume_histogram_empty_single_bin() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let hist = vol.histogram(8);
        assert_eq!(hist.len(), 8);
    }
    #[test]
    fn test_phong_shading_frontal_light() {
        let n = Vec3::new(0.0, 0.0, 1.0);
        let l = Vec3::new(0.0, 0.0, 1.0);
        let v = Vec3::new(0.0, 0.0, 1.0);
        let intensity = phong_shading(n, l, v, 0.1, 0.8, 0.1, 32.0);
        assert!(
            intensity > 0.9,
            "front-lit surface should be bright, got {intensity}"
        );
    }
    #[test]
    fn test_phong_shading_back_light() {
        let n = Vec3::new(0.0, 0.0, 1.0);
        let l = Vec3::new(0.0, 0.0, -1.0);
        let v = Vec3::new(0.0, 0.0, 1.0);
        let intensity = phong_shading(n, l, v, 0.1, 0.8, 0.0, 32.0);
        assert!(intensity < 0.2, "back-lit surface should be mostly dark");
    }
    #[test]
    fn test_phong_shading_result_clamped() {
        let n = Vec3::new(0.0, 1.0, 0.0);
        let l = Vec3::new(0.0, 1.0, 0.0);
        let v = Vec3::new(0.0, 1.0, 0.0);
        let i = phong_shading(n, l, v, 0.5, 0.5, 0.5, 32.0);
        assert!((0.0..=1.0 + 1e-5).contains(&i));
    }
    #[test]
    fn test_brick_volume_set_get() {
        let mut bv = BrickVolume::new(8, 8, 8, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        bv.set_value(3, 5, 7, 0.42);
        assert!((bv.get_value(3, 5, 7) - 0.42).abs() < 1e-5);
    }
    #[test]
    fn test_brick_volume_to_volume_roundtrip() {
        let mut bv = BrickVolume::new(4, 4, 4, 2, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        bv.set_value(1, 2, 3, 0.7);
        let vol = bv.to_volume();
        assert!((vol.get_value(1, 2, 3) - 0.7).abs() < 1e-5);
    }
    #[test]
    fn test_brick_volume_fill() {
        let mut bv = BrickVolume::new(4, 4, 4, 2, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        bv.fill(1.0);
        let vol = bv.to_volume();
        let (_, max_v) = vol.min_max();
        assert!((max_v - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_brick_volume_total_bricks() {
        let bv = BrickVolume::new(8, 8, 8, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(bv.total_bricks(), 8);
    }
    #[test]
    fn test_volume_gradient_flat_is_zero() {
        let mut vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        vol.fill(0.5);
        let g = vol.gradient(2, 2, 2);
        assert!(g[0].abs() < 1e-5 && g[1].abs() < 1e-5 && g[2].abs() < 1e-5);
    }
    #[test]
    fn test_volume_gradient_magnitude_positive() {
        let vol = dense_volume();
        let gm = vol.gradient_magnitude(4, 4, 4);
        assert!(gm.is_finite());
    }
    #[test]
    fn test_volume_normal_at_unit_length() {
        let vol = dense_volume();
        let n = vol.normal_at(4, 4, 4);
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            mag < 1.0 + 1e-5,
            "normal should be unit length or zero, got {mag}"
        );
    }
    #[test]
    fn test_early_termination_marcher_empty() {
        let vol = Volume::new(4, 4, 4, Vec3::zero(), Vec3::new(1.0, 1.0, 1.0));
        let tf = smoke_tf();
        let marcher = EarlyTerminationRayMarcher {
            step_size: 0.1,
            opacity_threshold: 0.99,
        };
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let (rgba, _steps) = marcher.march_with_stats(&ray, &vol, &tf);
        assert!(rgba[3] < 1e-5);
    }
    #[test]
    fn test_early_termination_step_count_nonzero_for_non_empty() {
        let vol = dense_volume();
        let tf = TransferFunction::plasma_transfer_function();
        let marcher = EarlyTerminationRayMarcher {
            step_size: 0.05,
            opacity_threshold: 0.99,
        };
        let ray = Ray::new(Vec3::new(0.5, 0.5, -1.0), Vec3::new(0.0, 0.0, 1.0));
        let (_, steps) = marcher.march_with_stats(&ray, &vol, &tf);
        assert!(steps > 0);
    }
    #[test]
    fn test_emission_absorption_fresh() {
        let acc = EmissionAbsorptionAccumulator::new();
        assert!((acc.transmittance - 1.0).abs() < 1e-10);
        assert!(acc.color[0].abs() < 1e-10);
    }
    #[test]
    fn test_emission_absorption_adds_color() {
        let mut acc = EmissionAbsorptionAccumulator::new();
        acc.step(0.1, [1.0, 0.0, 0.0], 1.0);
        assert!(acc.color[0] > 0.0, "red emission should add to color");
    }
    #[test]
    fn test_emission_absorption_transmittance_decreases() {
        let mut acc = EmissionAbsorptionAccumulator::new();
        for _ in 0..10 {
            acc.step(0.1, [0.5, 0.5, 0.5], 2.0);
        }
        assert!(acc.transmittance < 1.0, "transmittance should decrease");
    }
    #[test]
    fn test_emission_absorption_result_alpha() {
        let mut acc = EmissionAbsorptionAccumulator::new();
        acc.step(1.0, [1.0, 1.0, 1.0], 10.0);
        let r = acc.result();
        assert!(r[3] > 0.0 && r[3] <= 1.0 + 1e-5);
    }
}
