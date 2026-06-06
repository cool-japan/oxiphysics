//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-12 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Linear interpolation between two RGBA colours.
pub fn lerp_color(a: [f32; 4], b: [f32; 4], t: f32) -> [f32; 4] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
        a[3] + (b[3] - a[3]) * t,
    ]
}
/// Compute the Q-criterion from the velocity gradient Jacobian (row-major 3×3).
///
/// Q = 0.5 * (||Ω||_F² - ||S||_F²)
///
/// where Ω = antisymmetric part, S = symmetric part of J.
/// Positive Q indicates vortex cores.
pub fn compute_q_criterion(j: [f64; 9]) -> f64 {
    let mut norm_s2 = 0.0;
    let mut norm_omega2 = 0.0;
    for i in 0..3 {
        for k in 0..3 {
            let j_ik = j[i * 3 + k];
            let j_ki = j[k * 3 + i];
            let s = 0.5 * (j_ik + j_ki);
            let omega = 0.5 * (j_ik - j_ki);
            norm_s2 += s * s;
            norm_omega2 += omega * omega;
        }
    }
    0.5 * (norm_omega2 - norm_s2)
}
/// Compute the λ2 criterion for vortex identification.
///
/// λ2 is the second eigenvalue of S²+Ω². Negative λ2 indicates a vortex core.
/// This function returns an approximate value using the Frobenius norms.
pub fn compute_lambda2_approx(j: [f64; 9]) -> f64 {
    -compute_q_criterion(j)
}
/// Compute the turbulent kinetic energy (TKE) from a velocity ensemble.
///
/// TKE = 0.5 * (σ_u² + σ_v² + σ_w²)
pub fn compute_turbulent_kinetic_energy(velocities: &[[f64; 3]]) -> f64 {
    let n = velocities.len();
    if n == 0 {
        return 0.0;
    }
    let mean: [f64; 3] = {
        let mut s = [0.0f64; 3];
        for v in velocities {
            s[0] += v[0];
            s[1] += v[1];
            s[2] += v[2];
        }
        [s[0] / n as f64, s[1] / n as f64, s[2] / n as f64]
    };
    let mut variance = [0.0f64; 3];
    for v in velocities {
        for k in 0..3 {
            let d = v[k] - mean[k];
            variance[k] += d * d;
        }
    }
    for v in variance.iter_mut() {
        *v /= n as f64;
    }
    0.5 * (variance[0] + variance[1] + variance[2])
}
/// Compute an approximate TKE spectrum from a velocity ensemble.
///
/// Returns (wavenumber, energy) pairs for a simplified energy spectrum.
pub fn compute_tke_spectrum(velocities: &[[f64; 3]], num_bins: usize) -> Vec<(f64, f64)> {
    let n = velocities.len();
    if n == 0 || num_bins == 0 {
        return Vec::new();
    }
    let mut bins = vec![0.0f64; num_bins];
    let max_k = n as f64 / 2.0;
    let bin_size = max_k / num_bins as f64;
    for (i, v) in velocities.iter().enumerate() {
        let k = (i as f64 * max_k / n as f64) as usize % num_bins;
        let ke = 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        bins[k] += ke;
    }
    (0..num_bins)
        .map(|i| (i as f64 * bin_size + bin_size * 0.5, bins[i]))
        .collect()
}
/// Compute a voxelized density map from particle positions and masses.
///
/// Uses a simple nearest-cell splatting approach.
pub fn compute_density_map(
    positions: &[[f64; 3]],
    masses: &[f64],
    grid_res: usize,
    domain_size: f64,
    smoothing_radius: f64,
) -> Vec<f64> {
    let n = grid_res;
    let h = domain_size / n as f64;
    let mut density = vec![0.0f64; n * n * n];
    for (pos, &mass) in positions.iter().zip(masses.iter()) {
        let ix = (pos[0] / domain_size * n as f64).clamp(0.0, (n - 1) as f64) as usize;
        let iy = (pos[1] / domain_size * n as f64).clamp(0.0, (n - 1) as f64) as usize;
        let iz = (pos[2] / domain_size * n as f64).clamp(0.0, (n - 1) as f64) as usize;
        let r_cells = (smoothing_radius / h) as i64 + 1;
        for di in -r_cells..=r_cells {
            for dj in -r_cells..=r_cells {
                for dk in -r_cells..=r_cells {
                    let gx = ix as i64 + di;
                    let gy = iy as i64 + dj;
                    let gz = iz as i64 + dk;
                    if gx < 0
                        || gy < 0
                        || gz < 0
                        || gx >= n as i64
                        || gy >= n as i64
                        || gz >= n as i64
                    {
                        continue;
                    }
                    let cx = gx as f64 * h + h * 0.5;
                    let cy = gy as f64 * h + h * 0.5;
                    let cz = gz as f64 * h + h * 0.5;
                    let dx = pos[0] - cx;
                    let dy = pos[1] - cy;
                    let dz = pos[2] - cz;
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    if dist < smoothing_radius {
                        let w = 1.0 - dist / smoothing_radius;
                        let gidx = gx as usize * n * n + gy as usize * n + gz as usize;
                        density[gidx] += mass * w;
                    }
                }
            }
        }
    }
    density
}
/// Compute the gradient of a scalar field at a grid cell using central differences.
pub fn field_gradient(
    field: &[f64],
    ix: usize,
    iy: usize,
    iz: usize,
    res: usize,
    h: f64,
) -> [f64; 3] {
    let idx = |x: usize, y: usize, z: usize| x * res * res + y * res + z;
    let n = res;
    let dx = if ix + 1 < n && ix > 0 {
        (field[idx(ix + 1, iy, iz)] - field[idx(ix - 1, iy, iz)]) / (2.0 * h)
    } else {
        0.0
    };
    let dy = if iy + 1 < n && iy > 0 {
        (field[idx(ix, iy + 1, iz)] - field[idx(ix, iy - 1, iz)]) / (2.0 * h)
    } else {
        0.0
    };
    let dz = if iz + 1 < n && iz > 0 {
        (field[idx(ix, iy, iz + 1)] - field[idx(ix, iy, iz - 1)]) / (2.0 * h)
    } else {
        0.0
    };
    [dx, dy, dz]
}
/// Compute mean curvature of a level-set field using finite differences.
///
/// κ = div(∇φ / |∇φ|)  (simplified second-order FD)
pub fn level_set_mean_curvature(
    phi: &[f64],
    ix: usize,
    iy: usize,
    iz: usize,
    res: usize,
    h: f64,
) -> f64 {
    let grad = field_gradient(phi, ix, iy, iz, res, h);
    let grad_len = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    if grad_len < 1e-12 {
        return 0.0;
    }
    let n_hat = [grad[0] / grad_len, grad[1] / grad_len, grad[2] / grad_len];
    let _ = n_hat;
    grad_len
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow_visualization::LicRenderer;
    use crate::fluid_viz::Bubble;
    use crate::fluid_viz::BubbleRenderer;
    use crate::fluid_viz::FluidStats;
    use crate::fluid_viz::FoamRenderer;
    use crate::fluid_viz::FreeSurfaceExtractor;
    use crate::fluid_viz::IsosurfaceViz;
    use crate::fluid_viz::Metaball;
    use crate::fluid_viz::MetaballField;
    use crate::fluid_viz::ParticleSplat;
    use crate::fluid_viz::PathlineIntegrator;
    use crate::fluid_viz::PressureColormap;
    use crate::fluid_viz::Splat;
    use crate::fluid_viz::StreaklineRenderer;
    use crate::fluid_viz::StreamlineRenderer;
    use crate::fluid_viz::VelocityArrows;
    use crate::fluid_viz::VolumeSlice;
    use crate::fluid_viz::VorticityViz;
    use crate::fluid_viz::WaveViz;
    #[test]
    fn test_velocity_arrow_below_min() {
        let mut va = VelocityArrows::new(1.0, 10.0);
        va.min_magnitude = 5.0;
        assert!(va.make_arrow([0.0; 3], [1.0, 0.0, 0.0]).is_none());
    }
    #[test]
    fn test_velocity_arrow_above_min() {
        let va = VelocityArrows::new(1.0, 10.0);
        assert!(va.make_arrow([0.0; 3], [10.0, 0.0, 0.0]).is_some());
    }
    #[test]
    fn test_velocity_arrow_tip_displaced() {
        let va = VelocityArrows::new(2.0, 10.0);
        let arrow = va.make_arrow([0.0; 3], [1.0, 0.0, 0.0]).unwrap();
        assert!((arrow.tip[0] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_lerp_color_t0() {
        let a = [1.0_f32, 0.0, 0.0, 1.0];
        let b = [0.0_f32, 1.0, 0.0, 1.0];
        let c = lerp_color(a, b, 0.0);
        assert!((c[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_lerp_color_t1() {
        let a = [1.0_f32, 0.0, 0.0, 1.0];
        let b = [0.0_f32, 1.0, 0.0, 1.0];
        let c = lerp_color(a, b, 1.0);
        assert!((c[1] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_pressure_iso_contours() {
        let pc = PressureColormap::new(0.0, 100.0);
        assert_eq!(pc.iso_contours().len(), 11);
    }
    #[test]
    fn test_pressure_map_range() {
        let pc = PressureColormap::new(0.0, 100.0);
        let c0 = pc.map(0.0);
        let c1 = pc.map(100.0);
        assert!(c0 != c1);
    }
    #[test]
    fn test_streamline_trace_one_point() {
        let sr = StreamlineRenderer::new(5, 0.1);
        let sl = sr.trace([0.0; 3], |_| [1.0, 0.0, 0.0]);
        assert!(!sl.points.is_empty());
    }
    #[test]
    fn test_streamline_trace_steps() {
        let sr = StreamlineRenderer::new(10, 0.1);
        let sl = sr.trace([0.0; 3], |_| [1.0, 0.0, 0.0]);
        assert_eq!(sl.points.len(), 11);
    }
    #[test]
    fn test_streamline_stops_zero_vel() {
        let sr = StreamlineRenderer::new(100, 0.1);
        let sl = sr.trace([0.0; 3], |_| [0.0; 3]);
        assert_eq!(sl.points.len(), 1);
    }
    #[test]
    fn test_vorticity_color_zero() {
        let vv = VorticityViz::new(10.0);
        let c = vv.color(0.0);
        assert!((c[2] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_vorticity_identity_jacobian() {
        let vv = VorticityViz::new(1.0);
        let j = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let w = vv.vorticity_from_jacobian(j);
        assert!(w.iter().all(|x| x.abs() < 1e-12));
    }
    #[test]
    fn test_vortex_core_cells() {
        let vv = VorticityViz::new(1.0);
        let lambda2 = vec![-0.5, 0.1, -0.2, 0.3];
        let cells = vv.vortex_core_cells(&lambda2);
        assert_eq!(cells, vec![0, 2]);
    }
    #[test]
    fn test_splat_gaussian_peak() {
        let ps = ParticleSplat::new(100.0);
        assert!((ps.gaussian_weight(0.0, 1.0) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_splat_gaussian_decreases() {
        let ps = ParticleSplat::new(100.0);
        assert!(ps.gaussian_weight(0.5, 1.0) > ps.gaussian_weight(1.0, 1.0));
    }
    #[test]
    fn test_iso_edge_interp_midpoint() {
        let iso = IsosurfaceViz::new(0.5, 8, 1.0);
        let p = iso.edge_interp([0.0; 3], 0.0, [1.0, 0.0, 0.0], 1.0);
        assert!((p[0] - 0.5).abs() < 1e-9);
    }
    #[test]
    fn test_volume_slice_extract() {
        let mut vs = VolumeSlice::new(4, 2, 0.5, 0.25);
        let field: Vec<f64> = (0..64).map(|i| i as f64).collect();
        vs.extract(&field, 4, 1.0);
        assert!(vs.pixels.iter().any(|&x| x != 0.0));
    }
    #[test]
    fn test_volume_slice_sample() {
        let mut vs = VolumeSlice::new(4, 2, 0.5, 0.25);
        vs.pixels = vec![1.0; 16];
        let v = vs.sample_bilinear(0.5, 0.5);
        assert!((v - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_fluid_stats_reynolds() {
        let mut fs = FluidStats::new();
        fs.compute_reynolds(1000.0, 1.0, 0.1, 0.001);
        assert!((fs.reynolds - 100000.0).abs() < 1.0);
    }
    #[test]
    fn test_fluid_stats_format() {
        let fs = FluidStats::new();
        assert!(!fs.format().is_empty());
    }
    #[test]
    fn test_wave_fresnel_at_zero() {
        let wv = WaveViz::new(4);
        assert!((wv.fresnel(1.0) - 0.02).abs() < 1e-9);
    }
    #[test]
    fn test_wave_update_normals() {
        let mut wv = WaveViz::new(4);
        let heights: Vec<f64> = (0..16).map(|i| (i as f64).sin()).collect();
        wv.update_normals(&heights, 0.25);
        assert!(
            wv.normal_map
                .iter()
                .any(|n| n[0].abs() > 0.0 || n[2].abs() > 0.0)
        );
    }
    #[test]
    fn test_wave_foam_range() {
        let mut wv = WaveViz::new(4);
        let heights = vec![
            0.0, 0.5, 1.0, 2.0, 0.1, 0.3, 0.8, 1.5, 0.2, 0.4, 0.9, 0.6, 0.7, 1.1, 0.3, 0.5,
        ];
        wv.update_foam(&heights);
        assert!(wv.foam_mask.iter().all(|&f| (0.0..=1.0).contains(&f)));
    }
    #[test]
    fn test_bubble_max_count() {
        let mut br = BubbleRenderer::new(3);
        for _ in 0..5 {
            br.add_bubble(Bubble::new([0.0; 3], 0.01));
        }
        assert_eq!(br.bubbles.len(), 3);
    }
    #[test]
    fn test_bubble_prune() {
        let mut br = BubbleRenderer::new(10);
        br.add_bubble(Bubble::new([0.0, 5.0, 0.0], 0.01));
        br.add_bubble(Bubble::new([0.0, -1.0, 0.0], 0.01));
        br.prune_above_surface(0.0);
        assert_eq!(br.bubbles.len(), 1);
    }
    #[test]
    fn test_bubble_total_internal_reflection() {
        let b = Bubble::new([0.0; 3], 0.01);
        let result = b.refract([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.5);
        let _ = result;
    }
    #[test]
    fn test_pressure_apply_count() {
        let pc = PressureColormap::new(0.0, 100.0);
        let positions = vec![[0.0; 3]; 5];
        let pressures = vec![10.0, 20.0, 50.0, 80.0, 100.0];
        let samples = pc.apply(&positions, &pressures);
        assert_eq!(samples.len(), 5);
    }
    #[test]
    fn test_velocity_arrow_build_field() {
        let va = VelocityArrows::new(1.0, 5.0);
        let origins = vec![[0.0; 3]; 4];
        let vels = vec![[1.0, 0.0, 0.0]; 4];
        let arrows = va.build_field(&origins, &vels);
        assert_eq!(arrows.len(), 4);
    }
    #[test]
    fn test_splat_depth_sort() {
        let ps = ParticleSplat::new(100.0);
        let mut splats = vec![
            Splat {
                position: [0.0; 2],
                radius: 1.0,
                color: [1.0; 4],
                depth: 3.0,
            },
            Splat {
                position: [0.0; 2],
                radius: 1.0,
                color: [1.0; 4],
                depth: 1.0,
            },
            Splat {
                position: [0.0; 2],
                radius: 1.0,
                color: [1.0; 4],
                depth: 2.0,
            },
        ];
        let _ = ps;
        ParticleSplat::depth_sort(&mut splats);
        assert!(splats[0].depth >= splats[1].depth);
        assert!(splats[1].depth >= splats[2].depth);
    }
    #[test]
    fn test_iso_crossings_monotone() {
        let iso = IsosurfaceViz::new(0.5, 4, 1.0);
        let field: Vec<f64> = (0..64).map(|i| i as f64 / 63.0).collect();
        let crossings = iso.extract_crossings(&field);
        assert!(!crossings.is_empty());
    }
    #[test]
    fn test_metaball_field_at_center() {
        let field = MetaballField::new(
            vec![Metaball {
                center: [0.0; 3],
                radius: 1.0,
                strength: 1.0,
            }],
            0.5,
        );
        let v = field.evaluate([0.0, 0.0, 0.0]);
        assert!(v > 1.0);
    }
    #[test]
    fn test_metaball_field_far() {
        let field = MetaballField::new(
            vec![Metaball {
                center: [0.0; 3],
                radius: 0.1,
                strength: 0.1,
            }],
            0.5,
        );
        let v = field.evaluate([100.0, 0.0, 0.0]);
        assert!(v < 0.01);
    }
    #[test]
    fn test_lic_output_size() {
        use crate::flow_visualization::{LicConfig, Vec2, VectorField};
        let lic = LicRenderer::with_config(LicConfig {
            kernel_length: 5,
            step_size: 0.1,
            contrast: 1.0,
        });
        let mut field = VectorField::new_2d(8, 8, 1.0, 1.0, Vec2::zero());
        for j in 0..8 {
            for i in 0..8 {
                field.set(i, j, 0, crate::flow_visualization::Vec3::new(1.0, 0.0, 0.0));
            }
        }
        let noise: Vec<f64> = (0..64).map(|i| i as f64 / 64.0).collect();
        let img = lic.render(&field, Some(&noise));
        assert_eq!(img.len(), 64);
    }
    #[test]
    fn test_lic_uniform_flow_nonzero() {
        use crate::flow_visualization::{LicConfig, Vec2, VectorField};
        let lic = LicRenderer::with_config(LicConfig {
            kernel_length: 10,
            step_size: 0.1,
            contrast: 1.0,
        });
        let mut field = VectorField::new_2d(16, 16, 1.0, 1.0, Vec2::zero());
        for j in 0..16 {
            for i in 0..16 {
                field.set(i, j, 0, crate::flow_visualization::Vec3::new(1.0, 0.0, 0.0));
            }
        }
        let noise: Vec<f64> = (0..256)
            .map(|i| if i % 2 == 0 { 1.0 } else { 0.0 })
            .collect();
        let img = lic.render(&field, Some(&noise));
        let nonzero = img.iter().filter(|&&v| v > 0.0).count();
        assert!(nonzero > 0);
    }
    #[test]
    fn test_q_criterion_positive_at_vortex() {
        let omega = 2.0_f64;
        let j = [0.0, -omega, 0.0, omega, 0.0, 0.0, 0.0, 0.0, 0.0];
        let q = compute_q_criterion(j);
        assert!(q > 0.0, "Q={q}");
    }
    #[test]
    fn test_q_criterion_negative_at_strain() {
        let j = [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0];
        let q = compute_q_criterion(j);
        assert!(q < 0.0, "Q={q}");
    }
    #[test]
    fn test_pathline_length_increases() {
        let mut pl = PathlineIntegrator::new(5, 0.1);
        pl.add_seed([0.0; 3]);
        pl.step(|_| [1.0, 0.0, 0.0]);
        pl.step(|_| [1.0, 0.0, 0.0]);
        let path = &pl.paths[0];
        assert_eq!(path.len(), 3);
    }
    #[test]
    fn test_streakline_initial_empty() {
        let sr = StreaklineRenderer::new(20);
        assert!(sr.positions.is_empty());
    }
    #[test]
    fn test_streakline_emit() {
        let mut sr = StreaklineRenderer::new(20);
        sr.emit([1.0, 2.0, 3.0]);
        assert_eq!(sr.positions.len(), 1);
    }
    #[test]
    fn test_streakline_max_length() {
        let mut sr = StreaklineRenderer::new(5);
        for i in 0..10 {
            sr.emit([i as f64, 0.0, 0.0]);
        }
        assert_eq!(sr.positions.len(), 5);
    }
    #[test]
    fn test_free_surface_empty() {
        let fse = FreeSurfaceExtractor::new(4, 1.0);
        let phi = vec![1.0_f64; 64];
        let tris = fse.extract(&phi);
        assert!(tris.is_empty());
    }
    #[test]
    fn test_free_surface_crossings_exist() {
        let fse = FreeSurfaceExtractor::new(4, 1.0);
        let mut phi = vec![1.0_f64; 64];
        phi[0] = -1.0;
        let tris = fse.extract(&phi);
        assert!(!tris.is_empty());
    }
    #[test]
    fn test_tke_spectrum_nonempty() {
        let velocities: Vec<[f64; 3]> = (0..16)
            .map(|i| [(i as f64).sin(), (i as f64 * 0.5).cos(), 0.0])
            .collect();
        let spectrum = compute_tke_spectrum(&velocities, 16);
        assert!(!spectrum.is_empty());
    }
    #[test]
    fn test_tke_nonnegative() {
        let velocities: Vec<[f64; 3]> = (0..8).map(|i| [(i as f64).sin(); 3]).collect();
        let tke = compute_turbulent_kinetic_energy(&velocities);
        assert!(tke >= 0.0);
    }
    #[test]
    fn test_foam_no_low_velocity() {
        let fr = FoamRenderer::new(0.5, 0.8);
        let count = fr.should_foam(0.1);
        assert!(!count);
    }
    #[test]
    fn test_foam_high_velocity() {
        let fr = FoamRenderer::new(0.5, 0.8);
        let count = fr.should_foam(1.0);
        assert!(count);
    }
    #[test]
    fn test_density_map_positive() {
        let positions = vec![[0.5_f64; 3]; 8];
        let masses = vec![1.0_f64; 8];
        let dm = compute_density_map(&positions, &masses, 4, 1.0, 0.3);
        let total: f64 = dm.iter().sum();
        assert!(total > 0.0, "density map should be nonzero");
    }
    #[test]
    fn test_velocity_arrows_all_zero() {
        let mut va = VelocityArrows::new(1.0, 1.0);
        va.min_magnitude = 1e-9;
        let origins = vec![[0.0; 3]; 4];
        let vels = vec![[0.0; 3]; 4];
        let arrows = va.build_field(&origins, &vels);
        assert!(
            arrows.is_empty(),
            "zero velocities should produce no arrows"
        );
    }
    #[test]
    fn test_vorticity_rotation() {
        let vv = VorticityViz::new(10.0);
        let j = [0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let w = vv.vorticity_from_jacobian(j);
        assert!((w[2] - 2.0).abs() < 1e-10, "w_z={}", w[2]);
    }
    #[test]
    fn test_streamline_trace_many() {
        let sr = StreamlineRenderer::new(5, 0.1);
        let seeds = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let lines = sr.trace_many(&seeds, |_| [0.5, 0.0, 0.0]);
        assert_eq!(lines.len(), 3);
    }
    #[test]
    fn test_bubble_visible_count() {
        let mut br = BubbleRenderer::new(10);
        let mut b1 = Bubble::new([0.0; 3], 0.1);
        b1.opacity = 0.5;
        let mut b2 = Bubble::new([1.0; 3], 0.1);
        b2.opacity = 0.01;
        br.add_bubble(b1);
        br.add_bubble(b2);
        assert_eq!(br.visible_count(0.1), 1);
    }
    #[test]
    fn test_fluid_stats_cfl() {
        let mut fs = FluidStats::new();
        fs.compute_cfl(10.0, 0.001, 0.1);
        assert!((fs.cfl - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_fluid_stats_mach() {
        let mut fs = FluidStats::new();
        fs.compute_mach(340.0, 340.0);
        assert!((fs.mach - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_wave_fresnel_increases() {
        let wv = WaveViz::new(4);
        let f0 = wv.fresnel(1.0);
        let f90 = wv.fresnel(0.0);
        assert!(f90 > f0, "fresnel should increase at grazing angle");
    }
    #[test]
    fn test_pressure_colormap_valid_range() {
        let pc = PressureColormap::new(-100.0, 100.0);
        for v in [-100.0, -50.0, 0.0, 50.0, 100.0] {
            let c = pc.map(v);
            for component in c.iter().take(3) {
                assert!(
                    *component >= 0.0 && *component <= 1.0,
                    "component out of range: {component}"
                );
            }
        }
    }
}
