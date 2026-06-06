//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::primitives::Color;
use oxiphysics_core::math::Vec3;

use super::types::{
    ArrowGlyph, EnhancedLicParams, FlowTopologyPoint, FtleField, JobardLeferParams, LicParams,
    LicTexture, SeedStrategy, Streamline, StreamlineIntegrator, TopologyPointKind, VortexCore,
};

/// Compute the velocity-gradient tensor at a grid point using central finite
/// differences.
///
/// `h` is the grid spacing used for the central difference stencil.
pub(super) fn velocity_gradient(pos: Vec3, h: f64, field: &dyn Fn(Vec3) -> Vec3) -> [[f64; 3]; 3] {
    let two_h = 2.0 * h;
    let mut grad = [[0.0_f64; 3]; 3];
    for j in 0..3 {
        let mut fwd = pos;
        let mut bwd = pos;
        fwd[j] += h;
        bwd[j] -= h;
        let vf = field(fwd);
        let vb = field(bwd);
        for i in 0..3 {
            grad[i][j] = (vf[i] - vb[i]) / two_h;
        }
    }
    grad
}
/// Decompose velocity-gradient tensor `J` into symmetric (strain S) and
/// antisymmetric (rotation W) parts.
pub(super) fn decompose_gradient(j: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut s = [[0.0_f64; 3]; 3];
    let mut w = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            s[i][k] = 0.5 * (j[i][k] + j[k][i]);
            w[i][k] = 0.5 * (j[i][k] - j[k][i]);
        }
    }
    (s, w)
}
/// Compute the Q-criterion: Q = ½(|W|² – |S|²) where |·| is the Frobenius norm.
///
/// Positive Q indicates rotation-dominated regions (vortex cores).
pub fn q_criterion(j: &[[f64; 3]; 3]) -> f64 {
    let (s, w) = decompose_gradient(j);
    let frobenius_sq =
        |m: &[[f64; 3]; 3]| -> f64 { m.iter().flat_map(|row| row.iter()).map(|&x| x * x).sum() };
    0.5 * (frobenius_sq(&w) - frobenius_sq(&s))
}
/// Compute the three eigenvalues of the symmetric 3×3 matrix `m = S² + W²`
/// (used in the λ₂ criterion), returning them sorted in ascending order.
///
/// Uses Cardano's method for the characteristic polynomial of a symmetric 3×3 matrix.
pub fn lambda2_eigenvalues(j: &[[f64; 3]; 3]) -> [f64; 3] {
    use std::f64::consts::PI;
    let (s, w) = decompose_gradient(j);
    let mut a = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            let mut val = 0.0;
            for m in 0..3 {
                val += s[i][m] * s[m][k] + w[i][m] * w[m][k];
            }
            a[i][k] = val;
        }
    }
    let p1 = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    let trace = a[0][0] + a[1][1] + a[2][2];
    let q_val = (a[0][0] * a[1][1] + a[0][0] * a[2][2] + a[1][1] * a[2][2]
        - a[0][1] * a[0][1]
        - a[0][2] * a[0][2]
        - a[1][2] * a[1][2])
        / 3.0;
    let r_val = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[1][2])
        - a[0][1] * (a[0][1] * a[2][2] - a[1][2] * a[0][2])
        + a[0][2] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))
        / 2.0;
    if p1.abs() < 1e-24 {
        let mut eigs = [a[0][0], a[1][1], a[2][2]];
        eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        return eigs;
    }
    let p2 = trace / 3.0;
    let p = (p1 / 3.0).sqrt();
    let discr = (r_val - q_val * p2 - 2.0 * p2.powi(3)) / (2.0 * p * p * p);
    let discr_clamped = discr.clamp(-1.0, 1.0);
    let phi = discr_clamped.acos() / 3.0;
    let e1 = p2 + 2.0 * p * phi.cos();
    let e2 = p2 + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
    let e3 = p2 + 2.0 * p * (phi + 4.0 * PI / 3.0).cos();
    let mut eigs = [e1, e2, e3];
    eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}
/// Detect vortex cores from a 3-D grid of sample points.
///
/// A point is classified as a vortex core candidate when both:
///  - Q-criterion > `q_threshold`, **and**
///  - The second eigenvalue of S²+W² (λ₂) < `lambda2_threshold` (negative).
///
/// `h` is the finite-difference spacing used to estimate gradients.
pub fn detect_vortex_cores(
    sample_points: &[Vec3],
    velocity_field: &dyn Fn(Vec3) -> Vec3,
    h: f64,
    q_threshold: f64,
    lambda2_threshold: f64,
) -> Vec<VortexCore> {
    sample_points
        .iter()
        .filter_map(|&pos| {
            let j = velocity_gradient(pos, h, velocity_field);
            let q = q_criterion(&j);
            let eigs = lambda2_eigenvalues(&j);
            let l2 = eigs[1];
            if q > q_threshold && l2 < lambda2_threshold {
                Some(VortexCore {
                    position: pos,
                    q_criterion: q,
                    lambda2: l2,
                })
            } else {
                None
            }
        })
        .collect()
}
/// Compute a Line Integral Convolution texture for a 2-D velocity slice.
///
/// The velocity field is sampled in the XY plane at `z = params.z_plane`.
/// A white-noise input texture is convolved along local streamlines to produce
/// the characteristic striated LIC pattern.
///
/// `velocity_field` maps `[x, y, z]` to `[vx, vy, vz]`; only the x and y
/// components are used for the 2-D advection.
pub fn compute_lic(
    params: &LicParams,
    velocity_field: &dyn Fn([f64; 3]) -> [f64; 3],
) -> LicTexture {
    let w = params.width;
    let h = params.height;
    let kl = params.kernel_length;
    let noise: Vec<f64> = {
        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        (0..(w * h))
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (seed >> 33) as f64 / u32::MAX as f64
            })
            .collect()
    };
    let dx_world = (params.x_range[1] - params.x_range[0]) / w as f64;
    let dy_world = (params.y_range[1] - params.y_range[0]) / h as f64;
    let pixel_to_world = |col: f64, row: f64| -> [f64; 2] {
        [
            params.x_range[0] + (col + 0.5) * dx_world,
            params.y_range[0] + (row + 0.5) * dy_world,
        ]
    };
    let world_to_pixel = |x: f64, y: f64| -> (f64, f64) {
        (
            (x - params.x_range[0]) / dx_world - 0.5,
            (y - params.y_range[0]) / dy_world - 0.5,
        )
    };
    let sample_noise = |px: f64, py: f64| -> f64 {
        let col = px.floor() as isize;
        let row = py.floor() as isize;
        let fc = px - col as f64;
        let fr = py - row as f64;
        let clamp_col = |c: isize| c.clamp(0, w as isize - 1) as usize;
        let clamp_row = |r: isize| r.clamp(0, h as isize - 1) as usize;
        let n00 = noise[clamp_row(row) * w + clamp_col(col)];
        let n10 = noise[clamp_row(row) * w + clamp_col(col + 1)];
        let n01 = noise[clamp_row(row + 1) * w + clamp_col(col)];
        let n11 = noise[clamp_row(row + 1) * w + clamp_col(col + 1)];
        n00 * (1.0 - fc) * (1.0 - fr)
            + n10 * fc * (1.0 - fr)
            + n01 * (1.0 - fc) * fr
            + n11 * fc * fr
    };
    let mut output = vec![0.0_f64; w * h];
    for row in 0..h {
        for col in 0..w {
            let [wx, wy] = pixel_to_world(col as f64, row as f64);
            let mut sum = 0.0_f64;
            let mut count = 0;
            for sign in [-1.0_f64, 1.0] {
                let mut x = wx;
                let mut y = wy;
                let step = sign * dx_world.min(dy_world) * 0.5;
                for _ in 0..kl {
                    let v = velocity_field([x, y, params.z_plane]);
                    let vx = v[0];
                    let vy = v[1];
                    let vmag = (vx * vx + vy * vy).sqrt();
                    if vmag < 1e-12 {
                        break;
                    }
                    x += (vx / vmag) * step.abs();
                    y += (vy / vmag) * step.abs();
                    let (px, py) = world_to_pixel(x, y);
                    if px < 0.0 || py < 0.0 || px >= w as f64 || py >= h as f64 {
                        break;
                    }
                    sum += sample_noise(px, py);
                    count += 1;
                }
            }
            let (px0, py0) = world_to_pixel(wx, wy);
            sum += sample_noise(px0, py0);
            count += 1;
            output[row * w + col] = if count > 0 { sum / count as f64 } else { 0.0 };
        }
    }
    LicTexture {
        pixels: output,
        width: w,
        height: h,
    }
}
/// Trace a single streamline through a velocity field using forward Euler integration.
///
/// Starting from `seed`, the integration advances by `dt` each step,
/// up to `max_steps` or until the velocity becomes negligible.
pub fn trace_streamline(
    velocity_field: &dyn Fn(Vec3) -> Vec3,
    seed: Vec3,
    dt: f64,
    max_steps: usize,
) -> Streamline {
    let mut points = Vec::with_capacity(max_steps + 1);
    let mut magnitudes = Vec::with_capacity(max_steps + 1);
    let mut pos = seed;
    let vel0 = velocity_field(pos);
    points.push(pos);
    magnitudes.push(vel0.norm());
    for _ in 0..max_steps {
        let vel = velocity_field(pos);
        if vel.norm_squared() < 1e-20 {
            break;
        }
        pos += vel * dt;
        let v_next = velocity_field(pos);
        points.push(pos);
        magnitudes.push(v_next.norm());
    }
    Streamline::with_magnitudes(points, magnitudes, Color::white())
}
/// Trace multiple streamlines from a set of seed points.
pub fn trace_streamlines_grid(
    velocity_field: &dyn Fn(Vec3) -> Vec3,
    seeds: &[Vec3],
    dt: f64,
    max_steps: usize,
) -> Vec<Streamline> {
    seeds
        .iter()
        .map(|&seed| trace_streamline(velocity_field, seed, dt, max_steps))
        .collect()
}
/// Compute an enhanced LIC texture with multiple passes and optional
/// magnitude modulation.
pub fn compute_enhanced_lic(
    params: &EnhancedLicParams,
    velocity_field: &dyn Fn([f64; 3]) -> [f64; 3],
) -> LicTexture {
    let w = params.base.width;
    let h = params.base.height;
    let dx = (params.base.x_range[1] - params.base.x_range[0]) / w as f64;
    let dy = (params.base.y_range[1] - params.base.y_range[0]) / h as f64;
    let magnitudes: Vec<f64> = (0..h)
        .flat_map(|row| {
            (0..w).map(move |col| {
                let x = params.base.x_range[0] + (col as f64 + 0.5) * dx;
                let y = params.base.y_range[0] + (row as f64 + 0.5) * dy;
                let v = velocity_field([x, y, params.base.z_plane]);
                (v[0] * v[0] + v[1] * v[1]).sqrt()
            })
        })
        .collect();
    let max_mag = magnitudes.iter().cloned().fold(0.0_f64, f64::max);
    let mut tex = compute_lic(&params.base, velocity_field);
    for _ in 1..params.n_passes {
        for p in &mut tex.pixels {
            *p = 1.0 / (1.0 + (-6.0 * (*p - 0.5)).exp());
        }
    }
    if params.modulate_by_magnitude && max_mag > 1e-20 {
        for (p, &mag) in tex.pixels.iter_mut().zip(magnitudes.iter()) {
            *p *= 0.3 + 0.7 * (mag / max_mag);
        }
    }
    let min_v = tex.pixels.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_v = tex.pixels.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if (max_v - min_v) > 1e-12 {
        for p in &mut tex.pixels {
            *p = (*p - min_v) / (max_v - min_v);
        }
    }
    tex
}
/// Detect approximate 2-D flow topology critical points on a grid.
///
/// Scans cells of size `cell_size` across the given 2-D range, looking for
/// velocity sign changes that indicate critical points.
pub fn detect_topology_points(
    x_range: [f64; 2],
    y_range: [f64; 2],
    cell_size: f64,
    z_plane: f64,
    velocity_field: &dyn Fn([f64; 3]) -> [f64; 3],
    stagnation_threshold: f64,
) -> Vec<FlowTopologyPoint> {
    let mut points = Vec::new();
    let mut x = x_range[0];
    while x < x_range[1] {
        let mut y = y_range[0];
        while y < y_range[1] {
            let v = velocity_field([x + cell_size * 0.5, y + cell_size * 0.5, z_plane]);
            let vmag = (v[0] * v[0] + v[1] * v[1]).sqrt();
            if vmag < stagnation_threshold {
                let h = cell_size * 0.1;
                let cx = x + cell_size * 0.5;
                let cy = y + cell_size * 0.5;
                let dvx_dx = (velocity_field([cx + h, cy, z_plane])[0]
                    - velocity_field([cx - h, cy, z_plane])[0])
                    / (2.0 * h);
                let dvy_dy = (velocity_field([cx, cy + h, z_plane])[1]
                    - velocity_field([cx, cy - h, z_plane])[1])
                    / (2.0 * h);
                let det = dvx_dx * dvy_dy;
                let kind = if det < 0.0 {
                    TopologyPointKind::Saddle
                } else {
                    let dvx_dy = (velocity_field([cx, cy + h, z_plane])[0]
                        - velocity_field([cx, cy - h, z_plane])[0])
                        / (2.0 * h);
                    let dvy_dx = (velocity_field([cx + h, cy, z_plane])[1]
                        - velocity_field([cx - h, cy, z_plane])[1])
                        / (2.0 * h);
                    let curl = dvy_dx - dvx_dy;
                    if curl.abs() > 1e-6 {
                        TopologyPointKind::VortexCentre
                    } else {
                        TopologyPointKind::Stagnation
                    }
                };
                points.push(FlowTopologyPoint {
                    position: [cx, cy],
                    kind,
                    velocity_magnitude: vmag,
                });
            }
            y += cell_size;
        }
        x += cell_size;
    }
    points
}
/// Compute the forward FTLE field in a 2-D slice.
///
/// `integration_time` T is the advection duration.  `grid_spacing` controls
/// the spatial resolution of the output field.
pub fn compute_ftle_2d(
    x_range: [f64; 2],
    y_range: [f64; 2],
    z_plane: f64,
    width: usize,
    height: usize,
    integration_time: f64,
    dt: f64,
    velocity_field: &dyn Fn([f64; 3]) -> [f64; 3],
) -> FtleField {
    let dx = (x_range[1] - x_range[0]) / width as f64;
    let dy = (y_range[1] - y_range[0]) / height as f64;
    let n_steps = ((integration_time / dt).ceil() as usize).max(1);
    let advect = |mut x: f64, mut y: f64| -> (f64, f64) {
        for _ in 0..n_steps {
            let v = velocity_field([x, y, z_plane]);
            x += v[0] * dt;
            y += v[1] * dt;
        }
        (x, y)
    };
    let mut values = vec![0.0_f64; width * height];
    for row in 0..height {
        for col in 0..width {
            let x = x_range[0] + (col as f64 + 0.5) * dx;
            let y = y_range[0] + (row as f64 + 0.5) * dy;
            let eps = dx.min(dy) * 0.5;
            let (xpx, ypx) = advect(x + eps, y);
            let (xmx, ymx) = advect(x - eps, y);
            let (xpy, ypy) = advect(x, y + eps);
            let (xmy, ymy) = advect(x, y - eps);
            let dphi_dx = [(xpx - xmx) / (2.0 * eps), (ypx - ymx) / (2.0 * eps)];
            let dphi_dy = [(xpy - xmy) / (2.0 * eps), (ypy - ymy) / (2.0 * eps)];
            let c11 = dphi_dx[0] * dphi_dx[0] + dphi_dy[0] * dphi_dy[0];
            let c12 = dphi_dx[0] * dphi_dx[1] + dphi_dy[0] * dphi_dy[1];
            let c22 = dphi_dx[1] * dphi_dx[1] + dphi_dy[1] * dphi_dy[1];
            let trace = c11 + c22;
            let det = c11 * c22 - c12 * c12;
            let discriminant = ((trace * trace * 0.25 - det).max(0.0)).sqrt();
            let lambda_max = trace * 0.5 + discriminant;
            let ftle = if integration_time > 1e-20 && lambda_max > 1.0 {
                lambda_max.ln() / (2.0 * integration_time)
            } else {
                0.0
            };
            values[row * width + col] = ftle.max(0.0);
        }
    }
    FtleField {
        values,
        width,
        height,
        x_range,
        y_range,
    }
}
/// Generate seed points for streamline tracing according to the given strategy.
///
/// The bounding box is `[min, max]` (inclusive).
pub fn generate_seeds(min: [f64; 3], max: [f64; 3], strategy: SeedStrategy) -> Vec<Vec3> {
    match strategy {
        SeedStrategy::UniformGrid { nx, ny, nz } => {
            let nx = nx.max(1);
            let ny = ny.max(1);
            let nz = nz.max(1);
            let mut seeds = Vec::with_capacity(nx * ny * nz);
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let tx = if nx == 1 {
                            0.5
                        } else {
                            ix as f64 / (nx - 1) as f64
                        };
                        let ty = if ny == 1 {
                            0.5
                        } else {
                            iy as f64 / (ny - 1) as f64
                        };
                        let tz = if nz == 1 {
                            0.5
                        } else {
                            iz as f64 / (nz - 1) as f64
                        };
                        seeds.push(Vec3::new(
                            min[0] + (max[0] - min[0]) * tx,
                            min[1] + (max[1] - min[1]) * ty,
                            min[2] + (max[2] - min[2]) * tz,
                        ));
                    }
                }
            }
            seeds
        }
        SeedStrategy::Random { count, seed: rseed } => {
            let mut state = rseed.wrapping_add(12345);
            let mut seeds = Vec::with_capacity(count);
            let lcg = |s: u64| {
                s.wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407)
            };
            for _ in 0..count {
                state = lcg(state);
                let rx = ((state >> 33) as f64) / (u32::MAX as f64);
                state = lcg(state);
                let ry = ((state >> 33) as f64) / (u32::MAX as f64);
                state = lcg(state);
                let rz = ((state >> 33) as f64) / (u32::MAX as f64);
                seeds.push(Vec3::new(
                    min[0] + (max[0] - min[0]) * rx,
                    min[1] + (max[1] - min[1]) * ry,
                    min[2] + (max[2] - min[2]) * rz,
                ));
            }
            seeds
        }
    }
}
/// Map a velocity magnitude to a colour using a simple heat-map ramp
/// (blue → cyan → green → yellow → red).
///
/// `t` should be in `[0, 1]` (normalised magnitude).
pub fn velocity_color(t: f64) -> Color {
    let t = t.clamp(0.0, 1.0) as f32;
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (0.0_f32, s, 1.0 - s * 0.5)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (0.0, 0.5 + s * 0.5, 0.5 - s * 0.5)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (s, 1.0, 0.0)
    } else {
        let s = (t - 0.75) / 0.25;
        (1.0, 1.0 - s, 0.0)
    };
    Color { r, g, b, a: 1.0 }
}
/// Colour each segment of a streamline by its local velocity magnitude.
///
/// Returns one colour per point (same length as `sl.points`).
/// If `velocity_magnitudes` is empty, all colours are white.
pub fn color_by_velocity(sl: &Streamline, min_mag: f64, max_mag: f64) -> Vec<Color> {
    if sl.velocity_magnitudes.is_empty() {
        return vec![Color::white(); sl.points.len()];
    }
    let range = (max_mag - min_mag).max(1e-20);
    sl.velocity_magnitudes
        .iter()
        .map(|&m| velocity_color((m - min_mag) / range))
        .collect()
}
/// Generate arrow glyphs at regular arc-length intervals along a streamline.
///
/// - `interval`: arc-length spacing between arrows.
/// - `shaft_length`: length of each arrow shaft.
/// - `color`: colour applied to all arrows.
///
/// Returns an empty `Vec` if the streamline has fewer than 2 points.
pub fn arrow_glyphs_along_streamline(
    sl: &Streamline,
    interval: f64,
    shaft_length: f64,
    color: Color,
) -> Vec<ArrowGlyph> {
    if sl.points.len() < 2 {
        return Vec::new();
    }
    let mut glyphs = Vec::new();
    let mut accumulated = 0.0_f64;
    let mut next_arrow = interval * 0.5;
    for win in sl.points.windows(2) {
        let p0 = win[0];
        let p1 = win[1];
        let seg_len = (p1 - p0).norm();
        if seg_len < 1e-15 {
            continue;
        }
        while accumulated + seg_len >= next_arrow {
            let t = (next_arrow - accumulated) / seg_len;
            let pos = p0 + (p1 - p0) * t;
            let dir = (p1 - p0) / seg_len;
            glyphs.push(ArrowGlyph {
                origin: pos,
                direction: dir,
                shaft_length,
                color,
            });
            next_arrow += interval;
        }
        accumulated += seg_len;
    }
    glyphs
}
/// Compute evenly-spaced 2-D streamlines using the Jobard–Lefer algorithm.
///
/// The algorithm:
/// 1. Trace the first streamline from `initial_seed`.
/// 2. Generate candidate seeds at `d_sep` distance from existing streamlines.
/// 3. Accept a candidate if it is at least `d_test` from all existing lines.
/// 4. Repeat until no new candidates can be placed.
///
/// `velocity_field` maps `[x, y, z]` to `[vx, vy, vz]` (z component ignored for 2-D).
pub fn compute_jobard_lefer(
    params: &JobardLeferParams,
    velocity_field: &dyn Fn([f64; 3]) -> [f64; 3],
) -> Vec<Streamline> {
    let vf = |p: Vec3| -> Vec3 {
        let v = velocity_field([p.x, p.y, p.z]);
        Vec3::new(v[0], v[1], 0.0)
    };
    let [xmin, xmax, ymin, ymax] = params.domain;
    let trace_2d = |seed: Vec3| -> Vec<Vec3> {
        let mut pts_fwd = vec![seed];
        let mut pos = seed;
        for _ in 0..params.max_steps {
            let v = vf(pos);
            if v.norm() < 1e-12 {
                break;
            }
            let next = StreamlineIntegrator::RungeKutta4.step(pos, params.step_size, &vf);
            if next.x < xmin || next.x > xmax || next.y < ymin || next.y > ymax {
                break;
            }
            pos = next;
            pts_fwd.push(pos);
        }
        let mut pts_bwd = Vec::new();
        let mut pos = seed;
        for _ in 0..params.max_steps {
            let v = vf(pos);
            if v.norm() < 1e-12 {
                break;
            }
            let next = StreamlineIntegrator::RungeKutta4.step(pos, -params.step_size, &vf);
            if next.x < xmin || next.x > xmax || next.y < ymin || next.y > ymax {
                break;
            }
            pos = next;
            pts_bwd.push(pos);
        }
        pts_bwd.reverse();
        let mut all = pts_bwd;
        all.extend(pts_fwd);
        all
    };
    let min_dist_to_lines = |p: Vec3, lines: &[Streamline]| -> f64 {
        let mut min_d = f64::INFINITY;
        for sl in lines {
            for &q in &sl.points {
                let d = ((p.x - q.x).powi(2) + (p.y - q.y).powi(2)).sqrt();
                if d < min_d {
                    min_d = d;
                }
            }
        }
        min_d
    };
    let mut streamlines: Vec<Streamline> = Vec::new();
    let seed0 = Vec3::new(
        params.initial_seed[0],
        params.initial_seed[1],
        params.z_plane,
    );
    let pts0 = trace_2d(seed0);
    if pts0.len() >= 2 {
        streamlines.push(Streamline::new(pts0, Color::white()));
    }
    let d_sep = params.d_sep;
    let d_test = params.d_test;
    let mut candidate_queue: Vec<Vec3> = Vec::new();
    if let Some(first) = streamlines.first() {
        for &pt in &first.points {
            candidate_queue.push(Vec3::new(pt.x + d_sep, pt.y, pt.z));
            candidate_queue.push(Vec3::new(pt.x - d_sep, pt.y, pt.z));
            candidate_queue.push(Vec3::new(pt.x, pt.y + d_sep, pt.z));
            candidate_queue.push(Vec3::new(pt.x, pt.y - d_sep, pt.z));
        }
    }
    let max_lines = 200usize;
    let mut queue_idx = 0;
    while queue_idx < candidate_queue.len() && streamlines.len() < max_lines {
        let cand = candidate_queue[queue_idx];
        queue_idx += 1;
        if cand.x < xmin || cand.x > xmax || cand.y < ymin || cand.y > ymax {
            continue;
        }
        if min_dist_to_lines(cand, &streamlines) < d_test {
            continue;
        }
        let pts = trace_2d(cand);
        if pts.len() < 2 {
            continue;
        }
        let too_close = pts
            .iter()
            .step_by(5)
            .any(|&p| min_dist_to_lines(p, &streamlines) < d_test);
        if too_close {
            continue;
        }
        for &pt in pts.iter().step_by(4) {
            candidate_queue.push(Vec3::new(pt.x + d_sep, pt.y, pt.z));
            candidate_queue.push(Vec3::new(pt.x - d_sep, pt.y, pt.z));
            candidate_queue.push(Vec3::new(pt.x, pt.y + d_sep, pt.z));
            candidate_queue.push(Vec3::new(pt.x, pt.y - d_sep, pt.z));
        }
        streamlines.push(Streamline::new(pts, Color::white()));
    }
    streamlines
}
