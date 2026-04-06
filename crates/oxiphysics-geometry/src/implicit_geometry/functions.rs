//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    ContactPoint, IsoMeshResult, RayMarchResult, SdfCollisionResult, SdfGrid, Triangle,
};

#[inline]
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn len(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
#[inline]
pub(super) fn norm(v: [f64; 3]) -> [f64; 3] {
    let l = len(v).max(1e-300);
    [v[0] / l, v[1] / l, v[2] / l]
}
#[inline]
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn clamp01(v: f64) -> f64 {
    v.clamp(0.0, 1.0)
}
#[inline]
pub(super) fn clamp_f(v: f64, lo: f64, hi: f64) -> f64 {
    v.clamp(lo, hi)
}
#[inline]
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn len2(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}
#[inline]
pub(super) fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
/// Signed distance function trait.
///
/// Implementors return the signed distance from point `p` to the surface:
/// negative inside, zero on surface, positive outside.
pub trait Sdf: Send + Sync {
    /// Evaluate the signed distance at `p`.
    fn dist(&self, p: [f64; 3]) -> f64;
}
/// Polynomial smooth-union kernel (Quilez).
///
/// Returns a value ≤ min(a, b) with a smooth blend of width `k`.
pub fn sdf_smooth_union(a: f64, b: f64, k: f64) -> f64 {
    if k < 1e-15 {
        return a.min(b);
    }
    let h = (0.5 + 0.5 * (b - a) / k).clamp(0.0, 1.0);
    a * h + b * (1.0 - h) - k * h * (1.0 - h)
}
/// Polynomial smooth-intersection kernel.
pub fn sdf_smooth_intersection(a: f64, b: f64, k: f64) -> f64 {
    if k < 1e-15 {
        return a.max(b);
    }
    let h = (0.5 - 0.5 * (b - a) / k).clamp(0.0, 1.0);
    a * (1.0 - h) + b * h + k * h * (1.0 - h)
}
/// Polynomial smooth-difference kernel (subtract `b` from `a`).
pub fn sdf_smooth_difference(a: f64, b: f64, k: f64) -> f64 {
    sdf_smooth_intersection(a, -b, k)
}
/// Estimate the gradient of an SDF at `p` via central finite differences.
///
/// Returns the unnormalized gradient vector (norm ≈ 1 for exact SDFs).
pub fn sdf_gradient<S: Sdf>(sdf: &S, p: [f64; 3], eps: f64) -> [f64; 3] {
    let dx = sdf.dist([p[0] + eps, p[1], p[2]]) - sdf.dist([p[0] - eps, p[1], p[2]]);
    let dy = sdf.dist([p[0], p[1] + eps, p[2]]) - sdf.dist([p[0], p[1] - eps, p[2]]);
    let dz = sdf.dist([p[0], p[1], p[2] + eps]) - sdf.dist([p[0], p[1], p[2] - eps]);
    [dx / (2.0 * eps), dy / (2.0 * eps), dz / (2.0 * eps)]
}
/// Compute the unit outward surface normal at `p` by finite-difference gradient.
pub fn sdf_normal<S: Sdf>(sdf: &S, p: [f64; 3], eps: f64) -> [f64; 3] {
    norm(sdf_gradient(sdf, p, eps))
}
/// Mean curvature of an SDF at `p` estimated from the Laplacian of the SDF.
///
/// H ≈ −(1/2) · ∇²d / |∇d|  (signed; negative = concave toward outside).
pub fn sdf_mean_curvature<S: Sdf>(sdf: &S, p: [f64; 3], eps: f64) -> f64 {
    let d = sdf.dist(p);
    let dxp = sdf.dist([p[0] + eps, p[1], p[2]]);
    let dxm = sdf.dist([p[0] - eps, p[1], p[2]]);
    let dyp = sdf.dist([p[0], p[1] + eps, p[2]]);
    let dym = sdf.dist([p[0], p[1] - eps, p[2]]);
    let dzp = sdf.dist([p[0], p[1], p[2] + eps]);
    let dzm = sdf.dist([p[0], p[1], p[2] - eps]);
    let laplacian = (dxp + dxm + dyp + dym + dzp + dzm - 6.0 * d) / (eps * eps);
    let g = sdf_gradient(sdf, p, eps);
    let g_len = len(g).max(1e-30);
    -0.5 * laplacian / g_len
}
/// March a ray defined by `origin + t * direction` through an SDF.
///
/// * `max_steps` – step budget.
/// * `max_t`     – maximum travel distance.
/// * `surf_eps`  – surface-hit threshold.
///
/// Uses the sphere-tracing algorithm (Hart 1996).
pub fn ray_march<S: Sdf>(
    sdf: &S,
    origin: [f64; 3],
    direction: [f64; 3],
    max_steps: u32,
    max_t: f64,
    surf_eps: f64,
) -> RayMarchResult {
    let dir = norm(direction);
    let mut t = 0.0_f64;
    for step in 0..max_steps {
        let p = add(origin, scale(dir, t));
        let d = sdf.dist(p);
        if d.abs() < surf_eps {
            return RayMarchResult {
                hit: true,
                t,
                point: p,
                steps: step + 1,
            };
        }
        t += d.abs().max(surf_eps * 0.1);
        if t > max_t {
            break;
        }
    }
    RayMarchResult {
        hit: false,
        t,
        point: add(origin, scale(dir, t)),
        steps: max_steps,
    }
}
/// Ray march with over-relaxation for faster convergence.
///
/// Uses the relaxed sphere tracing of Keinert et al. (2014) with parameter `omega`.
/// A safe value is `omega = 1.2`.
pub fn ray_march_relaxed<S: Sdf>(
    sdf: &S,
    origin: [f64; 3],
    direction: [f64; 3],
    max_steps: u32,
    max_t: f64,
    surf_eps: f64,
    omega: f64,
) -> RayMarchResult {
    let dir = norm(direction);
    let mut t = 0.0_f64;
    let mut prev_d = f64::INFINITY;
    for step in 0..max_steps {
        let p = add(origin, scale(dir, t));
        let d = sdf.dist(p);
        if d.abs() < surf_eps {
            return RayMarchResult {
                hit: true,
                t,
                point: p,
                steps: step + 1,
            };
        }
        let step_length = if omega * d.abs() <= d.abs() + prev_d.abs() {
            omega * d.abs()
        } else {
            d.abs()
        };
        prev_d = d;
        t += step_length.max(surf_eps * 0.1);
        if t > max_t {
            break;
        }
    }
    RayMarchResult {
        hit: false,
        t,
        point: add(origin, scale(dir, t)),
        steps: max_steps,
    }
}
pub(super) const MC_EDGE_TABLE: [u16; 256] = [
    0x000, 0xfff, 0x00f, 0xff0, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0,
    0x00f, 0xfff, 0x000, 0x00f, 0xff0, 0x000, 0xfff, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0xff0, 0x00f,
    0xfff, 0x000, 0xf0f, 0x0f0, 0xf00, 0x0ff, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x000, 0xfff, 0x00f,
    0xff0, 0xf0f, 0x0f0, 0xf00, 0x0ff, 0xfff, 0x000, 0xff0, 0x00f, 0x0ff, 0xf00, 0x0f0, 0xf0f,
    0x00f, 0xff0, 0x000, 0xfff, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff, 0x000, 0xf00,
    0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff, 0x000, 0x000, 0xfff, 0x00f, 0xff0, 0x0f0, 0xf0f,
    0x0ff, 0xf00, 0xff0, 0x00f, 0xfff, 0x000, 0xf0f, 0x0f0, 0xf00, 0x0ff, 0x00f, 0xff0, 0x000,
    0xfff, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0x00f, 0xff0, 0x000, 0xfff,
    0xf0f, 0x0f0, 0xf00, 0x0ff, 0xfff, 0x000, 0xff0, 0x00f, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x000,
    0xfff, 0x00f, 0xff0, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff, 0x000, 0x000, 0xfff,
    0x00f, 0xff0, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff,
    0x000, 0x00f, 0xff0, 0x000, 0xfff, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0xff0, 0x00f, 0xfff, 0x000,
    0xf0f, 0x0f0, 0xf00, 0x0ff, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x000, 0xfff, 0x00f, 0xff0, 0xf0f,
    0x0f0, 0xf00, 0x0ff, 0xfff, 0x000, 0xff0, 0x00f, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0x00f, 0xff0,
    0x000, 0xfff, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff, 0x000, 0xf00, 0x0ff, 0xf0f,
    0x0f0, 0xff0, 0x00f, 0xfff, 0x000, 0x000, 0xfff, 0x00f, 0xff0, 0x0f0, 0xf0f, 0x0ff, 0xf00,
    0xff0, 0x00f, 0xfff, 0x000, 0xf0f, 0x0f0, 0xf00, 0x0ff, 0x00f, 0xff0, 0x000, 0xfff, 0x0ff,
    0xf00, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x0f0, 0xf0f, 0x00f, 0xff0, 0x000, 0xfff, 0xf0f, 0x0f0,
    0xf00, 0x0ff, 0xfff, 0x000, 0xff0, 0x00f, 0x0f0, 0xf0f, 0x0ff, 0xf00, 0x000, 0xfff, 0x00f,
    0xff0, 0xf00, 0x0ff, 0xf0f, 0x0f0, 0xff0, 0x00f, 0xfff, 0x000,
];
/// Helper: linearly interpolate between two cube vertices on an edge.
pub(super) fn interp_edge(p0: [f64; 3], d0: f64, p1: [f64; 3], d1: f64) -> [f64; 3] {
    if (d1 - d0).abs() < 1e-12 {
        return p0;
    }
    let t = d0 / (d0 - d1);
    [
        p0[0] + t * (p1[0] - p0[0]),
        p0[1] + t * (p1[1] - p0[1]),
        p0[2] + t * (p1[2] - p0[2]),
    ]
}
/// Extract the zero isosurface from an [`SdfGrid`] using Marching Cubes.
///
/// Returns an [`IsoMeshResult`] containing all output triangles.
/// This is an independent implementation distinct from the one in `implicit`.
pub fn extract_isosurface(grid: &SdfGrid) -> IsoMeshResult {
    let mut result = IsoMeshResult::default();
    if grid.nx < 2 || grid.ny < 2 || grid.nz < 2 {
        return result;
    }
    for iz in 0..grid.nz - 1 {
        for iy in 0..grid.ny - 1 {
            for ix in 0..grid.nx - 1 {
                let corners: [[usize; 3]; 8] = [
                    [ix, iy, iz],
                    [ix + 1, iy, iz],
                    [ix + 1, iy + 1, iz],
                    [ix, iy + 1, iz],
                    [ix, iy, iz + 1],
                    [ix + 1, iy, iz + 1],
                    [ix + 1, iy + 1, iz + 1],
                    [ix, iy + 1, iz + 1],
                ];
                let vals: [f64; 8] =
                    std::array::from_fn(|i| grid.get(corners[i][0], corners[i][1], corners[i][2]));
                let pts: [[f64; 3]; 8] = std::array::from_fn(|i| {
                    grid.cell_center(corners[i][0], corners[i][1], corners[i][2])
                });
                let mut cube_idx: usize = 0;
                for i in 0..8 {
                    if vals[i] < 0.0 {
                        cube_idx |= 1 << i;
                    }
                }
                if MC_EDGE_TABLE[cube_idx] == 0 {
                    continue;
                }
                let edge_mask = MC_EDGE_TABLE[cube_idx];
                let mut verts: [[f64; 3]; 12] = [[0.0; 3]; 12];
                let edges = [
                    (0, 1),
                    (1, 2),
                    (2, 3),
                    (3, 0),
                    (4, 5),
                    (5, 6),
                    (6, 7),
                    (7, 4),
                    (0, 4),
                    (1, 5),
                    (2, 6),
                    (3, 7),
                ];
                for (e, &(i, j)) in edges.iter().enumerate() {
                    if edge_mask & (1 << e) != 0 {
                        verts[e] = interp_edge(pts[i], vals[i], pts[j], vals[j]);
                    }
                }
                let active: Vec<usize> = (0..12).filter(|&e| edge_mask & (1 << e) != 0).collect();
                if active.len() >= 3 {
                    for i in 1..active.len() - 1 {
                        let tri = Triangle {
                            v0: verts[active[0]],
                            v1: verts[active[i]],
                            v2: verts[active[i + 1]],
                        };
                        result.triangles.push(tri);
                        result.vertex_count += 3;
                    }
                }
            }
        }
    }
    result
}
/// Test whether point `p` is inside the SDF and return contact data.
///
/// Returns `Some(result)` if the point is inside (`dist < 0`), else `None`.
pub fn sdf_point_query<S: Sdf>(sdf: &S, p: [f64; 3], eps: f64) -> Option<SdfCollisionResult> {
    let d = sdf.dist(p);
    if d >= 0.0 {
        return None;
    }
    let n = sdf_normal(sdf, p, eps);
    let contact = add(p, scale(n, -d));
    Some(SdfCollisionResult {
        depth: -d,
        normal: n,
        contact_point: contact,
    })
}
/// Simple swept-sphere vs SDF collision proxy.
///
/// Move sphere of `radius` along `velocity * dt` from `center`.
/// Returns the safe time-of-impact in \[0, dt\] and contact data if hit.
pub fn sdf_sphere_sweep<S: Sdf>(
    sdf: &S,
    center: [f64; 3],
    radius: f64,
    velocity: [f64; 3],
    dt: f64,
    eps: f64,
) -> Option<(f64, SdfCollisionResult)> {
    let speed = len(velocity);
    if speed < 1e-15 {
        return None;
    }
    let dir = scale(velocity, 1.0 / speed);
    let max_t = speed * dt;
    let result = ray_march_with_offset(sdf, center, dir, radius, 128, max_t, eps);
    if result.hit {
        let toi = result.t / speed;
        let n = sdf_normal(sdf, result.point, eps);
        Some((
            toi,
            SdfCollisionResult {
                depth: 0.0,
                normal: n,
                contact_point: result.point,
            },
        ))
    } else {
        None
    }
}
pub(super) fn ray_march_with_offset<S: Sdf>(
    sdf: &S,
    origin: [f64; 3],
    dir: [f64; 3],
    offset: f64,
    max_steps: u32,
    max_t: f64,
    eps: f64,
) -> RayMarchResult {
    let mut t = 0.0_f64;
    for step in 0..max_steps {
        let p = add(origin, scale(dir, t));
        let d = sdf.dist(p) - offset;
        if d.abs() < eps {
            return RayMarchResult {
                hit: true,
                t,
                point: p,
                steps: step + 1,
            };
        }
        t += d.abs().max(eps * 0.1);
        if t > max_t {
            break;
        }
    }
    RayMarchResult {
        hit: false,
        t,
        point: add(origin, scale(dir, t)),
        steps: max_steps,
    }
}
/// A simple hash-based value-noise function for SDF displacement.
///
/// Produces values in \[−1, 1\] from integer lattice positions.
pub(super) fn value_noise_lattice(ix: i32, iy: i32, iz: i32) -> f64 {
    let n = ix
        .wrapping_mul(1619)
        .wrapping_add(iy.wrapping_mul(31337))
        .wrapping_add(iz.wrapping_mul(6971))
        .wrapping_add(1376312589i32);
    let n = n ^ (n << 13);
    let n = n.wrapping_mul(
        n.wrapping_mul(n.wrapping_mul(15731).wrapping_add(789221))
            .wrapping_add(1376312589),
    );
    (n & 0x7fffffff) as f64 / 1073741824.0 - 1.0
}
/// Smooth step (fade) for noise interpolation.
#[inline]
pub(super) fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}
/// Sample trilinear value noise at continuous position `p`.
pub(super) fn value_noise(p: [f64; 3]) -> f64 {
    let ix = p[0].floor() as i32;
    let iy = p[1].floor() as i32;
    let iz = p[2].floor() as i32;
    let fx = fade(p[0] - ix as f64);
    let fy = fade(p[1] - iy as f64);
    let fz = fade(p[2] - iz as f64);
    let v000 = value_noise_lattice(ix, iy, iz);
    let v100 = value_noise_lattice(ix + 1, iy, iz);
    let v010 = value_noise_lattice(ix, iy + 1, iz);
    let v110 = value_noise_lattice(ix + 1, iy + 1, iz);
    let v001 = value_noise_lattice(ix, iy, iz + 1);
    let v101 = value_noise_lattice(ix + 1, iy, iz + 1);
    let v011 = value_noise_lattice(ix, iy + 1, iz + 1);
    let v111 = value_noise_lattice(ix + 1, iy + 1, iz + 1);
    let c00 = v000 * (1.0 - fx) + v100 * fx;
    let c10 = v010 * (1.0 - fx) + v110 * fx;
    let c01 = v001 * (1.0 - fx) + v101 * fx;
    let c11 = v011 * (1.0 - fx) + v111 * fx;
    let c0 = c00 * (1.0 - fy) + c10 * fy;
    let c1 = c01 * (1.0 - fy) + c11 * fy;
    c0 * (1.0 - fz) + c1 * fz
}
/// Fractal (octaved) value noise: fBm with `octaves` levels.
pub fn fbm_noise(p: [f64; 3], octaves: u32, lacunarity: f64, gain: f64) -> f64 {
    let mut value = 0.0_f64;
    let mut amplitude = 0.5_f64;
    let mut frequency = 1.0_f64;
    for _ in 0..octaves {
        value += amplitude * value_noise(scale(p, frequency));
        amplitude *= gain;
        frequency *= lacunarity;
    }
    value
}
/// Project point `p` to the nearest surface point of the SDF using
/// gradient descent on the absolute distance.
///
/// Returns the approximate surface point and the final signed distance.
pub fn sdf_closest_point<S: Sdf>(sdf: &S, p: [f64; 3], eps: f64, max_iter: u32) -> ([f64; 3], f64) {
    let mut q = p;
    for _ in 0..max_iter {
        let d = sdf.dist(q);
        if d.abs() < eps {
            break;
        }
        let g = sdf_gradient(sdf, q, eps);
        let step = d;
        q = sub(q, scale(g, step));
    }
    let d_final = sdf.dist(q);
    (q, d_final)
}
/// Estimate ambient occlusion at surface point `p` with normal `n` from an SDF.
///
/// Uses the analytic AO formula of Quilez (2016): sums SDF samples along the
/// normal direction weighted by 1/distance.
pub fn sdf_ambient_occlusion<S: Sdf>(
    sdf: &S,
    p: [f64; 3],
    n: [f64; 3],
    num_samples: u32,
    max_dist: f64,
) -> f64 {
    let mut occ = 0.0_f64;
    let mut scale_factor = 1.0_f64;
    for i in 1..=num_samples {
        let h = i as f64 / num_samples as f64 * max_dist;
        let sample_pt = add(p, scale(n, h));
        let d = sdf.dist(sample_pt);
        occ += (h - d) * scale_factor;
        scale_factor *= 0.5;
    }
    (1.0 - 2.0 * occ / max_dist).clamp(0.0, 1.0)
}
/// Estimate a soft-shadow factor for a point `p` toward a light at `light_pos`.
///
/// Returns a value in \[0, 1\]: 0 = fully shadowed, 1 = fully lit.
/// Uses Quilez's soft-shadow ray-march formula with `k` controlling penumbra sharpness.
pub fn sdf_soft_shadow<S: Sdf>(sdf: &S, p: [f64; 3], light_pos: [f64; 3], k: f64, eps: f64) -> f64 {
    let to_light = sub(light_pos, p);
    let dist_to_light = len(to_light);
    let dir = scale(to_light, 1.0 / dist_to_light);
    let mut t = 2.0 * eps;
    let mut shadow = 1.0_f64;
    let max_steps = 64u32;
    for _ in 0..max_steps {
        if t >= dist_to_light {
            break;
        }
        let q = add(p, scale(dir, t));
        let d = sdf.dist(q);
        if d < eps {
            return 0.0;
        }
        shadow = shadow.min(k * d / t);
        t += d;
    }
    shadow.clamp(0.0, 1.0)
}
/// Generate a contact manifold between a convex set of sample points and an SDF.
///
/// `points` are candidate contact points in world space.
/// Returns all contacts with `dist < tolerance`.
pub fn sdf_contact_manifold<S: Sdf>(
    sdf: &S,
    points: &[[f64; 3]],
    tolerance: f64,
    eps: f64,
) -> Vec<ContactPoint> {
    let mut contacts = Vec::new();
    for &p in points {
        let d = sdf.dist(p);
        if d < tolerance {
            let n = sdf_normal(sdf, p, eps);
            contacts.push(ContactPoint {
                position: p,
                normal: n,
                depth: (-d).max(0.0),
            });
        }
    }
    contacts
}
/// Count the number of grid cells that straddle the zero isosurface.
pub fn count_surface_cells(grid: &SdfGrid) -> usize {
    let mut count = 0usize;
    for iz in 0..grid.nz.saturating_sub(1) {
        for iy in 0..grid.ny.saturating_sub(1) {
            for ix in 0..grid.nx.saturating_sub(1) {
                let signs: [bool; 8] = [
                    grid.get(ix, iy, iz) < 0.0,
                    grid.get(ix + 1, iy, iz) < 0.0,
                    grid.get(ix + 1, iy + 1, iz) < 0.0,
                    grid.get(ix, iy + 1, iz) < 0.0,
                    grid.get(ix, iy, iz + 1) < 0.0,
                    grid.get(ix + 1, iy, iz + 1) < 0.0,
                    grid.get(ix + 1, iy + 1, iz + 1) < 0.0,
                    grid.get(ix, iy + 1, iz + 1) < 0.0,
                ];
                let all_in = signs.iter().all(|&s| s);
                let all_out = signs.iter().all(|&s| !s);
                if !all_in && !all_out {
                    count += 1;
                }
            }
        }
    }
    count
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::implicit_geometry::BoolOp;
    use crate::implicit_geometry::SdfBend;
    use crate::implicit_geometry::SdfBoundedProxy;
    use crate::implicit_geometry::SdfEllipsoid;
    use crate::implicit_geometry::SdfExtrude;
    use crate::implicit_geometry::SdfGyroid;
    use crate::implicit_geometry::SdfHexagonalPrism;
    use crate::implicit_geometry::SdfLineSegment;
    use crate::implicit_geometry::SdfNode;
    use crate::implicit_geometry::SdfNoiseDisplace;
    use crate::implicit_geometry::SdfRepeat;
    use crate::implicit_geometry::SdfRevolution;
    use crate::implicit_geometry::SdfScale;
    use crate::implicit_geometry::SdfShell;
    use crate::implicit_geometry::SdfTranslate;
    use crate::implicit_geometry::SdfTwist;
    use crate::implicit_geometry::types::SdfBox;
    use crate::implicit_geometry::types::SdfDifference;
    use crate::implicit_geometry::types::SdfIntersection;
    use crate::implicit_geometry::types::SdfOffset;
    use crate::implicit_geometry::types::SdfPlane;
    use crate::implicit_geometry::types::SdfSphere;
    use crate::implicit_geometry::types::SdfTorus;
    use crate::implicit_geometry::types::SdfUnion;
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_sphere_center_is_minus_radius() {
        let s = SdfSphere::new(2.0);
        approx_eq(s.dist([0.0, 0.0, 0.0]), -2.0, 1e-12);
    }
    #[test]
    fn test_sphere_on_surface() {
        let s = SdfSphere::new(1.0);
        assert!(approx_eq(s.dist([1.0, 0.0, 0.0]), 0.0, 1e-12));
    }
    #[test]
    fn test_sphere_outside() {
        let s = SdfSphere::new(1.0);
        assert!(s.dist([2.0, 0.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_plane_on_surface() {
        let p = SdfPlane::new([0.0, 1.0, 0.0], 0.0);
        assert!(approx_eq(p.dist([1.0, 0.0, 5.0]), 0.0, 1e-12));
    }
    #[test]
    fn test_plane_above() {
        let p = SdfPlane::new([0.0, 1.0, 0.0], 0.0);
        assert!(p.dist([0.0, 1.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_plane_below() {
        let p = SdfPlane::new([0.0, 1.0, 0.0], 0.0);
        assert!(p.dist([0.0, -1.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_ellipsoid_inside() {
        let e = SdfEllipsoid::new(2.0, 1.0, 1.5);
        assert!(e.dist([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_ellipsoid_surface_x() {
        let e = SdfEllipsoid::new(2.0, 1.0, 1.5);
        let d = e.dist([2.0, 0.0, 0.0]);
        assert!(approx_eq(d, 0.0, 1e-6), "d={d}");
    }
    #[test]
    fn test_box_inside() {
        let b = SdfBox::new(1.0, 1.0, 1.0);
        assert!(b.dist([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_box_corner() {
        let b = SdfBox::new(1.0, 1.0, 1.0);
        let d = b.dist([1.0, 1.0, 1.0]);
        assert!(approx_eq(d, 0.0, 1e-12));
    }
    #[test]
    fn test_box_outside_face() {
        let b = SdfBox::new(1.0, 2.0, 1.0);
        let d = b.dist([2.0, 0.0, 0.0]);
        assert!(approx_eq(d, 1.0, 1e-12));
    }
    #[test]
    fn test_torus_tube_surface() {
        let t = SdfTorus::new(3.0, 1.0);
        let d = t.dist([4.0, 0.0, 0.0]);
        assert!(approx_eq(d, 0.0, 1e-12), "d={d}");
    }
    #[test]
    fn test_torus_inside_tube() {
        let t = SdfTorus::new(3.0, 1.0);
        assert!(t.dist([3.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_line_segment_midpoint() {
        let s = SdfLineSegment::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.5);
        let d = s.dist([1.0, 0.0, 0.0]);
        assert!(approx_eq(d, -0.5, 1e-12));
    }
    #[test]
    fn test_line_segment_endpoint() {
        let s = SdfLineSegment::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.0);
        assert!(approx_eq(s.dist([0.0, 0.0, 0.0]), 0.0, 1e-12));
    }
    #[test]
    fn test_union_inside_either() {
        let a = SdfSphere::new(1.0);
        let b = SdfTranslate::new(SdfSphere::new(1.0), [3.0, 0.0, 0.0]);
        let u = SdfUnion::new(a, b);
        assert!(u.dist([0.0, 0.0, 0.0]) < 0.0);
        assert!(u.dist([3.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_intersection_inside_both() {
        let a = SdfSphere::new(2.0);
        let b = SdfBox::new(1.0, 1.0, 1.0);
        let i = SdfIntersection::new(a, b);
        assert!(i.dist([0.0, 0.0, 0.0]) < 0.0);
        assert!(i.dist([3.0, 0.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_difference_removes_inner() {
        let outer = SdfSphere::new(2.0);
        let inner = SdfSphere::new(1.0);
        let d = SdfDifference::new(outer, inner);
        assert!(d.dist([0.0, 0.0, 0.0]) > 0.0);
        assert!(d.dist([1.5, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_smooth_union_le_min() {
        let a = 0.5_f64;
        let b = 0.3_f64;
        let su = sdf_smooth_union(a, b, 0.2);
        assert!(su <= a.min(b) + 1e-12);
    }
    #[test]
    fn test_smooth_union_zero_k() {
        let a = 0.5_f64;
        let b = 0.3_f64;
        assert!(approx_eq(sdf_smooth_union(a, b, 0.0), a.min(b), 1e-12));
    }
    #[test]
    fn test_smooth_intersection_in_range() {
        let a = 0.5_f64;
        let b = 0.3_f64;
        let k = 0.2_f64;
        let si = sdf_smooth_intersection(a, b, k);
        assert!(si >= a.min(b) - k, "si={si} below lower bound");
        assert!(si <= a.max(b) + k, "si={si} above upper bound");
        assert!(approx_eq(
            sdf_smooth_intersection(a, b, 0.0),
            a.max(b),
            1e-12
        ));
    }
    #[test]
    fn test_sphere_gradient_unit_length() {
        let s = SdfSphere::new(1.0);
        let g = sdf_gradient(&s, [0.6, 0.5, 0.3], 1e-5);
        let l = len(g);
        assert!(approx_eq(l, 1.0, 0.01), "gradient length={l}");
    }
    #[test]
    fn test_sphere_normal_points_outward() {
        let s = SdfSphere::new(1.0);
        let p = [1.0, 0.0, 0.0];
        let n = sdf_normal(&s, p, 1e-5);
        assert!(n[0] > 0.9);
    }
    #[test]
    fn test_offset_grows_sphere() {
        let s = SdfSphere::new(1.0);
        let grown = SdfOffset::new(s, 0.5);
        assert!(grown.dist([1.0, 0.0, 0.0]) < 0.0);
        let d = grown.dist([1.5, 0.0, 0.0]);
        assert!(approx_eq(d, 0.0, 1e-12));
    }
    #[test]
    fn test_shell_makes_thin_surface() {
        let s = SdfSphere::new(1.0);
        let shell = SdfShell::new(s, 0.1);
        let d = shell.dist([1.0, 0.0, 0.0]);
        assert!(approx_eq(d, -0.05, 1e-12), "d={d}");
    }
    #[test]
    fn test_translate_moves_sphere() {
        let s = SdfTranslate::new(SdfSphere::new(1.0), [5.0, 0.0, 0.0]);
        assert!(s.dist([5.0, 0.0, 0.0]) < 0.0);
        assert!(s.dist([0.0, 0.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_scale_doubles_sphere() {
        let s = SdfScale::new(SdfSphere::new(1.0), 2.0);
        let d = s.dist([2.0, 0.0, 0.0]);
        assert!(approx_eq(d, 0.0, 1e-12));
    }
    #[test]
    fn test_ray_march_hits_sphere() {
        let s = SdfSphere::new(1.0);
        let res = ray_march(&s, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 128, 20.0, 1e-4);
        assert!(res.hit, "ray should hit sphere");
        assert!(approx_eq(res.t, 4.0, 0.01), "t={}", res.t);
    }
    #[test]
    fn test_ray_march_misses_sphere() {
        let s = SdfSphere::new(1.0);
        let res = ray_march(&s, [-5.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 128, 20.0, 1e-4);
        assert!(!res.hit);
    }
    #[test]
    fn test_ray_march_relaxed_hits() {
        let s = SdfSphere::new(1.0);
        let res = ray_march_relaxed(&s, [-5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 512, 20.0, 1e-4, 1.0);
        assert!(res.hit, "relaxed ray march should hit sphere");
    }
    #[test]
    fn test_grid_from_sphere_sdf() {
        let s = SdfSphere::new(1.0);
        let g = SdfGrid::from_sdf(&s, 10, 10, 10, [-2.0, -2.0, -2.0], 0.4);
        assert!(g.get(5, 5, 5) < 0.0);
    }
    #[test]
    fn test_grid_dilate_erode_inverse() {
        let s = SdfSphere::new(1.0);
        let mut g = SdfGrid::from_sdf(&s, 6, 6, 6, [-1.5, -1.5, -1.5], 0.5);
        let v_orig = g.get(3, 3, 3);
        g.dilate(0.3);
        g.erode(0.3);
        assert!(approx_eq(g.get(3, 3, 3), v_orig, 1e-12));
    }
    #[test]
    fn test_grid_interpolate_center() {
        let s = SdfSphere::new(1.0);
        let g = SdfGrid::from_sdf(&s, 20, 20, 20, [-2.0, -2.0, -2.0], 0.2);
        let interp = g.interpolate([0.0, 0.0, 0.0]);
        assert!(interp < -0.5, "interp={interp}");
    }
    #[test]
    fn test_extract_isosurface_has_triangles() {
        let s = SdfSphere::new(1.0);
        let g = SdfGrid::from_sdf(&s, 10, 10, 10, [-2.0, -2.0, -2.0], 0.4);
        let mesh = extract_isosurface(&g);
        assert!(
            !mesh.triangles.is_empty(),
            "isosurface should produce triangles"
        );
    }
    #[test]
    fn test_count_surface_cells() {
        let s = SdfSphere::new(1.0);
        let g = SdfGrid::from_sdf(&s, 10, 10, 10, [-2.0, -2.0, -2.0], 0.4);
        let cnt = count_surface_cells(&g);
        assert!(cnt > 0, "should have surface cells");
    }
    #[test]
    fn test_point_query_inside() {
        let s = SdfSphere::new(2.0);
        let res = sdf_point_query(&s, [0.0, 0.0, 0.0], 1e-5);
        assert!(res.is_some());
        assert!(res.unwrap().depth > 0.0);
    }
    #[test]
    fn test_point_query_outside() {
        let s = SdfSphere::new(1.0);
        let res = sdf_point_query(&s, [5.0, 0.0, 0.0], 1e-5);
        assert!(res.is_none());
    }
    #[test]
    fn test_contact_manifold_sphere() {
        let s = SdfSphere::new(2.0);
        let pts: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let contacts = sdf_contact_manifold(&s, &pts, 0.0, 1e-5);
        assert_eq!(contacts.len(), 2, "{} contacts", contacts.len());
    }
    #[test]
    fn test_noise_displace_varies() {
        let s = SdfSphere::new(1.0);
        let noisy = SdfNoiseDisplace::new(s, 3.0, 0.2, 4);
        let d1 = noisy.dist([1.0, 0.0, 0.0]);
        let d2 = noisy.dist([0.0, 1.0, 0.0]);
        let _ = (d1, d2);
    }
    #[test]
    fn test_fbm_noise_bounded() {
        for i in 0..20 {
            let v = fbm_noise(
                [i as f64 * 0.3, i as f64 * 0.7, i as f64 * 0.5],
                4,
                2.0,
                0.5,
            );
            assert!(v.abs() <= 2.0, "fbm out of range: {v}");
        }
    }
    #[test]
    fn test_sdf_node_union() {
        let a = SdfNode::leaf(SdfSphere::new(1.0));
        let b = SdfNode::leaf(SdfTranslate::new(SdfSphere::new(1.0), [4.0, 0.0, 0.0]));
        let u = SdfNode::combine(BoolOp::Union, a, b, 0.0);
        assert!(u.eval([0.0, 0.0, 0.0]) < 0.0);
        assert!(u.eval([4.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_sdf_node_smooth_union() {
        let a = SdfNode::leaf(SdfSphere::new(1.0));
        let b = SdfNode::leaf(SdfTranslate::new(SdfSphere::new(1.0), [1.5, 0.0, 0.0]));
        let su = SdfNode::combine(BoolOp::SmoothUnion, a, b, 0.5);
        let d = su.eval([0.75, 0.0, 0.0]);
        assert!(d < 0.0, "d={d}");
    }
    #[test]
    fn test_twist_identity_at_zero_strength() {
        let s = SdfBox::new(1.0, 1.0, 1.0);
        let t = SdfTwist::new(SdfBox::new(1.0, 1.0, 1.0), 0.0);
        let p = [0.5, 0.3, 0.2];
        assert!(approx_eq(s.dist(p), t.dist(p), 1e-12));
    }
    #[test]
    fn test_bend_identity_at_zero_strength() {
        let s = SdfBox::new(1.0, 1.0, 1.0);
        let b = SdfBend::new(SdfBox::new(1.0, 1.0, 1.0), 0.0);
        let p = [0.5, 0.3, 0.2];
        assert!(approx_eq(s.dist(p), b.dist(p), 1e-12));
    }
    #[test]
    fn test_repeat_periodicity() {
        let s = SdfRepeat::new(SdfSphere::new(0.4), 2.0, 2.0, 2.0);
        let d0 = s.dist([0.0, 0.0, 0.0]);
        let d1 = s.dist([2.0, 0.0, 0.0]);
        assert!(approx_eq(d0, d1, 1e-10), "d0={d0} d1={d1}");
    }
    #[test]
    fn test_gyroid_surface_exists() {
        let g = SdfGyroid::new(1.0, 0.05);
        let d0 = g.dist([0.0, 0.0, 0.0]);
        let d1 = g.dist([0.25, 0.0, 0.0]);
        let _ = (d0, d1);
        assert!(d0.is_finite());
        assert!(d1.is_finite());
    }
    #[test]
    fn test_bounded_proxy_outside_bounds() {
        let inner = SdfSphere::new(1.0);
        let proxy = SdfBoundedProxy::new(inner, [0.0, 0.0, 0.0], 1.5);
        let d = proxy.dist([10.0, 0.0, 0.0]);
        assert!(d > 0.0);
        assert!(approx_eq(d, 10.0 - 1.5, 0.01));
    }
    #[test]
    fn test_sphere_mean_curvature() {
        let s = SdfSphere::new(2.0);
        let p = [2.0, 0.0, 0.0];
        let h = sdf_mean_curvature(&s, p, 1e-4);
        assert!(h.is_finite(), "curvature should be finite: {h}");
    }
    #[test]
    fn test_revolution_sphere_equivalent() {
        let sphere_profile = |r: f64, y: f64| (r * r + y * y).sqrt() - 1.0;
        let rev = SdfRevolution::new(sphere_profile);
        let s = SdfSphere::new(1.0);
        let p = [0.6, 0.5, 0.3];
        assert!(approx_eq(rev.dist(p), s.dist(p), 1e-10));
    }
    #[test]
    fn test_extrude_circle_gives_cylinder() {
        let circle_profile = |xz: [f64; 2]| (xz[0] * xz[0] + xz[1] * xz[1]).sqrt() - 1.0;
        let cyl = SdfExtrude::new(circle_profile, 2.0);
        assert!(cyl.dist([0.0, 0.0, 0.0]) < 0.0);
        assert!(cyl.dist([2.0, 0.0, 0.0]) > 0.0);
        assert!(cyl.dist([0.0, 3.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_hex_prism_center_inside() {
        let h = SdfHexagonalPrism::new(1.0, 1.0);
        assert!(h.dist([0.0, 0.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_closest_point_sphere() {
        let s = SdfSphere::new(1.0);
        let (_q, d) = sdf_closest_point(&s, [0.5, 0.0, 0.0], 1e-5, 50);
        assert!(d.abs() < 1e-3, "final dist={d}");
    }
    #[test]
    fn test_ao_near_plane_unoccluded() {
        let s = SdfSphere::new(1.0);
        let n = [0.0, 1.0, 0.0];
        let p = [0.0, 1.0, 0.0];
        let ao = sdf_ambient_occlusion(&s, p, n, 5, 1.0);
        assert!(ao.is_finite());
        assert!((0.0..=1.0).contains(&ao));
    }
    #[test]
    fn test_soft_shadow_unobstructed() {
        let s = SdfSphere::new(0.5);
        let shadow = sdf_soft_shadow(&s, [-3.0, 0.0, 0.0], [3.0, 0.0, 0.0], 8.0, 1e-4);
        assert!(shadow < 1.0, "shadow={shadow}");
    }
}
