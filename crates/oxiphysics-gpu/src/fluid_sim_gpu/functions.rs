//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    FlipParticle, GpuBoundaryBox, LbmCellType, LbmD2Q9, MacGrid, SphConfig, SphKernels, SphParticle,
};

pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
pub(super) fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Compute SPH density for all particles (mock GPU kernel dispatch).
///
/// For each particle i, density_i = Σ_j m_j * W_poly6(|r_i - r_j|, h)
pub fn sph_compute_density(particles: &mut [SphParticle], config: &SphConfig) {
    let n = particles.len();
    let mut densities = vec![0.0f64; n];
    for i in 0..n {
        let mut rho = 0.0;
        for j in 0..n {
            let r_vec = sub3(particles[i].position, particles[j].position);
            let r = length3(r_vec);
            rho += particles[j].mass * SphKernels::poly6(r, config.h);
        }
        densities[i] = rho.max(1e-6);
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.density = densities[i];
    }
}
/// Compute SPH pressure from density (Tait equation).
pub fn sph_compute_pressure(particles: &mut [SphParticle], config: &SphConfig) {
    for p in particles.iter_mut() {
        let ratio = p.density / config.rest_density;
        p.pressure = config.pressure_k * (ratio.powi(7) - 1.0);
    }
}
/// Compute SPH forces: pressure gradient + viscosity + gravity + surface tension.
pub fn sph_compute_forces(particles: &mut [SphParticle], config: &SphConfig) {
    let n = particles.len();
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let mut f_pressure = [0.0f64; 3];
        let mut f_viscosity = [0.0f64; 3];
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_vec = sub3(particles[i].position, particles[j].position);
            let r = length3(r_vec);
            if r > config.h || r < 1e-15 {
                continue;
            }
            let pressure_factor = particles[i].pressure
                / (particles[i].density * particles[i].density)
                + particles[j].pressure / (particles[j].density * particles[j].density);
            let grad = SphKernels::spiky_grad(r_vec, r, config.h);
            f_pressure = add3(
                f_pressure,
                scale3(grad, -particles[j].mass * pressure_factor),
            );
            let v_diff = sub3(particles[j].velocity, particles[i].velocity);
            let lap = SphKernels::viscosity_laplacian(r, config.h);
            let visc_factor = config.viscosity * particles[j].mass * lap / particles[j].density;
            f_viscosity = add3(f_viscosity, scale3(v_diff, visc_factor));
        }
        let f_gravity = scale3(config.gravity, particles[i].density);
        forces[i] = add3(add3(f_pressure, f_viscosity), f_gravity);
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.force = forces[i];
    }
}
/// Integrate SPH particles using semi-implicit Euler.
pub fn sph_integrate(particles: &mut [SphParticle], config: &SphConfig) {
    for p in particles.iter_mut() {
        let accel = scale3(p.force, 1.0 / p.density.max(1e-6));
        p.velocity = add3(p.velocity, scale3(accel, config.dt));
        p.position = add3(p.position, scale3(p.velocity, config.dt));
    }
}
/// Full SPH step: density → pressure → forces → integrate.
pub fn sph_step(particles: &mut [SphParticle], config: &SphConfig) {
    sph_compute_density(particles, config);
    sph_compute_pressure(particles, config);
    sph_compute_forces(particles, config);
    sph_integrate(particles, config);
}
/// LBM D3Q27 velocity set index count.
pub const D3Q27_Q: usize = 27;
/// LBM D2Q9 velocity set index count.
pub const D2Q9_Q: usize = 9;
/// D2Q9 x-component velocity vectors.
pub const D2Q9_CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 y-component velocity vectors.
pub const D2Q9_CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// D2Q9 equilibrium weights.
pub const D2Q9_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D2Q9 opposite direction indices (for bounce-back).
pub const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/// Semi-Lagrangian advection of a scalar field `phi` on a MAC grid.
///
/// Uses RK2 (midpoint method) for back-tracing.
pub fn advect_scalar(phi: &[f64], grid: &MacGrid, dt: f64, gravity: [f64; 3]) -> Vec<f64> {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let dx = grid.dx;
    let mut phi_new = vec![0.0f64; nx * ny * nz];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let xc = (i as f64 + 0.5) * dx;
                let yc = (j as f64 + 0.5) * dx;
                let zc = (k as f64 + 0.5) * dx;
                let u0 = grid.interp_u(xc, yc, zc);
                let x_mid = xc - 0.5 * dt * (u0 + gravity[0]);
                let y_mid = yc - 0.5 * dt * gravity[1];
                let z_mid = zc - 0.5 * dt * gravity[2];
                let u_mid = grid.interp_u(x_mid, y_mid, z_mid);
                let x_back = xc - dt * (u_mid + gravity[0]);
                let y_back = yc - dt * gravity[1];
                let z_back = zc - dt * gravity[2];
                let ix = (x_back / dx - 0.5).floor() as isize;
                let iy = (y_back / dx - 0.5).floor() as isize;
                let iz = (z_back / dx - 0.5).floor() as isize;
                let ii = ix.clamp(0, nx as isize - 1) as usize;
                let jj = iy.clamp(0, ny as isize - 1) as usize;
                let kk = iz.clamp(0, nz as isize - 1) as usize;
                phi_new[k * nx * ny + j * nx + i] = phi[kk * nx * ny + jj * nx + ii];
            }
        }
    }
    phi_new
}
/// Compute curl (vorticity) of the velocity field at cell centers.
///
/// Returns a flat buffer of 3-vectors (ωx, ωy, ωz).
pub fn compute_vorticity(grid: &MacGrid) -> Vec<[f64; 3]> {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let inv2dx = 1.0 / (2.0 * grid.dx);
    let mut vorticity = vec![[0.0f64; 3]; nx * ny * nz];
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dwdy = (grid.get_w(i, j + 1, k) - grid.get_w(i, j - 1, k)) * inv2dx;
                let dvdz = (grid.get_v(i, j, k + 1) - grid.get_v(i, j, k - 1)) * inv2dx;
                let dudz = (grid.get_u(i, j, k + 1) - grid.get_u(i, j, k - 1)) * inv2dx;
                let dwdx = (grid.get_w(i + 1, j, k) - grid.get_w(i - 1, j, k)) * inv2dx;
                let dvdx = (grid.get_v(i + 1, j, k) - grid.get_v(i - 1, j, k)) * inv2dx;
                let dudy = (grid.get_u(i, j + 1, k) - grid.get_u(i, j - 1, k)) * inv2dx;
                let idx = k * nx * ny + j * nx + i;
                vorticity[idx] = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];
            }
        }
    }
    vorticity
}
/// Apply vorticity confinement forces to MAC grid velocities.
///
/// Vorticity confinement adds a force F = ε * (N × ω)
/// where N is the normalized gradient of |ω|.
pub fn vorticity_confinement(grid: &mut MacGrid, vorticity: &[[f64; 3]], epsilon: f64, dt: f64) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let inv2dx = 1.0 / (2.0 * grid.dx);
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let mag = |x: usize, y: usize, z: usize| {
                    let idx = z * nx * ny + y * nx + x;
                    length3(vorticity[idx])
                };
                let gx = (mag(i + 1, j, k) - mag(i - 1, j, k)) * inv2dx;
                let gy = (mag(i, j + 1, k) - mag(i, j - 1, k)) * inv2dx;
                let gz = (mag(i, j, k + 1) - mag(i, j, k - 1)) * inv2dx;
                let grad_len = (gx * gx + gy * gy + gz * gz).sqrt();
                if grad_len < 1e-12 {
                    continue;
                }
                let n_vec = [gx / grad_len, gy / grad_len, gz / grad_len];
                let idx = k * nx * ny + j * nx + i;
                let omega = vorticity[idx];
                let force = cross3(n_vec, omega);
                let force = scale3(force, epsilon);
                if i < nx {
                    let u = grid.get_u(i, j, k);
                    grid.set_u(i, j, k, u + dt * force[0]);
                }
                if j < ny {
                    let v = grid.get_v(i, j, k);
                    grid.set_v(i, j, k, v + dt * force[1]);
                }
                if k < nz {
                    let w = grid.get_w(i, j, k);
                    grid.set_w(i, j, k, w + dt * force[2]);
                }
            }
        }
    }
}
/// Surface tension force via Continuum Surface Force (CSF) method.
///
/// Requires a level-set function `phi` (negative inside fluid, positive outside).
/// Returns per-cell force vectors.
pub fn surface_tension_csf(
    phi: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    sigma: f64,
) -> Vec<[f64; 3]> {
    let inv_dx = 1.0 / dx;
    let inv2dx = 0.5 * inv_dx;
    let n = nx * ny * nz;
    let mut forces = vec![[0.0f64; 3]; n];
    let idx = |i: usize, j: usize, k: usize| k * nx * ny + j * nx + i;
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let c = idx(i, j, k);
                let gx = (phi[idx(i + 1, j, k)] - phi[idx(i - 1, j, k)]) * inv2dx;
                let gy = (phi[idx(i, j + 1, k)] - phi[idx(i, j - 1, k)]) * inv2dx;
                let gz = (phi[idx(i, j, k + 1)] - phi[idx(i, j, k - 1)]) * inv2dx;
                let grad_mag = (gx * gx + gy * gy + gz * gz).sqrt();
                if grad_mag < 1e-12 {
                    continue;
                }
                let nx_n = gx / grad_mag;
                let ny_n = gy / grad_mag;
                let nz_n = gz / grad_mag;
                let phi_xx = (phi[idx(i + 1, j, k)] - 2.0 * phi[c] + phi[idx(i - 1, j, k)])
                    * inv_dx
                    * inv_dx;
                let phi_yy = (phi[idx(i, j + 1, k)] - 2.0 * phi[c] + phi[idx(i, j - 1, k)])
                    * inv_dx
                    * inv_dx;
                let phi_zz = (phi[idx(i, j, k + 1)] - 2.0 * phi[c] + phi[idx(i, j, k - 1)])
                    * inv_dx
                    * inv_dx;
                let kappa = -(phi_xx + phi_yy + phi_zz) / grad_mag;
                let f_scale = sigma * kappa;
                forces[c] = [f_scale * nx_n, f_scale * ny_n, f_scale * nz_n];
            }
        }
    }
    forces
}
/// Particle-to-Grid (P2G) transfer: splat particle velocities onto MAC grid.
///
/// Uses trilinear weighting.
pub fn p2g_transfer(particles: &[FlipParticle], grid: &mut MacGrid) {
    let dx = grid.dx;
    let inv_dx = 1.0 / dx;
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let mut u_weight = vec![0.0f64; (nx + 1) * ny * nz];
    let mut v_weight = vec![0.0f64; nx * (ny + 1) * nz];
    let mut w_weight = vec![0.0f64; nx * ny * (nz + 1)];
    let mut u_num = vec![0.0f64; (nx + 1) * ny * nz];
    let mut v_num = vec![0.0f64; nx * (ny + 1) * nz];
    let mut w_num = vec![0.0f64; nx * ny * (nz + 1)];
    for p in particles {
        let [px, py, pz] = p.position;
        let iu = (px * inv_dx).floor() as isize;
        let ju = (py * inv_dx - 0.5).floor() as isize;
        let ku = (pz * inv_dx - 0.5).floor() as isize;
        let fxu = px * inv_dx - iu as f64;
        let fyu = py * inv_dx - 0.5 - ju as f64;
        let fzu = pz * inv_dx - 0.5 - ku as f64;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx_off in 0..2 {
                    let ni = (iu + dx_off as isize).clamp(0, nx as isize) as usize;
                    let nj = (ju + dy as isize).clamp(0, ny as isize - 1) as usize;
                    let nk = (ku + dz as isize).clamp(0, nz as isize - 1) as usize;
                    let wx = if dx_off == 0 { 1.0 - fxu } else { fxu };
                    let wy = if dy == 0 { 1.0 - fyu } else { fyu };
                    let wz = if dz == 0 { 1.0 - fzu } else { fzu };
                    let w = wx * wy * wz;
                    let uidx = nk * (nx + 1) * ny + nj * (nx + 1) + ni;
                    u_num[uidx] += w * p.velocity[0];
                    u_weight[uidx] += w;
                }
            }
        }
        let iv = (px * inv_dx - 0.5).floor() as isize;
        let jv = (py * inv_dx).floor() as isize;
        let kv = (pz * inv_dx - 0.5).floor() as isize;
        let fxv = px * inv_dx - 0.5 - iv as f64;
        let fyv = py * inv_dx - jv as f64;
        let fzv = pz * inv_dx - 0.5 - kv as f64;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx_off in 0..2 {
                    let ni = (iv + dx_off as isize).clamp(0, nx as isize - 1) as usize;
                    let nj = (jv + dy as isize).clamp(0, ny as isize) as usize;
                    let nk = (kv + dz as isize).clamp(0, nz as isize - 1) as usize;
                    let wx = if dx_off == 0 { 1.0 - fxv } else { fxv };
                    let wy = if dy == 0 { 1.0 - fyv } else { fyv };
                    let wz = if dz == 0 { 1.0 - fzv } else { fzv };
                    let wt = wx * wy * wz;
                    let vidx = nk * nx * (ny + 1) + nj * nx + ni;
                    v_num[vidx] += wt * p.velocity[1];
                    v_weight[vidx] += wt;
                }
            }
        }
        let iw = (px * inv_dx - 0.5).floor() as isize;
        let jw = (py * inv_dx - 0.5).floor() as isize;
        let kw = (pz * inv_dx).floor() as isize;
        let fxw = px * inv_dx - 0.5 - iw as f64;
        let fyw = py * inv_dx - 0.5 - jw as f64;
        let fzw = pz * inv_dx - kw as f64;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx_off in 0..2 {
                    let ni = (iw + dx_off as isize).clamp(0, nx as isize - 1) as usize;
                    let nj = (jw + dy as isize).clamp(0, ny as isize - 1) as usize;
                    let nk = (kw + dz as isize).clamp(0, nz as isize) as usize;
                    let wx = if dx_off == 0 { 1.0 - fxw } else { fxw };
                    let wy = if dy == 0 { 1.0 - fyw } else { fyw };
                    let wz = if dz == 0 { 1.0 - fzw } else { fzw };
                    let wt = wx * wy * wz;
                    let widx = nk * nx * ny + nj * nx + ni;
                    w_num[widx] += wt * p.velocity[2];
                    w_weight[widx] += wt;
                }
            }
        }
    }
    for i in 0..u_num.len() {
        grid.u[i] = if u_weight[i] > 1e-15 {
            u_num[i] / u_weight[i]
        } else {
            0.0
        };
    }
    for i in 0..v_num.len() {
        grid.v[i] = if v_weight[i] > 1e-15 {
            v_num[i] / v_weight[i]
        } else {
            0.0
        };
    }
    for i in 0..w_num.len() {
        grid.w[i] = if w_weight[i] > 1e-15 {
            w_num[i] / w_weight[i]
        } else {
            0.0
        };
    }
}
/// Grid-to-Particle (G2P) transfer: update particle velocities from MAC grid.
///
/// `flip_ratio` in \[0,1\]: 1.0 = pure FLIP, 0.0 = pure PIC.
pub fn g2p_transfer(
    particles: &mut [FlipParticle],
    grid_new: &MacGrid,
    grid_old: &MacGrid,
    flip_ratio: f64,
) {
    for p in particles.iter_mut() {
        let [px, py, pz] = p.position;
        let dx = grid_new.dx;
        let u_new = grid_new.interp_u(px, py, pz);
        let u_old = grid_old.interp_u(px, py, pz);
        let pic_vel = [u_new, 0.0, 0.0];
        let flip_delta = [u_new - u_old, 0.0, 0.0];
        let flip_vel = add3(p.velocity, flip_delta);
        p.velocity = [
            flip_ratio * flip_vel[0] + (1.0 - flip_ratio) * pic_vel[0],
            p.velocity[1],
            p.velocity[2],
        ];
        p.position = add3(p.position, scale3(p.velocity, dx));
    }
}
/// GPU-accelerated SPH density summation.
///
/// Simulates a GPU kernel where each "thread" computes the density contribution
/// from all neighbours within the smoothing radius h. In a real GPU
/// implementation this would be launched as one thread per particle pair with
/// atomic accumulation. Here we use a parallel-style double loop and write
/// results into a temporary buffer first to keep the interface identical to a
/// real GPU dispatch.
pub fn gpu_sph_density_parallel(particles: &mut [SphParticle], config: &SphConfig) {
    let n = particles.len();
    let mut densities = vec![0.0f64; n];
    for i in 0..n {
        let mut rho = 0.0;
        for j in 0..n {
            let r_vec = sub3(particles[i].position, particles[j].position);
            let r = length3(r_vec);
            rho += particles[j].mass * SphKernels::poly6(r, config.h);
        }
        densities[i] = rho.max(1e-6);
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.density = densities[i];
    }
}
/// GPU Jacobi pressure solver on a MAC grid.
///
/// Identical algorithm to `MacGrid::jacobi_pressure_solve` but expressed as a
/// function that could be launched as a compute shader (one thread per cell).
/// The indirection through a separate function makes the GPU dispatch boundary
/// explicit.
pub fn gpu_jacobi_pressure_solve(grid: &mut MacGrid, rho: f64, dt: f64, iterations: usize) {
    grid.compute_divergence();
    let scale = rho * grid.dx * grid.dx / dt;
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let mut p_new = grid.p.clone();
    for _ in 0..iterations {
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;
                    if grid.flags[idx] != 1 {
                        continue;
                    }
                    let mut nb_sum = 0.0;
                    let mut nb_cnt = 0u32;
                    macro_rules! nb {
                        ($ii:expr, $jj:expr, $kk:expr) => {{
                            nb_sum += grid.p[$kk * nx * ny + $jj * nx + $ii];
                            nb_cnt += 1;
                        }};
                    }
                    if i + 1 < nx {
                        nb!(i + 1, j, k);
                    }
                    if i > 0 {
                        nb!(i - 1, j, k);
                    }
                    if j + 1 < ny {
                        nb!(i, j + 1, k);
                    }
                    if j > 0 {
                        nb!(i, j - 1, k);
                    }
                    if k + 1 < nz {
                        nb!(i, j, k + 1);
                    }
                    if k > 0 {
                        nb!(i, j, k - 1);
                    }
                    if nb_cnt > 0 {
                        p_new[idx] = (nb_sum - scale * grid.div[idx]) / nb_cnt as f64;
                    }
                }
            }
        }
        grid.p.copy_from_slice(&p_new);
    }
}
/// GPU-style BGK collision step for the D2Q9 LBM grid.
///
/// In a real GPU implementation each cell maps to one thread.  The function
/// signature mirrors a compute-shader kernel: it takes the whole grid and
/// applies BGK in-place.
pub fn gpu_lbm_bgk_collide(lbm: &mut LbmD2Q9) {
    let nx = lbm.nx;
    let ny = lbm.ny;
    let inv_tau = lbm.inv_tau;
    let mut updates: Vec<(usize, usize, usize, f64)> = Vec::with_capacity(nx * ny * D2Q9_Q);
    for y in 0..ny {
        for x in 0..nx {
            if lbm.cell_type[y * nx + x] == LbmCellType::Solid {
                continue;
            }
            let rho = lbm.density(x, y);
            let [ux, uy] = lbm.velocity(x, y);
            for q in 0..D2Q9_Q {
                let f_eq = LbmD2Q9::f_equilibrium(rho, ux, uy, q);
                let f_old = lbm.get_f(x, y, q);
                let f_new = f_old - inv_tau * (f_old - f_eq);
                updates.push((x, y, q, f_new));
            }
        }
    }
    for (x, y, q, val) in updates {
        lbm.set_f(x, y, q, val);
    }
}
/// Expand a 10-bit integer into a 30-bit Morton code component (bit interleave).
pub fn morton_expand_bits(mut v: u32) -> u32 {
    v &= 0x000003ff;
    v = (v | (v << 16)) & 0xff0000ff;
    v = (v | (v << 8)) & 0x0300f00f;
    v = (v | (v << 4)) & 0x030c30c3;
    v = (v | (v << 2)) & 0x09249249;
    v
}
/// Encode 3D integer coordinates into a 30-bit Morton (Z-order curve) code.
pub fn morton_encode_3d(x: u32, y: u32, z: u32) -> u32 {
    morton_expand_bits(x) | (morton_expand_bits(y) << 1) | (morton_expand_bits(z) << 2)
}
/// Sort SPH particles by their Morton code for cache-friendly GPU access.
///
/// Positions are normalised into the unit cube `[0, domain_size]` and then
/// quantised to 10 bits per axis (1024 levels) before Morton encoding.
pub fn morton_sort_particles(particles: &mut Vec<SphParticle>, domain_size: [f64; 3]) {
    let bits = 1024u32;
    let mut keyed: Vec<(u32, SphParticle)> = particles
        .drain(..)
        .map(|p| {
            let xi = ((p.position[0] / domain_size[0]).clamp(0.0, 1.0) * (bits - 1) as f64) as u32;
            let yi = ((p.position[1] / domain_size[1]).clamp(0.0, 1.0) * (bits - 1) as f64) as u32;
            let zi = ((p.position[2] / domain_size[2]).clamp(0.0, 1.0) * (bits - 1) as f64) as u32;
            (morton_encode_3d(xi, yi, zi), p)
        })
        .collect();
    keyed.sort_by_key(|(code, _)| *code);
    *particles = keyed.into_iter().map(|(_, p)| p).collect();
}
/// GPU Euler integration: v += a*dt, x += v*dt.
///
/// Maps to one GPU thread per particle.
pub fn gpu_particle_integrate_euler(particles: &mut [SphParticle], dt: f64) {
    for p in particles.iter_mut() {
        let inv_rho = 1.0 / p.density.max(1e-6);
        for d in 0..3 {
            let accel = p.force[d] * inv_rho;
            p.velocity[d] += accel * dt;
            p.position[d] += p.velocity[d] * dt;
        }
    }
}
/// GPU Verlet integration using the previous time-step `dt_prev`.
///
/// x_new = x + v*dt + 0.5*a*dt^2 (Störmer–Verlet, velocity-explicit variant).
pub fn gpu_particle_integrate_verlet(particles: &mut [SphParticle], dt: f64, _dt_prev: f64) {
    for p in particles.iter_mut() {
        let inv_rho = 1.0 / p.density.max(1e-6);
        for d in 0..3 {
            let accel = p.force[d] * inv_rho;
            p.position[d] += p.velocity[d] * dt + 0.5 * accel * dt * dt;
            p.velocity[d] += accel * dt;
        }
    }
}
/// GPU boundary condition: clamp particles inside an AABB and reflect velocities.
///
/// One GPU thread per particle; no inter-thread communication needed.
pub fn gpu_apply_boundary_box(particles: &mut [SphParticle], bounds: &GpuBoundaryBox) {
    for p in particles.iter_mut() {
        for d in 0..3 {
            if p.position[d] < bounds.min[d] {
                p.position[d] = bounds.min[d];
                p.velocity[d] = p.velocity[d].abs() * bounds.restitution;
            }
            if p.position[d] > bounds.max[d] {
                p.position[d] = bounds.max[d];
                p.velocity[d] = -p.velocity[d].abs() * bounds.restitution;
            }
        }
    }
}
/// GPU parallel reduction for total kinetic energy.
///
/// In GPU hardware this would be a tree-reduction across threads in a work-group;
/// here we use an iterator fold which has the same semantics.
///
/// KE = Σ 0.5 * mass * |v|²
pub fn gpu_reduce_kinetic_energy(particles: &[SphParticle], mass: f64) -> f64 {
    particles.iter().fold(0.0, |acc, p| {
        let v2 = p.velocity[0] * p.velocity[0]
            + p.velocity[1] * p.velocity[1]
            + p.velocity[2] * p.velocity[2];
        acc + 0.5 * mass * v2
    })
}
/// GPU parallel reduction for total linear momentum.
///
/// Returns a 3-vector \[px, py, pz\] = Σ mass * v_i.
pub fn gpu_reduce_momentum(particles: &[SphParticle], mass: f64) -> [f64; 3] {
    particles.iter().fold([0.0f64; 3], |acc, p| {
        [
            acc[0] + mass * p.velocity[0],
            acc[1] + mass * p.velocity[1],
            acc[2] + mass * p.velocity[2],
        ]
    })
}
/// GPU semi-Lagrangian scalar field advection on a 2D regular grid.
///
/// Each cell is a separate GPU thread: back-traces the characteristic by dt
/// through the supplied velocity field `vel` (one 2-vector per cell) and
/// samples the field using bilinear interpolation.
pub fn gpu_advect_2d(
    field: &[f64],
    vel: &[[f64; 2]],
    nx: usize,
    ny: usize,
    dx: f64,
    dt: f64,
) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny];
    let inv_dx = 1.0 / dx;
    let sample = |x: f64, y: f64| -> f64 {
        let ix = (x * inv_dx - 0.5).floor() as isize;
        let iy = (y * inv_dx - 0.5).floor() as isize;
        let fx = x * inv_dx - 0.5 - ix as f64;
        let fy = y * inv_dx - 0.5 - iy as f64;
        let ci = |v: isize| v.clamp(0, nx as isize - 1) as usize;
        let cj = |v: isize| v.clamp(0, ny as isize - 1) as usize;
        let f00 = field[cj(iy) * nx + ci(ix)];
        let f10 = field[cj(iy) * nx + ci(ix + 1)];
        let f01 = field[cj(iy + 1) * nx + ci(ix)];
        let f11 = field[cj(iy + 1) * nx + ci(ix + 1)];
        let f0 = f00 + fx * (f10 - f00);
        let f1 = f01 + fx * (f11 - f01);
        f0 + fy * (f1 - f0)
    };
    for j in 0..ny {
        for i in 0..nx {
            let xc = (i as f64 + 0.5) * dx;
            let yc = (j as f64 + 0.5) * dx;
            let idx = j * nx + i;
            let [vx, vy] = vel[idx];
            let xb = xc - vx * dt;
            let yb = yc - vy * dt;
            out[idx] = sample(xb, yb);
        }
    }
    out
}
/// GPU Jacobi solver for the 2D pressure-Poisson equation on a staggered grid.
///
/// Solves ∇²p = rhs (passed as `div`) for `iterations` Jacobi sweeps.
/// The staggered MAC discretisation gives the standard 5-point stencil:
///   p\[i,j\] = (p\[i+1,j\] + p\[i-1,j\] + p\[i,j+1\] + p\[i,j-1\] - dx²*rhs\[i,j\]) / 4
pub fn gpu_pressure_poisson_jacobi_2d(
    pressure: &mut [f64],
    div: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    iterations: usize,
) {
    let dx2 = dx * dx;
    let mut p_new = pressure.to_vec();
    for _ in 0..iterations {
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let mut nb = 0.0;
                let mut cnt = 0u32;
                if i + 1 < nx {
                    nb += pressure[j * nx + i + 1];
                    cnt += 1;
                }
                if i > 0 {
                    nb += pressure[j * nx + i - 1];
                    cnt += 1;
                }
                if j + 1 < ny {
                    nb += pressure[(j + 1) * nx + i];
                    cnt += 1;
                }
                if j > 0 {
                    nb += pressure[(j - 1) * nx + i];
                    cnt += 1;
                }
                if cnt > 0 {
                    p_new[idx] = (nb - dx2 * div[idx]) / cnt as f64;
                }
            }
        }
        pressure.copy_from_slice(&p_new);
    }
}
